package org.broadinstitute.hellbender.dev.pipelines.bqsr;

import com.google.cloud.dataflow.sdk.coders.SerializableCoder;
import com.google.cloud.dataflow.sdk.options.PipelineOptions;
import com.google.cloud.dataflow.sdk.transforms.Combine;
import com.google.cloud.dataflow.sdk.transforms.DoFn;
import com.google.cloud.dataflow.sdk.transforms.ParDo;
import com.google.cloud.dataflow.sdk.transforms.join.CoGbkResult;
import com.google.cloud.dataflow.sdk.transforms.join.CoGroupByKey;
import com.google.cloud.dataflow.sdk.transforms.join.KeyedPCollectionTuple;
import com.google.cloud.dataflow.sdk.util.GcsUtil;
import com.google.cloud.dataflow.sdk.util.gcsfs.GcsPath;
import com.google.cloud.dataflow.sdk.values.KV;
import com.google.cloud.dataflow.sdk.values.PCollection;
import com.google.cloud.dataflow.sdk.values.TupleTag;
import com.google.common.base.Stopwatch;
import com.google.common.collect.ComparisonChain;
import com.google.common.collect.Iterables;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.util.Locatable;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.dev.tools.walkers.bqsr.BaseRecalibrationArgumentCollection;
import org.broadinstitute.hellbender.dev.tools.walkers.bqsr.BaseRecalibratorWorker;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.recalibration.RecalUtils;
import org.broadinstitute.hellbender.tools.recalibration.RecalibrationTables;
import org.broadinstitute.hellbender.tools.walkers.bqsr.RecalibrationEngine;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.dataflow.BucketUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.reference.ReferenceUtils;
import org.broadinstitute.hellbender.utils.test.BaseTest;

import java.io.*;
import java.nio.channels.Channels;
import java.security.GeneralSecurityException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.concurrent.TimeUnit;

/**
 * Base Quality Score Recalibration, phase 1
 * This contains the Dataflow-specific bits.
 */
public final class BaseRecalibratorDataflowUtils implements Serializable {

  private static final Logger logger = LogManager.getLogger(BaseRecalibratorDataflowUtils.class);
  // how many bases until we move on to the next shard
  static final int BLOCKSIZE = 1000000;
  private static final long serialVersionUID = 1L;
  public static final TupleTag<GATKRead> readTag = new TupleTag<>();
  public static final TupleTag<SimpleInterval> intervalTag = new TupleTag<>();
  public static final TupleTag<RecalibrationTables> tablesTag = new TupleTag<>();

  /**
   * Get a single BaseRecalOutput object that contains RecalibrationTables and more, output by phase 1 of BQSR.
   * <p>
   * The reference file (*.fasta) must also have a .dict and a .fasta.fai next to it.
   */
  public static PCollection<BaseRecalOutput> getRecalibrationOutput(SAMFileHeader readsHeader, PCollection<GATKRead> reads, String referenceFileName, BaseRecalibrationArgumentCollection toolArgs, PCollection<SimpleInterval> placesToIgnore) {
    PCollection<RecalibrationTables> recalibrationTables = getRecalibrationTables(readsHeader, reads, referenceFileName, toolArgs, placesToIgnore);
    return addQuantizationInfo(readsHeader, toolArgs, recalibrationTables);
  }

  /**
   * Get a single RecalibrationTables object, output by phase 1 of BQSR.
   * <p>
   * The reference file (*.fasta) must also have a .dict and a .fasta.fai next to it.
   */
  public static PCollection<RecalibrationTables> getRecalibrationTables(SAMFileHeader readsHeader, PCollection<GATKRead> reads, String referenceFileName, BaseRecalibrationArgumentCollection toolArgs, PCollection<SimpleInterval> placesToIgnore) {
    final ReadFilter readFilter = BaseRecalibratorWorker.readFilter(readsHeader);
    PCollection<GATKRead> filteredReads = reads.apply(new DataflowReadFilter(readFilter, readsHeader));
    PCollection<RecalibrationTables> stats = computeBlockStatistics(readsHeader, referenceFileName, toolArgs, groupByBlock(filteredReads, placesToIgnore));
    PCollection<RecalibrationTables> oneStat = aggregateStatistics(stats);
    return oneStat;
  }

  /**
   * Throws an exception if any of the files is missing.
   * (offlineauth can be null if the files are local)
   */
  public static void ensureReferenceIsReadable(final PipelineOptions popts, String filename) throws IOException, GeneralSecurityException {
    for (String fn : getRelatedFiles(filename)) {
      if (BucketUtils.isCloudStorageUrl(fn)) {
        // make sure we can access those files on Google Cloud Storage.
        GcsPath path = GcsPath.fromUri(fn);
        // this will throw if we can't get to the file.
        new GcsUtil.GcsUtilFactory().create(popts).fileSize(path);
      } else {
        // make sure we can access those files on the local filesystem.
        if (!new File(fn).canRead()) {
          if (!new File(fn).exists()) {
            throw new IOException("File not found: " + fn);
          } else {
            throw new IOException("File present, but not readable: " + fn);
          }
        }
      }
    }
  }

  // ---------------------------------------------------------------------------------------------------
  // non-public methods

  /*
  * addQuantizationInfo takes the computed RecalibrationTable and adds the QuantizationInfo and RequestedCovariates objects.
  * We call this triplet "BaseRecalOutput". It contains everything we need from phase 1 to continue onto phase 2 of BQSR.
  */
  private static PCollection<BaseRecalOutput> addQuantizationInfo(SAMFileHeader readsHeader, BaseRecalibrationArgumentCollection toolArgs, PCollection<RecalibrationTables> recal) {
    return recal.apply(ParDo
            .named("addQuantizationInfo")
            .of(new DoFn<RecalibrationTables, BaseRecalOutput>() {
              private static final long serialVersionUID = 1L;
              @Override
              public void processElement(ProcessContext c) {
                RecalibrationTables rt = c.element();
                BaseRecalibratorWorker baseRecalibratorWorker = BaseRecalibratorWorker.fromArgs(readsHeader, toolArgs);
                baseRecalibratorWorker.onTraversalStart(null);
                //BaseRecalOutput ret = new BaseRecalOutput(rt, baseRecalibratorWorker.getQuantizationInfo(rt), baseRecalibratorWorker.getRequestedCovariates());

                // Saving and loading back the report actually changes it. So we have to do it.
                // TODO: Figure out what it changes, and just do that instead of doing the whole rigamarole.
                try {
                    File temp = BaseTest.createTempFile("temp-recalibrationtable-",".tmp");
                    toolArgs.RAC.RECAL_TABLE = new PrintStream(temp);
                    RecalUtils.outputRecalibrationReport(toolArgs.RAC, baseRecalibratorWorker.getQuantizationInfo(rt), rt, baseRecalibratorWorker.getRequestedCovariates(), false);
                    BaseRecalOutput ret = new BaseRecalOutput(temp);
                    c.output(ret);
                } catch (FileNotFoundException e) {
                    throw new GATKException("can't find my own temporary file",e);
                }
              }
            })).setCoder(SerializableCoder.of(BaseRecalOutput.class));
  }

  private static String[] getRelatedFiles(String fastaFilename) {
    return new String[] { fastaFilename, ReferenceUtils.getFastaDictionaryFileName(fastaFilename), ReferenceUtils.getFastaIndexFileName(fastaFilename) };
  }

  // a key to group intervals by (aka sharding)
  private static String posKey(final Locatable loc) {
    return loc.getContig() + (loc.getStart() / BLOCKSIZE);
  }

  // a key to group intervals by (aka sharding)
  private static String posKey(final GATKRead r) {
    int start = r.getStart();
    return posKey(new SimpleInterval(r.getContig(),
            start,
            start));
  }

  // the shard that the last base of the interval lands in,
  // if different from the first base's. Null otherwise.
  private static String posKeyEnd(final Locatable loc) {
    if (loc.getStart()/BLOCKSIZE == (loc.getEnd())/BLOCKSIZE) return null;
    return loc.getContig() + (loc.getEnd() / BLOCKSIZE);
  }

  /**
   * returns union(reads,placestoIgnore).groupBy(x->posKey(x))
   */
  private static PCollection<KV<String, CoGbkResult>> groupByBlock(final PCollection<GATKRead> reads, final PCollection<SimpleInterval> placesToIgnore) {
    // shard reads
    PCollection<KV<String, GATKRead>> shardedReads = reads.apply(ParDo
            .named("shard reads")
            .of(
                    new DoFn<GATKRead, KV<String, GATKRead>>() {
                      private static final long serialVersionUID = 1L;

                      @Override
                      public void processElement( ProcessContext c ) {
                        GATKRead r = c.element();
                        c.output(KV.of(posKey(r), r));
                      }

                    }));
    // send ignores to every shard that overlaps with them (so, possibly more than one).
    // (for now we assume they can't span more than two)
    PCollection<KV<String, SimpleInterval>> shardedIgnore = placesToIgnore.apply(ParDo
            .named("shard known intervals")
            .of(
                    new DoFn<SimpleInterval, KV<String, SimpleInterval>>() {
                      private static final long serialVersionUID = 1L;
                      @Override
                      public void processElement( ProcessContext c ) {
                        SimpleInterval i = c.element();
                        String firstPos = posKey(i);
                        c.output(KV.of(firstPos, i));
                        String secondPos = posKeyEnd(i);
                        if ( null != secondPos ) {
                          c.output(KV.of(secondPos, i));
                        }
                      }
                    }));

    return KeyedPCollectionTuple.of(readTag, shardedReads)
        .and(intervalTag, shardedIgnore)
        .apply(CoGroupByKey.<String>create());
  }

  /**
   * For each input block, compute its statistics.
   * <p>
   * At a high level, just delegate to the CalibrationTablesBuilder.
   * More specifically, this downloads the reference files, sorts the known intervals for each shard,
   * and then delegates to CalibrationTablesBuilder. Then it outputs a log message about timings.
   */
  private static PCollection<RecalibrationTables> computeBlockStatistics(final SAMFileHeader readsHeader, String referenceFileName, final BaseRecalibrationArgumentCollection toolArgs, final PCollection<KV<String, CoGbkResult>> readsAndIgnores) {
    PCollection<RecalibrationTables> ret = readsAndIgnores.apply(ParDo
            .named("computeBlockStatistics")
            .of(new DoFn<KV<String, CoGbkResult>, RecalibrationTables>() {
              private static final long serialVersionUID = 1L;
              CalibrationTablesBuilder ct;
              Stopwatch timer;
              int nBlocks = 0;
              int nReads = 0;

              @Override
              public void startBundle(DoFn<KV<String, CoGbkResult>, RecalibrationTables>.Context c) throws Exception {
                timer = Stopwatch.createStarted();
                SAMFileHeader header = readsHeader;

                String localReference = referenceFileName;
                if (BucketUtils.isCloudStorageUrl(referenceFileName)) {
                  // the reference is on GCS, download all 3 files locally first.
                  boolean first = true;
                  for (String fname : getRelatedFiles(referenceFileName)) {
                    String localName = "reference";
                    int slash = fname.lastIndexOf('/');
                    if (slash >= 0) {
                      localName = fname.substring(slash + 1);
                    }
                    // download reference if necessary
                    if (new File(localName).exists()) {
                      if (first) localReference = localName;
                    } else {
                      try (
                              InputStream in = Channels.newInputStream(new GcsUtil.GcsUtilFactory().create(c.getPipelineOptions()).open(GcsPath.fromUri(fname)));
                              FileOutputStream fout = new FileOutputStream(localName)) {
                        final byte[] buf = new byte[1024 * 1024];
                        int count;
                        while ((count = in.read(buf)) > 0) {
                          fout.write(buf, 0, count);
                        }
                      }
                      if (first) localReference = localName;
                    }
                    first = false;
                  }
                  logger.info(String.format("Done downloading reference files (%s ms)\n", timer.elapsed(TimeUnit.MILLISECONDS)));
                } else if (BucketUtils.isHadoopUrl(referenceFileName)) {
                    // the reference is on HDFS, download all 3 files locally first.
                    boolean first = true;
                    for (String fname : getRelatedFiles(referenceFileName)) {
                        String localName = "reference";
                        int slash = fname.lastIndexOf('/');
                        if (slash >= 0) {
                            localName = fname.substring(slash + 1);
                        }
                        // download reference if necessary
                        if (new File(localName).exists()) {
                            if (first) localReference = localName;
                        } else {
                            try (
                                    InputStream in = BucketUtils.openFile(fname, null);
                                    FileOutputStream fout = new FileOutputStream(localName)) {
                                final byte[] buf = new byte[1024 * 1024];
                                int count;
                                while ((count = in.read(buf)) > 0) {
                                    fout.write(buf, 0, count);
                                }
                            }
                            if (first) localReference = localName;
                        }
                        first = false;
                    }
                    logger.info(String.format("Done downloading reference files (%s ms)\n", timer.elapsed(TimeUnit.MILLISECONDS)));
                }

                ct = new CalibrationTablesBuilder(header, localReference, toolArgs);
              }

              @Override
              public void processElement(ProcessContext c) throws Exception {
                nBlocks++;
                // get the reads
                KV<String, CoGbkResult> e = c.element();
                List<GATKRead> reads = new ArrayList<>();
                Iterable<GATKRead> readsIter = e.getValue().getAll(BaseRecalibratorDataflowUtils.readTag);
                Iterables.addAll(reads, readsIter);
                nReads += reads.size();
                // get the skip intervals
                List<SimpleInterval> skipIntervals = new ArrayList<>();
                Iterables.addAll(skipIntervals, e.getValue().getAll(BaseRecalibratorDataflowUtils.intervalTag));
                Collections.sort(skipIntervals, new Comparator<SimpleInterval>() {
                  @Override
                  public int compare(SimpleInterval o1, SimpleInterval o2) {
                    return ComparisonChain.start()
                            .compare(o1.getContig(), o2.getContig())
                            .compare(o1.getStart(), o2.getStart())
                            .result();
                  }
                });
                // update our statistics
                ct.add(reads, skipIntervals);
              }

              @Override
              public void finishBundle(DoFn<KV<String, CoGbkResult>, RecalibrationTables>.Context c) throws Exception {
                ct.done();
                c.output(ct.getRecalibrationTables());
                logger.info("Finishing a block statistics bundle. It took " + timer.elapsed(TimeUnit.MILLISECONDS) + " ms to process " + nBlocks + " blocks, " + nReads + " reads.");
            timer = null;
          }
        }));

    ret.setCoder(SerializableCoder.of(RecalibrationTables.class));
    return ret;
  }

  /**
   * Merge the statistics from each block. The resulting "collection" contains a single element, with the answer.
   */
  private static PCollection<RecalibrationTables> aggregateStatistics(final PCollection<RecalibrationTables> tables) {
    return tables
        // aggregate
        .apply(Combine.globally(new RecalibrationTablesMerger()))
        // call finalize on the result
        .apply(ParDo
            .named("finalizeRecalTables")
            .of(new DoFn<RecalibrationTables, RecalibrationTables>() {
              private static final long serialVersionUID = 1L;
              @Override
              public void processElement(ProcessContext c) throws Exception {
                RecalibrationTables tables = c.element();
                RecalibrationEngine.finalizeRecalibrationTables(tables);
                c.output(tables);
              }
            }));
  }

}


