package org.broadinstitute.hellbender.tools.dataflow;

import com.google.api.client.json.GenericJson;
import com.google.api.client.repackaged.com.google.common.annotations.VisibleForTesting;
import com.google.api.services.genomics.model.Read;
import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.coders.DefaultCoder;
import com.google.cloud.dataflow.sdk.coders.KvCoder;
import com.google.cloud.dataflow.sdk.coders.SerializableCoder;
import com.google.cloud.dataflow.sdk.coders.StringUtf8Coder;
import com.google.cloud.dataflow.sdk.io.TextIO;
import com.google.cloud.dataflow.sdk.transforms.*;
import com.google.cloud.dataflow.sdk.transforms.Partition.PartitionFn;
import com.google.cloud.dataflow.sdk.values.KV;
import com.google.cloud.dataflow.sdk.values.PCollection;
import com.google.cloud.dataflow.sdk.values.PCollectionList;
import com.google.cloud.dataflow.sdk.values.PCollectionView;
import com.google.cloud.genomics.dataflow.coders.GenericJsonCoder;
import com.google.common.collect.*;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.ArgumentCollection;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.IntervalArgumentCollection;
import org.broadinstitute.hellbender.cmdline.argumentcollections.OptionalIntervalArgumentCollection;
import org.broadinstitute.hellbender.cmdline.programgroups.DataFlowProgramGroup;
import org.broadinstitute.hellbender.engine.dataflow.DataflowCommandLineProgram;
import org.broadinstitute.hellbender.engine.dataflow.ReadsSource;
import org.broadinstitute.hellbender.tools.dataflow.transforms.MarkDuplicatesReadsKey;
import org.broadinstitute.hellbender.utils.GenomeLocSortedSet;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.dataflow.DataflowUtils;

import java.io.Serializable;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.StreamSupport;

import static org.broadinstitute.hellbender.tools.dataflow.GenomicsReadUtils.getAlignmentStart;
import static org.broadinstitute.hellbender.tools.dataflow.GenomicsReadUtils.getReferenceName;
import static org.broadinstitute.hellbender.tools.dataflow.GenomicsReadUtils.isPrimaryAlignment;

@CommandLineProgramProperties(
        usage="Marks duplicates on dataflow",
        usageShort="Mark Duplicates",
        programGroup = DataFlowProgramGroup.class)
public class MarkDuplicatesDataflow extends DataflowCommandLineProgram {

    //Bases below this quality will not be included in picking the best read from a set of duplicates.
    private static final int MIN_BASE_QUAL = 15;

    @Argument(doc="a prefix for the dataflow output files", shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, optional = false)
    protected String outputFile;

    @Argument(doc = "uri for the input bam, either a local file path or a gs:// bucket path",
            shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME, fullName = StandardArgumentDefinitions.INPUT_LONG_NAME,
            optional = false)
    protected String bam;

    @ArgumentCollection
    protected IntervalArgumentCollection intervalArgumentCollection = new OptionalIntervalArgumentCollection();

    @Override
    protected final void setupPipeline(final Pipeline pipeline) {
        final ReadsSource readsSource = new ReadsSource(bam, pipeline);
        final SAMFileHeader header = readsSource.getHeader();
        final SAMSequenceDictionary sequenceDictionary = header.getSequenceDictionary();
        final List<SimpleInterval> intervals = intervalArgumentCollection.intervalsSpecified() ? intervalArgumentCollection.getIntervals(sequenceDictionary):
                getAllIntervalsForReference(sequenceDictionary);

        final PCollectionView<SAMFileHeader> headerPcolView = pipeline.apply(Create.of(header)).apply(View.<SAMFileHeader>asSingleton());

        final PCollection<Read> preads = readsSource.getReadPCollection(intervals);
        final PCollectionList<Read> readsPartitioned = partitionReads(preads);

        final PCollection<Read> fragments = readsPartitioned.get(ReadsPartition.PRIMARY.ordinal());
        final PCollection<Read> fragmentsTransformed = transformFragments(headerPcolView, fragments);

        final PCollection<Read> pairs = readsPartitioned.get(ReadsPartition.PRIMARY.ordinal());
        final PCollection<Read> pairsTransformed = transformReads(headerPcolView, pairs);

        //no work on those
        final PCollection<Read> not_primary = readsPartitioned.get(ReadsPartition.NOT_PRIMARY.ordinal());

        final PCollection<Read> results = PCollectionList.of(fragmentsTransformed).and(pairsTransformed).and(not_primary).apply(Flatten.<Read>pCollections());

        final PCollection<String> pstrings = results.apply(DataflowUtils.convertToString());
        pstrings.apply(TextIO.Write.to(outputFile));
    }

    /**
     * (1) Reads are grouped by read group and read name.
     * (2) Reads are then paired together as follows:
     *   (a) The remaining reads (one per fragment key) are coordinate-sorted and paired
     *       consecutively together and emitted.
     *   (b) If a read is leftover, it is emitted, unmodified, as an unpaired end.
     * (3) Paired ends are grouped by a similar key as above but using both reads.
     * (4) Any unpaired end is emitted, unmodified.
     * (5) The remained paired ends are scored and all but the highest scoring are marked as
     *     duplicates. Both reads in the pair are emitted.
     */
    private PCollection<Read> transformReads(final PCollectionView<SAMFileHeader> headerPcolView, final PCollection<Read> pairs) {
        final PCollection<KV<String, Iterable<Read>>> pairsGroupedByKey = groupPairsByKey(pairs, headerPcolView);
        final PTransform<PCollection<? extends KV<String, Iterable<Read>>>, PCollection<KV<String, PairedEnds>>> markGroupedDuplicatePairs = markGroupedDuplicatePairs(headerPcolView);
        final PTransform<PCollection<? extends KV<String, Iterable<PairedEnds>>>, PCollection<Read>> markPairedEnds = markPairedEnds();

        return pairsGroupedByKey
                .apply(markGroupedDuplicatePairs)
                .setCoder(KvCoder.of(StringUtf8Coder.of(), GenericJsonCoder.of(PairedEnds.class)))
                .apply(GroupByKey.<String, PairedEnds>create())
                .apply(markPairedEnds);
    }

    private PTransform<PCollection<? extends KV<String, Iterable<PairedEnds>>>, PCollection<Read>> markPairedEnds() {
        return ParDo
                .named("mark paired ends")
                .of(new DoFn<KV<String, Iterable<PairedEnds>>, Read>() {
                    private static final long serialVersionUID = 1l;

                    @Override
                    public void processElement(final ProcessContext context) throws Exception {
                        final ImmutableListMultimap<Boolean, PairedEnds> paired = Multimaps.index(context.element().getValue(), pair -> pair.second() != null);

                        // As in Picard, unpaired ends left alone.
                        for (final PairedEnds pair : paired.get(false)) {
                            context.output(pair.first());
                        }

                        final List<PairedEnds> scored = Ordering.natural().reverse().onResultOf((PairedEnds pair) -> pair.score()).immutableSortedCopy(paired.get(true));
                        final PairedEnds best = Iterables.getFirst(scored, null);
                        if (best != null) {
                            context.output(best.first());
                            context.output(best.second());
                        }
                        for (final PairedEnds pair : Iterables.skip(scored, 1)) {
                            Read record = pair.first();
                            record.setDuplicateFragment(true);
                            context.output(record);
                            record = pair.second();
                            record.setDuplicateFragment(true);
                            context.output(record);
                        }
                    }
                });
    }

    private PTransform<PCollection<? extends KV<String, Iterable<Read>>>, PCollection<KV<String, PairedEnds>>> markGroupedDuplicatePairs(final PCollectionView<SAMFileHeader> headerPcolView) {
        final PTransform<PCollection<? extends KV<String, Iterable<Read>>>, PCollection<KV<String, PairedEnds>>> putIntoPairedEnds = ParDo
                .named("pair ends")
                .withSideInputs(Collections.singletonList(headerPcolView))
                .of(new DoFn<KV<String, Iterable<Read>>, KV<String, PairedEnds>>() {
                    @Override
                    public void processElement(final ProcessContext context) throws Exception {
                        final SAMFileHeader header = context.sideInput(headerPcolView);
                        final TreeSet<Read> sorted = Sets.newTreeSet(new CoordinateSorted(header));
                        for (final Read record : context.element().getValue()) {
                            sorted.add(record);
                        }
                        PairedEnds pair = null;
                        for (final Read record : sorted) {
                            if (pair == null) {
                                pair = PairedEnds.of(record);
                            } else {
                                pair.and(record);
                                context.output(KV.of(pair.key(header), pair));
                                pair = null;
                            }
                        }
                        if (pair != null) {
                            context.output(KV.of(pair.key(header), pair));
                        }
                    }
                });
        return putIntoPairedEnds;
    }

    @VisibleForTesting
    /**
     * Groups pairs by keys - keys are tuples of (library, contig, position, orientation).
     */
    PCollection<KV<String, Iterable<Read>>> groupPairsByKey(final PCollection<Read> pairs, final PCollectionView<SAMFileHeader> headerPcolView) {
        final PTransform<PCollection<? extends Read>, PCollection<KV<String, Read>>> makeKeysForPairs =
                ParDo.named("make keys for pairs")
                        .withSideInputs(Collections.singletonList(headerPcolView))
                        .of(new DoFn<Read, KV<String, Read>>() {
                            private static final long serialVersionUID = 1l;

                            @Override
                            public void processElement(final ProcessContext context) throws Exception {
                                final Read record = context.element();
                                if (GenomicsReadUtils.isPaired(record)) {
                                    final SAMFileHeader h = context.sideInput(headerPcolView);
                                    final String key = MarkDuplicatesReadsKey.keyForFragment(h, record);
                                    final KV<String, Read> kv = KV.of(key, record);
                                    context.output(kv);
                                }
                            }
                        });
        return pairs.apply(makeKeysForPairs).apply(GroupByKey.<String, Read>create());
    }

    //--------------


    private PCollectionList<Read> partitionReads(PCollection<Read> preads) {
        return preads.apply(
                    Partition.of(2, new PartitionFn<Read>() {
                        @Override
                        public int partitionFor(final Read read, final int n) {
                            return isPrimaryAlignment(read) ? ReadsPartition.PRIMARY.ordinal() : ReadsPartition.NOT_PRIMARY.ordinal();
                        }
                    }));
    }

    @VisibleForTesting
    /**
     * Takes the reads,
     * group them by library, contig, position and orientation,
     * within each group
     *   (a) if there are only fragments, mark all but the highest scoring as duplicates, or,
     *   (b) if at least one is marked as paired, mark all fragments as duplicates.
     *  Note: Emit only the fragments, as the paired reads are handled separately.
     */
    PCollection<Read> transformFragments(final PCollectionView<SAMFileHeader> headerPcolView, final PCollection<Read> fragments) {
        final PCollection<KV<String, Iterable<Read>>> readsGroupedByKey = groupReadsByKey(fragments, headerPcolView);
        final PTransform<PCollection<? extends KV<String, Iterable<Read>>>, PCollection<Read>> markGroupedDuplicateFragments = markGroupedDuplicateFragments();
        return readsGroupedByKey
                .apply(markGroupedDuplicateFragments);//no need to set up coder for Read (uses GenericJsonCoder)
    }

    private PTransform<PCollection<? extends KV<String, Iterable<Read>>>, PCollection<Read>> markGroupedDuplicateFragments() {
        return ParDo.named("mark dups")
                .of(new DoFn<KV<String, Iterable<Read>>, Read>() {
                    private static final long serialVersionUID = 1l;

                    @Override
                    public void processElement(final ProcessContext context) throws Exception {
                        final Map<Boolean, List<Read>> byPairing = StreamSupport.stream(context.element().getValue().spliterator(), false).collect(Collectors.partitioningBy(
                                read -> GenomicsReadUtils.isPaired(read)
                        ));
                        //
                        // The existence of any paired end means we mark all fragments as duplicates.
                        // Otherwise, mark all but the highest scoring fragment.
                        // Note the we emit only fragments from this mapper.
                        if (byPairing.get(true).isEmpty()) {
                            final List<Read> frags = Ordering.natural().reverse().onResultOf((Read read) -> score(read)).immutableSortedCopy(byPairing.get(false));
                            if (!frags.isEmpty()) {
                                context.output(Iterables.getFirst(frags, null));
                                for (final Read record : Iterables.skip(frags, 1)) {
                                    record.setDuplicateFragment(true);
                                    context.output(record);
                                }
                            }
                        } else {
                            for (final Read record : byPairing.get(false)) {
                                record.setDuplicateFragment(true);
                                context.output(record);
                            }
                        }
                    }
                });
    }

    //Note: copied from htsjdk.samtools.DuplicateScoringStrategy
    private static int score(final Read record) {
        if (record == null) {
            return 0;
        } else {
            return record.getAlignedQuality().stream().mapToInt(Integer::intValue).filter(i -> i >= MIN_BASE_QUAL).sum();
        }
    }

    @VisibleForTesting
    /**
     * Groups reads by keys - keys are tuples of (library, contig, position, orientation).
     */
    PCollection<KV<String, Iterable<Read>>> groupReadsByKey(final PCollection<Read> fragments, final PCollectionView<SAMFileHeader> headerPcolView) {
        final PTransform<PCollection<? extends Read>, PCollection<KV<String, Read>>> makeKeysForFragments =
                ParDo.named("make keys for reads")
                        .withSideInputs(Collections.singletonList(headerPcolView))
                        .of(new DoFn<Read, KV<String, Read>>() {
                            @Override
                            public void processElement(final ProcessContext context) throws Exception {
                                final Read record = context.element();
                                record.setDuplicateFragment(false);
                                final SAMFileHeader h = context.sideInput(headerPcolView);
                                final String key = MarkDuplicatesReadsKey.keyForFragment(h, record);
                                final KV<String, Read> kv = KV.of(key, record);
                                context.output(kv);
                            }
                        });
        return fragments.apply(makeKeysForFragments)
                .apply(GroupByKey.<String, Read>create());
    }

    private enum ReadsPartition {
        PRIMARY, NOT_PRIMARY
    }

    private List<SimpleInterval> getAllIntervalsForReference(final SAMSequenceDictionary sequenceDictionary) {
    return GenomeLocSortedSet.createSetFromSequenceDictionary(sequenceDictionary)
            .stream()
            .map(SimpleInterval::new)
            .collect(Collectors.toList());
    }

//    @DefaultCoder(GenericJsonCoder.class)
    public static final class PairedEnds extends GenericJson{
        private Read first, second;

        PairedEnds(final Read first) {
            this.first = first;
        }

        public static PairedEnds of(final Read first) {
            return new PairedEnds(first);
        }

        public PairedEnds and(final Read second) {
            if (second != null && GenomicsReadUtils.unclippedCoordinate(this.first) > GenomicsReadUtils.unclippedCoordinate(second)) {
                this.second = this.first;
                this.first = second;
            } else {
                this.second = second;
            }
            return this;
        }

        public String key(final SAMFileHeader header) {
            return MarkDuplicatesReadsKey.keyForPairedEnds(header, first, second);
        }

        public Read first() {
            return first;
        }

        public Read second() {
            return second;
        }

        public int score() {
            return MarkDuplicatesDataflow.score(first) + MarkDuplicatesDataflow.score(second);
        }

        @Override
        public PairedEnds clone() {
            return (PairedEnds) super.clone();
        }
    }

    // Necessary since the records that are passed around in this transform have empty headers
    // as it's not necessary to pass all that extra data around per record.
    public static final class CoordinateSorted implements Comparator<Read>, Serializable {
        private static final long serialVersionUID = 1l;
        private final SAMFileHeader header;

        public CoordinateSorted(final SAMFileHeader header) {
            this.header = header;
        }

        @Override
        public int compare(final Read lhs, final Read rhs) {
            final int res1 = Integer.compare(MarkDuplicatesReadsKey.index(header, getReferenceName(lhs)), MarkDuplicatesReadsKey.index(header, getReferenceName(rhs)));
            if (res1 != 0) return res1;

            final int res2 = Long.compare(getAlignmentStart(lhs), getAlignmentStart(rhs));
            if (res2 != 0) return res2;

            final int res3 = Boolean.compare(lhs.getDuplicateFragment(), rhs.getDuplicateFragment());
            if (res3 != 0) return res3;

            final int res4 = Boolean.compare(lhs.getFailedVendorQualityChecks(), rhs.getFailedVendorQualityChecks());
            if (res4 != 0) return res4;

            final int res5 = Integer.compare(lhs.getNumberReads(), rhs.getNumberReads());
            if (res5 != 0) return res5;

            final int res6 = Boolean.compare(lhs.getProperPlacement(), rhs.getProperPlacement());
            if (res6 != 0) return res6;

            final int res7 = Integer.compare(lhs.getReadNumber(), rhs.getReadNumber());
            if (res7 != 0) return res7;

            final int res8 = Boolean.compare(lhs.getSecondaryAlignment(), rhs.getSecondaryAlignment());
            if (res8 != 0) return res8;

            final int res9 = Boolean.compare(lhs.getSupplementaryAlignment(), rhs.getSupplementaryAlignment());
            if (res9 != 0) return res9;

            final int res10 = Integer.compare(GenomicsReadUtils.getMappingQuality(lhs), GenomicsReadUtils.getMappingQuality(rhs));
            if (res10 != 0) return res10;

            final int res11 = Integer.compare(MarkDuplicatesReadsKey.index(header, GenomicsReadUtils.getMateReferenceName(lhs)), MarkDuplicatesReadsKey.index(header, GenomicsReadUtils.getMateReferenceName(rhs)));
            if (res11 != 0) return res11;

            final int res12 = Long.compare(GenomicsReadUtils.getMateAlignmentStart(lhs), GenomicsReadUtils.getMateAlignmentStart(rhs));
            return res12 ;
        }
    }
}