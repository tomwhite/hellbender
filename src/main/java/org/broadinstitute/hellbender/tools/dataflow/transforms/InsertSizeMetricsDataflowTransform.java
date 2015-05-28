package org.broadinstitute.hellbender.tools.dataflow.transforms;


import com.google.api.services.genomics.model.Read;
import com.google.cloud.dataflow.sdk.coders.AvroCoder;
import com.google.cloud.dataflow.sdk.coders.DefaultCoder;
import com.google.cloud.dataflow.sdk.transforms.Combine;
import com.google.cloud.dataflow.sdk.transforms.Filter;
import com.google.cloud.dataflow.sdk.transforms.ParDo;
import com.google.cloud.dataflow.sdk.values.KV;
import com.google.cloud.dataflow.sdk.values.PCollection;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamPairUtil;
import htsjdk.samtools.metrics.MetricBase;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.Histogram;
import htsjdk.samtools.util.Log;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.ArgumentCollectionDefinition;
import org.broadinstitute.hellbender.engine.dataflow.DataFlowSAMFn;
import org.broadinstitute.hellbender.engine.dataflow.PTransformSAM;
import org.broadinstitute.hellbender.engine.dataflow.SAMSerializableFunction;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.metrics.MetricAccumulationLevel;
import org.broadinstitute.hellbender.tools.picard.analysis.InsertSizeMetrics;

import java.io.Serializable;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.EnumSet;
import java.util.List;
import java.util.Set;

public class InsertSizeMetricsDataflowTransform extends PTransformSAM<InsertSizeMetricsDataflowTransform.MetricsFileDataflow<InsertSizeMetrics,Integer>> {
    private final Arguments args;

    public static class Arguments implements ArgumentCollectionDefinition, Serializable {
        @Argument(doc = "Generate mean, sd and plots by trimming the data down to MEDIAN + DEVIATIONS*MEDIAN_ABSOLUTE_DEVIATION. " +
                "This is done because insert size data typically includes enough anomalous values from chimeras and other " +
                "artifacts to make the mean and sd grossly misleading regarding the real distribution.")
        public double DEVIATIONS = 10;

        @Argument(shortName = "W", doc = "Explicitly sets the Histogram width, overriding automatic truncation of Histogram tail. " +
                "Also, when calculating mean and standard deviation, only bins <= HISTOGRAM_WIDTH will be included.", optional = true)
        public Integer HISTOGRAM_WIDTH = null;

        @Argument(shortName = "M", doc = "When generating the Histogram, discard any data categories (out of FR, TANDEM, RF) that have fewer than this " +
                "percentage of overall reads. (Range: 0 to 0.5).")
        public float MINIMUM_PCT = 0.05f;

        @Argument(shortName = "LEVEL", doc = "The level(s) at which to accumulate metrics.  ")
        private Set<MetricAccumulationLevel> METRIC_ACCUMULATION_LEVEL = EnumSet.of(MetricAccumulationLevel.ALL_READS);

        @Override
        public void validate() {
            if (MINIMUM_PCT < 0 || MINIMUM_PCT > 0.5) {
                throw new UserException.BadArgumentValue("MINIMUM_PCT", "It must be between 0 and 0.5 so all data categories don't get discarded.");
            }
        }
    }

    private static final Log log = Log.getInstance(InsertSizeMetricsDataflowTransform.class);

    public InsertSizeMetricsDataflowTransform(Arguments args) {
        this.args = args;
    }

    @Override
    public PCollection<MetricsFileDataflow<InsertSizeMetrics, Integer>> apply(PCollection<Read> input) {

        PCollection<Read> filtered = input.apply(Filter.by(new SAMSerializableFunction<>(getHeader(), isSecondInMappedPair))).setName("Filter singletons and first of pair");

        PCollection<KV<Key, Integer>> kvPairs = filtered.apply(ParDo.of(new DataFlowSAMFn<KV<Key, Integer>>(getHeader()) {
            @Override
            protected void apply(SAMRecord read) {
                Integer metric = computeMetric(read);
                List<Key> keys = Key.getKeysForAllAggregationLevels(read, true,true,true,true);
                keys.stream().forEach(k -> output(KV.of(k,metric)));
            }
        })).setName("Calculate metric and key");

        Combine.CombineFn<Integer, DataflowHistogram<Integer>, DataflowHistogram<Integer>> combiner = new DataflowHistogrammer<>();

        PCollection<KV<Key,DataflowHistogram<Integer>>> histograms =   kvPairs.apply(Combine.<Key, Integer,DataflowHistogram<Integer>>perKey(combiner)).setName("Add reads to histograms");
        PCollection<MetricsFileDataflow<InsertSizeMetrics, Integer>> metricsFile = histograms.apply(Combine.globally(new CombineMetricsIntoFile(args.DEVIATIONS, args.HISTOGRAM_WIDTH)))
                //.setCoder(SerializableCoder.of((Class<MetricsFileDataflow<InsertSizeMetrics, Integer>>) new MetricsFileDataflow<InsertSizeMetrics, Integer>().getClass()))
                .setName("Add histograms and metrics to MetricsFile");

        return metricsFile;
    }


    public static class DataflowHistogrammer<K extends Comparable<K>> extends Combine.AccumulatingCombineFn<K, DataflowHistogram<K>, DataflowHistogram<K>>{
        @Override
        public DataflowHistogram<K> createAccumulator() {
            return new DataflowHistogram<>();
        }
    }


    public static class MetricsFileDataflow<BEAN extends MetricBase & Serializable , HKEY extends Comparable> extends MetricsFile<BEAN, HKEY> implements Serializable {
        @Override
        public String toString(){
            StringWriter writer = new StringWriter();
            write(writer);
            return writer.toString();
        }
    }

    public static class CombineMetricsIntoFile
            extends Combine.AccumulatingCombineFn<KV<Key,DataflowHistogram<Integer>>,CombineMetricsIntoFile,MetricsFileDataflow<InsertSizeMetrics,Integer>>
            implements Combine.AccumulatingCombineFn.Accumulator<KV<Key,DataflowHistogram<Integer>>,CombineMetricsIntoFile,MetricsFileDataflow<InsertSizeMetrics,Integer>> {

        private final MetricsFileDataflow<InsertSizeMetrics,Integer> metricsFile;
        private final double DEVIATIONS;
        private final Integer HISTOGRAM_WIDTH;  

        public CombineMetricsIntoFile(double deviations, Integer histogramWidth) {
            metricsFile = new MetricsFileDataflow<>();
            this.DEVIATIONS = deviations;
            this.HISTOGRAM_WIDTH = histogramWidth;
        }


        @Override
        public void addInput(KV<Key, DataflowHistogram<Integer>> input) {
            final DataflowHistogram<Integer> Histogram = input.getValue();
            final Key key = input.getKey();
            final SamPairUtil.PairOrientation pairOrientation = key.orientation;
            final double total = Histogram.getCount();

            // Only include a category if it has a sufficient percentage of the data in it
            if( true /*TODO total > totalInserts * args.MINIMUM_PCT */) {
                final InsertSizeMetrics metrics = new InsertSizeMetrics();
                metrics.SAMPLE = null;
                metrics.LIBRARY = null;
                metrics.READ_GROUP = null;
                metrics.PAIR_ORIENTATION = pairOrientation;
                metrics.READ_PAIRS = (long) total;
                metrics.MAX_INSERT_SIZE = (int) Histogram.getMax();
                metrics.MIN_INSERT_SIZE = (int) Histogram.getMin();
                metrics.MEDIAN_INSERT_SIZE = Histogram.getMedian();
                metrics.MEDIAN_ABSOLUTE_DEVIATION = Histogram.getMedianAbsoluteDeviation();

                final double median = Histogram.getMedian();
                double covered = 0;
                double low = median;
                double high = median;

                while (low >= Histogram.getMin() || high <= Histogram.getMax()) {
                    final htsjdk.samtools.util.Histogram<Integer>.Bin lowBin = Histogram.get((int) low);
                    if (lowBin != null) covered += lowBin.getValue();

                    if (low != high) {
                        final htsjdk.samtools.util.Histogram<Integer>.Bin highBin = Histogram.get((int) high);
                        if (highBin != null) covered += highBin.getValue();
                    }

                    final double percentCovered = covered / total;
                    final int distance = (int) (high - low) + 1;
                    if (percentCovered >= 0.1 && metrics.WIDTH_OF_10_PERCENT == 0)
                        metrics.WIDTH_OF_10_PERCENT = distance;
                    if (percentCovered >= 0.2 && metrics.WIDTH_OF_20_PERCENT == 0)
                        metrics.WIDTH_OF_20_PERCENT = distance;
                    if (percentCovered >= 0.3 && metrics.WIDTH_OF_30_PERCENT == 0)
                        metrics.WIDTH_OF_30_PERCENT = distance;
                    if (percentCovered >= 0.4 && metrics.WIDTH_OF_40_PERCENT == 0)
                        metrics.WIDTH_OF_40_PERCENT = distance;
                    if (percentCovered >= 0.5 && metrics.WIDTH_OF_50_PERCENT == 0)
                        metrics.WIDTH_OF_50_PERCENT = distance;
                    if (percentCovered >= 0.6 && metrics.WIDTH_OF_60_PERCENT == 0)
                        metrics.WIDTH_OF_60_PERCENT = distance;
                    if (percentCovered >= 0.7 && metrics.WIDTH_OF_70_PERCENT == 0)
                        metrics.WIDTH_OF_70_PERCENT = distance;
                    if (percentCovered >= 0.8 && metrics.WIDTH_OF_80_PERCENT == 0)
                        metrics.WIDTH_OF_80_PERCENT = distance;
                    if (percentCovered >= 0.9 && metrics.WIDTH_OF_90_PERCENT == 0)
                        metrics.WIDTH_OF_90_PERCENT = distance;
                    if (percentCovered >= 0.99 && metrics.WIDTH_OF_99_PERCENT == 0)
                        metrics.WIDTH_OF_99_PERCENT = distance;

                    --low;
                    ++high;
                }

                // Trim the Histogram down to get rid of outliers that would make the chart useless.
                final htsjdk.samtools.util.Histogram<Integer> trimmedHisto = Histogram; //alias it
                int actualWidth = inferHistogramWidth(HISTOGRAM_WIDTH, metrics, DEVIATIONS);

                trimmedHisto.trimByWidth(actualWidth);

                metrics.MEAN_INSERT_SIZE = trimmedHisto.getMean();
                metrics.STANDARD_DEVIATION = trimmedHisto.getStandardDeviation();

                metricsFile.addHistogram(trimmedHisto);
                metricsFile.addMetric(metrics);
            }
        }

        /**
         * If histogramWidth is null infer a value for it
         */
        private static int inferHistogramWidth(Integer histogramWidth, InsertSizeMetrics metrics, double deviations) {
            if (histogramWidth == null) {
                return (int) (metrics.MEDIAN_INSERT_SIZE + (deviations * metrics.MEDIAN_ABSOLUTE_DEVIATION));
            } else {
                return histogramWidth;
            }
        }

        @Override
        public void mergeAccumulator(CombineMetricsIntoFile other) {
            metricsFile.addAllMetrics(other.metricsFile.getMetrics());
            List<Histogram<Integer>> histograms = other.metricsFile.getAllHistograms();
            histograms.forEach(metricsFile::addHistogram);
        }


        @Override
        public MetricsFileDataflow<InsertSizeMetrics, Integer> extractOutput() {
            return metricsFile;
        }

        @Override
        public CombineMetricsIntoFile createAccumulator() {
            return new CombineMetricsIntoFile(this.DEVIATIONS, this.HISTOGRAM_WIDTH);
        }
    }


    ReadFilter isSecondInMappedPair = r -> r.getReadPairedFlag() &&
            !r.getReadUnmappedFlag() &&
            !r.getMateUnmappedFlag() &&
            !r.getFirstOfPairFlag() &&
            !r.isSecondaryOrSupplementary() &&
            !r.getDuplicateReadFlag() &&
            r.getInferredInsertSize() != 0;


    private Integer computeMetric(SAMRecord read) {
        return Math.abs(read.getInferredInsertSize());
    }

    @DefaultCoder(AvroCoder.class)
    public final static class InsertSizeGroup {
        private SamPairUtil.PairOrientation orientation;

        private InsertSizeGroup(){};

        public static InsertSizeGroup of(SAMRecord read){
            InsertSizeGroup group = new InsertSizeGroup();
            group.orientation =SamPairUtil.getPairOrientation(read);
            return group;
        }
    }


    @DefaultCoder(AvroCoder.class)
    public final static class Key {
        private SamPairUtil.PairOrientation orientation;
        private String readGroup;
        private String library;
        private String sample;

        public Key(){};

        public static Key of(final SAMRecord read,final boolean includeLibrary, final boolean includeReadGroup, final boolean includeSample){
            Key key = new Key();
            key.orientation = SamPairUtil.getPairOrientation(read);
            key.readGroup = includeReadGroup ? read.getReadGroup().getId() : null;
            key.library = includeLibrary ? read.getReadGroup().getLibrary() : null;
            key.sample = includeSample ? read.getReadGroup().getSample(): null;
            return key;
        }

        public static List<Key> getKeysForAllAggregationLevels(final SAMRecord read,final boolean includeAll, final boolean includeLibrary, final boolean includeReadGroup, final boolean includeSample){
            final List<Key> keys = new ArrayList<>();
            if(includeAll) {
                keys.add(Key.of(read, false, false, false));
            }
            if(includeLibrary){
                keys.add(Key.of(read,true, false, false));
            }
            if(includeReadGroup){
                keys.add(Key.of(read, true, true, false));
            }
            if(includeSample){
                keys.add(Key.of(read,true,true,true));
            }
            return keys;

        }

    }


}
