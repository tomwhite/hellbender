package org.broadinstitute.hellbender.tools.dataflow.transforms;


import com.google.api.services.genomics.model.Read;
import com.google.cloud.dataflow.sdk.transforms.Combine;
import com.google.cloud.dataflow.sdk.transforms.DoFn;
import com.google.cloud.dataflow.sdk.transforms.Filter;
import com.google.cloud.dataflow.sdk.values.PCollection;
import com.google.common.collect.*;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.metrics.MetricBase;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.Log;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.ArgumentCollectionDefinition;
import org.broadinstitute.hellbender.engine.dataflow.PTransformSAM;
import org.broadinstitute.hellbender.engine.dataflow.SAMSerializableFunction;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.metrics.MetricAccumulationLevel;
import org.broadinstitute.hellbender.tools.picard.analysis.InsertSizeMetrics;
import org.broadinstitute.hellbender.tools.picard.analysis.InsertSizeMetricsCollector;

import java.io.Serializable;
import java.util.EnumSet;
import java.util.Map;
import java.util.Set;

public class InsertSizeMetricsDataflowTransform extends PTransformSAM<InsertSizeMetricsCollector> {

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
        public void validate(){
            if (MINIMUM_PCT < 0 || MINIMUM_PCT > 0.5) {
                throw new UserException.BadArgumentValue("MINIMUM_PCT", "It must be between 0 and 0.5 so all data categories don't get discarded.");
            }
        }
    }
    private static final Log log = Log.getInstance(InsertSizeMetricsDataflowTransform.class);

    @Override
    public PCollection<MetricsFile<InsertSizeMetrics,Integer>>> apply(PCollection<Read> input) {
        input.apply(Filter.by(new SAMSerializableFunction<>(getHeaderString(), isMappedPair)))
            .apply(Combine.globally(new Combine.AccumulatingCombineFn<Integer, DataflowHistogram<Integer>, DataflowHistogram<Integer>>() {
                @Override
                public DataflowHistogram<Integer> createAccumulator() {
                    return new DataflowHistogram<Integer>();
                }
            }))
            .apply(new DoFn<DataflowHistogram<Integer>, MetricsFile<InsertSizeMetrics, Integer>>() {
                @Override
                public void processElement(ProcessContext c) throws Exception {
                    DataflowHistogram<Integer> histo = c.element();
                    MetricsFile file = new MetricsFileDataflow<InsertSizeMetrics, Integer>();
                    c.output(file);
                }
            });

    }


    public static class MetricsFileDataflow<BEAN extends MetricBase & Serializable, HKEY extends Comparable> extends MetricsFile<BEAN, HKEY> {
    }
//
//    public class HistogramAccumulator implements Combine.AccumulatingCombineFn.Accumulator<Integer, HistogramAccumulator , SortedMultiset<Integer>>{
//        public final Multiset<Integer> histogram = HashMultiset.create();
//
//
//        @Override
//        public void addInput(Integer input) {
//            histogram.add(input);
//        }
//
//        @Override
//        public void mergeAccumulator(HistogramAccumulator other) {
//            histogram.addAll(other.histogram);
//        }
//
//        @Override
//        public SortedMultiset<Integer> extractOutput() {
//            SortedMultiset<Integer> output = TreeMultiset.create();
//            output.addAll(histogram);
//            return output;
//        }
//    }


    ReadFilter isMappedPair = r -> r.getReadPairedFlag() &&
            !r.getReadUnmappedFlag() &&
            !r.getMateUnmappedFlag() &&
            !r.getFirstOfPairFlag() &&
            !r.isSecondaryOrSupplementary() &&
            !r.getDuplicateReadFlag() &&
            r.getInferredInsertSize() != 0;




}
