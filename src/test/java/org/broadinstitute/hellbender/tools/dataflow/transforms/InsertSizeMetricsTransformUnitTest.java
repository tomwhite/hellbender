package org.broadinstitute.hellbender.tools.dataflow.transforms;

import com.google.api.services.genomics.model.Read;
import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.coders.SerializableCoder;
import com.google.cloud.dataflow.sdk.testing.TestPipeline;
import com.google.cloud.dataflow.sdk.transforms.Combine;
import com.google.cloud.dataflow.sdk.values.PCollection;
import com.google.cloud.genomics.dataflow.utils.DataflowWorkarounds;
import com.google.common.collect.Lists;
import org.broadinstitute.hellbender.tools.picard.analysis.InsertSizeMetrics;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.dataflow.DataflowUtils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.annotations.Test;

import java.io.File;
import java.util.List;

public final class InsertSizeMetricsTransformUnitTest {

    @Test(groups = "dataflow")
    public void testInsertSizeMetricsTransform(){
        File bam = new File(BaseTest.publicTestDir, "org/broadinstitute/hellbender/tools/picard/analysis/CollectInsertSizeMetrics/insert_size_metrics_test.bam");
        Pipeline p = TestPipeline.create();
        DataflowWorkarounds.registerCoder(p,DataflowHistogram.class, SerializableCoder.of(DataflowHistogram.class) );
        DataflowWorkarounds.registerGenomicsCoders(p);
        List<SimpleInterval> intervals = Lists.newArrayList(new SimpleInterval("1", 1, 249250621));
        PCollection<Read> preads = DataflowUtils.getReadsFromLocalBams(p, intervals, Lists.newArrayList(bam));
        PCollection<InsertSizeMetricsDataflowTransform.MetricsFileDataflow<InsertSizeMetrics,Integer>> presult = preads.apply(new InsertSizeMetricsDataflowTransform(new InsertSizeMetricsDataflowTransform.Arguments()));

        p.run();

    }


    @Test(groups = "dataflow")
    public void testHistogramCombiner(){
        Combine.CombineFn combiner = new InsertSizeMetricsDataflowTransform.CombineMetricsIntoFile(10.0, null);
        List<DataflowHistogram<Integer>> histograms = Lists.newArrayList(new DataflowHistogram<Integer>(), new DataflowHistogram<Integer>(), new DataflowHistogram<Integer>());
        combiner.apply(histograms);
    }

}
