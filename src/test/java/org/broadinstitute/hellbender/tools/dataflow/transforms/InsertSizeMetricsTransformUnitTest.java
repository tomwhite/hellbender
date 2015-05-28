package org.broadinstitute.hellbender.tools.dataflow.transforms;

import com.google.api.services.genomics.model.Read;
import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.coders.SerializableCoder;
import com.google.cloud.dataflow.sdk.runners.DirectPipelineRunner;
import com.google.cloud.dataflow.sdk.testing.TestPipeline;
import com.google.cloud.dataflow.sdk.transforms.Combine;
import com.google.cloud.dataflow.sdk.transforms.DoFn;
import com.google.cloud.dataflow.sdk.util.SerializableUtils;
import com.google.cloud.dataflow.sdk.values.KV;
import com.google.cloud.dataflow.sdk.values.PCollection;
import com.google.cloud.genomics.dataflow.utils.DataflowWorkarounds;
import com.google.common.collect.Lists;
import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.engine.dataflow.ReadsSource;
import org.broadinstitute.hellbender.tools.picard.analysis.InsertSizeMetrics;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.dataflow.DataflowUtils;
import org.broadinstitute.hellbender.utils.read.ArtificialSAMUtils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.*;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public final class InsertSizeMetricsTransformUnitTest{

    @Test(groups = "dataflow")
    public void testInsertSizeMetricsTransform(){
        File bam = new File(BaseTest.publicTestDir, "org/broadinstitute/hellbender/tools/picard/analysis/CollectInsertSizeMetrics/insert_size_metrics_test.bam");
        Pipeline p = TestPipeline.create();
        p.getCoderRegistry().registerCoder(InsertSizeMetricsDataflowTransform.MetricsFileDataflow.class, SerializableCoder.of(InsertSizeMetricsDataflowTransform.MetricsFileDataflow.class));
        DataflowWorkarounds.registerCoder(p,DataflowHistogram.class, SerializableCoder.of(DataflowHistogram.class) );
        DataflowWorkarounds.registerGenomicsCoders(p);
        List<SimpleInterval> intervals = Lists.newArrayList(new SimpleInterval("1", 1, 249250621));

        ReadsSource source = new ReadsSource(bam.getAbsolutePath(), p);
        PCollection<Read> preads = source.getReadPCollection(intervals);


        InsertSizeMetricsDataflowTransform transform = new InsertSizeMetricsDataflowTransform(new InsertSizeMetricsDataflowTransform.Arguments());
        transform.setHeader(source.getHeader());

        PCollection<InsertSizeMetricsDataflowTransform.MetricsFileDataflow<InsertSizeMetrics,Integer>> presult = preads.apply(transform);
        //presult.apply(ParDo.of(new MetricsFileDataflowBooleanDoFn()));
        DirectPipelineRunner.EvaluationResults result = (DirectPipelineRunner.EvaluationResults)p.run();
        Assert.assertEquals(result.getPCollection(presult).get(0).toString(), "some string");

    }


    @Test
    public void testHistogrammer(){
        List<Integer> records = IntStream.rangeClosed(1,1000).boxed().collect(Collectors.toList());
        Combine.CombineFn<Integer, DataflowHistogram<Integer>, DataflowHistogram<Integer>> combiner = new InsertSizeMetricsDataflowTransform.DataflowHistogrammer<>();

        DataflowHistogram<Integer> result = combiner.apply(records);
        Assert.assertEquals(result.getCount(),1000.0);
        Assert.assertEquals(result.getMax(),1000.0);
        Assert.assertEquals(result.getMin(),1.0);
    }

    @Test(groups = "dataflow")
    public void testHistogramCombiner(){
        Combine.CombineFn<KV<InsertSizeMetricsDataflowTransform.Key, DataflowHistogram<Integer>>,?, InsertSizeMetricsDataflowTransform.MetricsFileDataflow<InsertSizeMetrics,Integer>> combiner = new InsertSizeMetricsDataflowTransform.CombineMetricsIntoFile(10.0, null);
        SAMRecord read1 = ArtificialSAMUtils.createPair(ArtificialSAMUtils.createArtificialSamHeader(), "Read1", 100, 4, 200, true, false).get(0);
        InsertSizeMetricsDataflowTransform.Key key = InsertSizeMetricsDataflowTransform.Key.of(read1, false, false, false);

        DataflowHistogram<Integer> h1= new DataflowHistogram<>();
        h1.addInput(10);
        h1.addInput(20);
        h1.addInput(10);

        DataflowHistogram<Integer> h2= new DataflowHistogram<>();
        h2.addInput(10);
        h2.addInput(20);
        h2.addInput(10);

        DataflowHistogram<Integer> h3= new DataflowHistogram<>();
        h3.addInput(10);
        h3.addInput(20);
        h3.addInput(10);

        List<DataflowHistogram<Integer>> histograms = Lists.newArrayList(h1, h2, h3);

        List<KV<InsertSizeMetricsDataflowTransform.Key,DataflowHistogram<Integer>>> keyedHistograms = histograms.stream().map(h -> KV.of(key, h)).collect(Collectors.toList());
        InsertSizeMetricsDataflowTransform.MetricsFileDataflow<InsertSizeMetrics, Integer> result = combiner.apply(keyedHistograms);
        Assert.assertEquals(result.getAllHistograms().size(), 3);
    }

    private static class MetricsFileDataflowBooleanDoFn extends DoFn<InsertSizeMetricsDataflowTransform.MetricsFileDataflow<InsertSizeMetrics, Integer>, Boolean> {

        @Override
        public void processElement(ProcessContext c) throws Exception {
            Assert.assertEquals(c.element().getAllHistograms().get(0).getCount(), "some string");
        }
    }

    public static class PrintLn<T> extends DoFn<T,String>{

        @Override
        public void processElement(ProcessContext c) throws Exception {
            String str = c.element().toString();
            System.out.println(str);
            c.output(str);
        }
    }

    @Test
    public void serializeMetricsFileTest(){
        InsertSizeMetricsDataflowTransform.MetricsFileDataflow<InsertSizeMetrics,Integer> metrics = new InsertSizeMetricsDataflowTransform.MetricsFileDataflow<>();
        metrics.addHistogram(new DataflowHistogram<>());

        InsertSizeMetricsDataflowTransform.MetricsFileDataflow<InsertSizeMetrics, Integer> newMetrics = SerializableUtils.ensureSerializableByCoder(SerializableCoder.of(InsertSizeMetricsDataflowTransform.MetricsFileDataflow.class), metrics, "error");
        Assert.assertEquals(newMetrics.getAllHistograms(),metrics.getAllHistograms());
    }

    @Test
    public void javaSerializeMetricsFileTest() throws IOException, ClassNotFoundException {
        final InsertSizeMetricsDataflowTransform.MetricsFileDataflow<InsertSizeMetrics,Integer> metrics = new InsertSizeMetricsDataflowTransform.MetricsFileDataflow<>();
        metrics.addHistogram(new DataflowHistogram<>());
        final ByteArrayOutputStream byteArrayStream = new ByteArrayOutputStream();
        final ObjectOutputStream out = new ObjectOutputStream(byteArrayStream);
        out.writeObject(metrics);
        out.close();
        final ObjectInputStream in = new ObjectInputStream(new ByteArrayInputStream(byteArrayStream.toByteArray()));
        final InsertSizeMetricsDataflowTransform.MetricsFileDataflow<InsertSizeMetrics,Integer> deserializedMetrics = (InsertSizeMetricsDataflowTransform.MetricsFileDataflow<InsertSizeMetrics,Integer>)in.readObject();
        in.close();
        Assert.assertEquals(deserializedMetrics.getAllHistograms(),metrics.getAllHistograms());
    }

}
