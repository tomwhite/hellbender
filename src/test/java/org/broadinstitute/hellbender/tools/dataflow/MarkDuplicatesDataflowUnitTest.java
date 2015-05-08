package org.broadinstitute.hellbender.tools.dataflow;

import com.google.cloud.dataflow.sdk.Pipeline;
import com.google.cloud.dataflow.sdk.testing.TestPipeline;
import com.google.cloud.genomics.dataflow.utils.DataflowWorkarounds;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.junit.Test;

public class MarkDuplicatesDataflowUnitTest extends BaseTest {

    @Test
    public void test1(){
        final MarkDuplicatesDataflow mddf = new MarkDuplicatesDataflow();

        Pipeline p = TestPipeline.create();
        DataflowWorkarounds.registerGenomicsCoders(p);
    }
}
