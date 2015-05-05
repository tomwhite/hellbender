package org.broadinstitute.hellbender.tools.dataflow.pipelines;

import org.broadinstitute.hellbender.cmdline.ArgumentCollection;
import org.broadinstitute.hellbender.engine.dataflow.PTransformSAM;
import org.broadinstitute.hellbender.tools.dataflow.transforms.InsertSizeMetricsDataflowTransform;

public class InsertSizeMetricsDataflow extends DataflowReadsPipeline {

    @ArgumentCollection
    InsertSizeMetricsDataflowTransform.Arguments arguments = new InsertSizeMetricsDataflowTransform.Arguments();


    @Override
    protected PTransformSAM<?> getTool() {
        arguments.validate();
        return new InsertSizeMetricsDataflowTransform(arguments);
    }
}
