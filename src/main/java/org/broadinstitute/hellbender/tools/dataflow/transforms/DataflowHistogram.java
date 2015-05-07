package org.broadinstitute.hellbender.tools.dataflow.transforms;

import com.google.cloud.dataflow.sdk.transforms.Combine;
import htsjdk.samtools.util.Histogram;

import java.io.Serializable;

public class DataflowHistogram<K extends Comparable<K>> extends Histogram<K> implements Serializable, Combine.AccumulatingCombineFn.Accumulator<K, DataflowHistogram<K>, Histogram<K>>{

    @Override
    public void addInput(K input) {
        this.increment(input);
    }

    @Override
    public void mergeAccumulator(DataflowHistogram<K> other) {
        this.addHistogram(other);
    }

    @Override
    public Histogram<K> extractOutput() {
        return this;
    }

}

