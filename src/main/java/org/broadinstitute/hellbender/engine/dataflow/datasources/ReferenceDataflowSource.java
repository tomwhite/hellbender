package org.broadinstitute.hellbender.engine.dataflow.datasources;

import com.google.cloud.dataflow.sdk.options.PipelineOptions;
import com.google.cloud.dataflow.sdk.transforms.SerializableFunction;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.dataflow.BucketUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;

import java.io.IOException;
import java.io.Serializable;

/**
 * Class to load a reference sequence from the Google Genomics API, or a file stored on HDFS or locally.
 *
 * This class needs to be mocked, so it cannot be declared final.
 */
public class ReferenceDataflowSource implements ReferenceSource, Serializable {
    private static final long serialVersionUID = 1L;

    private ReferenceSource referenceSource;
    private SerializableFunction<GATKRead, SimpleInterval> referenceWindowFunction;

    /**
     * @param pipelineOptions the pipeline options; must be GCSOptions if using the Google Genomics API
     * @param reference the name of the reference (if using the Google Genomics API), or a path to the reference file
     * @param referenceWindowFunction the custom reference window function used to map reads to desired reference bases
     */
    public ReferenceDataflowSource(final PipelineOptions pipelineOptions, final String reference,
                                   final SerializableFunction<GATKRead, SimpleInterval> referenceWindowFunction) {
        if (isFasta(reference)) {
            if (BucketUtils.isHadoopUrl(reference)) {
                referenceSource = new ReferenceHadoopSource(reference);
            } else {
                referenceSource = new ReferenceFileSource(reference);
            }
        } else { // use the Google Genomics API
            referenceSource = new RefAPISource(pipelineOptions, reference);
        }
        this.referenceWindowFunction = referenceWindowFunction;
    }

    private static boolean isFasta(String reference) {
        for (final String ext : ReferenceSequenceFileFactory.FASTA_EXTENSIONS) {
            if (reference.endsWith(ext)) {
                return true;
            }
        }
        return false;
    }

    /**
     * @return the custom reference window function used to map reads to desired reference bases
     */
    public SerializableFunction<GATKRead, SimpleInterval> getReferenceWindowFunction() {
        return referenceWindowFunction;
    }

    /**
     * Return reference bases for the given interval.
     * @param pipelineOptions the pipeline options; must be GCSOptions if using the Google Genomics API
     * @param interval the interval to return reference bases for
     * @return reference bases for the given interval
     */
    public ReferenceBases getReferenceBases(final PipelineOptions pipelineOptions, final SimpleInterval interval) throws IOException {
        return referenceSource.getReferenceBases(pipelineOptions, interval);
    }
}
