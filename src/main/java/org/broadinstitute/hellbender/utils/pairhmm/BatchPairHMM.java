package org.broadinstitute.hellbender.utils.pairhmm;

import org.broadinstitute.hellbender.utils.haplotype.Haplotype;

import java.util.List;

public interface BatchPairHMM {
    public void batchAdd(final List<Haplotype> haplotypes,
                         final byte[] readBases,
                         final byte[] readQuals,
                         final byte[] insertionGOP,
                         final byte[] deletionGOP,
                         final byte[] overallGCP);

    public double[] batchGetResult();
}
