package org.broadinstitute.hellbender.utils.genotyper;

import htsjdk.variant.variantcontext.Allele;

public interface AlleleList<A extends Allele> {

    public int alleleCount();

    public int alleleIndex(final A allele);

    public A alleleAt(final int index);
}
