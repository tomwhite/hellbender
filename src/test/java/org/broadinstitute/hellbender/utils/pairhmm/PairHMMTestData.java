package org.broadinstitute.hellbender.utils.pairhmm;

import htsjdk.samtools.SAMUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.text.XReadLines;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.zip.GZIPInputStream;

/**
 * Useful single class carrying test data for PairHMMs (for use in benchmarking and unit tests)
 */
public class PairHMMTestData {
    public final String ref;
    public final String nextRef;
    private final String read;
    public final byte[] baseQuals, insQuals, delQuals, gcp;
    public final double log10l;
    public final boolean newRead;

    PairHMMTestData(String ref, String nextRef, String read, byte[] baseQuals, byte[] insQuals, byte[] delQuals, byte[] gcp, double log10l, boolean newRead) {
        this.ref = ref;
        this.nextRef = nextRef;
        this.read = read;
        this.baseQuals = baseQuals;
        this.insQuals = insQuals;
        this.delQuals = delQuals;
        this.gcp = gcp;
        this.log10l = log10l;
        this.newRead = newRead;
    }

    PairHMMTestData(String ref, String nextRef, String read, final byte qual) {
        this.ref = ref;
        this.nextRef = nextRef;
        this.read = read;
        this.baseQuals = this.insQuals = this.delQuals =  Utils.dupBytes(qual, read.length());
        this.gcp =  Utils.dupBytes((byte) 10, read.length());
        this.log10l = -1;
        this.newRead = true;
    }

    public double runHMM(final PairHMM hmm) {
        hmm.initialize(getRead().length(), ref.length());
        return hmm.computeReadLikelihoodGivenHaplotypeLog10(ref.getBytes(), getRead().getBytes(),
                baseQuals, insQuals, delQuals, gcp, true, null);
    }

    @Override
    public String toString() {
        return "Info{" +
                "ref='" + ref + '\'' +
                ", nextRef=" + nextRef + '\'' +
                ", read='" + getRead() + '\'' +
                ", log10l=" + log10l + '\'' +
                ", newRead=" + newRead +
                '}';
    }

    public static double runHMMs(final PairHMM hmm, final List<PairHMMTestData> data, final boolean runSingly) {
        double result = 0;
        if ( runSingly ) {
            for ( final PairHMMTestData datum : data )
                result += datum.runHMM(hmm);
        } else {
            // running in batch mode
            final PairHMMTestData first = data.get(0);
            int maxHaplotypeLen = calcMaxHaplotypeLen(data);
            hmm.initialize(first.getRead().length(), maxHaplotypeLen);
            for ( final PairHMMTestData datum : data ) {
                result += hmm.computeReadLikelihoodGivenHaplotypeLog10(datum.ref.getBytes(), datum.getRead().getBytes(),
                        datum.baseQuals, datum.insQuals, datum.delQuals, datum.gcp, datum.newRead, datum.nextRef.getBytes());

            }
        }
        return result;
    }

    public static int calcMaxHaplotypeLen(final List<PairHMMTestData> data) {
        int maxHaplotypeLen = 0;
        for ( final PairHMMTestData datum : data )
            maxHaplotypeLen = Math.max(maxHaplotypeLen, datum.ref.length());
        return maxHaplotypeLen;
    }

    public static int calcMaxReadLen(final List<PairHMMTestData> data) {
        int maxReadLen = 0;
        for ( final PairHMMTestData datum : data )
            maxReadLen = Math.max(maxReadLen, datum.getRead().length());
        return maxReadLen;
    }

    public static Map<String, List<PairHMMTestData>> readLikelihoods(final File file) throws IOException {
        final Map<String, List<PairHMMTestData>> results = new LinkedHashMap<>();

        InputStream in = new FileInputStream(file);
        if ( file.getName().endsWith(".gz") ) {
            in = new GZIPInputStream(in);
        }

        String[] nextEntry;
        String[] thisEntry = null;
        for ( final String line : new XReadLines(in) ) {
            // peak at the next entry (to get the haplotype bases)
            nextEntry = line.split(" ");
            // process the current entry
            if (thisEntry != null) {
                final PairHMMTestData info = new PairHMMTestData(
                        thisEntry[0], nextEntry[0], thisEntry[1],
                        SAMUtils.fastqToPhred(thisEntry[2]),
                        SAMUtils.fastqToPhred(thisEntry[3]),
                        SAMUtils.fastqToPhred(thisEntry[4]),
                        SAMUtils.fastqToPhred(thisEntry[5]),
                        Double.parseDouble(thisEntry[6]),
                        ! results.containsKey(thisEntry[1]));

                if ( ! results.containsKey(info.read) )  {
                    results.put(info.read, new LinkedList<>());
                }
                final List<PairHMMTestData> byHap = results.get(info.read);
                byHap.add(info);
            }
            // update the current entry
            thisEntry = nextEntry;
        }
        // process the final entry
        final PairHMMTestData info = new PairHMMTestData(
                thisEntry[0], null, thisEntry[1],
                SAMUtils.fastqToPhred(thisEntry[2]),
                SAMUtils.fastqToPhred(thisEntry[3]),
                SAMUtils.fastqToPhred(thisEntry[4]),
                SAMUtils.fastqToPhred(thisEntry[5]),
                Double.parseDouble(thisEntry[6]),
                ! results.containsKey(thisEntry[1]));

        if ( ! results.containsKey(info.read) )  {
            results.put(info.read, new LinkedList<>());
        }
        final List<PairHMMTestData> byHap = results.get(info.read);
        byHap.add(info);

        return results;
    }


    /*
     * simplified likelihoods file reader that returns a list instead of a map
     *
     * readLikelihoods() method was reordering inputs, with the result that caching would be more efficient
     * This method simply returns a list of read/haplotype pairs in their original order, providing a more realistic caching scenario
     */
    public static List<PairHMMTestData> readLikelihoodsInOrder(final File file) throws IOException {
        final List<PairHMMTestData> results = new LinkedList<>();

        InputStream in = new FileInputStream(file);
        if ( file.getName().endsWith(".gz") ) {
            in = new GZIPInputStream(in);
        }

        String previousRead = null;
        String[] nextEntry;
        String[] thisEntry = null;
        for ( final String line : new XReadLines(in) ) {
            // peak at the next entry (to get the haplotype bases)
            nextEntry = line.split(" ");
            // process the current entry
            if (thisEntry != null) {
                final PairHMMTestData info = new PairHMMTestData(
                        thisEntry[0], nextEntry[0], thisEntry[1],
                        SAMUtils.fastqToPhred(thisEntry[2]),
                        SAMUtils.fastqToPhred(thisEntry[3]),
                        SAMUtils.fastqToPhred(thisEntry[4]),
                        SAMUtils.fastqToPhred(thisEntry[5]),
                        Double.parseDouble(thisEntry[6]),
                        !(thisEntry[1].equals(previousRead)));

                results.add(info);
                previousRead = info.getRead();
            }
            // update the current entry
            thisEntry = nextEntry;
        }
        // process the final entry
        final PairHMMTestData info = new PairHMMTestData(
                thisEntry[0], null, thisEntry[1],
                SAMUtils.fastqToPhred(thisEntry[2]),
                SAMUtils.fastqToPhred(thisEntry[3]),
                SAMUtils.fastqToPhred(thisEntry[4]),
                SAMUtils.fastqToPhred(thisEntry[5]),
                Double.parseDouble(thisEntry[6]),
                !(thisEntry[1].equals(previousRead)));

        results.add(info);

        return results;
    }

    public String getRead() {
        return read;
    }
}
