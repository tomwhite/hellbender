package org.broadinstitute.hellbender.tools.dataflow.pipelines;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.text.XReadLines;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Objects;


public class MarkDuplicatesDataflowIntegrationTest extends CommandLineProgramTest{

    @DataProvider(name = "md")
    public Object[][] md(){
        return new Object[][]{
                {new File(getTestDataDir(),"example.chr1.1-1K.unmarkedDups.bam"), 90},  //90 total reads, 70 mapped reads, 6 dups
/*
90 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 secondary
0 + 0 supplimentary
6 + 0 duplicates
77 + 0 mapped (85.56%:nan%)
90 + 0 paired in sequencing
44 + 0 read1
46 + 0 read2
58 + 0 properly paired (64.44%:nan%)
64 + 0 with itself and mate mapped
13 + 0 singletons (14.44%:nan%)
6 + 0 with mate mapped to a different chr
3 + 0 with mate mapped to a different chr (mapQ>=5)
*/
//                {new File(getTestDataDir(),"flag_stat.bam"), 11},
//                {new File(getTestDataDir(),"dataflow/count_bases.bam"), 7},
//                {new File(getTestDataDir(),"BQSR/HiSeq.1mb.1RG.2k_lines.bam"), 7},   //NOTE this blows up htsjdk.samtools.SAMException: Illegal MD pattern: 0C0T0T1A61A34 for read 20GAVAAXX100126:8:5:8975:183748 with CIGAR 5M4D96M

        };
    }

    @Test(groups = "dataflow", dataProvider = "md")
    public void testMarkDuplicatesDataflowIntegrationTestLocal(File input, long expected) throws IOException {

        ArgumentsBuilder args = new ArgumentsBuilder();
        args.add("--"+StandardArgumentDefinitions.INPUT_LONG_NAME); args.add(input.getPath());
        args.add("--" + StandardArgumentDefinitions.OUTPUT_LONG_NAME);
//        args.add("--" + "secret client-secrets.json");
        File placeHolder = createTempFile("markdups", ".txt");
        args.add(placeHolder.getPath());

        runCommandLine(args.getArgsArray());
        File outputFile = findDataflowOutput(placeHolder);

        Assert.assertTrue(outputFile.exists());
        Assert.assertEquals(new XReadLines(outputFile).readLines().size(), expected);
    }
}