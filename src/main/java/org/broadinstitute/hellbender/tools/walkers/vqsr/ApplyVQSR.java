package org.broadinstitute.hellbender.tools.walkers.vqsr;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.*;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.VariantProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.UserException;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collector;
import java.util.stream.Collectors;

/**
 * Applies cuts to the input vcf file (by adding filter lines) to achieve the desired novel truth sensitivity levels which were specified during VariantRecalibration
 *
 * <p>
 * Using the tranche file generated by the previous step the ApplyVQSR walker looks at each variant's VQSLOD value
 * and decides which tranche it falls in. Variants in tranches that fall below the specified truth sensitivity filter level
 * have their filter field annotated with its tranche level. This will result in a call set that simultaneously is filtered
 * to the desired level but also has the information necessary to pull out more variants for a higher sensitivity but a
 * slightly lower quality level.
 *
 * <h3>Input</h3>
 * <p>
 * The input raw variants to be recalibrated.
 * <p>
 * The recalibration table file in VCF format that was generated by the VariantRecalibrator walker.
 * <p>
 * The tranches file that was generated by the VariantRecalibrator walker.
 *
 * <h3>Output</h3>
 * <p>
 * A recalibrated VCF file in which each variant is annotated with its VQSLOD and filtered if the score is below the desired quality level.
 *
 * <h3>Examples</h3>
 * <pre>
 * java -Xmx3g -jar GenomeAnalysisTK.jar \
 *   -T ApplyVQSR \
 *   -R reference/human_g1k_v37.fasta \
 *   -input NA12878.HiSeq.WGS.bwa.cleaned.raw.subset.b37.vcf \
 *   --ts_filter_level 99.0 \
 *   -tranchesFile path/to/output.tranches \
 *   -recalFile path/to/output.recal \
 *   -mode SNP \
 *   -o path/to/output.recalibrated.filtered.vcf
 * </pre>
 *
 */
@CommandLineProgramProperties(
        usage = "Applies the VQSR table to the input variant file.",
        usageShort = "ApplyVQSR",
        programGroup = VariantProgramGroup.class
)
public final class ApplyVQSR extends VariantWalker {

    private static final Logger logger = LogManager.getLogger(ApplyVQSR.class);

    public static final String VQS_LOD_KEY = "VQSLOD"; // Log odds ratio of being a true variant versus being false under the trained gaussian mixture model
    public static final String CULPRIT_KEY = "culprit"; // The annotation which was the worst performing in the Gaussian mixture model, likely the reason why the variant was filtered out
    public static final String NEGATIVE_LABEL_KEY = "NEGATIVE_TRAIN_SITE"; // this variant was used in the negative training set
    public static final String POSITIVE_LABEL_KEY = "POSITIVE_TRAIN_SITE"; // this variant was used in the positive training set

    public static final String LOW_VQSLOD_FILTER_NAME = "LOW_VQSLOD";
    private final double DEFAULT_VQSLOD_CUTOFF = 0.0;

    /////////////////////////////
    // Inputs
    /////////////////////////////

    @Argument(fullName="recal_file", shortName="recalFile", doc="The input recal file used by ApplyVQSR", optional = false)
    protected FeatureInput<VariantContext> recal;

    @Argument(fullName="tranches_file", shortName="tranchesFile", doc="The input tranches file describing where to cut the data", optional = true)
    protected File TRANCHES_FILE;

    /////////////////////////////
    // Outputs
    /////////////////////////////
    @Argument( shortName = "vcfOut", fullName = "vcfOut" , doc="The output filtered and recalibrated VCF file in which each variant is annotated with its VQSLOD value")
    private File vcfOut = null;

    /////////////////////////////
    // Command Line Arguments
    /////////////////////////////
    @Argument(fullName="ts_filter_level", shortName="ts_filter_level", doc="The truth sensitivity level at which to start filtering", optional = true)
    protected Double TS_FILTER_LEVEL = null;

    @Argument(fullName="lodCutoff", shortName="lodCutoff", doc="The VQSLOD score below which to start filtering", optional = true)
    protected Double VQSLOD_CUTOFF = null;

    /**
     * For this to work properly, the -ignoreFilter argument should also be applied to the VariantRecalibration command.
     */
    @Argument(fullName="ignore_filter", shortName="ignoreFilter", doc="If specified, the recalibration will be applied to variants marked as filtered by the specified filter name in the input VCF file", optional = true)
    private List<String> IGNORE_INPUT_FILTERS = null;

    @Argument(fullName="ignore_all_filters", shortName="ignoreAllFilters", doc="If specified, the variant recalibrator will ignore all input filters. Useful to rerun the VQSR from a filtered output file.", optional = true)
    private Boolean IGNORE_ALL_FILTERS = false;

    @Argument(fullName="excludeFiltered", shortName="ef", doc="Don't output filtered loci after applying the recalibration", optional = true)
    protected Boolean EXCLUDE_FILTERED = false;

    @Argument(fullName = "mode", shortName = "mode", doc = "Recalibration mode to employ: 1.) SNP for recalibrating only SNPs (emitting indels untouched in the output VCF); 2.) INDEL for indels; and 3.) BOTH for recalibrating both SNPs and indels simultaneously.", optional = true)
    public VariantRecalibratorArgumentCollection.Mode MODE = VariantRecalibratorArgumentCollection.Mode.SNP;

    /////////////////////////////
    // Private Member Variables
    /////////////////////////////
    final private List<Tranche> tranches = new ArrayList<>();
    final private Set<String> ignoreInputFilterSet = new TreeSet<>();

    private VariantContextWriter vcfWriter;


    //---------------------------------------------------------------------------------------------------------------
    //
    // initialize
    //
    //---------------------------------------------------------------------------------------------------------------

    @Override
    public void onTraversalStart() {
        tranches.addAll(initializeTrancheSet());
        ignoreInputFilterSet.addAll(makeListOfInputsToIgnore());

        VariantContextWriterBuilder builder = new VariantContextWriterBuilder()
                      .unsetBuffering();

        vcfWriter = builder
                      .setOutputFile(vcfOut)
                      .setOutputFileType(VariantContextWriterBuilder.OutputType.VCF)
                      .setOptions(VariantContextWriterBuilder.NO_OPTIONS)
                      .build();
        vcfWriter.writeHeader(makeVCFHeaderForOutput());
    }

    /**
     * Collector that additionally reverses the list.
     */
    public static <T> Collector<T, ?, List<T>> toListReversed() {
        return Collectors.collectingAndThen(Collectors.toList(), l -> {Collections.reverse(l); return l;} );
    }

    private List<Tranche> initializeTrancheSet() {
        if( TS_FILTER_LEVEL == null) return Collections.emptyList();

        try {
            return Tranche.readTranches(TRANCHES_FILE).stream()
                            .filter(t -> (t.targetTruthSensitivity >= TS_FILTER_LEVEL))
                            .collect(toListReversed());         // tranches must be ordered from best (lowest truth sensitivity) to worst (highest truth sensitivity)
        } catch (IOException e) {
            throw new UserException.CouldNotReadInputFile(TRANCHES_FILE, e);
        }
    }

    private VCFHeader makeVCFHeaderForOutput() {
        VCFHeader inputHeader = getHeaderForVariants();
        final Set<VCFHeaderLine> headerLines = makeNewVcfHeaderLines(inputHeader);
        final ArrayList<String> sampleNames = inputHeader.getSampleNamesInOrder();
        return new VCFHeader(headerLines, sampleNames);
    }

    private List<String> makeListOfInputsToIgnore() {
        return IGNORE_INPUT_FILTERS == null ? Collections.emptyList() : IGNORE_INPUT_FILTERS;
    }

    private Set<VCFHeaderLine> makeNewVcfHeaderLines(VCFHeader inputHeader) {
        final Set<VCFHeaderLine> headerLines = new HashSet<>();
        headerLines.addAll(inputHeader.getMetaDataInSortedOrder());
        headerLines.addAll(makeVQSRStandardHeaderLines());

        if( TS_FILTER_LEVEL != null ) {
            // if the user specifies both ts_filter_level and lodCutoff then throw a user error
            if( VQSLOD_CUTOFF != null ) {
                throw new UserException.CommandLineException("Arguments --ts_filter_level and --lodCutoff are mutually exclusive. Please only specify one option.");
            }

            if( tranches.size() >= 2 ) {
                for( int i = 0; i < tranches.size() - 1; i++ ) {
                    final Tranche t = tranches.get(i);
                    headerLines.add(new VCFFilterHeaderLine(t.name, String.format("Truth sensitivity tranche level for " + t.model.toString() + " model at VQS Lod: " + t.minVQSLod + " <= x < " + tranches.get(i + 1).minVQSLod)));
                }
            }
            if( tranches.size() >= 1 ) {
                headerLines.add(new VCFFilterHeaderLine(tranches.get(0).name + "+", String.format("Truth sensitivity tranche level for " + tranches.get(0).model.toString() + " model at VQS Lod < " + tranches.get(0).minVQSLod)));
            } else {
                throw new UserException.BadInput("No tranches were found in the file or were above the truth sensitivity filter level " + TS_FILTER_LEVEL);
            }

            logger.info("Keeping all variants in tranche " + tranches.get(tranches.size()-1));
        } else {
            if( VQSLOD_CUTOFF == null ) {
                VQSLOD_CUTOFF = DEFAULT_VQSLOD_CUTOFF;
            }
            headerLines.add(new VCFFilterHeaderLine(LOW_VQSLOD_FILTER_NAME, "VQSLOD < " + VQSLOD_CUTOFF));
            logger.info("Keeping all variants with VQSLOD >= " + VQSLOD_CUTOFF);
        }
        return headerLines;
    }

    public static Set<VCFHeaderLine> makeVQSRStandardHeaderLines() {
        final Set<VCFHeaderLine> res = new HashSet<>();
        res.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.END_KEY));
        res.add(new VCFInfoHeaderLine(VQS_LOD_KEY, 1, VCFHeaderLineType.Float, "Log odds ratio of being a true variant versus being false under the trained gaussian mixture model"));
        res.add(new VCFInfoHeaderLine(CULPRIT_KEY, 1, VCFHeaderLineType.String, "The annotation which was the worst performing in the Gaussian mixture model, likely the reason why the variant was filtered out"));
        res.add(new VCFInfoHeaderLine(POSITIVE_LABEL_KEY, 1, VCFHeaderLineType.Flag, "This variant was used to build the positive training set of good variants"));
        res.add(new VCFInfoHeaderLine(NEGATIVE_LABEL_KEY, 1, VCFHeaderLineType.Flag, "This variant was used to build the negative training set of bad variants"));
        return res;
    }

    //---------------------------------------------------------------------------------------------------------------
    //
    // map
    //
    //---------------------------------------------------------------------------------------------------------------

    @Override
    public void apply(VariantContext vc, ReadsContext readsContext, ReferenceContext ref, FeatureContext featureContext) {
        final List<VariantContext> recals =  featureContext.getValues(recal);

        if( VQSRUtils.checkVariationClass( vc, MODE ) && (IGNORE_ALL_FILTERS || vc.isNotFiltered() || ignoreInputFilterSet.containsAll(vc.getFilters())) ) {

            final VariantContext recalDatum = getMatchingRecalVC(vc, recals);
            if( recalDatum == null ) {
                throw new UserException.BadInput("Encountered input variant which isn't found in the input recal file. Please make sure VariantRecalibrator and ApplyVQSR were run on the same set of input variants. First seen at: " + vc );
            }

            final String lodString = recalDatum.getAttributeAsString(VQS_LOD_KEY, null);
            if( lodString == null ) {
                throw new UserException.BadInput("Encountered a malformed record in the input recal file. There is no lod for the record at: " + vc );
            }
            final double lod;
            try {
                lod = Double.valueOf(lodString);
            } catch (NumberFormatException e) {
                throw new UserException.BadInput("Encountered a malformed record in the input recal file. The lod is unreadable for the record at: " + vc );
            }

            VariantContextBuilder builder = new VariantContextBuilder(vc);

            // Annotate the new record with its VQSLOD and the worst performing annotation
            builder.attribute(VQS_LOD_KEY, lod);
            builder.attribute(CULPRIT_KEY, recalDatum.getAttribute(CULPRIT_KEY));
            if ( recalDatum.hasAttribute(POSITIVE_LABEL_KEY)) {
                builder.attribute(POSITIVE_LABEL_KEY, true);
            }
            if ( recalDatum.hasAttribute(NEGATIVE_LABEL_KEY)) {
                builder.attribute(NEGATIVE_LABEL_KEY, true);
            }

            final String filterString = generateFilterString(lod);

            if( filterString.equals(VCFConstants.PASSES_FILTERS_v4) ) {
                builder.passFilters();
            } else {
                builder.filters(filterString);
            }

            final VariantContext outputVC = builder.make();
            if( !EXCLUDE_FILTERED || outputVC.isNotFiltered() ) {
                vcfWriter.add( outputVC );
            }
        } else { // valid VC but not compatible with this mode, so just emit the variant untouched
            vcfWriter.add(vc);
        }
    }

    /**
     * Generate the VCF filter string for this record based on the provided lod score
     * @param lod non-null double
     * @return the String to use as the VCF filter field
     */
    protected String generateFilterString( final double lod ) {
        String filterString = null;
        if( TS_FILTER_LEVEL != null ) {
            for( int i = tranches.size() - 1; i >= 0; i-- ) {
                final Tranche tranche = tranches.get(i);
                if( lod >= tranche.minVQSLod ) {
                    if( i == tranches.size() - 1 ) {
                        filterString = VCFConstants.PASSES_FILTERS_v4;
                    } else {
                        filterString = tranche.name;
                    }
                    break;
                }
            }

            if( filterString == null ) {
                filterString = tranches.get(0).name+"+";
            }
        } else {
            filterString = ( lod < VQSLOD_CUTOFF ? LOW_VQSLOD_FILTER_NAME : VCFConstants.PASSES_FILTERS_v4 );
        }

        return filterString;
    }

    private static VariantContext getMatchingRecalVC(final VariantContext target, final List<VariantContext> recalVCs) {
        for( final VariantContext recalVC : recalVCs ) {
            if ( target.getEnd() == recalVC.getEnd() ) {
                return recalVC;
            }
        }

        return null;
    }

}

