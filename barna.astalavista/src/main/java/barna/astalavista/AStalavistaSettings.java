package barna.astalavista;

import barna.commons.log.Log;
import barna.commons.parameters.*;
import barna.model.ASEvent;
import barna.model.Graph;
import barna.model.Species;
import barna.model.Transcript;
//import barna.model.gff.GTFschema;
import barna.model.constants.Constants;
import barna.model.splicegraph.SplicingGraph;

import java.io.File;
import java.io.OutputStream;
import java.util.EnumSet;
import java.util.Iterator;

/**
 * Created with IntelliJ IDEA.
 * User: micha
 * Date: 6/18/12
 * Time: 1:14 PM
 * To change this template use File | Settings | File Templates.
 */
public class AStalavistaSettings extends ParameterSchema /*GTFschema*/ {

    // TODO {SplicingGraph.class, null, SJextractor.class, AttributeExtractor.class};	// null= LaVista.class


    /**
     * Print parameters and descriptions.
     */
    public static final Parameter<Boolean> HELP = Parameters.booleanParameter("HELP",
            "print parameters and descriptions", false, null).longOption("printParameters").shortOption('h');

    /**
     * Path to the reference annotation.
     */
    public static final Parameter<File> REF_FILE = Parameters.fileParameter("INPUT",
            "path to the GTF reference annotation",
            null, new ParameterValidator() {
        @Override
        public void validate(ParameterSchema schema, Parameter parameter) throws ParameterException {
            File refFile = (File) schema.get(parameter);
            if (refFile == null) {
                throw new ParameterException("Reference annotation cannot be null!");
            }
            if (!refFile.exists()) {
                throw new ParameterException("The reference annotation " + refFile.getAbsolutePath() + " could not be found!");
            }
        }
    }, null).longOption("in").shortOption('i');

    /**
     * Filters on the input elements
     */
    public static enum InputOptions {
        /* consider canonical splice sites only */
        CSS,
        /* acceptable introns */
        IOK
    };

    /**
     * Parameter collecting criteria for elements of the input to be considered
     */
    public static final Parameter<EnumSet<OutputOptions>> INOPTIONS = Parameters.enumSetParameter(
            "INOPTIONS",
            "Toggle criteria for elements of the input",
            EnumSet.noneOf(OutputOptions.class),
            OutputOptions.class,
            null).longOption("io");


    /**
     * File with GeneID parameters / splice site profiles.
     */
    public static final Parameter<File> GENEID_PARAM = Parameters.fileParameter("GENEID_PARAM",
            "name and path of a file with the GeneID models for splice sites",
            null, null, null).longOption("gparam");


    /**
     * File with variant information (vcf format).
     */
    public static final Parameter<File> VARIANTS = Parameters.fileParameter("VARIANTS",
            "name and path of a file with the variant information (vcf)",
            null, new ParameterValidator() {

            @Override
            public void validate(ParameterSchema schema, Parameter parameter) throws ParameterException {
                File vcf = (File) schema.get(parameter);

                if (vcf == null || (!vcf.exists())) {
                    throw new ParameterException("VCF file not valid: "+ vcf== null? "null": vcf.getAbsolutePath());
                }
            }
    }).longOption("vcf");

    /**
     * Keyword to redirect program output to standard out stream.
     */
    public static final String STDOUT= "stdout";

    /**
     * Path to the GTF output annotation.
     */
    public static final Parameter<File> OUTPUT = Parameters.fileParameter("OUTPUT",
            "keyword '"+STDOUT+"' for standard output, or a path to the GTF output file",
            null, null, null).longOption("out").shortOption('o');


    /**
     * Path to the directory with the genomic sequences,
     * i.e., one fasta file per chromosome/scaffold/contig
     * with a file name corresponding to the identifiers of
     * the first column in the GTF annotation.
     */
    public static final Parameter<File> GENOME = Parameters.fileParameter("GENOME",
                    "path to the directory with the genomic sequences,\n" +
                    "i.e., one fasta file per chromosome/scaffold/contig\n" +
                    "with a file name corresponding to the identifiers of\n" +
                    "the first column in the GTF annotation", null, new ParameterValidator() {
        @Override
        public void validate(ParameterSchema schema, Parameter parameter) throws ParameterException {
            File genomeFile = (File) schema.get(parameter);
            if (genomeFile == null) {
                throw new ParameterException("You have to specify a genome directory");
            }
        }
    }, null).longOption("genome").shortOption('g');

    /**
     * The temporary directory.
     */
    public static final Parameter<File> TMP_DIR = Parameters.fileParameter("TMP_DIR", "The temporary directory", new File(System.getProperty("java.io.tmpdir")), new ParameterValidator() {
        public void validate(ParameterSchema schema, Parameter parameter) throws ParameterException {
            File file = (File) schema.get(parameter);
            if (file == null) {
                schema.set(parameter, new File(System.getProperty(Constants.PROPERTY_TMPDIR)));
            }
            if (!file.exists()) {
                throw new ParameterException("The temporary directory " + file.getAbsolutePath()
                        + " could not be found!");
            }
            if (!file.canWrite()) {
                throw new ParameterException("The temporary directory " + file.getAbsolutePath()
                        + " is not writable!");
            }

        }
    }).longOption("tmp");


    /**
     * Dimension of the AS events to be extracted,
     * retrieves 'complete' events <TM> for parameter
     * values < 2.
     */
    public static final Parameter<Integer> DIMENSION = Parameters.intParameter("DIMENSION",
            "Dimension of the AS events to be extracted, retrieves 'complete' events <TM>\n" +
            "for parameter values < 2",
            2, new ParameterValidator() {
        @Override
        public void validate(ParameterSchema schema, Parameter parameter) throws ParameterException {
            int v= (Integer) schema.get(parameter);
            if (v< 2)
                EventExtractor.n= (-1); // complete events
            else
                EventExtractor.n= 2;
        }
    }).longOption("dim").shortOption('k');


    /**
     * Different types of variation found in exon-intron structures
     * of transcripts:
     * <ul><li>AS= alternative splicing when comprising
     * at least one alternative splice site.
     * Types of alternative splicing can either be &quot;internal&quot;
     * and delimited by two common sites, or &quot;external&quot;
     * comprising at least one alternative splice site in addition
     * to alternative 5'- or 3' transcript structures.</li>
     * extending transcript structures by additional (splice) sites
     * to the 5'- or the 3'-end</li>
     * <li>VS= variable sites is any other form of sites that differ
     * between overlapping transcript structures</li>
     * @see barna.model.ASEvent#getType()
     * @see <a href="http://www.ploscompbiol.org/article/info%3Adoi%2F10.1371%2Fjournal.pcbi.1000147">
     *     http://www.ploscompbiol.org/article/info%3Adoi%2F10.1371%2Fjournal.pcbi.1000147</a><li>DS= additional splicing that are flanked by a common site and
     */
    public static enum EventTypes {
        /** external AS events */
        ASE,
        /** internal AS events */
        ASI,
        /** adDitional Splicing events */
        DSP,
        /** Variable Site events */
        VST,
    };


    /**
     * Parameter for the list of event types that is to be considered.
     */
    public static final Parameter<EnumSet<EventTypes>> EVENTS = Parameters.enumSetParameter(
            "EVENTS",
            " Type of events that is considered",
            EnumSet.of(EventTypes.ASI),
            EventTypes.class,
            null).longOption("events").shortOption('e');


    /**
     * Level of intron confidence, below which introns are trusted without checks.
     * The default is to trust all introns (i.e., ic= 255). Introns are assigned a
     * confidency class:
     * <ul>
     * <li>0 for 'RefSeq' appears in the source field of the annotation</li>
     * <li>1 for 'mRNA' appears in the source field of the annotation</li>
     * <li>2 for 'EST' appears in the source field of the annotation</li>
     * </ul>
     * All introns in transcripts of confidence level > threshold are discarded.
     */
    public static final Parameter<Integer> INTRON_CONFIDENCE = Parameters.intParameter("INTRON_CONFIDENCE",
            "Confidence level for introns in the annotation",
            Transcript.ID_SRC_MOST_INCONFIDENT, new ParameterValidator() {
        @Override
        public void validate(ParameterSchema schema, Parameter parameter) throws ParameterException {
            int x= schema.get(INTRON_CONFIDENCE);
            if (x> Byte.MAX_VALUE|| x< 0)
                throw new IllegalArgumentException("Invalid confidence level "+ x);
        }
    }).longOption("ic");

    /**
     * Level of confidence for edges (i.e., annotated transcription starts/poly-adenylation sites).
     * The default is to trust no annotated edge and to extend overlapping first/last exons of a
     * transcript to their most extreme position:
     * <ul>
     * <li>0 if 'RefSeq' appears in the source field of the annotation</li>
     * <li>1 if 'mRNA' appears in the source field of the annotation</li>
     * <li>2 if 'EST' appears in the source field of the annotation</li>
     * <li>3 if if none of the above applies</li>
     * </ul>
     * All transcript edges of confidence level > edgeConfidence will be extended in case the
     * annotation shows another exon with the same adjacent splice site and an earlier/later
     * start/end.
     */
    public static final Parameter<Integer> EDGE_CONFIDENCE = Parameters.intParameter("EDGE_CONFIDENCE",
            "Transcript edge confidence level",
            Transcript.ID_SRC_MOST_INCONFIDENT, new ParameterValidator() {
        @Override
        public void validate(ParameterSchema schema, Parameter parameter) throws ParameterException {

            int x= schema.get(EDGE_CONFIDENCE);
            if (x> Byte.MAX_VALUE|| x< 0)
                throw new IllegalArgumentException("Invalid confidence level "+ x);
        }
    }).longOption("ec");

    /**
     * Consider only canonical splice sites during graph construction.
     */
    public static final Parameter<Boolean> CANONICAL = Parameters.booleanParameter("CANONICAL",
            "consider only canonical sites", false, null).longOption("css");

    /**
     * Flags to control output options for events:
     * <ul>
     *     <li>CP3: predict 3'-complete</li>
     *     <li>FLT: output flank type, i.e. 'constitutive' or 'alternative'</li>
     *     <li>NMD: predict NMD</li>
     *     <li>SEQ: output sequences of flanking splice sites</li>
     * </ul>
     */
    public static enum OutputOptions {
        /* predict 3'-complete */
        CP3,
        /* output flank type, ie 'constitutive' or 'alternative' */
        FLT,
        /* predict NMD */
        NMD,
        /* output splice site sequences of event flanks */
        SEQ
    };

    /**
     * Parameter collecting flags for output options
     */
    public static final Parameter<EnumSet<OutputOptions>> OUTOPTIONS = Parameters.enumSetParameter(
            "OUTOPTIONS",
            "Toggle optional attributes to be output",
            EnumSet.noneOf(OutputOptions.class),
            OutputOptions.class,
            null).longOption("oo");



    /**
     * Checks whether a folder with genomic sequences is necessary in order
     * to complete all tasks specified with <code>this</code> parameter set.
     * @return <code>true</code> if the genomic sequence is required given the
     * current parameters, <code>false</code> otherwise
     */
    public boolean requiresGenomicSequence() {
        if (get(AStalavistaSettings.OUTOPTIONS).contains(OutputOptions.SEQ))
            return true;
        if (get(AStalavistaSettings.INOPTIONS).contains(InputOptions.CSS)
                || get(AStalavistaSettings.INOPTIONS).contains(InputOptions.IOK))
            return true;

        return false;
    }

    @Override
    public void write(OutputStream out) {
        super.write(out);
/*        if (isPrintParameters()) {
            FluxSimulatorSettings settings = new FluxSimulatorSettings();
            settings.write(System.out);
            return false;

            System.err.println("Here is a list of the options I understand:\n");
            System.err.println("-i, --input <input file>");
            System.err.println("This is a bit important, I cannot work without an input annotation. I want a GTF file with " +
                    "transcript annotations (exon features, with a mandatory optional attribute named \'transcript_id\') " +
                    "IN THE SAME COLUMN (i.e., if the transcript identifier of the 1st line is in column #10, it has to be in " +
                    "all lines of the file in column #10. The rest of the file should comply with the standard as specified " +
                    "at http://mblab.wustl.edu/GTF2.html.\n"+
                    "There may also be CDS features, but they become only interesting when checking for additional things " +
                    "as NMD probability etc.."+
                    "\n");

            // ...

            System.err.println("AStalavista.");
            System.exit(0);
*/
    }



    /**
     * Score splice site sequences.
     * @deprecated replace by OutputSite.SITESCORES
     */
    public static final Parameter<File> SCORE_SITES = Parameters.fileParameter("SCORE_SITES",
            "score splice site sequences", null, null, null).longOption("scoresites");


}
