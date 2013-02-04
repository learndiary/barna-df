package barna.astalavista;

import barna.commons.parameters.*;
import barna.model.Transcript;
//import barna.model.gff.GTFschema;
import barna.model.constants.Constants;

import java.io.File;
import java.io.OutputStream;
import java.util.EnumSet;

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
    public static final Parameter<File> IN_FILE = Parameters.fileParameter("IN_FILE",
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


    /*
     * Parameter specifying a collection of 'source' tags that are read from the input,
     * possibly other flag for non-/coding transcripts to be considered
     *
    public static final Parameter<EnumSet<EventOptions>> IN_OPTIONS = Parameters.enumSetParameter(
            "IN_OPTIONS",
            "Toggle criteria for elements of the input",
            EnumSet.noneOf(EventOptions.class),
            EventOptions.class,
            null).longOption("io");
     */

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
     * @deprecated to be refactored to IN_OPTIONS
     */
    public static final Parameter<Integer> INTRON_CONFIDENCE = Parameters.intParameter("INTRON_CONFIDENCE",
            "Confidence level for introns in the annotation",
            Transcript.ID_SRC_MOST_INCONFIDENT, new ParameterValidator() {
        @Override
        public void validate(ParameterSchema schema, Parameter parameter) throws ParameterException {
            int x = schema.get(INTRON_CONFIDENCE);
            if (x > Byte.MAX_VALUE || x < 0)
                throw new IllegalArgumentException("Invalid confidence level " + x);
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
     * @deprecated to be refactored to IN_OPTIONS
     */
    public static final Parameter<Integer> EDGE_CONFIDENCE = Parameters.intParameter("EDGE_CONFIDENCE",
            "Transcript edge confidence level",
            Transcript.ID_SRC_MOST_INCONFIDENT, new ParameterValidator() {
        @Override
        public void validate(ParameterSchema schema, Parameter parameter) throws ParameterException {

            int x = schema.get(EDGE_CONFIDENCE);
            if (x > Byte.MAX_VALUE || x < 0)
                throw new IllegalArgumentException("Invalid confidence level " + x);
        }
    }).longOption("ec");

    /**
     * File with GeneID parameters / splice site profiles.
     */
    public static final Parameter<File> GENE_ID = Parameters.fileParameter("GENE_ID",
            "name and path of a file with the GeneID models for splice sites",
            null, null, null).longOption("gid").shortOption('g');


    /**
     * File with variant information (vcf format).
     */
    public static final Parameter<File> VARIANT_FILE = Parameters.fileParameter("VARIANT_FILE",
            "name and path of a file with the variant information (vcf)",
            null, new ParameterValidator() {

            @Override
            public void validate(ParameterSchema schema, Parameter parameter) throws ParameterException {
                File vcf = (File) schema.get(parameter);

                if (vcf == null || (!vcf.exists())) {
                    throw new ParameterException("VCF file not valid: "+ vcf== null? "null": vcf.getAbsolutePath());
                }
            }
    }).longOption("vcf").shortOption('v');

    /**
     * Path to the GTF output annotation.
     */
    public static final Parameter<File> EVENTS_FILE = Parameters.fileParameter("EVENTS_FILE",
            "a path to the GTF output file for events",
            null, new ParameterValidator() {
        @Override
        public void validate(ParameterSchema schema, Parameter parameter) throws ParameterException {
            File f = (File) schema.get(parameter);
            if (f == null|| !f.getParentFile().canWrite()) {
                throw new ParameterException("Invalid output file "+ f.getAbsolutePath());
            }
        }
    }).longOption("eo").shortOption('o');

    /**
     * Path to the GTF output annotation.
     */
    public static final Parameter<File> SITES_FILE = Parameters.fileParameter("SITES_FILE",
            "a path to the VCF output file for sites",
            null, new ParameterValidator() {
        @Override
        public void validate(ParameterSchema schema, Parameter parameter) throws ParameterException {
            File f = (File) schema.get(parameter);
            if (f == null|| !f.getParentFile().canWrite()) {
                throw new ParameterException("Invalid output file "+ f.getAbsolutePath());
            }
        }
    }).longOption("so").shortOption('f');

    /**
     * Path to the directory with the genomic sequences,
     * i.e., one fasta file per chromosome/scaffold/contig
     * with a file name corresponding to the identifiers of
     * the first column in the GTF annotation.
     */
    public static final Parameter<File> CHR_SEQ = Parameters.fileParameter("CHR_SEQ",
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
    }, null).longOption("chr").shortOption('c');

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
    public static final Parameter<Integer> EVENTS_DIMENSION = Parameters.intParameter("EVENTS_DIMENSION",
            "Dimension of the AS events to be extracted, retrieves 'complete' events <TM>\n" +
            "for parameter values < 2",
            2, new ParameterValidator() {
        @Override
        public void validate(ParameterSchema schema, Parameter parameter) throws ParameterException {
            int v= (Integer) schema.get(parameter);
        }
    }).longOption("ed").shortOption('d');


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
            "Type of events that are considered",
            EnumSet.of(EventTypes.ASI),
            EventTypes.class,
            null).longOption("ev").shortOption('e');

    /**
     * Flags to control output options for events:
     * <ul>
     *     <li>CP3: predict 3'-complete</li>
     *     <li>FLT: output flank type, i.e. 'constitutive' or 'alternative'</li>
     *     <li>NMD: predict NMD</li>
     *     <li>SEQ: output sequences of flanking splice sites</li>
     *     <li>IOK: consider only events with acceptable introns</li>
     * </ul>
     */
    public static enum EventOptions {
        /* predict 3'-complete */
        CP3,
        /* consider only introns with canonical splice sites */
        CSS,
        /* output flank type, ie 'constitutive' or 'alternative' */
        FLT,
        /* acceptable introns */
        IOK,
        /* predict NMD */
        NMD,
        /* output splice site sequences of event flanks */
        SEQ,
    };

    /**
     * Parameter collecting flags for event output options
     */
    public static final Parameter<EnumSet<EventOptions>> EVENTS_OPT = Parameters.enumSetParameter(
            "EVENTS_OPT",
            "Toggle optional attributes to be output",
            EnumSet.noneOf(EventOptions.class),
            EventOptions.class,
            null).longOption("ep").shortOption('p');



    /**
     * Enumeration of the different types for splice sites.
     */
    public static enum SiteTypes {
        /** Splice Site Donor */
        SSD,
        /** Splice Site Acceptor */
        SSA,
        /** Transcription Start Site */
        TSS,
        /** Cleavage Site */
        CLV,
        /** Soft Start */
        SST,
        /** Soft End */
        SND,
        /** Start Codon */
        AUG,
        /** Stop Codon */
        STP,
    };

    /**
     * Parameter for the list of site types that is output.
     */
    public static final Parameter<EnumSet<SiteTypes>> SITES = Parameters.enumSetParameter(
            "SITES",
            "Types of sites that are output",
            EnumSet.noneOf(SiteTypes.class),
            SiteTypes.class,
            null).longOption("ss").shortOption('s');

    /**
     * Flags to control output options for sites:
     * <ul>
     *     <li>SSS: compute scores for splice sites</li>
     * </ul>
     */
    public static enum SiteOptions {
        /** compute splice site score */
        SSS,
    };

    /**
     * Parameter collecting flags for site output options
     */
    public static final Parameter<EnumSet<SiteOptions>> SITES_OPT = Parameters.enumSetParameter(
            "SITES_OPT",
            "Toggle optional site attributes to be output",
            EnumSet.noneOf(SiteOptions.class),
            SiteOptions.class,
            null).longOption("sp").shortOption('t');

    /**
     * Checks whether a folder with genomic sequences is necessary in order
     * to complete all tasks specified with <code>this</code> parameter set.
     * @return <code>true</code> if the genomic sequence is required given the
     * current parameters, <code>false</code> otherwise
     */
    public boolean requiresGenomicSequence() {
        if (get(AStalavistaSettings.EVENTS_OPT).contains(EventOptions.SEQ))
            return true;
        if (get(AStalavistaSettings.SITES_OPT).contains(EventOptions.CSS)
                || get(AStalavistaSettings.EVENTS_OPT).contains(EventOptions.IOK))
            return true;

        return false;
    }

}
