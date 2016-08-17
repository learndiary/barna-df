package barna.astalavista;

import barna.commons.parameters.*;
import barna.model.Transcript;

import java.io.File;
import java.util.EnumSet;

/**
 * Created with IntelliJ IDEA.
 * User: micha
 * Date: 2/5/13
 * Time: 2:56 PM
 */
public class AStaSettings extends AStalavistaSettings {

    /**
     * Level of intron confidence, below which introns are trusted without checks.
     * The default is to trust all introns (i.e., ic= 255). Introns are assigned a
     * confidency class:
     * <ul>
     * <li>0 for 'RefSeq' appears in the source field of the annotation</li>
     * <li>1 for 'mRNA' appears in the source field of the annotation</li>
     * <li>2 for 'EST' appears in the source field of the annotation</li>
     * </ul>
     * All introns in transcripts of confidence level &gt; threshold are discarded.
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
     * All transcript edges of confidence level &gt; edgeConfidence will be extended in case the
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
     * Path to the GTF output file for events.
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
     * Dimension of the AS events to be extracted,
     * retrieves 'complete' events &trade; for parameter
     * values &lt; 2.
     */
    public static final Parameter<Integer> EVENTS_DIMENSION = Parameters.intParameter("EVENTS_DIMENSION",
            "Dimension of the AS events to be extracted, retrieves 'complete' events <TM>\n" +
                    "for parameter values < 2",
            2, new ParameterValidator() {
        @Override
        public void validate(ParameterSchema schema, Parameter parameter) throws ParameterException {
            int v = (Integer) schema.get(parameter);
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
     * to alternative 5&#8217;- or 3&#8217; transcript structures.</li>
     * <li>extending transcript structures by additional (splice) sites
     * to the 5&#8217;- or the 3&#8217;-end</li>
     * <li>DS= additional splicing that are flanked by a common site and</li>
     * <li>VS= variable sites is any other form of sites that differ
     * between overlapping transcript structures</li>
     * </ul>
     * @see barna.model.ASEvent#getType()
     * @see <a href="http://www.ploscompbiol.org/article/info%3Adoi%2F10.1371%2Fjournal.pcbi.1000147">http://www.ploscompbiol.org/article/info%3Adoi%2F10.1371%2Fjournal.pcbi.1000147</a>
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
    public static final Parameter<EnumSet<EventOptions>> EVENTS_ATR = Parameters.enumSetParameter(
            "EVENTS_ATR",
            "Toggle optional attributes to be output",
            EnumSet.noneOf(EventOptions.class),
            EventOptions.class,
            null).longOption("ea").shortOption('a');

    /**
     * Checks whether a folder with genomic sequences is necessary in order
     * to complete all tasks specified with <code>this</code> parameter set.
     * @return <code>true</code> if the genomic sequence is required given the
     * current parameters, <code>false</code> otherwise
     */
    public boolean requiresGenomicSequence() {
        if (super.requiresGenomicSequence()
            || get(AStaSettings.EVENTS_ATR).contains(EventOptions.CSS)
            || get(AStaSettings.EVENTS_ATR).contains(EventOptions.SEQ)
            || get(AStaSettings.EVENTS_ATR).contains(EventOptions.IOK))
            return true;
        // otherwise
        return false;
    }

}
