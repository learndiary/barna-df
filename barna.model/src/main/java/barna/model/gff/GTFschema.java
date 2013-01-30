package barna.model.gff;

import barna.commons.parameters.Parameter;
import barna.commons.parameters.ParameterSchema;
import barna.commons.parameters.Parameters;

import java.util.EnumSet;

/**
 * Class to reuse I/O features across tools, e.g.
 * in Astalavista, Flux, etc.
 *
 */
public class GTFschema extends ParameterSchema {
    /**
     * Parameter values to specify attributes for outputting loci.
     */
    public static enum OutputLocus {
        LOCUSGEOMETRY  // for locus lines, geometrical attributes are written
    };

    /**
     * Parameter values to specify attributes for outputting sites.
     */
    public static enum OutputSite {
        SITESEQ,
        SITESCORES
    };

    /**
     * Parameter values to specify attributes for outputting events.
     */
    public static enum OutputEvent {
        ASEVENT,
        DSEVENT,
        VSEVENT,
        SITESCORES
    }

    /**
     * Parameter for specifying the features and attributes output as loci.
     * If locus features are not to be output, the parameter is set to
     * <code>EnumSet.noneOf(OutputLocus.class)</code> (default).
     */
    public static final Parameter<EnumSet<OutputLocus>> OUTPUT_LOCUS = Parameters.enumSetParameter(
            "OUTPUT_LOCUS",
            "features and attributes to output for loci",
            EnumSet.noneOf(OutputLocus.class),
            OutputLocus.class,
            null);

    /**
     * Parameter for specifying the features and attributes output as sites.
     * If site features are not to be output, the parameter is set to
     * <code>EnumSet.noneOf(OutputSite.class)</code> (default).
     */
    public static final Parameter<EnumSet<OutputSite>> OUTPUT_SITE = Parameters.enumSetParameter(
            "OUTPUT_SITE",
            "features and attributes to output for sites",
            EnumSet.noneOf(OutputSite.class),
            OutputSite.class,
            null);

    /**
     * Parameter for specifying which features and attributes are output as events.
     * If event features are not to be output, the parameter is set to
     * <code>EnumSet.noneOf(OutputEvent.class)</code>. Default is to
     * output AS events.
     */
    public static final Parameter<EnumSet<OutputEvent>> OUTPUT_EVENT = Parameters.enumSetParameter(
            "OUTPUT_EVENT",
            "features and attributes to output for events",
            EnumSet.of(OutputEvent.ASEVENT),
            OutputEvent.class,
            null).shortOption('e');

}
