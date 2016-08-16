package barna.scorer;

import barna.astalavista.AStalavistaSettings;
import barna.commons.parameters.*;

import java.io.File;
import java.util.EnumSet;

/**
 * Created with IntelliJ IDEA.
 * User: micha
 * Date: 2/5/13
 * Time: 2:51 PM
 */
public class ScorerSettings extends AStalavistaSettings {

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
    }).longOption("vcf"); // .shortOption('v'); // blocked by version

    /**
     * Path to the VCF output file.
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
            EnumSet.of(SiteOptions.SSS),
            SiteOptions.class,
            null).longOption("sp"); // .shortOption('t'); // blocked by tool

    /**
     * Checks whether a folder with genomic sequences is necessary in order
     * to complete all tasks specified with <code>this</code> parameter set.
     * @return <code>true</code> if the genomic sequence is required given the
     * current parameters, <code>false</code> otherwise
     */
    @Override
    public boolean requiresGenomicSequence() {
        if (super.requiresGenomicSequence()
            || get(ScorerSettings.SITES_OPT).contains(SiteOptions.SSS))
            return true;
        // otherwise
        return false;
    }


}
