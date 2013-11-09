package barna.io.sam;

/**
 * Wrapper class for SAM file format specific constants.
 * User: micha
 */
public class SAMConstants {

    /**
     * Number of reported alignments that contains the query in the current record.
     * Although optional, the <code>NM</code> tag should be present according to the
     * specification.
     * Example: <code>NM:i:2</code> corresponds to 2 mismatches.
     */
    public static String SAM_OPTION_NH= "NH";

    /**
     * <code>XT:A:U</code> identifies unique mappings.
     * Note that tags starting with `X', `Y' and `Z' or tags containing lowercase letters
     * in either position are reserved for local use and will not be formally dened in any future version of
     * this specication.
     */
    public static String SAM_OPTION_XT= "XT";

    /**
     * Edit distance to the reference, including ambiguous bases but excluding clipping.
     */
    public static String SAM_OPTION_NM= "NM";


    /**
     * <code>XS:A:+/-</code> identifies the strand of the mapping.
     * Note that tags starting with `X', `Y' and `Z' or tags containing lowercase letters
     * in either position are reserved for local use and will not be formally dened in any future version of
     * this specication.
     */
    public static String SAM_OPTION_XS= "XS";

}
