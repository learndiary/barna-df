package barna.geneid;

import java.text.DecimalFormat;

/**
 * Created with IntelliJ IDEA.
 * User: micha
 * Date: 9/8/12
 * Time: 6:07 PM
 * To change this template use File | Settings | File Templates.
 */
public class GeneIDconstants {

    /**
     * Infinity: score functions
     */
    public static final double INF= Double.POSITIVE_INFINITY;   // 1.7976931348623157E+308;

    /**
     * Array range in C: 0..N-1
     */
    public static final int COFFSET= 1;

    /**
     * Forward strand
     */
    public static final int FORWARD= 0;

    /**
     * Reverse strand
     */
    public static final int REVERSE= 1;

    /**
     * Number of different coding frames.
     */
    public static final int FRAMES= 3;

    /**
     * Dictionary definitions (hash)
     */
    public static final int MAXTYPE= 50;

    /**
     * Maximum length for strings (mess)
     */
    public static final int MAXSTRING= 600;

    /**
     * Infinity: positions in sequence
     */
    public static final int INFI= Integer.MAX_VALUE; // 2147483647;

    /**
     * Infinity: word in Gene model
     */
    public static final String SINFI= "Infinity";

    /* Mark rules up as blocking in Gene model  */
    public static final int BLOCK= 1;
    public static final int NONBLOCK= 0;

    /* Maximum number of isochores              */
    public static final int MAXISOCHORES= 4;

    /* Name of default parameter file           */
    public static final String PARAMETERFILE= "/human.070704.param.param";   // "param.default";

    public static final DecimalFormat DF_92F= new DecimalFormat("#########.##");
    public static final DecimalFormat DF_93F= new DecimalFormat("#########.###");

    /**
     * Maximum profile dimension (length)
     */
    public static final int PROFILEDIM= 100;

    /**
     * ACCEPTOR_CONTEXT 25
     */
    public static final int ACCEPTOR_CONTEXT= 50;
    public static final int PPT_ACC_MAXDIST= 14;

    /**
     * used to be 15.
     */
    public static final int U12BP_PENALTY_SCALING_FACTOR= 6;

    /**
     * used to be 15
     */
    public static final int U2BP_PENALTY_SCALING_FACTOR= 0;

    /**
     * Minimum distance between branch point
     * and acceptor site
     */
    public static final int MIN_U12BPACC_DIST= 7;
    public static final int MIN_U2BPACC_DIST= 15;
    public static final int OPT_U12BP_DIST= 12;
    public static final int OPT_U2BP_DIST= 25;

    public static final String sACC= "Acceptor";
    public static final String sDON= "Donor";
    public static final String sSTA= "Start";
    public static final String sSTO= "Stop";
    public static final String sPOL= "PolyA";
    public static final String sTSS= "TSS";
    public static final String sTES= "TES";
    public static final String sPPT= "PolyPyrimidineTract";
    public static final String sBP= "BranchPoint";

    /* Intron Types                           */
    public static final String sU2type= "U2";
    public static final String sU12type= "U12";

    /* Header profiles                          */
    public static final String sMarkov= "Markov_order";
    public static final String sprofilePolyA= "PolyA_Signal_profile";
    public static final String sprofilePPT= "Poly_Pyrimidine_Tract_profile";
    public static final String sprofileBP=  "Branch_point_profile";
    public static final String sprofileACC= "Acceptor_profile";
    public static final String sprofileU12BP=  "U12_Branch_point_profile";
    public static final String sprofileU12gtagACC= "U12gtag_Acceptor_profile";
    public static final String sprofileU12atacACC= "U12atac_Acceptor_profile";
    public static final String sprofileDON= "Donor_profile";
    public static final String sprofileU2gcagDON= "U2gcag_Donor_profile";
    public static final String sprofileU12gtagDON= "U12gtag_Donor_profile";
    public static final String sprofileU12atacDON= "U12atac_Donor_profile";
    public static final String sprofileU2gtaDON= "U2gta_Donor_profile";
    public static final String sprofileU2gtgDON= "U2gtg_Donor_profile";
    public static final String sprofileU2gtyDON= "U2gty_Donor_profile";

    /* U12 splice site (sum d + a) and exon (sum d + a) score thresholds    */
    public static final String sU12_SPLICE_SCORE_THRESH= "U12_Splice_Score_Threshold";
    public static final String sU12_EXON_SCORE_THRESH= "U12_Exon_Score_Threshold";

    /* U12 acceptor splice site initial threshold factor (added to final threshold given in param file) */
    public static final int U12ACC_CUTOFF_FACTOR= -7;

    /* Header exon weights                          */
    public static final String sExon_weights= "Exon_weights";
    public static final String sU12_EXON_WEIGHT= "U12_Exon_weight";

    /* Header evidence factor and weight */
    public static final String sEVIDENCEF= "Evidence_Factor";
    public static final String sEVIDENCEW= "Evidence_Exon_Weight";
    public static final String sBKGD_SUBTRACT_FLANK_LENGTH= "BKGD_SUBTRACT_FLANK_LENGTH";
    public static final String sNUMISO= "number_of_isochores";

    /* Recursive splice site thresholds */
    public static final int RDT= 4; /*donor*/
    public static final int RAT= 4; /*acceptor*/
    public static final String sRSSMARKOVSCORE= "RSS_Markov_Score";
    public static final String sRSS_DONOR_SCORE_CUTOFF= "RSS_Donor_Score_Cutoff";
    public static final String sRSS_ACCEPTOR_SCORE_CUTOFF= "RSS_Acceptor_Score_Cutoff";

    /* Field group in gff: Not grouped exons    */
    public static final String NOGROUP= "NON_GROUPED";
    public static final int NOVALUE= 0;

    /* Maximum oligo (word) length (Markov)     */
    public static final int OLIGOLENGTH= 10;

    /* Length of every processed fragment       */
    public static final int LENGTHSi= 220000;

}
