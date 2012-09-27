package barna.geneid;

/**
 * Created with IntelliJ IDEA.
 * User: micha
 * Date: 9/8/12
 * Time: 7:03 PM
 * To change this template use File | Settings | File Templates.
 */
public class GeneIDsettings {

    /* geneid setup flags  (from geneid.c) */

    /* sites to print */
    int SFP=0, SDP=0, SAP=0, STP=0;

    /* exons to print */
    int EFP=0, EIP=0, ETP=0, EXP=0, ESP=0, EOP = 0;

    /* introns to print */
    int PRINTINT = 0;

    /* Partial or full prediction engine */
    int GENAMIC = 1, GENEID = 1;

    /* Only forward or reverse prediction engine */
    int FWD=1, RVS=1;

    /* switch ORF prediction on */
    int scanORF = 0;

    /* Input annotations or homology to protein information/reads to UTR prediction */
    int EVD = 0, SRP = 0;
    static int UTR=0;

    /* Output formats */
    int GFF = 0, GFF3 = 0, X10 = 0, XML = 0, cDNA = 0, PSEQ = 0, tDNA = 0;

    /* Verbose flag (memory/processing information) */
    int BEG=0, VRB=0;

    /* Score for regions not-supported by protein homology */
    static int NO_SCORE;

    /* Force single prediction: 1 gene */
    static int SGE= 0;

    /* Length of flank around exons to subtract background RNA-seq signal */
    static int BKGD_SUBTRACT_FLANK_LENGTH = 0;



    /* User defined lower limit */
    long LOW=0;

    /* User defined upper limit */
    long HI=0;

    /* Millions of reads mapped */
    float MRM= 15f;

    /* Optional Predicted Gene Prefix */
    String genePrefix= "";


    /* Generic maximum values: sites, exons and backup elements */
    long NUMSITES,NUMEXONS,MAXBACKUPSITES,MAXBACKUPEXONS,NUMU12SITES,NUMU12EXONS,NUMU12U12EXONS;

    /* Accounting time and results */
    long[] m;   //account *m;



}
