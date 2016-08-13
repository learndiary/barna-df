package barna.geneid;

/**
 * Created with IntelliJ IDEA.
 * User: micha
 * Date: 9/8/12
 * Time: 1:39 PM
 * To change this template use File | Settings | File Templates.
 */
public class GParam {

    /**
     * Lazily compute the number of nucleotides <b>upstream</b> of a <u>donor</u> site
     * that are required to compute its score according to the profile.
     * @param profile a donor profile
     * @return number of nucleotides to be read upstream of the splice site
     */
    public static int getDonorFlank5(Profile profile) {
        //return profile.getOffset()+ profile.getOrder();
        return profile.getOffset()+ profile.getOrder()+ 1; // 1st order -> di-nucleotides (+1)
    }

    /**
     * Lazily compute the number of nucleotides <b>downstream</b> of a <u>donor</u> site
     * that are required to compute its score according to the profile.
     * @param profile a donor profile
     * @return number of nucleotides to be read downstream of the splice site
     */
    public static int getDonorFlank3(Profile profile) {
        //return profile.getDimension()- profile.getOffset()- profile.getOrder()- 2;
        return profile.getDimension()- getDonorFlank5(profile)- 2+ profile.getOrder();  // (-2) donor di-nucleotide
    }

    /**
     * Lazily compute the number of nucleotides <b>upstream</b> of an <u>acceptor</u> site
     * that are required to compute its score according to the profile.
     * @param profile an acceptor profile
     * @return number of nucleotides to be read upstream of the splice site
     */
    public static int getAcceptorFlank5(Profile profile) {
        return profile.getOffset()- 1;
    }

    /**
     * Lazily compute the number of nucleotides <b>downstream</b> of an <u>acceptor</u> site
     * that are required to compute its score according to the profile.
     * @param profile an acceptor profile
     * @return number of nucleotides to be read downstream of the splice site
     */
    public static int getAcceptorFlank3(Profile profile) {
        return profile.getDimension()- getAcceptorFlank5(profile)- 2+ profile.getOrder();
    }

    public static final int MAXENTRY= 97;
    public static final int FRAMES= 3;

    public Profile getU2gcagDonorProfile() {
        return U2gcagDonorProfile;
    }

    public Profile getU2gtaDonorProfile() {
        return U2gtaDonorProfile;
    }

    public Profile getU2gtgDonorProfile() {
        return U2gtgDonorProfile;
    }

    public Profile getU2gtyDonorProfile() {
        return U2gtyDonorProfile;
    }

    public Profile getU12gtagAcceptorProfile() {
        return U12gtagAcceptorProfile;
    }

    public Profile getU12gtagDonorProfile() {
        return U12gtagDonorProfile;
    }

    public Profile getU12atacAcceptorProfile() {
        return U12atacAcceptorProfile;
    }

    public Profile getU12atacDonorProfile() {
        return U12atacDonorProfile;
    }


    public static class ParamExons {
        float siteFactor;

        float exonFactor;
        float OligoCutoff;

        float HSPFactor;

        float ExonWeight;
/*   float U12atacExonWeight; */
/*   float U12gtagExonWeight; */
        float ExonCutoff;
    }


    int leftValue;
    int rightValue;

    Profile PolyASignalProfile;
    Profile StartProfile;
    Profile AcceptorProfile;
    Profile PolyPTractProfile;
    Profile BranchPointProfile;

    public Profile getAcceptorProfile() {
        return AcceptorProfile;
    }

    public Profile getDonorProfile() {
        return DonorProfile;
    }

    Profile DonorProfile;
    protected Profile U2gcagDonorProfile;
    protected Profile U2gtaDonorProfile;
    protected Profile U2gtgDonorProfile;
    protected Profile U2gtyDonorProfile;
    protected Profile U12gtagAcceptorProfile;
    Profile U12BranchPointProfile;
    protected Profile U12gtagDonorProfile;
    protected Profile U12atacAcceptorProfile;
    protected Profile U12atacDonorProfile;
    Profile StopProfile;

    float[][] OligoLogsIni= new float[3][];
    float[][] OligoLogsTran= new float[3][];

    long OligoDim;
    long OligoDim_1;
    int  OligoLength;

    ParamExons Initial= new ParamExons();
    ParamExons Internal= new ParamExons();
    ParamExons Terminal= new ParamExons();
    ParamExons Single= new ParamExons();
    ParamExons utr= new ParamExons();

    float[][] OligoDistIni= new float[FRAMES][];
    float[][] OligoDistTran= new float[FRAMES][];

    int MaxDonors;

    Dictionary D;
    int[]  nc;
    int[]  ne;
    long[] md;
    long[] Md;
    int[][] UC= new int[MAXENTRY][MAXENTRY];
    int[][] DE= new int[MAXENTRY][MAXENTRY];
    int[] block= new int[MAXENTRY];
    int nclass;

    /* Detection of recursive splice sites */
    float RSSMARKOVSCORE = 0;
    float RSSDON = GeneIDconstants.RDT;
    float RSSACC = GeneIDconstants.RAT;

    /* Increase/decrease exon weight value (exon score) */
    float EvidenceEW = 0;
    float EvidenceFactor = 1;
    float U12EW = 0;

    public float getU12SpliceScoreThresh() {
        return u12SpliceScoreThresh;
    }

    float u12SpliceScoreThresh = Float.NaN;
    float U12_EXON_SCORE_THRESH = -1000;
    float EW = GeneIDconstants.NOVALUE;

    /* Detection of PolyPTracts in Acceptors */
    static int PPT=0;

    /* Detection of BranchPoints in Acceptors */
    static int BP=0;

    /* Detection of recursive splice sites */
    int RSS=0;

    /* Detection of U12 introns */
    int U12=0;

    /* Detection of U12gtag sites (acceptor uses BranchPoint)*/
    int U12GTAG=0;

    /* Detection of U12atac sites (acceptor uses BranchPoint)*/
    int U12ATAC=0;

    /* Detection of U2gcag sites */
    int U2GCAG=0;

    /* Detection of U2gta donor sites */
    int U2GTA=0;

    /* Detection of U2gtg donor sites */
    int U2GTG=0;

    /* Detection of U2gty donor sites */
    int U2GTY=0;

    /* Detection of PolyA Signal */
    int PAS=0;

    /* Splice classes: the number of compatible splice site combinations used in genamic for joining exons */
    short SPLICECLASSES = 1;

    /**
     * from RequestMemory.c function RequestMemoryParams()
     */
    public GParam() {
        /* 0. Main structure: gparam */

        /* 1. Profiles for signals */

        /* 2. Markov model: initial and transition values */
        int OligoDim = (int) Math.pow(4,  GeneIDconstants.OLIGOLENGTH);
        OligoLogsIni[0]= new float[OligoDim];
        OligoLogsIni[1]= new float[OligoDim];
        OligoLogsIni[2]= new float[OligoDim];
        OligoLogsTran[0]= new float[OligoDim];
        OligoLogsTran[1]= new float[OligoDim];
        OligoLogsTran[2]= new float[OligoDim];

        /* 3. Markov temporary data structures to compute every split: LENGTHSi */
        OligoDistIni[0]= new float[GeneIDconstants.LENGTHSi];
        OligoDistIni[1]= new float[GeneIDconstants.LENGTHSi];
        OligoDistIni[2]= new float[GeneIDconstants.LENGTHSi];
        OligoDistTran[0]= new float[GeneIDconstants.LENGTHSi];
        OligoDistTran[1]= new float[GeneIDconstants.LENGTHSi];
        OligoDistTran[2]= new float[GeneIDconstants.LENGTHSi];

        /* 3. Markov temporary data structures to compute every split: LENGTHSi */
        OligoDistIni[0]= new float[GeneIDconstants.LENGTHSi];
        OligoDistIni[1]= new float[GeneIDconstants.LENGTHSi];
        OligoDistIni[2]= new float[GeneIDconstants.LENGTHSi];
        OligoDistTran[0]= new float[GeneIDconstants.LENGTHSi];
        OligoDistTran[1]= new float[GeneIDconstants.LENGTHSi];
        OligoDistTran[2]= new float[GeneIDconstants.LENGTHSi];

        /* 4. Exons score parameters */


        /* Allocating space for global parameters (gene model) */
        nc= new int[MAXENTRY];
        ne= new int[MAXENTRY];
        md= new long[MAXENTRY];
        Md= new long[MAXENTRY];

    }

}
