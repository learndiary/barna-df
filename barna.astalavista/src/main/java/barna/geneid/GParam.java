package barna.geneid;

/**
 * Created with IntelliJ IDEA.
 * User: micha
 * Date: 9/8/12
 * Time: 1:39 PM
 * To change this template use File | Settings | File Templates.
 */
public class GParam {

    public static final int MAXENTRY= 97;
    public static final int FRAMES= 3;


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
    Profile DonorProfile;
    Profile U2gcagDonorProfile;
    Profile U2gtaDonorProfile;
    Profile U2gtgDonorProfile;
    Profile U2gtyDonorProfile;
    Profile U12gtagAcceptorProfile;
    Profile U12BranchPointProfile;
    Profile U12gtagDonorProfile;
    Profile U12atacAcceptorProfile;
    Profile U12atacDonorProfile;
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
    float U12_SPLICE_SCORE_THRESH = -1000;
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
