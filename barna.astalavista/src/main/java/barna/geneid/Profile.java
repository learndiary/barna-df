package barna.geneid;

import barna.astalavista.AStalavista;
import barna.commons.log.Log;
import barna.scorer.ScorerSettings;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.text.DecimalFormat;

/**
 * Created with IntelliJ IDEA.
 * User: micha
 * Date: 9/8/12
 * Time: 10:58 AM
 * To change this template use File | Settings | File Templates.
 */
public class Profile {

    /**
     * Length of the profile as by <code>(offset + 2)</code> plus
     * the number of nucleotides that are considered after the
     * di-nucleotide.
     */
    int    dimension;

    public int getAcc_context() {
        return acc_context;
    }

    public float getAfactor() {
        return afactor;
    }

    public float getBfactor() {
        return bfactor;
    }

    public float getCutoff() {
        return cutoff;
    }

    public int getDimension() {
        return dimension;
    }

    public int getDimensionTrans() {
        return dimensionTrans;
    }

    public int getDist() {
        return dist;
    }

    public int getOffset() {
        return offset;
    }

    public int getOpt_dist() {
        return opt_dist;
    }

    public int getOrder() {
        return order;
    }

    public float getPenalty_factor() {
        return penalty_factor;
    }

    public float[][] getTransitionValues() {
        return transitionValues;
    }

    /**
     * Position of the di-nucleotide in the profile, i.e.,
     * number of exonic positions before the "GT" of donors
     * and number of intronic positions before the "AG" of
     * acceptors.
     */
    int    offset;

    float  cutoff;
    int    order;
    float  afactor= 0f;
    float  bfactor= 1f;
    int    acc_context= GeneIDconstants.ACCEPTOR_CONTEXT;
    int    dist= GeneIDconstants.MIN_U12BPACC_DIST;
    int    opt_dist= GeneIDconstants.OPT_U12BP_DIST;
    float    penalty_factor= (float) GeneIDconstants.U12BP_PENALTY_SCALING_FACTOR;

    int dimensionTrans; // reduced from long
    float[][]  transitionValues= new float[GeneIDconstants.PROFILEDIM][];


    public Profile(int order, int dimension) {

        this.order= order;
        this.dimension= dimension;

        /* Transition probabilities in every position of the PWA */
        dimensionTrans = (int) Math.pow(5, this.order+ 1);

        for(int i= 0; i < this.dimension; i++)
            transitionValues[i] = new float[dimensionTrans];

    }

    /**
     * @deprecated obviously does the same as readLine()
     * @param rootFile
     * @return
     * @throws Exception
     */
    static String readHeader(BufferedReader rootFile) throws Exception {
        return readLine(rootFile);
    }

    /**
    * Read numeric values or headers of values:
    * read one line skipping comments and empty lines.
    */
    static String readLine(BufferedReader rootFile) throws Exception {

        String res= null;
        try {
            if (rootFile.ready())
                res= rootFile.readLine();
            while((res!= null)&& (res.startsWith("#")|| res.trim().length()== 0))
                res= rootFile.readLine();
        } catch (Exception e) {
            res= null;
        }

        if (res == null)
            throw new RuntimeException("Parameter file: unexpected end of reading");
        return res;
    }

    /**
     * Read information useful to predict Start codons, donors, stop codons
     *
     * @param rootFile
s     * @param signal sequence to score
     * @param H 0 for not reading header, 1 for reading the header
     */
    static Profile readProfile(BufferedReader rootFile, String signal, int H) throws Exception {

        // Definition parameters: Length, offset, cutoff and order (Markov chain)
        String line= null;
        if (H==1)
            line= readHeader(rootFile);

        line= readLine(rootFile);
        String[] ss= line.split("\\s");
        if (ss.length< 4) {
            String mess= "Wrong format: Definition parameters in "+ signal+" profile";
            throw new RuntimeException(mess);
        }


        int dimension= Integer.parseInt(ss[0]);
        int order= Integer.parseInt(ss[3]);

        /* Memory to allocate the data with these fixed dimensions */
        Profile p= new Profile(order, dimension);

        p.offset= Integer.parseInt(ss[1]);
        p.cutoff= Float.parseFloat(ss[2]);

        if (ss.length>= 6) {
            p.afactor= Float.parseFloat(ss[4]);
            p.bfactor= Float.parseFloat(ss[5]);
        }

        if (ss.length>= 10) {
            p.acc_context= Integer.parseInt(ss[6]);
            p.dist= Integer.parseInt(ss[7]);
            p.opt_dist= Integer.parseInt(ss[8]);
            p.penalty_factor= Float.parseFloat(ss[9]);
        }

        /* Useful to check everything is OK */
        Log.info("Reading... "+ signal+ ": "+ p.toString());

        /* Prepare and read values */
        p.setProfile(rootFile, signal);

        return p;
    }


    void setProfile(BufferedReader rootFile, String signal) throws Exception {

        // TODO read models of different order in a systematic way

        /* According to the order of Markov chain, select a different method */
        /* Position weight array: transition probabilities in every position */
        switch(order) {

            case 0:
                /* 5 combinations / pos */
                for(int i= 0; i < dimension; i++) {
                    readACGT(i, 0, rootFile, signal);
                }
                break;

            case 1:
                /* 25 combinations / pos */
                for(int i= 0; i < dimension; i++) {

                    /* Reading AX,CX,GX and TX: creating XN */
                    int jlen= (dimensionTrans- 5);
                    for(int j= 0; j < dimensionTrans- 5; j= j+5) {
                        readACGT(i, j, rootFile, signal);
                    }

                    /* Creating NA,NC,NG,NT transition values */
                    for(int x= 0; x < 4; x++) {
                        transitionValues[i][jlen+ x] =
                                (transitionValues[i][x] +
                                        transitionValues[i][5+x] +
                                        transitionValues[i][10+x] +
                                        transitionValues[i][15+x]) / 4;
                    }

                    /* Creating the value NN */
                    transitionValues[i][(dimensionTrans- 5)+ 4] =
                            (transitionValues[i][jlen] +
                                    transitionValues[i][jlen+ 1] +
                                    transitionValues[i][jlen+ 2] +
                                    transitionValues[i][jlen+ 3]) / 4;
                }
                break;

            case 2:
                /* 125 combinations / pos: there are "dimension" positions */
                for(int i= 0; i < dimension; i++) {

                    /* Reading AXX,CXX,GXX,TXX and creating ANX,CNX,GNX,TNX */
                    int jlen= dimensionTrans- 25;
                    for(int j= 0; j < jlen; j= j+ 25) {

                        /* 20 = 5^order - 5 */
                        int ylen= 20;
                        for(int y= 0; y < ylen; y= y+ 5) {
                            readACGT(i, j+ y, rootFile, signal);
                        }

                        /* Creating XNA,XNC,XNG,XNT */
                        for(int x=0; x < 4; x++) {
                            transitionValues[i][j+ ylen+ x] =
                                    (transitionValues[i][j+x] +
                                            transitionValues[i][j+5+x] +
                                            transitionValues[i][j+10+x] +
                                            transitionValues[i][j+15+x]) / 4;
                        }

                        /* Creating the value XNN */
                        transitionValues[i][j+ 4+ ylen] =
                                (transitionValues[i][j+ ylen] +
                                        transitionValues[i][j+ 1+ ylen] +
                                        transitionValues[i][j+ 2+ ylen] +
                                        transitionValues[i][j+ 3+ ylen]) / 4;
                    }

                    /* Creating NAX,NCX,NGX,NTX (j=100)*/
                    int ylen= 20;
                    for(int y=0; y < ylen; y=y+5) {

                        /* Creating NAX... */
                        int xlen= 4;
                        for(int x= 0; x < xlen; x++) {
                            transitionValues[i][jlen+ y+ x] =
                                    (transitionValues[i][y+ x] +
                                            transitionValues[i][25+ y+ x] +
                                            transitionValues[i][50+ y+ x] +
                                            transitionValues[i][75+ y+ x]) / 4;
                        }

                        /* Computing NAN (x=4)*/
                        transitionValues[i][jlen+ y+ xlen] =
                                (transitionValues[i][jlen+ y] +
                                        transitionValues[i][jlen+ y+ 1] +
                                        transitionValues[i][jlen+ y+ 2] +
                                        transitionValues[i][jlen+ y+ 3]) / 4;
                    }

                    /* Creating NNX (j=100,y=20)*/
                    int xlen= 4;
                    for(int x=0; x < xlen; x++) {
                        transitionValues[i][jlen+ ylen+ x] =
                                (transitionValues[i][jlen+ x] +
                                        transitionValues[i][jlen+ 5+ x] +
                                        transitionValues[i][jlen+ 10+ x] +
                                        transitionValues[i][jlen+ 15+ x]) /4;
                    }

                    /* Finally, creating NNN (j=100,y=20,x=4) */
                    transitionValues[i][jlen+ ylen+ xlen] =
                            (transitionValues[i][jlen+ ylen] +
                                    transitionValues[i][jlen+ ylen+ 1] +
                                    transitionValues[i][jlen+ ylen+ 2] +
                                    transitionValues[i][jlen+ ylen+ 3]) /4;
                }
                break;

            default:
                throw new RuntimeException("Sorry, Markov order higher than 2 not implemented yet");
                //break;
        }
    }

    private void readACGT(int i, int j, BufferedReader rootFile, String signal) throws Exception {

        /* Reading A,C,G and T: creating N transition */
        /* Reading AX,CX,GX and TX: creating XN */
        /* Reading AXX,CXX,GXX,TXX and creating ANX,CNX,GNX,TNX */
        for (int x = 0; x < 4; x++) {
            String line= readLine(rootFile);
            String[] ss= line.split(" ");
            try {
                // if ((sscanf(line,"%*d %*s %f", &(p.transitionValues[i][j])))!=1)
                transitionValues[i][j+ x]= Float.parseFloat(ss[2]);
            } catch (Exception e) {
                String mess= "Wrong format: Transition values in "+ signal+ " profile";
                throw new RuntimeException(mess, e);
            }
        }

        /* Generating the N value */
        /* Generating the XN value every iteration: AN,CN,GN,TN */
        /* Generating the XXN value (5th value): XAN,XCN,XGN,XTN */
        transitionValues[i][j+ 4] =
                (transitionValues[i][0] +
                        transitionValues[i][1] +
                        transitionValues[i][2] +
                        transitionValues[i][3]) / 4;

    }

    @Override
    public String toString() {

        String mess= dimension+ "\t"+ offset+ "\t"+ order+ "\t"+ dimensionTrans
                    + "\t" + new DecimalFormat("#####.##").format(cutoff);
        return mess;
    }

    /**
     * Read information useful to predict Acceptor splice sites
     */
    static void readProfileSpliceSites(BufferedReader rootFile, GParam gp) throws Exception {

        Profile p;
        int u12bp=0;
        int u12gtagAcc=0;
        int u12atacAcc=0;
        int u12gtagDon=0;
        int u12atacDon=0;

        /* A. Optional profiles: U12GTAG, U12ATAC Donor, acceptor and branch points
and U2 branch points and Poly Pyrimidine Tract */

        String line= readHeader(rootFile);
        String header= line.trim();
        if (header.length()== 0) {
            String mess= "Wrong format: header in optional profile for splice sites";
            throw new RuntimeException(mess);
        }

        while(!header.equalsIgnoreCase(GeneIDconstants.sprofileACC)) {

            /* Read optional profiles: sprofilePPT,sprofileBP,sprofileU12BP,sprofileU12gtagACC,sprofileU12atacACC */

            if (header.equalsIgnoreCase(GeneIDconstants.sprofileU12BP)) {
                u12bp++;

                /* Reading the U12 Branch Point profile */
                gp.U12BranchPointProfile= readProfile(rootFile, GeneIDconstants.sBP, 0);

            } else {

                if (header.equalsIgnoreCase(GeneIDconstants.sprofileU12gtagACC)) {
                    u12gtagAcc++;

                    /* Reading the U12gtag acceptor profile */
                    gp.U12gtagAcceptorProfile= readProfile(rootFile, GeneIDconstants.sACC, 0);

                } else {

                    if (header.equalsIgnoreCase(GeneIDconstants.sprofileU12atacACC)) {
                        u12atacAcc++;

                        /* Reading the U12atac acceptor profile */
                        gp.U12atacAcceptorProfile= readProfile(rootFile, GeneIDconstants.sACC, 0);

                    } else {

                        if (header.equalsIgnoreCase(GeneIDconstants.sprofilePPT)) {

                            /* Switch on the acceptor prediction using PPTs */
                            gp.PPT++;

                            /* Reading the Poly Pyrimidine Tract profile */
                            gp.PolyPTractProfile= readProfile(rootFile, GeneIDconstants.sPPT, 0);

                        } else {

                            if (header.equalsIgnoreCase(GeneIDconstants.sprofileBP)) {
                                /* Switch on the acceptor prediction using BPs */
                                gp.BP++;

                                /* Reading the Branch Point profile */
                                gp.BranchPointProfile= readProfile(rootFile, GeneIDconstants.sBP, 0);

                            } else {
                                String mess= "Wrong format: profile name "+ header+" \n" +
                                        "\tis not admitted for acceptors [only " +
                                        GeneIDconstants.sprofileACC+ ", "+ GeneIDconstants.sprofilePPT+ ", "+ GeneIDconstants.sprofileBP+ ", "+
                                        GeneIDconstants.sprofileU12BP+ ", "+ GeneIDconstants.sprofileU12gtagACC+ " or "+ GeneIDconstants.sprofileU12atacACC+ "]";
                                throw new RuntimeException(mess);
                            }
                        }
                    }
                }
            }

            /* Next profile for acceptor site */
            line= readHeader(rootFile);
            header= line.trim();
            if (header.length()== 0) {
                String mess= "Wrong format: header in optional profile for acceptor sites";
                throw new RuntimeException(mess);
            }
        }

        /* Reading the Acceptor site profile */
        gp.AcceptorProfile= readProfile(rootFile, GeneIDconstants.sACC, 0);

        line= readHeader(rootFile);
        header= line.trim();
        if (header.length()== 0) {
            String mess= "Wrong format: header in optional profile for splice sites";
            throw new RuntimeException(mess);
        }

        /* Read optional profiles: sprofileU12gtagDON,sprofileU12atacDON */

        while(!header.equalsIgnoreCase(GeneIDconstants.sprofileDON)) {

            if (header.equalsIgnoreCase(GeneIDconstants.sprofileU12gtagDON)) {

                u12gtagDon++;

                /* Reading the U12gtag donor profile */
                gp.U12gtagDonorProfile= readProfile(rootFile, GeneIDconstants.sDON, 0);

            } else {

                if (header.equalsIgnoreCase(GeneIDconstants.sprofileU12atacDON)) {

                    u12atacDon++;

                    /* Reading the U12atac donor profile */
                    gp.U12atacDonorProfile= readProfile(rootFile, GeneIDconstants.sDON, 0);

                } else {

                    if (header.equalsIgnoreCase(GeneIDconstants.sprofileU2gcagDON)) {

                        gp.U2GCAG++;

                        /* Reading the U2gcag donor profile */
                        gp.U2gcagDonorProfile= readProfile(rootFile, GeneIDconstants.sDON, 0);

                    } else {

                        if (header.equalsIgnoreCase(GeneIDconstants.sprofileU2gtaDON)) {

                            gp.U2GTA++;

                            /* Reading the U2gta donor profile */
                            gp.U2gtaDonorProfile= readProfile(rootFile, GeneIDconstants.sDON, 0);

                        } else {

                            if (header.equalsIgnoreCase(GeneIDconstants.sprofileU2gtgDON)) {

                                gp.U2GTG++;

                                /* Reading the U2gtg donor profile */
                                gp.U2gtgDonorProfile= readProfile(rootFile, GeneIDconstants.sDON, 0);

                            } else {

                                if (header.equalsIgnoreCase(GeneIDconstants.sprofileU2gtyDON)) {

                                    gp.U2GTY++;

                                    /* Reading the U2gty donor profile */
                                    gp.U2gtyDonorProfile= readProfile(rootFile, GeneIDconstants.sDON, 0);

                                } else {

                                    String mess= "Wrong format: profile name "+ header+ " \n" +
                                            "\tis not admitted for donors " +
                                            "[only "+ GeneIDconstants.sprofileDON+ ", "+ GeneIDconstants.sprofileU12gtagDON+ ", "+ GeneIDconstants.sprofileU12atacDON+ ", "+
                                            GeneIDconstants.sprofileU2gcagDON+ ", "+ GeneIDconstants.sprofileU2gtaDON+ ", "+ GeneIDconstants.sprofileU2gtgDON+ " or "+
                                            GeneIDconstants.sprofileU2gtyDON+ "]";
                                    throw new RuntimeException(mess);
                                }
                            }
                        }
                    }
                }
            }
            /* Next profile for donor site */
            line= readHeader(rootFile);
            header= line.trim();
            if (header.length()== 0) {
                String mess= "Wrong format: header in optional profile for donor sites";
                throw new RuntimeException(mess);
            }
        }

        /* Reading the Donor site profile */
        gp.DonorProfile= readProfile(rootFile, GeneIDconstants.sDON, 0);

        /* Switch on the site prediction of U12gtag and U12atac introns */
        if (u12bp!= 0 && u12gtagAcc!= 0 && u12gtagDon!= 0){ gp.U12GTAG++;gp.SPLICECLASSES++;}
        if (u12bp!= 0 && u12atacAcc!= 0 && u12atacDon!= 0){ gp.U12ATAC++;gp.SPLICECLASSES++;}

    }

    /**
     * Read information about signal and exon prediction in one isochore
     * - isochores are specific DNA regions according to the G+C content -
     */
    static void readIsochore(BufferedReader rootFile, GParam gp) throws Exception {

        /* 1. read boundaries of isochores */
        String line= readHeader(rootFile);
        line= readLine(rootFile);
        String[] ss= line.split("\\s");
        try {
            gp.leftValue= Integer.parseInt(ss[0]);
            gp.rightValue= Integer.parseInt(ss[1]);
            if (ss.length> 2)
                throw new ArrayIndexOutOfBoundsException("Array too long "+ ss.length+ ", expected 2!");
        } catch (Exception e) {
            throw new RuntimeException("Wrong format: isochore boundaries (G+C percent)", e);
        }
        Log.message("Isochores boundaries(min/max percentage): "+ gp.leftValue+ ","+ gp.rightValue);

        /* 2. read cutoff (final score) to accept one predicted exon */
        line= readHeader(rootFile);
        line= readLine(rootFile);
        ss= line.split("\\s");
        try {
            gp.Initial.ExonCutoff= Float.parseFloat(ss[0]);
            gp.Internal.ExonCutoff= Float.parseFloat(ss[1]);
            gp.Terminal.ExonCutoff= Float.parseFloat(ss[2]);
            gp.Single.ExonCutoff= Float.parseFloat(ss[3]);
            if (ss.length> 4)
                gp.utr.ExonCutoff= Float.parseFloat(ss[4]);
        } catch (Exception e) {
            throw new RuntimeException("Wrong format: exon score cutoffs (number/type)", e);
        }
        Log.message("Exon cutoffs: \t" +
                GeneIDconstants.DF_93F.format(gp.Initial.ExonCutoff)+ "\t"+
                GeneIDconstants.DF_93F.format(gp.Internal.ExonCutoff)+ "\t"+
                GeneIDconstants.DF_93F.format(gp.Terminal.ExonCutoff)+ "\t"+
                GeneIDconstants.DF_93F.format(gp.Single.ExonCutoff)+ "\t"+
                GeneIDconstants.DF_93F.format(gp.utr.ExonCutoff));

        /* 3. read cutoff (potential coding score) to accept one predicted exon */
        line= readHeader(rootFile);
        line= readLine(rootFile);
        ss= line.split("\\s");
        try {
            gp.Initial.OligoCutoff= Float.parseFloat(ss[0]);
            gp.Internal.OligoCutoff= Float.parseFloat(ss[1]);
            gp.Terminal.OligoCutoff= Float.parseFloat(ss[2]);
            gp.Single.OligoCutoff= Float.parseFloat(ss[3]);
            if (ss.length> 4)
                throw new ArrayIndexOutOfBoundsException("Array too long "+ ss.length+ ", expected 4!");
        } catch (Exception e) {
            throw new RuntimeException("Wrong format: potential coding score cutoffs (number/type)", e);
        }
        Log.message("Oligo cutoffs: \t" +
                GeneIDconstants.DF_93F.format(gp.Initial.OligoCutoff)+ "\t"+
                GeneIDconstants.DF_93F.format(gp.Internal.OligoCutoff)+ "\t"+
                GeneIDconstants.DF_93F.format(gp.Terminal.OligoCutoff)+ "\t"+
                GeneIDconstants.DF_93F.format(gp.Single.OligoCutoff));

        /* 4. Weight of signals in final exon score */
        line= readHeader(rootFile);
        line= readLine(rootFile);
        try {
            gp.Initial.siteFactor= Float.parseFloat(ss[0]);
            gp.Internal.siteFactor= Float.parseFloat(ss[1]);
            gp.Terminal.siteFactor= Float.parseFloat(ss[2]);
            gp.Single.siteFactor= Float.parseFloat(ss[3]);
            if (ss.length> 4)
                gp.utr.siteFactor= Float.parseFloat(ss[4]);
        } catch (Exception e) {
            throw new RuntimeException("Wrong format: weight of signal scores (number/type)", e);
        }
        Log.message("Site factors: \t" +
                GeneIDconstants.DF_93F.format(gp.Initial.siteFactor)+ "\t"+
                GeneIDconstants.DF_93F.format(gp.Internal.siteFactor)+ "\t"+
                GeneIDconstants.DF_93F.format(gp.Terminal.siteFactor)+ "\t"+
                GeneIDconstants.DF_93F.format(gp.Single.siteFactor)+ "\t"+
                GeneIDconstants.DF_93F.format(gp.utr.siteFactor));

        /* 5. Weight of coding potential in final exon score */
        line= readHeader(rootFile);
        line= readLine(rootFile);
        try {
            gp.Initial.exonFactor= Float.parseFloat(ss[0]);
            gp.Internal.exonFactor= Float.parseFloat(ss[1]);
            gp.Terminal.exonFactor= Float.parseFloat(ss[2]);
            gp.Single.exonFactor= Float.parseFloat(ss[3]);
            if (ss.length> 4)
                throw new ArrayIndexOutOfBoundsException("Array too long "+ ss.length+ ", expected 4!");
        } catch (Exception e) {
            throw new RuntimeException("Wrong format: weight of coding potential scores (number/type)", e);
        }
        Log.message("Exon factors: \t" +
                GeneIDconstants.DF_93F.format(gp.Initial.exonFactor)+ "\t"+
                GeneIDconstants.DF_93F.format(gp.Internal.exonFactor)+ "\t"+
                GeneIDconstants.DF_93F.format(gp.Terminal.exonFactor)+ "\t"+
                GeneIDconstants.DF_93F.format(gp.Single.exonFactor));

        /* 6. Weight of homology information in final exon score */
        line= readHeader(rootFile);
        line= readLine(rootFile);
        try {
            gp.Initial.HSPFactor= Float.parseFloat(ss[0]);
            gp.Internal.HSPFactor= Float.parseFloat(ss[1]);
            gp.Terminal.HSPFactor= Float.parseFloat(ss[2]);
            gp.Single.HSPFactor= Float.parseFloat(ss[3]);
            if (ss.length> 4)
                gp.utr.HSPFactor= Float.parseFloat(ss[4]);
        } catch (Exception e) {
            throw new RuntimeException("Wrong format: weight of homology scores (number/type)", e);
        }
        Log.message("HSP factors: \t" +
                GeneIDconstants.DF_93F.format(gp.Initial.HSPFactor)+ "\t"+
                GeneIDconstants.DF_93F.format(gp.Internal.HSPFactor)+ "\t"+
                GeneIDconstants.DF_93F.format(gp.Terminal.HSPFactor)+ "\t"+
                GeneIDconstants.DF_93F.format(gp.Single.HSPFactor)+ "\t"+
                GeneIDconstants.DF_93F.format(gp.utr.HSPFactor));

        /* 7. read weights to correct the score of exons after general cutoff */
        String header= readHeaderCheck(rootFile,
                "Wrong format: header for exon weights and optional U12 score threshold");

        while((!header.equalsIgnoreCase(GeneIDconstants.sExon_weights))&& (!header.equals("Exon_weigths"))) {
            /* 1. Read RSSMARKOVSCORE for markov score to assign non-exonic recursively spliced elements */
            if(header.equalsIgnoreCase(GeneIDconstants.sRSSMARKOVSCORE)) {
                gp.RSSMARKOVSCORE= readLineParseFloat(rootFile, GeneIDconstants.sRSSMARKOVSCORE);
            /* 1. Read Evidence Exon Weight */
            } else if(header.equalsIgnoreCase(GeneIDconstants.sEVIDENCEW)) {
                gp.EvidenceEW= readLineParseFloat(rootFile, GeneIDconstants.sEVIDENCEW);
            /* 1. Read Evidence Exon Factor */
            } else if(header.equalsIgnoreCase(GeneIDconstants.sEVIDENCEF)) {
                gp.EvidenceFactor= readLineParseFloat(rootFile, GeneIDconstants.sEVIDENCEF);
            /* 1. Read RSS_Donor_Score_Cutoff */
            } else if(header.equalsIgnoreCase(GeneIDconstants.sRSS_DONOR_SCORE_CUTOFF)) {
                gp.RSSDON= readLineParseFloat(rootFile, GeneIDconstants.sRSS_DONOR_SCORE_CUTOFF);
            /* 1. Read RSSMARKOVSCORE for markov score to assign non-exonic recursively spliced elements */
            } else if(header.equalsIgnoreCase(GeneIDconstants.sRSS_ACCEPTOR_SCORE_CUTOFF)) {
                gp.RSSACC= readLineParseFloat(rootFile, GeneIDconstants.sRSS_ACCEPTOR_SCORE_CUTOFF);
            /* 1. Read u12SpliceScoreThresh for sum of U12 donor and acceptor splice scores */
            } else if(header.equalsIgnoreCase(GeneIDconstants.sU12_SPLICE_SCORE_THRESH)) {
                gp.u12SpliceScoreThresh = readLineParseFloat(rootFile, GeneIDconstants.sU12_SPLICE_SCORE_THRESH);
            /* 1. Read U12_EXON_SCORE_THRESH for sum of U12 donor and acceptor exon scores */
            } else if(header.equalsIgnoreCase(GeneIDconstants.sU12_EXON_SCORE_THRESH)) {
                gp.U12_EXON_SCORE_THRESH= readLineParseFloat(rootFile, GeneIDconstants.sU12_EXON_SCORE_THRESH);
            /* 1. Read U12_EXON_WEIGHT, an additional exon weight that applies to exons flanking U12 introns */
            } else if(header.equalsIgnoreCase(GeneIDconstants.sU12_EXON_WEIGHT)) {
                gp.U12EW= readLineParseFloat(rootFile, GeneIDconstants.sU12_EXON_WEIGHT);
            }
            header= readHeaderCheck(rootFile, "Wrong format: header for exon weights");
        }
        line= readLine(rootFile);
        ss= line.split("\\s");
        try {
            gp.Initial.ExonWeight= Float.parseFloat(ss[0]);
            gp.Internal.ExonWeight= Float.parseFloat(ss[1]);
            gp.Terminal.ExonWeight= Float.parseFloat(ss[2]);
            gp.Single.ExonWeight= Float.parseFloat(ss[3]);
            if (ss.length> 4)
                gp.utr.ExonWeight= Float.parseFloat(ss[4]);
        } catch (Exception e) {
            throw new RuntimeException("Wrong format: exon weight values (number/type)", e);
        }
        Log.info("Exon weights: \t" +
                GeneIDconstants.DF_93F.format(gp.Initial.ExonWeight)+ "\t"+
                GeneIDconstants.DF_93F.format(gp.Internal.ExonWeight)+ "\t"+
                GeneIDconstants.DF_93F.format(gp.Terminal.ExonWeight)+ "\t"+
                GeneIDconstants.DF_93F.format(gp.Single.ExonWeight));

        /* 8. Read splice site profiles */
        /* (a).start codon profile */
        gp.StartProfile= readProfile(rootFile, GeneIDconstants.sSTA,1);

        /* (b).acceptor and donor site profiles */
        readProfileSpliceSites(rootFile, gp);

        /* (c).donor site profile */
        /* ReadProfile(RootFile, gp.DonorProfile , sDON,1); */

        /* (d).stop codon profile */
        gp.StopProfile= readProfile(rootFile, GeneIDconstants.sSTO,1);

        /* 9. read coding potential log-likelihood values (Markov chains) */
        header= readHeaderCheck(rootFile, "Wrong format: header "+ line);

        while ((!header.equalsIgnoreCase(GeneIDconstants.sMarkov))&& (!header.equalsIgnoreCase("Markov_oligo_logs_file"))) {

            /* printMess(header); */

            if (header.equalsIgnoreCase(GeneIDconstants.sprofilePolyA)) {
                gp.PAS++;
                Log.info("Reading PolyA Signal Profile");
                /* Reading the U2gty donor profile */
                gp.PolyASignalProfile= readProfile(rootFile, GeneIDconstants.sPOL, 0);
            }

            header= readHeaderCheck(rootFile, "Wrong format: header "+ line);
        }
        /* Next profile for Markov order */
        line= readLine(rootFile);
        try {
            gp.OligoLength= Integer.parseInt(line.trim());
        } catch (Exception e) {
            throw new RuntimeException("Wrong format: oligonucleotide length", e);
        }
        Log.message("Oligonucleotide (word) length: "+ gp.OligoLength);

        /* (a). Initial probability matrix */
        Log.message("Reading Markov Initial likelihood matrix");

        /* Computing the right number of initial values to read */
        gp.OligoDim= (int) Math.pow(4, gp.OligoLength);
        Log.message("Used oligo array size: "+ (gp.OligoDim * 3));

        line= readHeader(rootFile);
        for(int j = 0; j < gp.OligoDim * 3; j++) {
            line= readLine(rootFile);
            ss= line.split("\\s");
            try {
                // sscanf(line, "%*s %d %d %f", &i, &f, &lscore))!=3
                int i= Integer.parseInt(ss[1]);
                int f= Integer.parseInt(ss[2]);
                float lscore= Float.parseFloat(ss[3]);
                gp.OligoLogsIni[f][i]=lscore;
            } catch (Exception e) {
                throw new RuntimeException("Wrong format/nunber (%s): Initial Markov value" +
                        "\n"+ line, e);
            }
        }

        /* (b). Transition probability matrix */
        Log.message("Reading Markov Transition likelihood matrix");

        double OligoLength_1= gp.OligoLength + 1;
        gp.OligoDim_1= (int) Math.pow(4, OligoLength_1);

        Log.message("Used oligo array size: " + (gp.OligoDim_1 * 3));

        line= readHeader(rootFile);
        for(int j = 0; j < gp.OligoDim_1 * 3; j++) {
            line= readLine(rootFile);
            ss= line.split("\\s");
            try {
                // sscanf(line, "%*s %d %d %f", &i, &f, &lscore))!=3
                int i= Integer.parseInt(ss[1]);
                int f= Integer.parseInt(ss[2]);
                float lscore= Float.parseFloat(ss[3]);
                gp.OligoLogsTran[f][i]=lscore;
            } catch (Exception e) {
                throw new RuntimeException("Wrong format/number (%s): Transition Markov value" +
                        "\n"+ line, e);
            }
        }

        /* 10. read maximum number of donors per acceptor site (BuildExons) */
        line= readHeader(rootFile);
        line= readLine(rootFile);
        try {
            gp.MaxDonors= Integer.parseInt(line.trim());
        } catch (Exception e) {
            throw new RuntimeException("Bad format: MaxDonors value", e);
        }
        Log.message("Maximum donors by acceptor = "+ gp.MaxDonors);

    }

    /**
     * Reads a header line and checks for a non-empty header.
     * @param rootFile file with the site models
     * @param mess error message if no header was found
     * @return header identifier
     * @throws Exception for I/O errors and empty headers
     */
    private static String readHeaderCheck(BufferedReader rootFile, String mess) throws Exception {
        String line= readHeader(rootFile);
        String header= line.trim();
        if (header.length()== 0)
            throw new RuntimeException(mess);
        return header;
    }

    /**
     * Reads a line and to be parsed into a single float value.
     * @param rootFile file with the site models
     * @param parID name of the parameter that is to be parsed from the line
     * @return float value that has been parsed from the next line
     * @throws Exception I/O or parsing errors
     */
    private static float readLineParseFloat(BufferedReader rootFile, String parID) throws Exception {
        String line= readLine(rootFile);
        float value= Float.NaN;
        try {
            value= Float.parseFloat(line);
            Log.message(parID + ":\t" +
                    GeneIDconstants.DF_92F.format(value));
            return value;
        } catch (Exception e) {
            throw new RuntimeException("Wrong format: "+ parID+ " value scores (number/type)");
        }

    }


    /*
     * Read the input of statistics data model
     * @param name name of an explicitly provided geneID parameter file
     * @param settings the geneID settings
     * @return the geneID parameters
     */
    public static GParam[] readParam (String name, GeneIDsettings settings) throws Exception {

        /* 0. Select parameters filename for reading it */
        /* Filename must be: option P, env.var GENEID or default (none) */
        String Geneid= System.getProperty("GENEID");

        /* rootfile will be the parameter file handle descriptor */
        BufferedReader rootFile= null;

        /* (a) Using -P option */
        //the original GeneID code uses this branch
        //if (!GeneIDconstants.PARAMETERFILE.equals(name)) {
        if (name!= null) {
            Log.message("Loading parameter file "+ name);
            try {
                rootFile= new BufferedReader(new FileReader(name));
            } catch (Exception e) {
                Log.error("Parameter file (-P) can not be open to read: "+ e.getMessage());
                e.printStackTrace();
            }

        /* (b) Using GENEID environment var */
        } else if (Geneid!= null) {
            Log.message("Loading parameter file from GENEID env. var "+ Geneid);
            try {
                rootFile= new BufferedReader(new FileReader(Geneid));
            } catch (Exception e) {
                Log.error("Parameter file (GENEID env.var) can not be open to read");
            }

        /* (c) Using default parameter file */
        } else {
            Log.message("Loading default parameter file (human) "+ GeneIDconstants.PARAMETERFILE);
            try {
                // cannot be loaded directly by a file reader (i.e., dynamic resource)
                //rootFile= new BufferedReader(
                //        new FileReader(new File(Profile.class.getResource(GeneIDconstants.PARAMETERFILE).getFile())
                //));
                rootFile= new BufferedReader(new InputStreamReader(
                        Profile.class.getResourceAsStream(GeneIDconstants.PARAMETERFILE)));
            } catch (Exception e) {
                Log.error("Parameter file (default) can not be open to read");
            }
        }
        if (rootFile== null)
            throw new RuntimeException("No parameter file!");


        /* 1. Read NO_SCORE penalty for nucleotides not supported by homology */
        String line= readHeader(rootFile);
        GeneIDsettings.NO_SCORE= (int) readLineParseFloat(rootFile, "NO_SCORE");

        /* 2. Read the number of isochores */
        String header= readHeaderCheck(rootFile, "Wrong format: header for number of isochores");
        while(!header.equalsIgnoreCase(GeneIDconstants.sNUMISO)) {

            /* BKGD_SUBTRACT_FLANK_LENGTH */
            if(header.equalsIgnoreCase(GeneIDconstants.sBKGD_SUBTRACT_FLANK_LENGTH)) {
                GeneIDsettings.BKGD_SUBTRACT_FLANK_LENGTH= (int) readLineParseFloat(rootFile, GeneIDconstants.sBKGD_SUBTRACT_FLANK_LENGTH);
            }

            line= readHeaderCheck(rootFile, "Wrong format: header for number of isochores");
        }

        int nIsochores= (int) readLineParseFloat(rootFile, "Number of isochores");
        if (nIsochores > GeneIDconstants.MAXISOCHORES || nIsochores<= 0)
            Log.error("Wrong value: number of isochores(MAXISOCHORES)");

        /* 3. Reading every one of the isochores */
        GParam[] isochores= new GParam[nIsochores];
        for(int i= 0; i< nIsochores; i++) {
            Log.message("Reading isochore "+ (i+1));
            isochores[i]= new GParam();
            readIsochore(rootFile, isochores[i]);
        }

        /* 4. Reading the GeneModel */
        line= readHeader(rootFile);

        /* Ready to update dictionary of exon types */
        isochores[0].D= new Dictionary();
        Log.message("Dictionary ready to acquire information");

        if (GeneIDsettings.SGE!= 0) {
            Log.message("Using an internal Gene Model");
            isochores[0].nclass = isochores[0].D.forceGeneModel(
                    isochores[0].nc,
                    isochores[0].ne,
                    isochores[0].UC,
                    isochores[0].DE,
                    isochores[0].md,
                    isochores[0].Md,
                    isochores[0].block);

            Log.message(isochores[0].nclass+ " Gene Model rules have been read and saved");

        } else {

            Log.message("Reading Gene Model rules");
            isochores[0].nclass = isochores[0].D.readGeneModel(
                    rootFile,
                    isochores[0].nc,
                    isochores[0].ne,
                    isochores[0].UC,
                    isochores[0].DE,
                    isochores[0].md,
                    isochores[0].Md,
                    isochores[0].block);

            Log.message(isochores[0].nclass+ " Gene Model rules have been read and saved");
        }

        /* Replication of gene model information for each isochore */
        Dictionary.shareGeneModel(isochores, nIsochores);

        return isochores;
    }

}

