package barna.astafunk.HMM;

/**
 * @author Vitor Lima Coelho
 * @version 1.2
 *
 */


import barna.astafunk.utils.FunkSettings;

/**
 * This class describes a typical profile HMM used by HMMER3.x.
 */
public class ProfileHMM {

    /* Constants */

    /**
     * Index of m &rarr; m transition values
     * @see #stateTransitionMatrix
     */
    public static final int MM = 0; // m -> m

    /**
     * Index of m &rarr; i transition values
     * @see #stateTransitionMatrix
     */
    public static final int MI = 1; // m -> i

    /**
     * Index of m &rarr; d transition values
     * @see #stateTransitionMatrix
     */
    public static final int MD = 2; // m -> d

    /**
     * Index of i &rarr; m transition values
     * @see #stateTransitionMatrix
     */
    public static final int IM = 3; // i -> m

    /**
     * Index of i &rarr; i transition values
     * @see #stateTransitionMatrix
     */
    public static final int II = 4; // i -> i

    /**
     * Index of d &rarr; m transition values
     * @see #stateTransitionMatrix
     */
    public static final int DM = 5; // d -> m

    /**
     * Index of d &rarr; d transition values
     * @see #stateTransitionMatrix
     */
    public static final int DD = 6; // d -> d

    /**
     * Index of N (N-term) state values
     * @see #xsc
     */
    public static final int N = 0;

    /**
     * Index of E (End) state values
     * @see #xsc
     */
    public static final int E = 1;

    /**
     * Index of C (C-term) state values
     * @see #xsc
     */
    public static final int C = 2;

    /**
     * Index of J (Join) state values
     * @see #xsc
     */
    public static final int J = 3;

    /**
     * Index of MOVE condition
     * @see #xsc
     */
    public static final int MOVE = 0;

    /**
     * Index of LOOP condition
     * @see #xsc
     */
    public static final int LOOP = 1;


    /**
     * The natural logarithm of 2.
     */

    public static final double LOG2= Math.log(2);


    /* End constants*/

   /* Begin profile HMM meta-info

    /**
     * The unique identifier for the save file format version
     */
    private String version;

    /**
     * The model name.
     */
    private String name;

    /**
     * The accession number.
     */
    private String acc;

    /**
     * The description line.
     */
    private String description;

    /**
     * The number of match states in the model.
     */
    private int length;

    /**
     * The symbol alphabet type.
     */
    private String alphabet;

    /**
     * The reference annotation flag.
     */
    private String rf;

    /**
     * The consensus structure annotation flag.
     */
    private String cs;

    /**
     * The map annotation flag.
     */
    private String map;

    /**
     * Date building of HMM
     */
    private String date;

    /**
     * The number of sequences that the HMM was trained on.
     */
    private int nseq;

    /**
     * The effective total number of sequences determined
     * by hmmbuild during sequence weighting, for combining
     * observed counts with Dirichlet prior information
     * in parameterizing the model.
     */
    private String effn;

    /**
     * Training alignment checksum.
     */
    private String cksum;

    /**
     * Pfam gathering thresholds GA1.
     */
    private double ga1;

    /**
     * Pfam gathering thresholds GA2 or Theta.
     */
    private double ga2;

    /**
     * Pfam trusted cutoffs TC1 and TC2.
     */
    private String tc;

    /**
     * Pfam noise cutoffs NC1 and NC2.
     */
    private String nc;

    /**
     * Statistical parameters needed for E-value calculations. See STATS parameter in HMMER documentation.
     */
    private String stats_local_msv;

    /**
     * Statistical parameters needed for E-value calculations. See STATS parameter in HMMER documentation.
     */
    private String stats_local_viterbi;

    /**
     * Statistical parameters needed for E-value calculations. See STATS parameter in HMMER documentation.
     */
    private String stats_local_forward;


    /* End profile HMM meta-info */

     /* Begin delta factor variables */

    /**
     * Maximal gap open cost of inserting.
     */

    public double ALPHA;

    /**
     * Maximal gap extend cost of inserting
     */
    public double BETA;

    /**
     * Max inserting probability
     */
    public double MAX_INSERT_PROB;


    /**
     * Index of AA with max inserting probability
     */

    public int INDEX_MAX_PROB;

    /**
     * Delta factor
     */

    private int DELTA_FACTOR;

    /* End delta factor variables */

    /**
     * The matrix with match emission probability distributions.
     */
    private double [][] matchEmissionMatrix;

    /**
     * The matrix with insert emission probability distributions.
     */
    private double [][] insertEmissionMatrix;

    /**
     * The matrix with state transition probabilities.
     */
    private double [][] stateTransitionMatrix;

    /**
     * these are the model’s overall average match state emission probabilities,
     * which are used as a background residue composition in the “filter null” model
     */

    private double [] compo;

    /**
     * Scores for extended states {E,J,N,C},
     * either for LOOP or MOVE condition.
     */
    public double[][] xsc= null;

    /**
     * expected # of J's: 0 or 1, uni vs. multihit
     * see hmmer.h
     */
    public int nj;


    /**
     * Array with the needed score to hit a domain from a state of the model.
     * For example, if we are in state 2, position 2 of this vector says that
     * need optimumScoreArray[2] - GA2 (bit score) to hit a domain;
     */
    private double[] optimumScoreArray;


    /**
     * The ProfileHMM constructor. It initializes the distributions and probabilities matrices.
     * @param length An integer number of match states in the model.
     * @param multiHit indicating whether a model can hit the sequence multiple times (= loop)
     */
    protected ProfileHMM(int length, boolean multiHit) {

        this.length = length;

        // E state loop/move probabilities:
        // nonzero for MOVE allows loops/multihits
        // N,C,J transitions are set later by length config
        // see modelconfig.c
        xsc= new double[4][]; // N,E,C,J state
        for (int i = 0; i < xsc.length; i++)
            xsc[i]= new double[2];  // MOVE and LOOP
        if (multiHit) {
            xsc[E][MOVE] = -LOG2;
            xsc[E][LOOP] = -LOG2;
            nj = 1; // number of J's
        } else {
              xsc[E][MOVE] = 0.0f;
            xsc[E][LOOP] = -FunkSettings.INF;
            nj = 0; // number of J's
        }
    }

    /**
     * Creates an instance of the given length,
     * without allowing for multiple hits
     * @param length An integer number of match states in the model.
     */
    public ProfileHMM(int length) {
        this(length, false);
    }

    /**
     * @return Name of profile HMM.
     */
    public String getName() {
        return name;
    }

    /**
     * @param name Name of profile HMM to set.
     */
    public void setName(String name) {
        this.name = name;
    }

    /**
     * @return Accession ID
     */
    public String getAcc() {
        return acc;
    }

    /**
     * @param acc Accession ID to set.
     */
    public void setAcc(String acc) {
        this.acc = acc;
    }

    /**
     * @return Length of profile HMM
     */
    public int getLength() {
        return length;
    }

    /**
     * @param length Model length (number of states) to set.
     */
    public void setLength(int length) {
        this.length = length;
    }

    /**
     * @return Description of profile HMM.
     */
    public String getDescription() {
        return description;
    }

    /**
     * @param description Description of profile HMM to set.
     */
    public void setDescription(String description) {
        this.description = description;
    }

    /**
     * @return Matching emission matrix.
     */
    public double[][] getMatchEmissionMatrix() {
        return matchEmissionMatrix;
    }

    /**
     * @param matchEmissionMatrix Match emission matrix to set.
     */
    public void setMatchEmissionMatrix(double[][] matchEmissionMatrix) {
        this.matchEmissionMatrix = matchEmissionMatrix;
    }

    /**
     * Method to return inserting emission matrix.
     * @return Inserting emission matrix.
     */
    public double[][] getInsertEmissionMatrix() {
        return insertEmissionMatrix;
    }

    /**
     * Method to set inserting emission matrix.
     * @param insertEmissionMatrix Inserting emission matrix.
     */
    public void setInsertEmissionMatrix(double[][] insertEmissionMatrix) {
        this.insertEmissionMatrix = insertEmissionMatrix;
    }

    /**
     * Method to return state transition matrix.
     * @return State transition matrix
     */
    public double[][] getStateTransitionMatrix() {
        return stateTransitionMatrix;
    }

    /**
     * Method to set state transition matrix.
     * @param stateTransitionMatrix State transition matrix
     */
    public void setStateTransitionMatrix(double[][] stateTransitionMatrix) {
        this.stateTransitionMatrix = stateTransitionMatrix;
    }

    /**
     * Method to output attributes of profile HMM.
     */
    public void show(){
        System.out.println("NAME " + getName());
        System.out.println("ACC" + getAcc());
        System.out.println("DESC " + getDescription());
        System.out.println("LENG " + getLength());

        System.out.println("Match Emission Matrix ");

        for(int i = 0; i < length + 1; i++){

            System.out.println(i);

            for (int j = 0; i < 20; j++){

                System.out.println(matchEmissionMatrix[i][j] + " ");

            }

            System.out.println();

        }

        System.out.println("Insert Emission Matrix ");

        for(int i = 0; i < length + 1; i++){

            System.out.println(i);

            for (int j = 0; i < 20; j++){

                System.out.println(insertEmissionMatrix[i][j]);

            }

            System.out.println();

        }

        System.out.println("State Transition Matrix ");

        for(int i = 0; i < length + 1; i++){

            System.out.println(i);

            for (int j = 0; i < 7; j++){

                System.out.println(stateTransitionMatrix[i][j]);

            }

            System.out.println();
        }
    }

    /**
     * Method to print matching emission matrix.
     */
    public void printMatchEmissionMatrix(){

        for(int i = 0; i < length+1; i++){
            for(int j = 0; j < 20; j++){

                System.out.print(matchEmissionMatrix[i][j] + " ");
            }
            System.out.println();
        }

    }

    /**
     * Method to print inserting emission matrix.
     */
    public void printInsertEmissionMatrix(){

        for(int i = 0; i < length+1; i++){
            for(int j = 0; j < 20; j++){

                System.out.print(insertEmissionMatrix[i][j] + " ");
            }
            System.out.println();
        }

    }

    /**
     * Method to print state transition matrix.
     */
    public void printStateTransitionMatrix(){

        for(int i = 0; i < length+1; i++){
            for(int j = 0; j < 7; j++){

                System.out.print(stateTransitionMatrix[i][j] + " ");
            }
            System.out.println();
        }

    }

    /**
     * @return Gathering threshold 1 (sequence).
     */
    public double getGa1() {
        return ga1;
    }

    /**
     * @param ga1 Gathering threshold 1 (sequence) to set.
     */
    public void setGa1(double ga1) {
        this.ga1 = ga1;
    }

    /**
     * @return  Gathering threshold 2 (domain).
     */
    public double getGa2() {
        return ga2;
    }

    /**
     * @param ga2 Gathering threshold 2 (domain) to set.
     */
    public void setGa2(double ga2) {
        this.ga2 = ga2;
    }

    /**
     * @return String with TC (trusted cutoffs) parameters TC1 e TC2
     */
    public String getTc() {
        return tc;
    }

    /**
     * @param tc String with TC (trusted cutoffs TC1 e TC2) parameters to set.
     */
    public void setTc(String tc) {
        this.tc = tc;
    }

    /**
     * @return String with NC (noise cutoffs NC1 and NC2) parameters
     */
    public String getNc() {
        return nc;
    }

    /**
     * @param nc String with NC (noise cutoffs NC1 and NC2) parameters to set.
     */
    public void setNc(String nc) {
        this.nc = nc;
    }

    /**
     * @return COMPO array.
     */
    public double[] getCompo() {
        return compo;
    }

    /**
     * @param compo COMPO array to set.
     */
    public void setCompo(double[] compo) {
        this.compo = compo;
    }

    /**
     * @return Alphabet array.
     */

    public char[] getAlphabet(){

        char [] d;

        if(alphabet.toLowerCase().equals("dna")){

            d = new char[4];

            d[0] = 'A';
            d[1] = 'C';
            d[2] = 'G';
            d[3] = 'T';

        }else if(alphabet.toLowerCase().equals("rna")){

            d = new char[4];

            d[0] = 'A';
            d[1] = 'C';
            d[2] = 'G';
            d[3] = 'U';

        }else {

            d = new char[20];

            d[0] = 'A';
            d[1] = 'C';
            d[2] = 'D';
            d[3] = 'E';
            d[4] = 'F';
            d[5] = 'G';
            d[6] = 'H';
            d[7] = 'I';
            d[8] = 'K';
            d[9] = 'L';
            d[10] = 'M';
            d[11] = 'N';
            d[12] = 'P';
            d[13] = 'Q';
            d[14] = 'R';
            d[15] = 'S';
            d[16] = 'T';
            d[17] = 'V';
            d[18] = 'W';
            d[19] = 'Y';

        }

        return d;

    }

    /**
     * Method to return null model. A null model is used to calculate HMMER log odds scores.
     * The null model states the expected background
     * occurrence frequencies of the 20 amino acids or the 4 nucleotide bases.
     * @return Null model array.
     */

    public double[] getNullModel(){

        double [] d;

        if(alphabet.toLowerCase().equals("dna") || alphabet.toLowerCase().equals("rna")){

            d = new double[4];

            double neglog = 1.38629;

            d[0] = neglog;
            d[1] = neglog;
            d[2] = neglog;
            d[3] = neglog;

        }else {

            d = new double[20];

            // negative natural log of probabilities: from source code of HMMER 3.1

            /*d[0] = 2.58336;    // A
            d[1] = 2.93692;    //# C
            d[2] = 2.93201;    //# D
            d[3] = 2.76139;    //# E
            d[4] = 3.20001;    //# F
            d[5] = 2.68168;    //# G
            d[6] = 3.79843;    //# H
            d[7] = 2.85973;    //# I
            d[8] = 2.82349;    //# K
            d[9] = 2.37087;    //# L
            d[10] = 3.74782;   //# M
            d[11] = 3.0946;   //# N
            d[12] = 3.0106;   //# P
            d[13] = 3.21312;   //# Q
            d[14] = 2.96476;   //# R
            d[15] = 2.62812;   //# S
            d[16] = 2.85677;   //# T
            d[17] =  2.7295;  //# V
            d[18] = 4.38099;   //# W
            d[19] = 3.44249;   //# Y*/

            // Probabilities: from source code of HMMER 3.1

            d[0] =0.075520; // A
            d[1] =0.016973; // C
            d[2] =0.053029; // D
            d[3] =0.063204; //E
            d[4] =0.040762; //# F
            d[5] =0.068448 ;//# G
            d[6] =0.022406 ;//# H
            d[7] =0.057284 ;//# I
            d[8] =0.059398 ;//# K
            d[9] =0.093399 ;//# L
            d[10] = 0.023569; //# M
            d[11] = 0.045293 ;//# N
            d[12] = 0.049262 ;//# P
            d[13] = 0.040231 ;//# Q
            d[14] = 0.051573 ;//# R
            d[15] = 0.072214 ;//# S
            d[16] =  0.057454; //# T
            d[17] =  0.065252; //# V
            d[18] = 0.012513 ;//# W
            d[19] = 0.031985 ;//# Y
        }
        return d;

    }

    /**
     *  Method to set alphabet.
     * @param alphabet Symbol alphabet type
     */
    public void setAlphabet(String alphabet) {
        this.alphabet = alphabet;
    }

    /**
     * Purpose:   Given a model already configured for scoring, in some
     *            particular algorithm mode; reset the expected length
     *            distribution of the profile for a new mean of <code>L</code>.
     *
     *            This doesn't affect the length distribution of the null
     *            model. That must also be reset, using p7_bg_SetLength().
     *
     *            We want this routine to run as fast as possible, because
     *            the caller needs to dynamically reconfigure the model
     *            for the length of each target sequence in a database
     *            search. The profile has precalculated gm&rarr;nj,
     *            the number of times the J state is expected to be used,
     *            based on the E state loop transition in the current
     *            configuration.
     * see modelconfig.c
     * @param L length of the unannotated (target) sequence
     */
    public void reconfigLength(int L) {
          // Configure N,J,C transitions so they bear L/(2+nj) of the total
          // unannotated sequence length L.
        double pmove = (2.0f + (double) nj) / ((double) L + 2.0f + nj); // 2/(L+2) for sw; 3/(L+3) for fs
        double ploop = 1.0f - pmove;
        xsc[N][LOOP] =  xsc[C][LOOP] = xsc[J][LOOP] = Math.log(ploop);
        xsc[N][MOVE] =  xsc[C][MOVE] = xsc[J][MOVE] = Math.log(pmove);
    }

    /**
     * @return Variable ALPHA to set.
     */
    public double getALPHA() {
        return ALPHA;
    }

    /**
     * @param ALPHA Variable ALPHA to set.
     */
    public void setALPHA(double ALPHA) {
        this.ALPHA = ALPHA;
    }

    /**
     * @return Variable BETA.
     */
    public double getBETA() {
        return BETA;
    }

    /**
     * @param BETA Variable BETA to set.
     */
    public void setBETA(double BETA) {
        this.BETA = BETA;
    }


    /**
     * @param delta Delta factor to set.
     */
    public void setDeltaFator(int delta) {

        this.DELTA_FACTOR = delta;
    }


    /**
     * @return Delta factor.
     */
    public int getDeltaFactor() {
        return DELTA_FACTOR;
    }

    /**
     * @return Optimum Score Array.
     */
    public double[] getOptimumScoreArray() {
        return optimumScoreArray;
    }

    /**
     * @return Optimum Suffix Score on interval [x,y] of profile consensus.
     */
    public double getOptimumScoreArray(int i, int j, int sequenceLength) {
        if(i == length)
            return 0;
        int temp = i + sequenceLength - j;
        if(temp > length)
            temp = length;
        return optimumScoreArray[temp] - optimumScoreArray[i];
    }

    /**
     * @param optimumScoreArray Optimum score array to set.
     */
    public void setOptimumScoreArray(double[] optimumScoreArray) {

        this.optimumScoreArray = optimumScoreArray;
    }
}
