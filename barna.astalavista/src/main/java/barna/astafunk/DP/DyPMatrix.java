package barna.astafunk.DP;

import barna.astafunk.Astafunk;
import barna.astafunk.HMM.ProfileHMM;
import barna.astafunk.Tsearch;
import barna.astafunk.utils.FunkSettings;
import barna.astafunk.utils.Utils;
import barna.model.ASEvent;
import barna.model.DirectedRegion;
import barna.model.Transcript;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * @author  Michael Sammeth and Vitor Coelho
 * @version 3.0
 * This class implements a dynamic programming matrix for profile HMM search.
 */

/**
 * This class describes a dynamic programming matrix
 * to search profile HMM against an amino acid sequence by Viterbi algorithm (see Durbin 1998).
 */
public class DyPMatrix {

    /**
     * Horizontal or inserting score matrix.
     */
    private double [][] dpmatrix_h;

    /**
     * Vertical or deleting score matrix.
     */
    private double [][] dpmatrix_v;

    /**
     * Diagonal or matching score matrix.
     */
    private double [][] dpmatrix_d;

    /**
     * The score matrix.
     */
    private DyPCell [][] dpmatrix_score;

    /**
     * N state, N-terminal unaligned sequence state.
     * Emits on transition with K emission probabilities
     * for K symbols in the alphabet.
     */
    private double[] seqstate_n;

    /**
     * Begin state (for entering main model). Non-emitter.
     */
    private double[] seqstate_b;

    /**
     *  End state (for exiting main model). Non-emitter.
     */
    private double[] seqstate_e;

    /**
     * C-terminal unaligned sequence state.
     * Emits on transition with K emission probabilities
     * for K symbols in the alphabet.
     */
    private double[] seqstate_c;

    /**
     * Joining segment unaligned sequence state.
     * Emits on transition with K emission probabilities
     * for K symbols in the alphabet.
     */
    private double[] seqstate_j;

    /**
     * Input sequence. Amino acid sequence of transcript.
     */
    private String sequence;

    /**
     * Sequence length.
     */
    private int sequenceLength;

    /**
     * Number of states of a profile HMM.
     */
    private int lengthHMM;

    /**
     * List of possible hits
     */
    private List<DyPCell> tempHits;

    /**
     * Profile HMM.
     */
    private ProfileHMM hmm;

    /**
     * Constructor of DyPMatrix object.
     * @param sequence Amino acid sequence of transcript.
     * @param hmm profile HMM.
     */
    public DyPMatrix(String sequence, ProfileHMM hmm){

        this.sequence = sequence;

        this.sequenceLength = sequence.length();

        this.hmm = hmm;

        this.lengthHMM = hmm.getLength();

        this.dpmatrix_h = new double[lengthHMM+1][sequenceLength+1];

        this.dpmatrix_d = new double[lengthHMM+1][sequenceLength+1];

        this.dpmatrix_v = new double[lengthHMM+1][sequenceLength+1];

        this.dpmatrix_score = new DyPCell[lengthHMM+1][sequenceLength+1];

        // states for entering, leaving and looping dp matrix
        seqstate_n= new double[sequenceLength+ 1];
        seqstate_b= new double[sequenceLength+ 1];
        seqstate_e= new double[sequenceLength+ 1];
        seqstate_c= new double[sequenceLength+ 1];
        seqstate_j= new double[sequenceLength+ 1];

    }

    /**
     * Redefine DP matrix for a new sequence
     */
    public void initializeNewSequence(String newSequence){
        this.sequence = newSequence;
        this.sequenceLength = newSequence.length();

        this.dpmatrix_h = new double[lengthHMM+1][sequenceLength+1];

        this.dpmatrix_d = new double[lengthHMM+1][sequenceLength+1];

        this.dpmatrix_v = new double[lengthHMM+1][sequenceLength+1];

        this.dpmatrix_score = new DyPCell[lengthHMM+1][sequenceLength+1];

        // states for entering, leaving and looping dp matrix
        seqstate_n= new double[sequenceLength+ 1];
        seqstate_b= new double[sequenceLength+ 1];
        seqstate_e= new double[sequenceLength+ 1];
        seqstate_c= new double[sequenceLength+ 1];
        seqstate_j= new double[sequenceLength+ 1];
    }

    /**
     * Method to initialize dynamic programming score matrix.
     */
    public void initializeScoreMatrix(){
        for(int i = 0; i < lengthHMM+1; i++){
            for (int j = 0; j < sequenceLength+1; j++){
                DyPCell c = new DyPCell(i,j);
                dpmatrix_score[i][j] = c;
            }
        }
    }

    /**
     * Main DP recursion with Branch-and-bound approach and glocal search.
     * see generic_fwdback.c and generic_viterbi.c
     */
    public void alignSequenceBnB(){

        initializeScoreMatrix();

        tempHits = new ArrayList<DyPCell>(lengthHMM);

        double [][] stateMatrix = hmm.getStateTransitionMatrix(); // state transition
        double [][] matchMatrix = hmm.getMatchEmissionMatrix();   // match emission
        double [][] insertMatrix = hmm.getInsertEmissionMatrix(); // insert emission

        String alphabet = new String(hmm.getAlphabet()); // char[] -> String

        double horizontalScore1;
        double horizontalScore2;
        double diagonalScore1;
        double diagonalScore2;
        double diagonalScore3;
        double diagonalScore4;
        double verticalScore1;
        double verticalScore2;

        // reconfig model for sequence length
        hmm.reconfigLength(sequenceLength);

        double laScore= -FunkSettings.INF; //Glocal
        double p1 = 0.997151; // null transition probability (hmmer3.1 code)

        // S->N, p=1
        seqstate_n[0]= 0;

        // S->N->B, no N-tail
        seqstate_b[0]= hmm.xsc[ProfileHMM.N][ProfileHMM.MOVE];

        // block end, c-term, join here
        seqstate_e[0]= seqstate_c[0]= seqstate_j[0]= -FunkSettings.INF;

        // init/block extension of first column => require extended state transition
        for(int i = 0; i < lengthHMM + 1; i++) {
            dpmatrix_d[i][0] =  dpmatrix_h[i][0] = dpmatrix_v[i][0] = -FunkSettings.INF;
            dpmatrix_score[i][0].setScore(-FunkSettings.INF);
        }

        char [] aminoacid = sequence.toCharArray();
        double [] nullProbArray = hmm.getNullModel();

        for(int j = 1; j < sequenceLength+ 1; j++){
            for(int i = 1; i < lengthHMM + 1; i++){

                // init/block extension of first line => require extended state transition
                dpmatrix_d[0][j] = dpmatrix_h[0][j] = dpmatrix_v[0][j] = -FunkSettings.INF;
                dpmatrix_score[0][j].setScore(-FunkSettings.INF);
                seqstate_e[j]= -FunkSettings.INF;

                char alpha = aminoacid[j-1];
                int AA = alphabet.indexOf(alpha); // index of amino acid
                double nullProb = nullProbArray[AA]*p1;

                if(Double.isNaN(dpmatrix_score[i-1][j].getScore()) && Double.isNaN(dpmatrix_score[i-1][j-1].getScore())
                        && Double.isNaN(dpmatrix_score[i][j-1].getScore())){
                    dpmatrix_score[i][j].setScore(Double.NaN);
                    dpmatrix_v[i][j] = Double.NaN;
                    dpmatrix_h[i][j] = Double.NaN;
                    dpmatrix_d[i][j] = Double.NaN;
                    continue;
                }

                // Diagonal scores
                if(Double.isNaN(dpmatrix_d[i-1][j-1])){
                    diagonalScore1 = Double.NaN;
                }
                else{
                    diagonalScore1 = dpmatrix_d[i-1][j-1]
                            + Utils.prob2LogOdd(stateMatrix[i - 1][ProfileHMM.MM]); // m -> m
                }

                if(Double.isNaN(dpmatrix_h[i-1][j-1])){
                    diagonalScore2 = Double.NaN;
                }
                else{
                    diagonalScore2 = dpmatrix_h[i-1][j-1]
                        + Utils.prob2LogOdd(stateMatrix[i - 1][ProfileHMM.IM]);  // i -> m
                }


                if(Double.isNaN(dpmatrix_v[i-1][j-1])){
                    diagonalScore3 = Double.NaN;
                }
                else{
                    diagonalScore3 = dpmatrix_v[i-1][j-1]
                            +  Utils.prob2LogOdd(stateMatrix[i - 1][ProfileHMM.DM]); // d -> m

                }
                // BM transition value stored before first node in the place of MM
                if(Double.isNaN(seqstate_b[j-1])){
                    diagonalScore4 = Double.NaN;
                }
                else{
                    // We should have a array of B->Mk transitions              
                   if (i > 1)
                       diagonalScore4= -FunkSettings.INF; // block B -> M_x transitions for all x>1
                   else
                       diagonalScore4 = seqstate_b[j - 1] + Utils.prob2LogOdd(stateMatrix[0][ProfileHMM.MM]); // b -> m

                }

                double maxDiagonal = Utils.maxWithNaN(diagonalScore1, diagonalScore2, diagonalScore3, diagonalScore4);
                dpmatrix_d[i][j] = Utils.prob2LogOdd(matchMatrix[i][AA], nullProb) + maxDiagonal;
                seqstate_e[j] = Utils.maxWithNaN(seqstate_e[j], dpmatrix_d[i][j] + laScore);


                // Horizontal scores => Insertion
                if(Double.isNaN(dpmatrix_d[i][j-1])){
                    horizontalScore1 = Double.NaN;
                }
                else{
                    horizontalScore1 =  dpmatrix_d[i][j-1] + Utils.prob2LogOdd(stateMatrix[i][ProfileHMM.MI]); // m -> i
                }


                if(Double.isNaN(dpmatrix_h[i][j-1])){
                    horizontalScore2 = Double.NaN;
                }
                else{
                    horizontalScore2 = dpmatrix_h[i][j-1] + Utils.prob2LogOdd(stateMatrix[i][ProfileHMM.II]);  // i -> i
                }

                double maxHorizontal = Utils.maxWithNaN(horizontalScore1, horizontalScore2);

                if (i== lengthHMM){
                    dpmatrix_h[i][j] = -FunkSettings.INF;    // last row, force to go to E
                }
                else{
                    dpmatrix_h[i][j] = Utils.prob2LogOdd(insertMatrix[i][AA], nullProb) +
                            maxHorizontal;
                }


                // Vertical scores => Deletion
                if(Double.isNaN(dpmatrix_v[i-1][j])){
                    verticalScore1 = Double.NaN;
                }
                else{
                    verticalScore1 = dpmatrix_v[i-1][j] + Utils.prob2LogOdd(stateMatrix[i - 1][ProfileHMM.DD]);    // d -> d
                }

                if(Double.isNaN(dpmatrix_d[i-1][j])){
                    verticalScore2 = Double.NaN;
                }
                else{
                    verticalScore2 = dpmatrix_d[i-1][j] + Utils.prob2LogOdd(stateMatrix[i - 1][ProfileHMM.MD]);   // m -> d
                }


                double maxVertical = Utils.maxWithNaN(verticalScore1, verticalScore2);
                dpmatrix_v[i][j] = maxVertical;


                if(!Double.isNaN(dpmatrix_v[i][j]) && dpmatrix_v[i][j]  // d->m->m->m->...
                        + hmm.getOptimumScoreArray(i, j, sequenceLength)  - Utils.prob2LogOdd(sequence.length()) < hmm.getGa2()){
                    dpmatrix_v[i][j] = Double.NaN;

                }

                if(!Double.isNaN(dpmatrix_h[i][j]) && dpmatrix_h[i][j]  // I->m->m->m->...
                        + hmm.getOptimumScoreArray(i, j, sequenceLength) - Utils.prob2LogOdd(sequence.length()) < hmm.getGa2()){
                    dpmatrix_h[i][j] = Double.NaN;

                }

                if(!Double.isNaN(dpmatrix_d[i][j]) && dpmatrix_d[i][j]  // M->m->m->m->...
                        + hmm.getOptimumScoreArray(i, j, sequenceLength) - Utils.prob2LogOdd(sequence.length()) < hmm.getGa2()){
                    dpmatrix_d[i][j] = Double.NaN;

                }

                double maxx = Utils.maxWithNaN(dpmatrix_v[i][j], dpmatrix_h[i][j], dpmatrix_d[i][j]);
                dpmatrix_score[i][j].setScore(maxx);

                if (maxx == dpmatrix_h[i][j]){
                    dpmatrix_score[i][j].setPrevCell(dpmatrix_score[i][j-1]);

                }else if (maxx== dpmatrix_v[i][j]){
                    dpmatrix_score[i][j].setPrevCell(dpmatrix_score[i-1][j]);

                }else if (maxDiagonal != diagonalScore4){   // no B->M transition
                        dpmatrix_score[i][j].setPrevCell(dpmatrix_score[i-1][j-1]);
                }

                // E (End) state update, only relevant for isGlocal alignment                    
                seqstate_e[j] = Utils.maxWithNaN(seqstate_e[j], dpmatrix_d[i][j],
                            dpmatrix_v[i][j]) + (i == lengthHMM ? laScore : 0);   // last row or isGlocal alignment

                // update Join, C-term, N-term, Begin states
                // J state
                seqstate_j[j]= Utils.maxWithNaN(seqstate_j[j- 1]+ hmm.xsc[ProfileHMM.J][ProfileHMM.LOOP],
                        seqstate_e[j]+ hmm.xsc[ProfileHMM.E][ProfileHMM.LOOP]);

                // C state
                seqstate_c[j] = Utils.maxWithNaN(seqstate_c[j-1] + hmm.xsc[ProfileHMM.C][ProfileHMM.LOOP],
                        seqstate_e[j] + hmm.xsc[ProfileHMM.E][ProfileHMM.MOVE]);
                
                // N state
                seqstate_n[j] = seqstate_n[j-1] + hmm.xsc[ProfileHMM.N][ProfileHMM.LOOP];

                // B state
                seqstate_b[j] = Utils.maxWithNaN(seqstate_n[j] + hmm.xsc[ProfileHMM.N][ProfileHMM.MOVE],seqstate_j[j] +
                        hmm.xsc[ProfileHMM.J][ProfileHMM.MOVE]);


                if(i == lengthHMM && !Double.isNaN(dpmatrix_score[i][j].getScore())){
                    if(dpmatrix_score[i][j].getScore() - Utils.prob2LogOdd(sequence.length()) >= hmm.getGa2()){
                        tempHits.add(dpmatrix_score[i][j]);
                    }
                }
            } // end iterating model            
        }
    }

    /**
     * Main DP recursion of local search (without BnB)
     * see generic_fwdback.c and generic_viterbi.c
     */
    public void alignSequenceLocal(){

        initializeScoreMatrix();

        tempHits = new ArrayList<DyPCell>(lengthHMM);

        double [][] stateMatrix = hmm.getStateTransitionMatrix(); // state transition
        double [][] matchMatrix = hmm.getMatchEmissionMatrix();   // match emission
        double [][] insertMatrix = hmm.getInsertEmissionMatrix(); // insert emission

        String alphabet = new String(hmm.getAlphabet()); // char[] -> String

        double horizontalScore1;
        double horizontalScore2;
        double diagonalScore1;
        double diagonalScore2;
        double diagonalScore3;
        double diagonalScore4;
        double verticalScore1;
        double verticalScore2;

        // reconfig model for sequence length
        hmm.reconfigLength(sequenceLength);

        // flag for isGlocal alignment

        double laScore= 0; // isGlocal is always false

        double p1 = 0.997151;


        // Init
        seqstate_n[0]= 0;   // S->N, p=1
        seqstate_b[0]= hmm.xsc[ProfileHMM.N][ProfileHMM.MOVE]; // S->N->B, no N-tail
        seqstate_e[0]= seqstate_c[0]= seqstate_j[0]= -FunkSettings.INF;  // block end, c-term, join here

        dpmatrix_score[0][0].setScore(-FunkSettings.INF);

        // init/block extension of first column => require extended state transition
        for(int i = 0; i < lengthHMM + 1; i++) {
            dpmatrix_d[i][0] =  dpmatrix_h[i][0] = dpmatrix_v[i][0] = -FunkSettings.INF;
            dpmatrix_score[i][0].setScore(-FunkSettings.INF);
        }


        char [] aminoacid = sequence.toCharArray();

        double [] nullProbArray = hmm.getNullModel();

        for(int i = 1; i < lengthHMM + 1; i++){

            for(int j = 1; j < sequenceLength+ 1; j++){

                // init/block extension of first line => require extended state transition

                dpmatrix_d[0][j] = dpmatrix_h[0][j] = dpmatrix_v[0][j] = -FunkSettings.INF;

                dpmatrix_score[0][j].setScore(-FunkSettings.INF);

                seqstate_e[j]= -FunkSettings.INF;

                char alpha = aminoacid[j-1];

                int AA = alphabet.indexOf(alpha); // index of amino acid

                double nullProb = nullProbArray[AA]*p1;

                // Diagonal scores
                diagonalScore1 = dpmatrix_d[i-1][j-1]
                            + Utils.prob2LogOdd(stateMatrix[i - 1][ProfileHMM.MM]); // m -> m

                diagonalScore2 = dpmatrix_h[i-1][j-1]
                            + Utils.prob2LogOdd(stateMatrix[i - 1][ProfileHMM.IM]);  // i -> m

                    diagonalScore3 = dpmatrix_v[i-1][j-1]
                            +  Utils.prob2LogOdd(stateMatrix[i - 1][ProfileHMM.DM]); // d -> m

                // BM transition value stored before first node in the place of MM
                diagonalScore4 = seqstate_b[j - 1] + Utils.prob2LogOdd(stateMatrix[0][ProfileHMM.MM]); // b -> m

                double maxDiagonal = Utils.maxWithNaN(diagonalScore1, diagonalScore2, diagonalScore3, diagonalScore4);

                dpmatrix_d[i][j] = Utils.prob2LogOdd(matchMatrix[i][AA], nullProb) + maxDiagonal;
                seqstate_e[j] = Utils.maxWithNaN(seqstate_e[j], dpmatrix_d[i][j] + laScore);

                // Horizontal scores => Insertion
                horizontalScore1 =  dpmatrix_d[i][j-1] + Utils.prob2LogOdd(stateMatrix[i][ProfileHMM.MI]); // m -> i
                horizontalScore2 = dpmatrix_h[i][j-1] + Utils.prob2LogOdd(stateMatrix[i][ProfileHMM.II]);  // i -> i

                double maxHorizontal = Utils.maxWithNaN(horizontalScore1, horizontalScore2);

                if (i== lengthHMM){
                    dpmatrix_h[i][j] = -FunkSettings.INF;    // last row, force to go to E
                }
                else{
                    dpmatrix_h[i][j] = Utils.prob2LogOdd(insertMatrix[i][AA], nullProb) +
                            maxHorizontal;
                }

                // Vertical scores => Deletion
                verticalScore1 = dpmatrix_v[i-1][j] + Utils.prob2LogOdd(stateMatrix[i - 1][ProfileHMM.DD]);    // d -> d
                verticalScore2 = dpmatrix_d[i-1][j] + Utils.prob2LogOdd(stateMatrix[i - 1][ProfileHMM.MD]);   // m -> d

                double maxVertical = Utils.maxWithNaN(verticalScore1, verticalScore2);
                dpmatrix_v[i][j] = maxVertical;

                double maxx = Utils.maxWithNaN(dpmatrix_v[i][j], dpmatrix_h[i][j], dpmatrix_d[i][j]);
                dpmatrix_score[i][j].setScore(maxx);

                if (maxx == dpmatrix_h[i][j]){
                    dpmatrix_score[i][j].setPrevCell(dpmatrix_score[i][j-1]);

                }else if (maxx== dpmatrix_v[i][j]){
                    dpmatrix_score[i][j].setPrevCell(dpmatrix_score[i-1][j]);

                }else if (maxDiagonal != diagonalScore4){   // no B->M transition
                    dpmatrix_score[i][j].setPrevCell(dpmatrix_score[i-1][j-1]);
                }

                // E (End) state update, only relevant for isGlocal alignment
                seqstate_e[j] = Utils.maxWithNaN(seqstate_e[j], dpmatrix_d[i][j],
                        dpmatrix_v[i][j]) + (i == lengthHMM ? laScore : 0);   // last row or isGlocal alignment

                // update Join, C-term, N-term, Begin states
                // J state

                seqstate_j[j]= Utils.maxWithNaN(seqstate_j[j- 1]+ hmm.xsc[ProfileHMM.J][ProfileHMM.LOOP],
                        seqstate_e[j]+ hmm.xsc[ProfileHMM.E][ProfileHMM.LOOP]);

                // C state
                seqstate_c[j] = Utils.maxWithNaN(seqstate_c[j-1] + hmm.xsc[ProfileHMM.C][ProfileHMM.LOOP],
                        seqstate_e[j] + hmm.xsc[ProfileHMM.E][ProfileHMM.MOVE]);
                // N state
                seqstate_n[j] = seqstate_n[j-1] + hmm.xsc[ProfileHMM.N][ProfileHMM.LOOP];

                // B state
                seqstate_b[j] = Utils.maxWithNaN(seqstate_n[j] + hmm.xsc[ProfileHMM.N][ProfileHMM.MOVE],seqstate_j[j] +
                        hmm.xsc[ProfileHMM.J][ProfileHMM.MOVE]);

                if(!Double.isNaN(dpmatrix_score[i][j].getScore())){
                    if(dpmatrix_score[i][j].getScore() - Utils.prob2LogOdd(sequence.length()) >= hmm.getGa2()){
                        tempHits.add(dpmatrix_score[i][j]);
                    }
                }
            } // end iterating model
        }
    }
    /**
     * Iterates list of possible hits (hitList) of dynamic programming matrix and get possible hits.
     * @param hitList List of possible non-overlapped hits.
     * @param downPosition A integer N that is the position in the sequence where starts the alignment.
     * @param t Reference transcript of a variant.
     * @param groupID String of concatenated gene/transcript IDs representing a variant.
     * @param eventRegions List of coordinates of AS regions.
     * @param firstSource First source of variant.
     * @param lastSink Last sink of variant.
     * @return Hit List.
     */
    public List<Hit> getHits(List<Hit> hitList, int downPosition, Transcript t, String groupID, DirectedRegion[] eventRegions,
                             int firstSource, int lastSink){

        if(hitList == null)
            hitList = new ArrayList<Hit>();

        for(DyPCell c: tempHits){
            if(!Double.isNaN(c.getScore()))
                hitList = scorePath(hitList, downPosition, c, t, groupID, eventRegions, firstSource, lastSink);
        }
        return hitList;
    }

    public List<Hit> getConstitutiveHits(List<Hit> hitList, Transcript t, List<ASEvent> events){

        if(hitList == null)
            hitList = new ArrayList<Hit>();

        for(DyPCell c: tempHits){
            if(!Double.isNaN(c.getScore()))
                hitList = scoreConstitutiveHitPath(hitList, c, t, events);
        }
        return hitList;
    }

    /**
     * Walk along the path of dynamic matrix, retrieving the previous cells of alignment to reconstruct
     * the track.
     * @param hitList List of possible non-overlapped hits.
     * @param downPosition A integer N that is the position in the sequence where starts the aligment.
     * @param c The cell with possible end alignment position
     * @param t Reference transcript of a variant
     * @param groupID String of concatenated gene/transcript IDs representing a variant
     * @param eventRegions List of coordinates of AS regions
     * @param firstSource First source of variant
     * @param lastSink Last sink of variant
     * @return List of possible hits
     */
    public List<Hit> scorePath(List<Hit> hitList, int downPosition, DyPCell c, Transcript t,
                                              String groupID, DirectedRegion[] eventRegions, int firstSource, int lastSink){

        // Insert a correction factor of sequence length
        double score = c.getScore() - Utils.prob2LogOdd(sequenceLength);

        if(score < hmm.getGa2()){ // test if it is a valid score
            return hitList;
        }

        int row = c.getRow();
        int endModel = c.getRow();
        int col = c.getCol();
        int endAlignment = c.getCol();

        while ((c = dpmatrix_score[row][col].getPrevCell())!= null ) {
            row = c.getRow();
            col = c.getCol();
        }

        int realStartAlignment;
        int realEndAlignment;

        realStartAlignment = downPosition + col;
        realEndAlignment = downPosition + endAlignment;

        Hit newHit = new Hit(hmm.getName(), hmm.getAcc(), hmm.getDescription(), score, realStartAlignment, realEndAlignment,
                row, endModel, lengthHMM, t.getTranscriptID());

        newHit.setGenomicInfo(t.getGene().getChromosome(),
                t.getGene().getGeneID(),
                Utils.getCDSPosition(t, realStartAlignment),
                Utils.getCDSPosition(t, realEndAlignment));

        newHit.setGroupID(groupID);

        newHit.setSource(t.isForward() ? firstSource : -1 * firstSource);
        newHit.setSink(t.isForward()? lastSink:-1*lastSink);

        if(Astafunk.isAllTranscriptHits()){
            selectAllTranscriptHits(hitList,newHit,t,eventRegions);
        }if(Astafunk.isBestGeneHits()){
            selectBestGeneHits(hitList, newHit,t,eventRegions);
        }else{
            selectBestTranscriptHits(hitList,newHit,t,eventRegions);
        }
        return hitList;
    }

    private List<Hit> scoreConstitutiveHitPath(List<Hit> hitList, DyPCell c, Transcript t, List<ASEvent> events) {
        // Insert a correction factor of sequence length
        double score = c.getScore() - Utils.prob2LogOdd(sequenceLength);

        if(score < hmm.getGa2()){ // test if it is a valid score
            return hitList;
        }

        int row = c.getRow();
        int endModel = c.getRow();
        int col = c.getCol();
        int endAlignment = c.getCol();

        while ((c = dpmatrix_score[row][col].getPrevCell())!= null ) {
            row = c.getRow();
            col = c.getCol();
        }

        int realStartAlignment;
        int realEndAlignment;

        realStartAlignment = col;
        realEndAlignment = endAlignment;

        Hit newHit = new Hit(hmm.getName(), hmm.getAcc(), hmm.getDescription(), score, realStartAlignment, realEndAlignment,
                row, endModel, lengthHMM, t.getTranscriptID());

        newHit.setGenomicInfo(t.getGene().getChromosome(),
                t.getGene().getGeneID(),
                Utils.getCDSPosition(t, realStartAlignment),
                Utils.getCDSPosition(t, realEndAlignment));


        hitList = selectNonAsHits(hitList, newHit, t, events);

        return hitList;
    }

    /**
     * Select constitutive hits, i.e., hits outside AS regions
     * @param hitList List of hits
     * @param newHit New hit to be added
     * @param t Transcript
     * @param events List of events of AS genes or null if it is an non-AS gene
     * @return List of constitutive hits
     */
    private List<Hit> selectNonAsHits(List<Hit> hitList, Hit newHit, Transcript t, List<ASEvent> events) {

        int startHit, endHit;

        //Because Event Region are given in format [x,y] where Math.abs(x) < Math.abs(y)
        // but Hit Regions are also given in format [y,x] where Math.abs(x) > Math.abs(y) (backward)

        if(newHit.getGenomicStartAlignment() > 0){
            startHit = newHit.getGenomicStartAlignment();
            endHit = newHit.getGenomicEndAlignment();
        }else{
            startHit = newHit.getGenomicEndAlignment();
            endHit = newHit.getGenomicStartAlignment();
        }

        /**
         * 1 - Check if there is a event that contains the transcript t
         *      1.1 - If an event contains the transcript, check if the prediction falls into AS region
         *      1.2 - If
         *          1.2.1 - NO, continue;
         *          1.2.2 - YES, returns;
         */

        if(events!=null) {
            for (ASEvent event : events) {
                if (Tsearch.containsTranscript(t.getTranscriptID(), event)) {

                    DirectedRegion dr = event.getRegionEvent();

                    int dirStart = Math.abs(dr.getStart());
                    int dirEnd = Math.abs(dr.getEnd());

                    int hitStart = Math.abs(startHit);
                    int hitEnd = Math.abs(endHit);

                    if (dirStart <= hitEnd && dirEnd >= hitStart) {
                        return hitList;
                    }
                }
            }
        }

        HitComparator hc = new HitComparator();
        Collections.sort(hitList, hc); // sort by score

        List<Hit> hitAuxList = new ArrayList<Hit>();

        for(Hit h: hitList) {

            // has overlapping
            if(newHit.getGenomicStartAlignment() <= h.getGenomicEndAlignment() &&
                    newHit.getGenomicEndAlignment() >= h.getGenomicStartAlignment()) {

                // CDS to AA position of h1
                // in amino acids
                int posA = t.getPrevCDSPosition(h.getGenomicStartAlignment()) / 3;
                int posB = t.getPrevCDSPosition(h.getGenomicEndAlignment()) / 3;

                // CDS to AA position of h2
                // in amino acids
                int posX = t.getPrevCDSPosition(newHit.getGenomicStartAlignment()) / 3;
                int posY = t.getPrevCDSPosition(newHit.getGenomicEndAlignment()) / 3;

                int overlappingHitRegion = posB - posA + 1; // default value 100% of overlapping

                if(posX <= posA && posA <= posY){
                    overlappingHitRegion = posY - posA + 1;
                }else if(posX<=posB && posB <= posY){
                    overlappingHitRegion = posB - posX + 1;
                }

                float diff = (float) overlappingHitRegion / (posB - posA + 1);

                // if hit overlapping > threshold.
                if (diff > Astafunk.getOverlappingThreshold()
                        && h.getAcc().equals(newHit.getAcc())) {
                    // has overlapping
                    if(h.getScore() >= newHit.getScore()){
                        return hitList;
                    }else{
                        hitAuxList.add(h);
                    }
                }
            }
        }
        hitList.removeAll(hitAuxList);
        hitList.add(newHit);
        hitAuxList.clear();
        return hitList;
    }

    /**
     * Print score matrix.
     */
    public void getScoreMatrix(){
            for(DyPCell[] row: dpmatrix_score){
                for(DyPCell c: row){
                    System.out.print(c.getScore() + "\t");
                }
                System.out.println();
            }
    }

    /**
     *
     * @return Amino acid string.
     */
    public String getSequence() {
        return sequence;
    }

    /**
     *
     * @param sequence Amino acid string to set.
     */
    public void setSequence(String sequence) {
        this.sequence = sequence;
    }

    /**
     * Select the best hits of a gene (non-overlapped alignments between all transcripts).
     * @param hitList List of Hits.
     * @param newHit New hit to be added.
     * @param t Target transcript.
     * @param eventRegions Event regions of gene of the target transcript.
     * @return List of best hits for a AS gene
     */
    public List<Hit> selectBestGeneHits(List<Hit> hitList, Hit newHit, Transcript t, DirectedRegion[] eventRegions){

        boolean insert = true;
        int startHit, endHit;

        //Because Event Region are given in format [x,y] where Math.abs(x) < Math.abs(y)
        // but Hit Regions are also given in format [y,x] where Math.abs(x) > Math.abs(y) (backward)

        if(newHit.getGenomicStartAlignment() > 0){
            startHit = newHit.getGenomicStartAlignment();
            endHit = newHit.getGenomicEndAlignment();
        }else{
            startHit = newHit.getGenomicEndAlignment();
            endHit = newHit.getGenomicStartAlignment();
        }

        List<Hit> hitAuxList = new ArrayList<Hit>();

        if(eventRegions != null){
            for(DirectedRegion dr: eventRegions){

                int dirStart = Math.abs(dr.getStart());
                int dirEnd = Math.abs(dr.getEnd());

                int hitStart = Math.abs(startHit);
                int hitEnd = Math.abs(endHit);

                if(dirStart > hitEnd || dirEnd < hitStart){
                    insert = false;
                }else{
                    insert = true;
                    break;
                }
            }
        }

        if(!insert)
            return hitList;

        HitComparator hc = new HitComparator();
        Collections.sort(hitList, hc); // sort by score

        for(Hit h: hitList) {

            // has overlapping
            if(newHit.getGenomicStartAlignment() <= h.getGenomicEndAlignment() &&
                    newHit.getGenomicEndAlignment() >= h.getGenomicStartAlignment()) {

                // CDS to AA position of h1
                // in amino acids
                int posA = t.getPrevCDSPosition(h.getGenomicStartAlignment()) / 3;
                int posB = t.getPrevCDSPosition(h.getGenomicEndAlignment()) / 3;

                // CDS to AA position of h2
                // in amino acids
                int posX = t.getPrevCDSPosition(newHit.getGenomicStartAlignment()) / 3;
                int posY = t.getPrevCDSPosition(newHit.getGenomicEndAlignment()) / 3;

                int overlappingHitRegion = posB - posA + 1; // default value 100% of overlapping

                if(posX <= posA && posA <= posY){
                    overlappingHitRegion = posY - posA + 1;
                }else if(posX<=posB && posB <= posY){
                    overlappingHitRegion = posB - posX + 1;
                }

                float diff = (float) overlappingHitRegion / (posB - posA + 1);

                // if hit overlapping > threshold.
                if (diff > Astafunk.getOverlappingThreshold()) {
                    // has overlapping
                    if(h.getScore() >= newHit.getScore()){
                        return hitList;
                    }else{
                        hitAuxList.add(h);
                    }
                }
            }
        }

        hitList.removeAll(hitAuxList);
        hitList.add(newHit);
        hitAuxList.clear();
        return hitList;
    }

    /**
     * Select all hits of each transcripts of a AS gene (overlapped alignments from differents HMMs)
     * If the hits are from different domains or transcripts, there is not overlapping.
     * @param hitList List of Hits.
     * @param newHit New hit to be added.
     * @param t Target transcript.
     * @param eventRegions Event regions of gene of the target transcript.
     * @return List of all hits for a transcript
     */
    public List<Hit> selectAllTranscriptHits(List<Hit> hitList, Hit newHit, Transcript t, DirectedRegion[] eventRegions){

        boolean insert = true;
        int startHit, endHit;

        //Because Event Region are given in format [x,y] where Math.abs(x) < Math.abs(y)
        // but Hit Regions are also given in format [y,x] where Math.abs(x) > Math.abs(y) (backward)

        if(newHit.getGenomicStartAlignment() > 0){
            startHit = newHit.getGenomicStartAlignment();
            endHit = newHit.getGenomicEndAlignment();
        }else{
            startHit = newHit.getGenomicEndAlignment();
            endHit = newHit.getGenomicStartAlignment();
        }

        List<Hit> hitAuxList = new ArrayList<Hit>();

        if(eventRegions != null){
            for(DirectedRegion dr: eventRegions){

                int dirStart = Math.abs(dr.getStart());
                int dirEnd = Math.abs(dr.getEnd());

                int hitStart = Math.abs(startHit);
                int hitEnd = Math.abs(endHit);

                if(dirStart > hitEnd || dirEnd < hitStart){
                    insert = false;
                }else{
                    insert = true;
                    break;
                }
            }
        }

        if(!insert)
            return hitList;

        HitComparator hc = new HitComparator();
        Collections.sort(hitList, hc); // sort by score

        for(Hit h: hitList) {

            // has overlapping
            if(newHit.getGenomicStartAlignment() <= h.getGenomicEndAlignment() &&
                    newHit.getGenomicEndAlignment() >= h.getGenomicStartAlignment()) {



                // CDS to AA position of h1
                // in amino acids
                int posA = t.getPrevCDSPosition(h.getGenomicStartAlignment()) / 3;
                int posB = t.getPrevCDSPosition(h.getGenomicEndAlignment()) / 3;

                // CDS to AA position of h2
                // in amino acids
                int posX = t.getPrevCDSPosition(newHit.getGenomicStartAlignment()) / 3;
                int posY = t.getPrevCDSPosition(newHit.getGenomicEndAlignment()) / 3;

                int overlappingHitRegion = posB - posA + 1; // default value 100% of overlapping

                if(posX <= posA && posA <= posY){
                    overlappingHitRegion = posY - posA + 1;
                }else if(posX<=posB && posB <= posY){
                    overlappingHitRegion = posB - posX + 1;
                }

                float diff = (float) overlappingHitRegion / (posB - posA + 1);

                // if hit overlapping > threshold AND same domain AND same transcript
                // AND new hit score < older hit score are the conditions to reject a new hit
                if (diff > Astafunk.getOverlappingThreshold()
                        && h.getAcc().equals(newHit.getAcc())
                        && h.getTranscriptID().equals(newHit.getTranscriptID())) {
                    // has overlapping
                    if(h.getScore() >= newHit.getScore()){
                        return hitList;
                    }else{
                        hitAuxList.add(h);
                    }
                }
            }
        }

        hitList.removeAll(hitAuxList);
        hitList.add(newHit);
        hitAuxList.clear();
        return hitList;
    }

    /**
     * Select best hits of each transcripts of a AS gene (non-overlapped alignments; deafault selection)
     * If the hits are from different domains or transcripts, there is not overlapping.
     * @param hitList List of Hits.
     * @param newHit New hit to be added.
     * @param t Target transcript.
     * @param eventRegions Event regions of gene of the target transcript.
     * @return List of best hits for a transcript
     */
    public List<Hit> selectBestTranscriptHits(List<Hit> hitList, Hit newHit, Transcript t, DirectedRegion[] eventRegions){

        boolean insert = true;
        int startHit, endHit;

        //Because Event Region are given in format [x,y] where Math.abs(x) < Math.abs(y)
        // but Hit Regions are also given in format [y,x] where Math.abs(x) > Math.abs(y) (backward)

        if(newHit.getGenomicStartAlignment() > 0){
            startHit = newHit.getGenomicStartAlignment();
            endHit = newHit.getGenomicEndAlignment();
        }else{
            startHit = newHit.getGenomicEndAlignment();
            endHit = newHit.getGenomicStartAlignment();
        }

        List<Hit> hitAuxList = new ArrayList<Hit>();

        if(eventRegions != null){
            for(DirectedRegion dr: eventRegions){

                int dirStart = Math.abs(dr.getStart());
                int dirEnd = Math.abs(dr.getEnd());

                int hitStart = Math.abs(startHit);
                int hitEnd = Math.abs(endHit);

                if(dirStart > hitEnd || dirEnd < hitStart){
                    insert = false;
                }else{
                    insert = true;
                    break;
                }
            }
        }

        if(!insert)
            return hitList;

        HitComparator hc = new HitComparator();
        Collections.sort(hitList, hc); // sort by score

        for(Hit h: hitList) {

            // has overlapping
            if(newHit.getGenomicStartAlignment() <= h.getGenomicEndAlignment() &&
                    newHit.getGenomicEndAlignment() >= h.getGenomicStartAlignment()) {



                // CDS to AA position of h1
                // in amino acids
                int posA = t.getPrevCDSPosition(h.getGenomicStartAlignment()) / 3;
                int posB = t.getPrevCDSPosition(h.getGenomicEndAlignment()) / 3;

                // CDS to AA position of h2
                // in amino acids
                int posX = t.getPrevCDSPosition(newHit.getGenomicStartAlignment()) / 3;
                int posY = t.getPrevCDSPosition(newHit.getGenomicEndAlignment()) / 3;

                int overlappingHitRegion = posB - posA + 1; // default value 100% of overlapping

                if(posX <= posA && posA <= posY){
                    overlappingHitRegion = posY - posA + 1;
                }else if(posX<=posB && posB <= posY){
                    overlappingHitRegion = posB - posX + 1;
                }

                float diff = (float) overlappingHitRegion / (posB - posA + 1);

                // if hit overlapping > threshold AND same domain AND same transcript
                // AND new hit score < older hit score are the conditions to reject a new hit
                if (diff > Astafunk.getOverlappingThreshold()
                        && h.getTranscriptID().equals(newHit.getTranscriptID())) {
                    // has overlapping
                    if(h.getScore() >= newHit.getScore()){
                        return hitList;
                    }else{
                        hitAuxList.add(h);
                    }
                }
            }
        }
        hitList.removeAll(hitAuxList);
        hitList.add(newHit);
        hitAuxList.clear();
        return hitList;
    }

}