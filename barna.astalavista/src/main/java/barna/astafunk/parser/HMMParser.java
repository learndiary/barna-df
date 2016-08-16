/**
 * @author Vitor Lima Coelho
 * @version 1.1
 */
package barna.astafunk.parser;

import barna.commons.log.Log;
import barna.astafunk.HMM.ProfileHMM;
import barna.astafunk.utils.Utils;
import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

/**
 * This class describes a profile HMM parser.
 */
public class HMMParser{

    /**
     * An BufferRead attribute. Reads text from a character-input stream, buffering characters so as to
     * provide for the efficient reading of characters, arrays, and lines.
     * @see java.io.BufferedReader
     */
    BufferedReader stdin;

    /**
     * Initialize the BufferRead stdin attribute
     * @param path profile HMM file.
     * @throws FileNotFoundException
     */
    public HMMParser(String path) throws FileNotFoundException {
        this.stdin = new BufferedReader(new FileReader(path), 256);
    }

    /**
     * This method parses a profile HMM and returns the attribute ID values.
     * @param attributeID An attribute ID from profile HMM file format (e.g., ACC, DESC, GA1, etc).
     * @return A string with the value(s) from attribute ID.
     */
    public String parseAttribute(String attributeID) throws IOException {

        // A line of input file.
        String auxLine;
        // Returned values string
        String attribute = "";

        String attributeAux = attributeID;
        if(attributeID.equals("GA1") || attributeID.equals("GA2")){
            attributeAux = "GA";
        }

        while((auxLine = stdin.readLine()) != null){
            String [] line = auxLine.split("\\s++");


            if(line[0].equals(attributeAux)){
                for(int i = 1; i < line.length; i++){
                    if(attributeID.equals("LENG"))
                        attribute += line[i];
                    else if (attributeID.equals("GA1"))
                        attribute = line[1];
                    else if (attributeID.equals("GA2"))
                        attribute = line[2].replace(";","");
                    else attribute += line[i] + " ";
                }
                break;
            }else if((attributeAux.equals("ACC") && (line[0].equals("DESC") || line[0].equals("LENG")))
                    || (attributeAux.equals("DESC") && line[0].equals("LENG"))){
                break;
            }else if(attributeAux.equals("GA") && line[0].equals("HMM")){
                attribute = "0";
                break;
            }
        }
        return attribute;
    }

    private static double parseValue(String s) {

        return Utils.naturalNegLog2Prob(s);

        // TODO debug
//        if (s.equals("*"))
//            return Double.NEGATIVE_INFINITY;
//        return -Double.parseDouble(s);
    }

    /**
     * Read HMM file and parse match emission probability distributions lines.
     * @param hmm Target profile HMM.
     * @return Match emission probability distributions matrix.
     * @throws IOException
     */
    private ProfileHMM parseModel(ProfileHMM hmm) throws IOException {

        String buffer;
        String [] tokens;

        int alphabetLength = hmm.getAlphabet().length;
        int modelLength = hmm.getLength();

        double [][] insertMatrix = new double[modelLength+1][alphabetLength];
        double [][] matchMatrix = new double[modelLength+1][alphabetLength];
        double [][] stateMatrix = new double[modelLength+1][7];

        double []compo = new double[alphabetLength];
        buffer = stdin.readLine().trim();

        while(!buffer.startsWith("COMPO")){
            buffer = stdin.readLine().trim();
        }

        tokens = buffer.trim().split("\\s++");

        for(int x = 0; x < alphabetLength; x++)
            compo[x] = parseValue(tokens[x + 1]);

        buffer = stdin.readLine();
        // parse emission line 0
        tokens = buffer.trim().split("\\s++");

        for(int x = 0; x < alphabetLength; x++){
            insertMatrix[0][x] = parseValue(tokens[x]);
        }

        // parse 0 state transitions
        buffer = stdin.readLine();
        tokens = buffer.trim().split("\\s++");

        for(int x = 0; x < 7; x++)
            stateMatrix[0][x] = parseValue(tokens[x]);

        double minAlpha = 0;
        double minBeta = 0;
        double maxInsertProb = 0;

        for(int k = 1; k < modelLength+1; k++){
            // match emission of line 1
            buffer= stdin.readLine();
            tokens = buffer.trim().split("\\s++");

            for (int x = 0; x < alphabetLength; x++){
                matchMatrix[k][x] = parseValue(tokens[x + 1]);
            }
            // Insert emissions of line 2
            buffer= stdin.readLine().trim();
            tokens = buffer.split("\\s++");


            // get the max probability of insertion

            for (int x = 0; x < alphabetLength; x++){
                double ins_prob = parseValue(tokens[x]);
                insertMatrix[k][x] = ins_prob;

                if (ins_prob > maxInsertProb){
                    maxInsertProb = ins_prob;
                    hmm.INDEX_MAX_PROB = x;

                }
            }
            // transitions of line 3
            buffer= stdin.readLine().trim();
            tokens = buffer.split("\\s++");

            for (int x = 0; x < 7; x++){
                double value = parseValue(tokens[x]);
                stateMatrix[k][x] = value;
                // Calculate alpha and beta
                // Here, we are working with probabilities.
                // So, if I have the highest probability, I have the lowest cost in neg nat log (what we want).
                if(x == 1){
                    if(value >minAlpha)
                        minAlpha = value;
                }else if(x == 4)
                    if(value > minBeta)
                        minBeta = value;
            }
        }
        stdin.readLine();
        hmm.ALPHA = minAlpha;
        hmm.BETA = minBeta;
        hmm.MAX_INSERT_PROB = maxInsertProb;
        hmm.setMatchEmissionMatrix(matchMatrix);
        hmm.setInsertEmissionMatrix(insertMatrix);
        hmm.setStateTransitionMatrix(stateMatrix);
        calcDeltaFactor(hmm);

        return hmm;
    }

    /**
     * Parse header and model of profile HMM file. The attributes are parsed in sequence like in file.
     * @return A profile HMM object
     * @throws IOException
     */
    public ProfileHMM parse() throws IOException {
        /**
         * In this section, we parsing the HMM file.
         */
        //Reading profile HMM file.
        // Parsing attributes
        String st_name = parseAttribute("NAME");
        String st_acc = parseAttribute("ACC").trim(); // ACC is not mandatory, so it will not work for Pfam-B, e.g.
        String st_desc = parseAttribute("DESC"); // DESC is not mandatory!
        String st_length = parseAttribute("LENG");
        String st_alpha = parseAttribute("ALPH");

        //Per-domain threshold.
        String st_ga2 = parseAttribute("GA2"); // GA is not mandatory, so it will not work for Pfam-B, e.g.

        //String st_tc = parseAttribute("TC"); // TC is not mandatory
        //String st_nc = parseAttribute("NC"); // NC is not mandatory
        int length = Integer.parseInt(st_length);

        //Creating a profile HMM object
        ProfileHMM phmm = new ProfileHMM(length);

        // Setting the attributes
        phmm.setName(st_name);
        phmm.setAcc(st_acc);
        phmm.setDescription(st_desc);

        phmm.setGa2(Double.parseDouble(st_ga2));
        //phmm.setTc(st_tc);
        //phmm.setNc(st_nc);
        phmm.setLength(length);
        phmm.setAlphabet(st_alpha);

        return parseModel(phmm);
    }

    /**
     * Calculate delta factor (Alignment Window restriction).
     * @param hmm profile HMM
     */

    private static void calcDeltaFactor(ProfileHMM hmm) {
        int INDEX_AA = 0;
        int PROB  = 1;

        double [] optimumArrayScore = new double[hmm.getLength()+1];
        double [][] stateMatrix = hmm.getStateTransitionMatrix(); // state transition
        double [][] matchMatrix = hmm.getMatchEmissionMatrix();   // match emission
        double [] nullProbArray = hmm.getNullModel();
        double omega =  0;
        optimumArrayScore[0] = omega;
        double score;

        for(int i = 1; i < hmm.getLength() + 1 ; i++){
            double [] index_prob = calcMinIndexProb(matchMatrix, i);

            int AA = (int)(index_prob[INDEX_AA]); // index of amino acid

            score = omega + Utils.prob2LogOdd(stateMatrix[i - 1][ProfileHMM.MM]);

            omega = Utils.prob2LogOdd(index_prob[PROB], nullProbArray[AA] * 0.997151) + score;

            optimumArrayScore[i] = omega;

        } // end iterating model

        hmm.setOptimumScoreArray(optimumArrayScore);

        double alpha = Utils.prob2LogOdd(hmm.ALPHA)
                + Utils.prob2LogOdd(hmm.MAX_INSERT_PROB); // ALPHA is the highest probability (lowest cost)

        double beta = Utils.prob2LogOdd(hmm.BETA)
                + Utils.prob2LogOdd(hmm.MAX_INSERT_PROB); // minimum insert score

        double ga2 = hmm.getGa2();

        double delta = ((omega  - Utils.prob2LogOdd(hmm.getLength()) - ga2
                - Math.abs(alpha))/Math.abs(beta)) + (((hmm.getLength()*3)-1))/3;

        delta = Math.floor(delta+1);

        hmm.setDeltaFator((int) delta);
    }

    /**
     * Find index in match matrix with lowest cost (highest prob)
     * @param matchMatrix Matrix with match probability distribution
     * @param row State or row of matrix
     * @return Index of amino acid in matrix
     */
    public static double[] calcMinIndexProb(double[][] matchMatrix, int row){

        double [] index_prob = new double[2];
        double minProb = 0;

        for(int j = 0; j < 20; j++){
            if(matchMatrix[row][j] > minProb){
                index_prob[0] = j;
                index_prob[1] = matchMatrix[row][j];
                minProb = matchMatrix[row][j];
            }
        }
        return index_prob;
    }


    public String nextModelLine() throws IOException {
        return stdin.readLine();
    }

    /**
     * Parse all models in a HMM file.
     * @return List of all profile hmm in HMM datafile
     * @throws IOException
     */
    public List<ProfileHMM> parseAll() throws IOException {

        List<ProfileHMM> list = new ArrayList<ProfileHMM>();
        String continueInFile = "s";

        Log.message("\nReading profile HMM dataset...");

        int i = 0;
        //Utils.update(i, 20000);

        int totalEstimative = 18000;

        Utils.initBar();

        Utils.updateBar(i, totalEstimative);

        while(continueInFile!=null){
            // Profile HMM ready
            ProfileHMM phmm = parse();
            list.add(phmm);
            continueInFile = nextModelLine();
            i++;
            Utils.updateBar(i, totalEstimative);
        }

        Utils.updateBar(totalEstimative, totalEstimative);


        Log.message(" done.\n");
        stdin.close();
        return list;
    }

    /**
     * Parse all models in a HMM file.
     * @return List of all profile hmm in HMM datafile
     * @throws IOException
     */
    public HashMap<String, ProfileHMM> parseAll2() throws IOException {

        HashMap<String, ProfileHMM> hashmap = new HashMap<String, ProfileHMM>();

        String continueInFile = "s";

        Log.message("\nReading profile HMM dataset...");

       // int i = 0;
        //Utils.update(i, 20000);

       // int totalEstimative = 18000;

      //  Utils.initBar();

      //  Utils.updateBar(i, totalEstimative);

        while(continueInFile!=null){
            // Profile HMM ready
            ProfileHMM phmm = parse();
            hashmap.put(phmm.getAcc(), phmm);
            continueInFile = nextModelLine();
     //       i++;
    //        Utils.updateBar(i, totalEstimative);
        }
    //    Utils.updateBar(totalEstimative, totalEstimative);
        Log.message(" done.\n");
        stdin.close();
        return hashmap;
    }
}