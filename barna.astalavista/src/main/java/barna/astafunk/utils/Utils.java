package barna.astafunk.utils;


import barna.astafunk.DP.Hit;
import barna.model.DirectedRegion;
import barna.model.Transcript;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * @version *
 * @autor vitorlc
 * @since 17/05/14
 */

/**
 * Class with several useful methods.
 */
public class Utils {


    public static int finishedThread;

    /**
     * Print the labels of the hit output.
     */
    public static void printHeaderHitList(){
        System.out.println("#" +"\t" + "chr" + "\t" + "gene_cluster_name" + "\t" + "name_hmm" + "\t" + "acc" + "\t"
                + "description" + "\t" + "bitscore" + "\t" + "start_seq" + "\t" + "end_seq" + "\t"
                + "start_genomic" + "\t" + "end_genomic"
                +"\t" + "first_source" + "\t" + "last_sink" + "\t" + "start_model" + "\t" + "end_model" + "\t" + "length_model"
                + "\t" + "sequence" + "\t" + "variants");
    }

    /**
     * Convert a number in log2 base
     * @param n A number 10-base
     * @return Log2(n)
     */
    public static double prob2LogOdd(double n){
        return prob2LogOdd(n, 1);
    }

    /**
     * Calculate log2(e/n)
     * @param e Number 10-base
     * @param n A number 10-base
     * @return log2 of (e divided by n).
     */
    public static double prob2LogOdd(double e, double n){

        if(e/n == 0)
            return -FunkSettings.INF;

        return Math.log(e/n)/Math.log(2);
    }

    /**
     * Converts a natural negative log number into probability.
     * @param number A number.
     * @return Probability
     */
    public static double naturalNegLog2Prob(String number){
        if(number.equals("*"))
            return 0;
        else
            return Math.exp((-1)*Double.parseDouble(number));
    }

    /**
     * Returns the greater of two {@code double} values.
     * If an argument is NaN, the other will be the greater.
     * @param d1 double number 1
     * @param d2 double number 2
     * @return the larger of {@code d1} and {@code d2}.
     */
    public static double maxWithNaN(double d1, double d2){
        double max;
        if(!Double.isNaN(d1) && !Double.isNaN(d2)){
            max = Math.max(d1,d2);
        }else if(!Double.isNaN(d1) && Double.isNaN(d2)){
            max = d1;
        }else if(!Double.isNaN(d2) && Double.isNaN(d1)){
            max = d2;
        }else
            max = Double.NaN;
        return max;
    }

    /**
     * Returns the greater of three {@code double} values.
     * If an argument is NaN, the other will be the greater.
     * @param d1 one number
     * @param d2 another number
     * @param d3 another number
     * @return the larger of {@code d1}, {@code d2} and {@code 3}.
     */
    public static double maxWithNaN(double d1, double d2, double d3){
        return maxWithNaN(maxWithNaN(d1, d2), d3);
    }

    /**
     * Returns the greater of two {@code double} values.
     * If an argument is NaN, the other will be the greater.
     * @param d1 a number
     * @param d2 another number
     * @param d3 another number
     * @param d4 another number
     * @return the larger of {@code d1}, {@code d2}, {@code d3} and {@code d4}.
     */
    public static double maxWithNaN(double d1, double d2, double d3, double d4){
        return maxWithNaN(maxWithNaN(d1, d2, d3), d4);
    }

    /**
     * Get a valid position on proteomic coordinates.
     * @param pos [1...|CDS Sequence|]
     * @return proteomic coordinate.
     */
    public static int getCDSPosition(Transcript transcript, int pos){

        DirectedRegion[] reg = transcript.getCDSRegions();

        Arrays.sort(reg, new DirectedRegion.DirectedPositionComparator());

        int x;
        int dist = 0;

        for(x = 0; dist <= pos*3 && x < reg.length; x++ ){
            dist+=reg[x].getLength();
        }

        if(x > 0){
            --x;
            dist-= reg[x].getLength();
        }

        int genomicPos;
        if(transcript.isForward())
            genomicPos= reg[x].getStart() + (pos*3 - dist);
        else
            genomicPos = reg[x].getEnd() + (pos*3 - dist);
        return genomicPos;
    }

    /**
     * Print error message given by code {@code error}.
     * @param t Transcript.
     * @param groupID Variant group ID of transcript {@code t}
     * @param error Error code
     * @return Error message.
     */
    public static String printErrorOutput(Transcript t, String groupID, int error) {

        String info = null;

        switch (error){
            case 1: info="There is not CDS event.";
                break;
            case 2: info="The ORF is not valid or annotated.";
                break;
            case 4: info ="There are not valid or annotated ORF reference of this gene.";
                break;
        }
        return "ERROR" + "\t" +
                t.getChromosome() + "\t"+
                t.getGene().getGeneID() + "\t"+
                info + "\t" +
                t.getTranscriptID() + "\t" +
                groupID + "\t";
    }

    /**
     * Print the hit output.
     * @param hitHashMap Hash map of hits. Key 0 are the list of all hits found on the whole sequence.
     *                   Key 1 are the list of hits found on AS region.
     */
    public static void printDefaultHit(HashMap<Integer, List<Hit>> hitHashMap) {
        for(Map.Entry<Integer, List<Hit>> entry: hitHashMap.entrySet()){
            if(entry.getKey().equals(1)){
                for(Hit hit : entry.getValue()) {
                    System.out.println("EHIT" + "\t" + hit.toString());
                }
            }else{
                for(Hit hit : entry.getValue()) {
                    System.out.println("HIT" + "\t" + hit.toString());
                }
            }
        }
    }

    /**
     * Print the hit output.
     * @param hitList List of hits. Key 0 are the list of all hits found on the whole sequence.
     *                   Key 1 are the list of hits found on AS region.
     */
    public static void printDefaultHit(List<Hit> hitList) {
                for(Hit hit : hitList) {
                    System.out.println("EHIT" + "\t" + hit.toString());
                }
    }

    /**
     * Is astafunk running?
     * @param done Arbitrary value
     * @param total random value
     */
    public static void update(int done, int total) {
        char[] workchars = {'|', '/', '-', '\\'};

        System.err.print("\r"+workchars[done % workchars.length]);

        if (done == total) {
            System.out.flush();
        }
    }

    private static StringBuilder progress;

    /**
     * called whenever the progress bar needs to be updated.
     * that is whenever progress was made.
     *
     * @param done an int representing the work done so far
     * @param total an int representing the total work
     */
    public static void updateBar(int done, int total) {

        if (progress == null)
            initBar();
        char[] workchars = {'|', '/', '-', '\\'};
        String format = "\r%3d%% %s %c";

        int percent = (++done * 100) / total;

        if(done > total)
            percent = (total * 100) / total;

        int extrachars = (percent / 2) - progress.length();

        while (extrachars-- > 0) {
            progress.append('#');
        }

        System.err.printf(format, percent, progress,
                workchars[done % workchars.length]);

        if (done == total) {
            System.out.flush();
            System.out.println();
            initBar();
        }
    }

    public static void initBar() {

        finishedThread = 0;
        progress = new StringBuilder(60);
    }


}