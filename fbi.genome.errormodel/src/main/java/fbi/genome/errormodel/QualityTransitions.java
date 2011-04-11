package fbi.genome.errormodel;

/**
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class QualityTransitions {

    /**
     * Transition matrix
     */
    private long [][][] transitions;
    /**
     * Number of reads
     */
    private int[][] reads;

    /**
     * The distribution for position 0
     */
    private int[] initialDistribution;

    /**
     * Number of quality states
     */
    private int size;

    /**
     * The read length
     */
    private int readLength;

    /**
     * The number of reads added
     */
    private int numReads;

    QualityTransitions(int size, int readLength) {
        this.size = size;
        this.readLength = readLength;
        transitions = new long[readLength][size][size];
        reads = new int[readLength][size];
        initialDistribution = new int[size];
    }

    void addRead(Read r){
        numReads++;
        int[] readQualities = r.getQualities();

        // position 0 distribution
        initialDistribution[readQualities[0]]++;

        // positions 1..n transitions
        for (int i = 1; i < r.getLength(); i++) {
                int q0 = readQualities[i-1];
                int q1 = readQualities[i];

                transitions[i][q0][q1]++;
                reads[i][q0]++;
        }
    }


    public int getQuality(int position, int lastQualityValue, double random){
        // if start, sample from initial position 0 distribution
        if(position == 0) return getInitialQuality(random);

        // otherwise sample from transitions
        double sum = 0;
        for (int i = 0; i < size; i++) {

            double numberOfReadsFromQ0 = (double) reads[position][lastQualityValue];
            if(numberOfReadsFromQ0 == 0) {
                numberOfReadsFromQ0 = 1d/size; // equally distributed ?
            }
            sum += transitions[position][lastQualityValue][i]/ numberOfReadsFromQ0;
            if(sum >= random) return i;
        }
        return -1;
    }

    private int getInitialQuality(double r){
        double sum = 0;
        for (int i = 0; i < size; i++) {
            sum += initialDistribution[i]/(double)numReads;
            if(sum >= r) return i;
        }
        return -1;
    }

}
