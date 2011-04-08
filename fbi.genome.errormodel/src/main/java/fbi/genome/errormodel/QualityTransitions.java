package fbi.genome.errormodel;

/**
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
class QualityTransitions {

    /**
     * Transition matrix
     */
    private long [][][] transitions;
    /**
     * Number of reads
     */
    private int[][] reads;
    /**
     * Number of quality states
     */
    private int size;
    /**
     * The read length
     */
    private int readLength;

    public QualityTransitions(int size, int readLength) {
        this.size = size;
        this.readLength = readLength;
        transitions = new long[readLength][size][size];
        reads = new int[readLength][size];
    }

    void addRead(Read r){
        for (int i = 0; i < r.getLength()-1; i++) {
            int q0 = r.getQualities()[i];
            int q1 = r.getQualities()[i+1];
            transitions[i][q0][q1]++;
            reads[i][q0]++;
        }
    }


    public double getTransition(int position, int q1, int q2){
        double sum = (double) reads[position][q1];
        if(sum == 0) return 1d/size; // equally distributed ?
        return transitions[position][q1][q2]/ sum;
    }

    public int getNext(int position, int current, double r){
        long[] cs = transitions[position][current];
        double sum = 0;
        for (int i = 0; i < size; i++) {
            sum += getTransition(position,current, i);
            if(sum >= r) return i;
        }
        return -1;
    }


}
