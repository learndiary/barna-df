package fbi.genome.errormodel;

/**
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
class QualityTransitions {

    /**
     * Transition matrix
     */
    private long [][] transitions;
    /**
     * Number of reads
     */
    private int[] reads;
    /**
     * Number of quality states
     */
    private int size;

    public QualityTransitions(int size) {
        this.size = size;
        transitions = new long[size][size];
        reads = new int[size];
    }

    void addRead(Read r){
        for (int i = 0; i < r.getLength()-1; i++) {
            int q0 = r.getQualities()[i];
            int q1 = r.getQualities()[i+1];
            transitions[q0][q1]++;
            reads[q0]++;
        }
    }


    public double getTransition(int q1, int q2){
        double sum = (double) reads[q1];
        if(sum == 0) return 1d/size; // equally distributed ?
        return transitions[q1][q2]/ sum;
    }

    public int getNext(int current, double r){
        long[] cs = transitions[current];
        double sum = 0;
        for (int i = 0; i < size; i++) {
            sum += getTransition(current, i);
            if(sum >= r) return i;
        }
        return -1;
    }


}
