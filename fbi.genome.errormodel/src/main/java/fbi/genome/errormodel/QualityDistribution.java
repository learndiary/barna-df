package fbi.genome.errormodel;

/**
 * Distribution of all quality values, independent of the position
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
class QualityDistribution extends Distribution{
    /**
     * If true, we trim the end of reads to cu away all qualities scores == 2
     */
    private static final boolean TRIM_ILLUMINA_BAD_ENDS = false;
    public QualityDistribution(int size) {
        super(size);
    }

    public void addRead(Read read){
        // trim ends

        int finalPositon =read.getLength()-1;
        int[] q = read.getQualities();
        if(TRIM_ILLUMINA_BAD_ENDS){
            for(;finalPositon > 0;finalPositon--){
                if(q[finalPositon] != 2) break;
            }
        }
        reads+=finalPositon+1;

        //for (int i = 0; i < q.length; i++) {
        for (int i = 0; i < finalPositon; i++) {
            values[q[i]]++;
        }
    }
}
