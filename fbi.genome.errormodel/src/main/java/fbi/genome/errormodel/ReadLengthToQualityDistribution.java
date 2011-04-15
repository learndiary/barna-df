package fbi.genome.errormodel;

/**
 * Distribution of the quality per read position
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
class ReadLengthToQualityDistribution extends Distribution{

    public ReadLengthToQualityDistribution(int size) {
        super(size);
    }

    public void addRead(Read read){
        reads++;
        int[] q = read.getQualities();
        for (int i = 0; i < size; i++) {
            values[i] += q[i];
        }
    }
}
