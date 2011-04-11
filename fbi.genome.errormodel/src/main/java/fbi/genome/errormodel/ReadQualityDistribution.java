package fbi.genome.errormodel;

/**
 * Distribution of the average read quality
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
class ReadQualityDistribution extends Distribution{

    public ReadQualityDistribution(int size) {
        super(size);
    }

    public void addRead(Read read){
        reads++;
        int[] q = read.getQualities();
        int s = 0;
        for (int i = 0; i < q.length; i++) {
            s += q[i];
        }
        int avg = s/read.getLength();
        values[avg]++;
    }
}
