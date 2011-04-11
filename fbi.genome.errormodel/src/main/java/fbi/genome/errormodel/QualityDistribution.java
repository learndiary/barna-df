package fbi.genome.errormodel;

/**
 * Distribution of all quality values, independent of the position
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
class QualityDistribution extends Distribution{

    public QualityDistribution(int size) {
        super(size);
    }

    public void addRead(Read read){
        reads+=read.getLength();
        int[] q = read.getQualities();
        for (int i = 0; i < q.length; i++) {
            values[q[i]]++;
        }
    }
}
