package fbi.genome.sequencing.rnaseq.simulation.fragmentation;

import fbi.commons.ByteArrayCharSequence;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

/**
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class FragmentNebulization implements FragmentProcessor{
    // 2.85 - <0.5%
    private static final double CUT_OFF_GAUSSIAN_VAL = 2.85f;

    /**
     * Number of Iterations for recursive nebulization.
     */
    private int nebuRecursionDepth = -1;
    //private double nebuC;
    private double lambda = 0;
    private double M = 0;

    private int[] index1;

    private Random rndBreak = new Random();
    private Random rndBP = new Random();
    //private double thold;
    private double maxLen;
    //private RandomDataImpl randomData;

    public FragmentNebulization(final double lambda, final double m, final double thold, double maxLen) {
        if(maxLen <= 0) throw new IllegalArgumentException("Max length <= 0 not permitted !");
        //randomData = new RandomDataImpl();
        //this.thold = thold;
        this.maxLen = maxLen;
        //double nebuC = lambda * (1.5d - Math.pow(-Math.log(0.5d), 1d / M));
        //double nebuC = 0d;


        // expected recursion depth
        // p_t<= 0.1
        // tmax= ceil( log2((maxlen-C)/(lambda*(-ln(0.9)^(1/M))) )
        // e.g. len= 1500 -> ceil(0.53), len= 15k -> ceil(4.37)
        // maxLen= 10000; lambda= 900; nebuC= 486;
        this.nebuRecursionDepth = (int) Math.ceil(
                //1d - (1d/Math.exp(maxLen/lambda))
                //Math.log10((maxLen - nebuC) / (lambda * Math.pow(-Math.log(1d - thold), 1d / M))) / Math.log10(2)
                maxLen / lambda

        );
        //this.nebuC = nebuC;
        this.lambda = lambda;
        M = m;

    }

    @Override
    public List<Fragment> process(final ByteArrayCharSequence id, final ByteArrayCharSequence cs, int start, final int end, final int len) {
        if (nebuRecursionDepth < 1) {
            List<Fragment> ll = new ArrayList<Fragment>();
            ll.add(new Fragment(id, start, end));
            return ll;
        }

        int recDepth = nebuRecursionDepth;

        if (index1 == null) {
            index1 = new int[(int) Math.pow(2, recDepth)];
        }
        Arrays.fill(index1, -1);
            index1[0] = len;
        int fragmentNb = 1;

        // now break it!
        for (int i = 0; i < recDepth; ++i) {
            for (int j = 0; index1[j] > 0; ++j) {

                // breakpoint location
                // N(length/2,sigma)= (N(0,1)*sigma*(length/2)+ length/2
                int L = index1[j];

                //int bp= (int) randomData.nextGaussian(len/ 2d, len/ 4d);
                //int bp= (int) ((rndBP.nextGaussian()* sigma
                //		* ((L-1)/2d))+ (L-1)/2d);	// bp index [0;L[
                int bp = (int) nextGaussianDouble(rndBP, 0, L - 1);

                // breaking probability (pb)
                // pb= 1- exp^(-((x-C)/lambda)^M)
                int minL = (bp < (L - bp) ? bp + 1 : L - bp - 1);
                //double pb = minL < nebuC ? 0d : 1d - Math.exp(-Math.pow((minL - nebuC) / lambda, M));
                double pb = 1d - (1d/Math.exp(Math.pow(minL / lambda, M)));
                double r = rndBreak.nextDouble();
                if (r > pb) {
                    continue;
                }

                // fragment j breaks
                int rest = fragmentNb - j - 1;
                if (rest > 0) {
                    assert (j + 2 + rest <= index1.length);
                    System.arraycopy(index1, j + 1, index1, j + 2, rest);
                }
                assert (bp + 1 > 0 && L - bp - 1 > 0);
                index1[j] = bp + 1;
                index1[j + 1] = L - bp - 1;    // L- (bp+1)
                ++fragmentNb;
                ++j;
            }

        }

        // write result to disk
        List<Fragment> fragments = new ArrayList<Fragment>();
        for (int j = 0; index1[j] > 0; ++j) {
            int s = start;
            int e = s + index1[j];
            Fragment fragment = new Fragment(id, s, e);
            fragments.add(fragment);
            //cs.replace(0, start);
            //start += index1[j];
            //cs.replace(1, start - 1);
            //fragments.add(cs.toString());
        }
        return fragments;
    }


    /**
     * min<= r <= max
     *
     * @param random
     * @param min
     * @param max
     * @return
     */
    private static strictfp double nextGaussianDouble(Random random, double min, double max) {
        double rdm = 3d;    // gaussian value, stddev 1
        while (rdm < -CUT_OFF_GAUSSIAN_VAL || rdm > CUT_OFF_GAUSSIAN_VAL) {
            rdm = random.nextGaussian();
        }
        double mid = min + (max - min) / 2f;
        double realValue = mid + (rdm * (max - mid) / CUT_OFF_GAUSSIAN_VAL);    // 0..CUT_OFF_GAUSSIAN = mid..max
        assert (realValue >= min && realValue <= max);
        return realValue;
    }

    @Override
    public String getName() {
        return "Nebulization";
    }

    @Override
    public String getConfiguration() {
        StringBuffer b = new StringBuffer();
        b.append("\t\t").append("Lambda: ").append(this.lambda).append("\n");
        //b.append("\t\t").append("Threshold: ").append(this.thold).append("\n");
        b.append("\t\t").append("M: ").append(M).append("\n");
        //b.append("\t\t").append("C: ").append(nebuC).append("\n");
        b.append("\t\t").append("Max Length: ").append(maxLen).append("\n");
        b.append("\t\t").append("Recursions: ").append(nebuRecursionDepth).append("\n");
        return b.toString();

    }

    @Override
    public String done() {
        return null;
    }

}
