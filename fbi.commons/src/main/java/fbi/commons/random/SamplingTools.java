package fbi.commons.random;

import fbi.commons.Log;
import org.apache.commons.math.random.RandomDataImpl;

import java.util.Random;

/**
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class SamplingTools {
    public static double sampleLinear(double rank, double par1, double par2) {
        return par1 * rank + par2;
    }

    public static double samplePareto(RandomDataImpl rnd, double par1) {
        double d = 0;
        while (d == 0) {
            d = rnd.nextUniform(0, 1);
        }
        double val = (1 / Math.pow(d, (1d / par1)));
        return val;
    }

    /**
     * Generates pseudo-random numbers under a linear distribution
     * with a<>0.
     *
     * @param u
     * @param a
     * @param len
     * @param more5
     * @return
     */
    public static int sampleInversionLinear(double u, double a, int len, boolean more5) {
//		if (a== 0) {
//			System.err.println("[OOPS] linear function failed a==0");
//			return 0;
//		}
//
//		boolean invert= false;
//		if (a< 0) {
//			a= -a;
//			invert= true;
//		}
//
        double b = 0;//len* a;
        //u*= len;

//		double r= -(b/a)+ Math.sqrt(Math.pow(b/a, 2)+ (2*u/a));
//		double max= -(b/a)+ Math.sqrt(Math.pow(b/a, 2)+ (2/a));

        double r = Math.sqrt(2 * u / a);
        double max = Math.sqrt(2 / a);
        r /= max;
        r *= len;
        if (more5) {
            r = len - r;
        }

//		if (true/*new Double(r).equals(Double.NaN)*/) {
//			double x1= Math.pow(b/a, 2);
//			double x2= 2*u/a;
//			double x3= Math.sqrt(x1+ x2);
//			double x4= -(b/a);
//			double x5= x4+ x3;
//			System.currentTimeMillis();
//		}


        return (int) Math.round(r);
    }

    public static double sampleLinearUnderTrpt(double u, int len, double raiseAlong) {
        double a = raiseAlong;
        double b = (a > 0) ? 0 : len * a;
        double r = sampleInversionLinear(u, a, b);
        return r;
    }

    public static double sampleMPVexpr(double rank, double par1, double par2) {
//		double val= par2/ Math.pow(rank, par1)
//			* Math.exp(- (Math.pow(rank, 2)/ 122000000)
//					- (rank/7000));
        double val = par2 / Math.pow(rank, par1);    // par1= 0,6  par2= (41627d/ 2731598d)
        return val;
    }

    public static double sampleMPVdecay(double rank, double par1) {
//		double val= Math.exp(- (Math.pow(rank, 2)/ Math.pow(par1, 2))
//				- (rank/par1));
        double val = Math.exp(-(Math.pow(rank, 2) / 122000000)
                - (rank / 7000));

        return val;
    }


    private class FragmentA implements Comparable {
        int start = 0, end = 0;

        public FragmentA(int a, int b) {
            start = a;
            end = b;
        }

        public int length() {
            return end - start + 1;
        }

        //@Override
        public int compareTo(Object o) {
            int len2 = ((FragmentA) o).length();
            return length() - len2;
        }

        @Override
        public String toString() {
            return "[" + start + "," + end + "]";
        }
    }

    /**
     * Generates pseudo-random numbers under a linear distribution
     * with a<>0.
     *
     * @param u
     * @param a
     * @param b
     * @return
     */
    // f^{-1}(x)= (y-b)/a
    public static double sampleInversionLinear(double u, double a, double b) {
        if (a == 0) {
            Log.error("[OOPS] linear function failed a==0");
            return u;
        }

        boolean invert = false;
        if (a < 0) {
            a = -a;
            invert = true;
        }

        double r = -(b / a) + Math.sqrt(Math.pow(b / a, 2) + (2 * u / a));

        if (true/*new Double(r).equals(Double.NaN)*/) {
            double x1 = Math.pow(b / a, 2);
            double x2 = 2 * u / a;
            double x3 = Math.sqrt(x1 + x2);
            double x4 = -(b / a);
            double x5 = x4 + x3;
        }


        return r;
    }

    /**
     * Returns a gaussian value sampled from a normal (gaussian)
     * distribution with >90% of the function covering the range
     * between the specified boundaries. Subsequently, the mean
     * of the distribution is at <code>min+ ((max-min)/ 2)</code>
     * and the standard deviation is <code>(max- mean)/ 2.85</code>
     * since <0.5 of a gaussian distribution exceed 2.85.
     */
    private static final double CUT_OFF_GAUSSIAN_VAL = 2.85f;

    public static double nextGaussianBetween(Random random, double min, double max) {
        double rdm = 3d;    // gaussian value, stddev 1
        while (rdm < -CUT_OFF_GAUSSIAN_VAL || rdm > CUT_OFF_GAUSSIAN_VAL) {
            rdm = random.nextGaussian();
        }
        double mid = ((double) min) + (max - min) / 2f;
        double realValue = mid + (rdm * (max - mid) / CUT_OFF_GAUSSIAN_VAL);
        return realValue;
    }
}
