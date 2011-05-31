package fbi.genome.sequencing.rnaseq.simulation.fragmentation;

import fbi.commons.ByteArrayCharSequence;
import fbi.genome.sequencing.rnaseq.simulation.PWM;

import java.io.File;
import java.util.*;

/**
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class FragmentEnzymatic implements FragmentProcessor{
    private Map<CharSequence, double[]> mapWeightSense = null;
    private Map<CharSequence, double[]> mapWeightAsense = null;
    private PWM pwmSense;
    private PWM pwmAsense;
    Random rnd1 = new Random();
    Random rnd2 = new Random();
    Random rnd3 = new Random();
    private int leftFlank;
    private int rightFlank;


    public FragmentEnzymatic(File ezMotif, Map<CharSequence, CharSequence> mapTxSeq, int leftFlank, final int rightFlank) throws Exception {
        this.leftFlank = leftFlank;
        this.rightFlank = rightFlank;
        if(ezMotif == null) throw new NullPointerException("Null Motif file for enzymatic fragmentation not permitted !");
        pwmSense = PWM.create(ezMotif);
        mapWeightSense = Fragmenter.getMapWeight(mapTxSeq, mapTxSeq, pwmSense);
        pwmSense.invert();
        pwmAsense = pwmSense;
        mapWeightAsense = Fragmenter.getMapWeight(mapTxSeq, mapTxSeq, pwmAsense);
    }

    @Override
    public List<Fragment> process(final ByteArrayCharSequence id, final ByteArrayCharSequence cs, final int start, final int end, final int len) {
        assert (mapWeightSense != null && mapWeightAsense != null);
        double[] wsense = mapWeightSense.get(id), wasense = mapWeightAsense.get(id);

        // short-hand value to get number of enzyme attacks as function of length
        // could be a parameter
        int n = (int) Math.sqrt(len);
        int k = 0;
        int[] pos = new int[n];

        // get breakpoints
        double norm = Fragmenter.toCDF(wsense, start, end, leftFlank, rightFlank);
        for (int i = 0; i < n; i++) {

            // localize potential breakpoint
            boolean sense = rnd1.nextFloat() < 0.5;
            double[] a = (sense ? wsense : wasense);
            int bp = Arrays.binarySearch(a, start, end + 1, rnd2.nextDouble());
            bp = (bp >= 0 ? bp : -(bp + 1));
            if (bp >= end) {
                continue;
            }

            // lookup in array with bp positions
            double pb = (bp == 0 ? a[bp] : a[bp] - a[bp - 1]);    // breaking probability
            // ok, the pwm recognizes the base where it is cut
            // but has no concept about 5'/3'. If cutting always
            // consistently 5' of 0-position, the bp in sense
            // marks the start of the next fragment, when motif
            // is recognized in anti-sense, we got to add 1.
            // (e.g., palindrome ACGT cut at G is found at
            // position 1166 in sense, and position 1165 in
            // anti-sense)
            bp = (sense ? bp : bp + 1);    // set bp consistently to 1st nt of next fragment
            int idx = Arrays.binarySearch(pos, 0, k, bp);
            if (idx >= 0) {
                continue;
            } else {
                idx = -(idx + 1);
            }

            // breaking probability
            pb /= (sense ? pwmSense.getMaximumP() : pwmSense.getMaximumP()); // should be the same if only transposed
            if (rnd3.nextDouble() <= pb) {
                System.arraycopy(pos, idx, pos, idx + 1, (k - idx) + 1);
                pos[idx] = bp;
                ++k;
            }
        }
        Fragmenter.toPDF(wsense, start, end, norm, leftFlank, rightFlank);

        // break
        int last = start;
        List<Fragment> fragments = new ArrayList<Fragment>();
        for (int i = 0; i < k; i++) {
            int now = start + pos[i] - 1;    // end
            int nuLen = now - last;
            if (nuLen == 0) {
                continue;
            }
            cs.replace(0, last);
            cs.replace(1, now);
            fragments.add(new Fragment(id, last, now));
            //fragments.add(cs.toString());
            //rw.writeLine(cs, fos);
            last = now + 1;    // start
        }
        int now = end;
        int nuLen = now - last + 1;
        //updateMedian(nuLen);

        cs.replace(0, last);
        cs.replace(1, now);
        //rw.writeLine(cs, fos);
        //fragments.add(cs.toString());
        fragments.add(new Fragment(id, last, now));
        return fragments;
    }

    @Override
    public String getName() {
        return "Enzymatic Digestion";
    }

    @Override
    public String getConfiguration() {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public String done() {
        return null;
    }
}
