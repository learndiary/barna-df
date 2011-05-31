package fbi.genome.sequencing.rnaseq.simulation.fragmentation;

import fbi.commons.ByteArrayCharSequence;
import fbi.genome.sequencing.rnaseq.simulation.distributions.AbstractDistribution;

import java.util.Arrays;
import java.util.List;
import java.util.Random;

/**
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class FragmentFilterRejection implements FragmentProcessor {
    private AbstractDistribution[] d;
    private boolean probDistr;
    private Random rndGel = new Random();

    public FragmentFilterRejection(AbstractDistribution[] d, boolean probDistr) {
        this.d = d;
        this.probDistr = probDistr;
    }

    @Override
    public List<Fragment> process(final ByteArrayCharSequence id, final ByteArrayCharSequence cs, final int start, final int end, final int len) {
        // get (possibly cumulative) probability for length being in result
        double plen = 0d;
        for (int i = 0; i < d.length; i++) {
            double p = (probDistr ? d[i].getP(len) : d[i].getRelFreq(len));
            plen += d[i].getWeight() * p;
        }

        // Bernoulli trial
        if (plen > 1 || rndGel.nextDouble() < plen) {
            return Arrays.asList(new Fragment(id, start, end));
        } else {
            return null;
        }
    }


    @Override
    public String getName() {
        return "Segregating cDNA (Rejection Filter)";
    }

    @Override
    public String getConfiguration() {
        return null;
    }

    @Override
    public String done() {
        return null;
    }

}
