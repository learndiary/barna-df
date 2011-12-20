/*
 * This file is part of the Flux Library.
 *
 * The code of the Flux Library may be freely distributed and modified under the terms of the
 * European Union Public Licence (EUPL) published on the web site <http://www.osor.eu/eupl/european-union-public-licence-eupl-v.1.1>.
 * Copyright for the code is held jointly by the individual authors, who should be listed
 * in @author doc comments. According to Article 5 and Article 11 of the EUPL, publications that
 * include results produced by the Flux Library are liable to reference the Work,
 * see the Flux Library homepage <http://flux.sammeth.net> for more information.
 */

package barna.genome.sequencing.rnaseq.simulation.fragmentation;

import barna.commons.ByteArrayCharSequence;
import barna.genome.sequencing.rnaseq.simulation.distributions.AbstractDistribution;

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
        return "Segregating cDNA ("+(!probDistr ? "Acceptance" : "Rejection")+")";
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
