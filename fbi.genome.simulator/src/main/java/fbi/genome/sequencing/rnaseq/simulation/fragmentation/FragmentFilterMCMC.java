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

package fbi.genome.sequencing.rnaseq.simulation.fragmentation;

import fbi.commons.ByteArrayCharSequence;
import fbi.genome.sequencing.rnaseq.simulation.distributions.AbstractDistribution;

import java.util.Arrays;
import java.util.List;
import java.util.Random;

/**
 * see <code>http://personal.strath.ac.uk/gary.koop/extra_material_on_metropolis.pdf</code>,
 * <code>http://www.ps.uci.edu/~markm/numerical_methods/Metropolis%96Hastings%20algorithm%20-%20Wikipedia,%20the%20free%20encyclopedia.pdf</code>
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class FragmentFilterMCMC implements FragmentProcessor {
    /**
     * Last length
     */
    private int lastLen = -1;
    /**
     * Last propability
     */
    private double lastP = -1;
    /**
     * Random generator
     */
    private Random rndGel = new Random();
    private AbstractDistribution dGenerate;
    private AbstractDistribution[] dProposal;


    public FragmentFilterMCMC(final AbstractDistribution dGenerate, final AbstractDistribution[] dProposal) {
        this.dGenerate = dGenerate;
        this.dProposal = dProposal;
    }

    @Override
    public List<Fragment> process(final ByteArrayCharSequence id, final ByteArrayCharSequence cs, final int start, final int end, final int len) {
        // first value always accepted to init algorithm (but not output)
        double p = 0d;
        for (int i = 0; i < dProposal.length; i++) {
            p += dProposal[i].getP(len);
        }
        if (lastLen < 0) {
            lastLen = len;
            lastP = p;
            return null;
        }

        // Metropolis/Hastings/Ema
        double a1 = p / lastP;
        double a2 = dGenerate.getP(lastLen, len) / dGenerate.getP(len, lastLen);
        double a = a1 * a2;

        // accept
        if (a >= 1 || rndGel.nextDouble() <= a) {
            lastLen = len;
            lastP = p;
            return Arrays.asList(new Fragment(id, start, end));
        } else {
            return null;
        }
    }

    @Override
    public String getName() {
        return "Segregating cDNA (MCMC Filter)";
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
