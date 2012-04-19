/*
 * Copyright (c) 2010, Micha Sammeth
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * The names of its contributors may be not used to endorse or promote
 *       products derived from this software without specific prior written
 *       permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL MICHA SAMMETH BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

package barna.flux.simulator.fragmentation;

import barna.commons.ByteArrayCharSequence;
import barna.flux.simulator.distributions.AbstractDistribution;

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

        // Metropolis/Hastings
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
