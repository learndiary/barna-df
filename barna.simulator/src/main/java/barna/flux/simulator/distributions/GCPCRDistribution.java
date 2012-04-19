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

package barna.flux.simulator.distributions;

import barna.commons.log.Log;

/**
 * Create a set of {@link PCRDistribution} instances.
 *
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class GCPCRDistribution {
    /**
     * The number of bins
     */
    private int bins;
    /**
     * The number of generations
     */
    private int generations;
    /**
     * The set of distributions
     */
    private PCRDistribution[] distributions;

    /**
     * Return the number of generations used for this distribution
     *
     * @return generations the number of generations
     */
    public int getGenerations() {
        return generations;
    }

    /**
     * Get the next sample value
     *
     * @param random a random value
     * @param gc_probability the gc probability
     * @return sample the next sample from the PCRDistribution that is linked to the given pc probability
     */
    public double getNext(double random, double gc_probability){
        if(gc_probability < 0 || gc_probability > 1){
            Log.error("PCR/GC probability " + gc_probability );
            gc_probability = Math.max(0.0, Math.min(1.0, gc_probability));
        }
        int bin = (int) (Math.ceil(bins * gc_probability)-1);
        bin = Math.max(0, bin); // in case gc_p is 0
        return distributions[bin].getNext(random);
    }

    /**
     * Print the distribution for the given p to the log
     *
     * @param p the p that selects the bin
     */
    public void printBin(double p){
        if(p < 0 || p > 1){
            Log.error("PCR/GC probability " + p );
            p = Math.max(0.0, Math.min(1.0, p));
        }
        int bin = (int) (Math.ceil(bins * p)-1);
        bin = Math.max(0, bin); // in case gc_p is 0

        for(int i = 0; i< distributions[bin].v.length; i++){
            Log.println(i + " " + distributions[bin].v[i] + "");
        }
    }

    /**
     * Create a new instance
     *
     * @param bins the number of bins
     * @param generations the number of generations
     * @return dist the distribution
     */
    public static GCPCRDistribution create(int bins, int generations){
        GCPCRDistribution dist = new GCPCRDistribution();
        dist.bins = bins;
        dist.generations = generations;

        double splitter = 1.0/(double)bins;
        dist.distributions = new PCRDistribution[bins];

        double current = splitter;
        for (int i = 0; i < bins; i++) {
            //Log.info("PCR", "Creating distribution ("+generations+" generations) for p<" + StringUtils.fprint(current, 2));
            dist.distributions[i] = PCRDistribution.create(generations, current);
            current+=splitter;
        }
        return dist;
    }

}
