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
