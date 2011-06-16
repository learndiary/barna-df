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

package fbi.genome.sequencing.rnaseq.simulation.distributions;

/**
 * Generic wrapper for empiric or analytic distributions.
 *
 * @author micha
 */
public abstract class AbstractDistribution {

    public abstract double getP(double x);

    public abstract double getRelFreq(double x);

    public abstract double getP(double x, double mean);

    public abstract double getMean();

    /**
     * The weight of the distribution in composites.
     */
    double weight = Double.NaN;

    /**
     * The mean of the distribution.
     */
    double mean = Double.NaN;

    public double getWeight() {
        return weight;
    }

    public void setWeight(double weight) {
        this.weight = weight;
    }


}
