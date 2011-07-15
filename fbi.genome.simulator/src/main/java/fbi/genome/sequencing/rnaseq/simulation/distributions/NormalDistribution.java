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

public class NormalDistribution extends AbstractDistribution {

    double mean = Double.NaN, sd = Double.NaN;
    double sdSquare = Double.NaN;

    public NormalDistribution(double mean) {
        this.mean = mean;
    }

    public NormalDistribution(double mean, double sd) {
        this(mean);
        this.sd = sd;
        this.sdSquare = Math.pow(sd, 2d);
    }

    public double getSd() {
        return sd;
    }

    public double getP(double x) {

        return getP(x, getMean());
    }

    public double getP(double x, double mean) {

        double a1 = Math.pow(x - mean, 2d);
        double p = Math.exp(-a1 / (2d * sdSquare)) / Math.sqrt(2 * Math.PI * sdSquare);

        return p;
    }

    public double getRelFreq(double x) {
        return getRelFreq(x, getMean());
    }

    public double getRelFreq(double x, double mean) {

        return (getP(x, mean) / getMax(mean));
    }

    private double getMax(double mean) {
        return getP(mean, mean);
    }

    public double getMean() {
        return mean;
    }

    @Override
    public String toString() {
        return "Normal-Distribution mean: " + getMean() + " SD: " + getSd();
    }


}
