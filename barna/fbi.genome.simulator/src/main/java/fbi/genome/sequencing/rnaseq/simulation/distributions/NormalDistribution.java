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

    double mean = Double.NaN;
    double sd = Double.NaN;
    private boolean pdf;
    double sdSquare = Double.NaN;

    public NormalDistribution(double mean) {
        this(mean, Double.NaN);
    }

    public NormalDistribution(double mean, double sd) {
        this(mean, sd, false);
    }

    public NormalDistribution(double mean, double sd, boolean pdf) {
        this.mean = mean;
        this.sd = sd;
        this.pdf = pdf;
        this.sdSquare = Math.pow(sd, 2d);
    }

    public double getSd() {
        return sd;
    }

    public double getP(double x) {
        if(pdf) return pdf(x, getMean());
        return getP(x, getMean());
    }

    public double getP(double x, double mean) {
        //if(pdf) return pdf(x, mean);
        if(pdf) return Phi(x);
        double a1 = Math.pow(x - mean, 2d);
        double p = Math.exp(-a1 / (2d * sdSquare)) / Math.sqrt(2 * Math.PI * sdSquare);

        return p;
    }

 // Compute z such that Phi(z) = y via bisection search
    public static double PhiInverse(double y) {
        return PhiInverse(y, .00000001, -8, 8);
    }

    // bisection search
    private static double PhiInverse(double y, double delta, double lo, double hi) {
        double mid = lo + (hi - lo) / 2;
        if (hi - lo < delta) return mid;
        if (Phi(mid) > y) return PhiInverse(y, delta, lo, mid);
        else              return PhiInverse(y, delta, mid, hi);
    }

    // return phi(x) = standard Gaussian pdf
    public static double phi(double x) {
        return Math.exp(-x*x / 2) / Math.sqrt(2 * Math.PI);
    }

    // return phi(x, mu, signma) = Gaussian pdf with mean mu and stddev sigma
    public static double phi(double x, double mu, double sigma) {
        return phi((x - mu) / sigma) / sigma;
    }

    // return Phi(z) = standard Gaussian cdf using Taylor approximation
    public static double Phi(double z) {
        if (z < -8.0) return 0.0;
        if (z >  8.0) return 1.0;
        double sum = 0.0, term = z;
        for (int i = 3; sum + term != sum; i += 2) {
            sum  = sum + term;
            term = term * z * z / i;
        }
        return 0.5 + sum * phi(z);
    }


    public double pdf(double x){
        return pdf(x, getMean());
    }

    public double pdf(double x, double mean){
        double a = 1.0 / (Math.sqrt(2.0 * Math.PI * sdSquare));
        double b = Math.exp(-(Math.pow((x-mean),2))/(2.0*sdSquare));
        return a*b;
    }

    public double getRelFreq(double x) {
        return getRelFreq(x, getMean());
    }

    public double getRelFreq(double x, double mean) {

        return (getP(x, mean) / getMax(mean));
    }

    public double getMax(double mean) {
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
