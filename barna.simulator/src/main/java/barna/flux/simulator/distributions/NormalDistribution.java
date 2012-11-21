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

public class NormalDistribution extends AbstractDistribution {

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

    /**
     * Convert this to an empirical distribution
     *
     * @param bins number of bins for the empirical distribution
     * @return emp the empirical distribution
     */
    public EmpiricalDistribution toEmpirical(int bins){
        double min = 0;
        double max = 3.0 * getMean();
        double[] a = new double[bins+2];
        for(int i = 1; i<a.length-1; i++){
            double x = (i / (double) (a.length-2)) * max;
            a[i]  = getP(x);
        }
        EmpiricalDistribution empDist = new EmpiricalDistribution(a, min, max, getMean());
        empDist.sum = 1.0;
        return empDist;

    }


    @Override
    public String toString() {
        return "Normal-Distribution mean: " + getMean() + " SD: " + getSd();
    }


}
