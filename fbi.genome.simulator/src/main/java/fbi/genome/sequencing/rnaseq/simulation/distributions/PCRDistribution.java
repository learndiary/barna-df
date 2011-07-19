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

import fbi.commons.Log;
import fbi.commons.StringUtils;
import org.apache.commons.math.special.Gamma;

import java.util.Arrays;

/**
 * PCR amplification distribution
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class PCRDistribution {
    /**
     * Gamma cache for creation
     */
    private static double[] gammaCache;
    /**
     * Number of generations
     */
    private int generations;
    /**
     * Probability that a molecule is amplified
     */
    private double p;
    /**
     * Distribution
     */
    private double[] v;

    public PCRDistribution(final int generations, final double p, final double[] v) {
        this.generations = generations;
        this.p = p;
        this.v = v;
    }

    public int getNext(double random){
        double s = 0;
        for (int i = 0; i < v.length; i++) {
            s += v[i];
            if(s >= random){
                return i;
            }
        }
        return v.length-1;
    }



    private static void initCache(int generations){
        gammaCache = new double[(int) Math.pow(2, generations)];
        Arrays.fill(gammaCache, -1);
    }

    private static double logbico(double n, double k){
        return factln(n) - factln(k) - factln(n-k);
    }

    private static double factln(final double n) {
        if(n >= gammaCache.length) return Gamma.logGamma(n + 1.0);
        if(gammaCache[((int) n)] < 0){
            gammaCache[((int) n)] = Gamma.logGamma(n + 1.0);
        }
        return gammaCache[((int) n)];
    }

    private static double pbico(double n, double k, double p){
        if (n<k) return 0;
        double v = logbico(n,k) + (Math.log(p) * k) + (Math.log(1 - p) * (n - k));
        return Math.exp(v);
    }

    private static double psd(int n, double p, double[] v){
        int n2 = (int) Math.floor(n/2.0);
        double sum = 0;
        for (int k = 0; k <= n2; k++) {
            double ps_1 = 0;
            int nk = (n-k)-1;
            if(nk < v.length){
                ps_1 = v[nk];
            }
            double tmp = pbico(n-k, k, p);
            double np = ps_1 * tmp;
            sum += np;
        }
        return sum;
    }

    public static PCRDistribution create(int generations, double p){

        Log.progressStart("Creating PCR Distributions (p=" + p + ", generations=" + generations + ")");
        initCache(generations);
        double[] v = null;
        for(int j=1; j<= generations; j++ ){
            Log.progress(j, generations);
            if(j == 1){
                // initial probs
                v = new double[2];
                v[0] = 1.0-p;
                v[1] = p;
            }else{
                double[] nv = new double[(int) Math.pow(2,j)];
                for (int i = 0; i < nv.length; i++) {
                    nv[i] = psd(i+1, p, v);
                }
                v = nv;
            }
        }
        Log.progressFinish(StringUtils.OK, true);
        // reset cache
        gammaCache = null;

        return new PCRDistribution(generations, p, v);
    }

    public static void main(String[] args) {
        double p = 0.2;
        int g = 15;
        PCRDistribution pcrDistribution = create(g, p);

        for(int i=0; i<pcrDistribution.v.length;i++){
            System.out.println(i + " " + pcrDistribution.v[i]);
        }

    }
}
