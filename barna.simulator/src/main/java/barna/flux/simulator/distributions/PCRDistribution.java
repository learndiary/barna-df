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
import barna.commons.utils.StringUtils;
import org.apache.commons.math.special.Gamma;

import java.util.Arrays;

/**
 * PCR amplification distribution for a given duplication probability
 * and a given number of generations.
 * <p>
 *     To compute a new distribution, use the static {@link #create(int, double)} method.
 * </p>
 *
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class PCRDistribution {
    /**
     * Gamma cache for creation
     */
    private double[] gammaCache;
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
    double[] v;



    /**
     * Creaet a new instance with the given number of generations, a given duplication
     * probability and the probability distribution
     *
     * @param generations the number of generations
     * @param p the duplication probability
     * @param v the probability distribution
     */
    public PCRDistribution(final int generations, final double p, final double[] v) {
        this.generations = generations;
        this.p = p;
        this.v = v;
    }

    /**
     * INTERNAL: create an empty distributtion
     *
     * @param generations
     * @param p
     */
    protected PCRDistribution(final int generations, final double p) {
        this(generations, p, null);
    }

    /**
     * Get the next sample from the probability distribution
     *
     * @param random a random value {@code 0<=r<=1.0}
     * @return sample a sample from the distribution
     */
    public int getNext(double random){
        double s = 0;
        int zeroIndex = -1;
        for (int i = 0; i < v.length; i++) {
            s += v[i];
            if(s >= random){
                return i;
            }
            if(v[i] == 0 && zeroIndex < 0) zeroIndex = i;
            else if(v[i] > 0) zeroIndex = -1;
        }
        return zeroIndex > 0 ? zeroIndex-1 : v.length-1;
    }


    /**
     * Helper method to initialize the gamma cache for a given number of generations
     *
     * @param generations the number of generations
     */
    private void initCache(int generations){
        gammaCache = new double[(int) Math.pow(2, generations)];
        Arrays.fill(gammaCache, -1);
    }

    private double logbico(double n, double k){
        return factln(n) - factln(k) - factln(n-k);
    }

    private double factln(final double n) {
        if(n >= gammaCache.length) return Gamma.logGamma(n + 1.0);
        if(gammaCache[((int) n)] < 0){
            gammaCache[((int) n)] = Gamma.logGamma(n + 1.0);
        }
        return gammaCache[((int) n)];
    }

    private double pbico(double n, double k, double p){
        if (n<k) return 0;
        double v = logbico(n,k) + (Math.log(p) * k) + (Math.log(1 - p) * (n - k));
        return Math.exp(v);
    }


    private double psd(int n, double p, double[] v){
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

    /**
     * Compute a new distribution. Note that the running time grows
     * exponentially with the number of generations
     *
     * @param generations the number of generations
     * @param p the duplication probability
     * @return distribution the distribution
     */
    public static PCRDistribution create(int generations, double p){
        Log.progressStart("Creating PCR Distributions (p=" + StringUtils.fprint(p, 2) + ", generations=" + generations + ")");
        PCRDistribution dist = new PCRDistribution(generations, p);
        dist.initCache(generations);
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
                    nv[i] = dist.psd(i+1, p, v);
                }
                v = nv;
            }
        }
        Log.progressFinish(StringUtils.OK, true);
        dist.v = v;
        // reset cache
        dist.gammaCache = null;
        return dist;
    }
}
