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

package barna.model.commons;

import java.util.Arrays;

public class Coverage {

    /**
     * Array for coverage
     */
    int[] coverage= null; 

    /**
     * Length of the actual array (underlying data array can be longer)
     */
    int length= -1;
    
    /**
     * Initializes underlying array, re-use if already instantiated.
     * @param length the length that has to be provided by the unerlying array
     */
    public Coverage(int length) {
    	reset(length);
	}

    /**
     * Resets/resizes underlying array to meet the required length.
     * @param length the length that is required for the coverage profile
     */
    public void reset(int length) {
		this.length= length;
        if (coverage== null|| coverage.length< length) 
        	coverage= new int[length];
        else
        	Arrays.fill(coverage, 0, length, 0);
    }
    
    /**
     * Increments the array at the given position
     * @param pos position to be incremented by <code>1</code>
     */
    public void increment(int pos) {
    	if (pos< 0 || pos>= length)
    		//throw new RuntimeException(pos+ " out of range [0;"+ length+ "[");
    		return; // lazily allow under- and overflows
    	++coverage[pos];
    }

    /**
     * Computes the mean coverage, with or without considering 0-values. 
     * @param excludeZero flag whether positions with 0-counts 
     * should be excluded
     * @return
     */
    public double getMean(boolean excludeZero) {
        // exclude 0-positions
        double avgCov= 0d;
        int cnt= 0;
        for (int i = 0; i < length; i++) {
        	if (excludeZero&& coverage[i]== 0)
        		continue;
        	++cnt;
			avgCov+= coverage[i];
        }
        avgCov/= cnt;

        return avgCov;
    }
    
    /**
     * Computes the chi-square metrics of the coverage profile.
     * @param excludeZero flag whether positions with 0-counts 
     * should be excluded
     * @return the chi-square value of the current profile
     */
    public long getChiSquare(boolean excludeZero) {
    	
        double avgCov= getMean(excludeZero);
        double x2= 0d;
        for (int i = 0; i < length; i++) { 
        	if (excludeZero&& coverage[i]== 0)
        		continue;
			x2+= (coverage[i]- avgCov)* (coverage[i]- avgCov);
        }
		x2/= avgCov;
		
		return Math.round(x2);
    }

    /**
     * Computes the coefficient of variation (CV) of the current profile.
     * Uses Anscombe's residuals to transform assumed Poisson-distribution
     * to something close to normally distributed.
     * @param excludeZero flag whether positions with 0-counts 
     * should be excluded
     * @return the coefficient of variation computed as described
     */
    public double getCV(boolean excludeZero) {
    	double avgCov= getMean(excludeZero);
		double mean= 0, min= Double.MAX_VALUE;
		int cnt= 0;
		for (int i = 0; i < length; i++) {
			if (excludeZero&& coverage[i]== 0)
				continue;
			++cnt;
			double a= coverage[i];
			// Anscombe residuals [Hansen et al. 2010]
			a= (3d/ 2d)* (Math.pow(coverage[i], 2d/3d)- Math.pow(avgCov, 2d/3d))/ Math.pow(avgCov, 1d/6d);
			mean+= a;
			if (a< min)
				min= a;
		}
		mean/= cnt;
		mean+= 2* Math.abs(min);
		double cv= 0;
		for (int i = 0; i < length; i++) {
			if (excludeZero&& coverage[i]== 0)
				continue;
			double a= coverage[i];
			a= (3d/ 2d)* (Math.pow(coverage[i], 2d/3d)- Math.pow(avgCov, 2d/3d))/ Math.pow(avgCov, 1d/6d);
			a+= 2* Math.abs(min);
			cv+= (a- mean)* (a- mean);
		}
		cv/= cnt;
		cv= Math.sqrt(cv);	// sdev
		cv/= mean;
    	
		return cv;
    }

    /**
     * Computes the fraction of transcript that is covered.
     * @return proportion of transcript covered
     */
	public float getFractionCovered() {
		int cnt= 0;
		for (int i = 0; i < length; i++) 
			if (coverage[i]> 0)
				++cnt;
		float f= cnt/ (float) length;
		return f;
	}
}
