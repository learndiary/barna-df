package fbi.genome.model.commons;

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
}
