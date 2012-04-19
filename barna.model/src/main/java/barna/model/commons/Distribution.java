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
import java.util.HashMap;

public class Distribution {

	double[] arrayD= null;	// sorted
	int[] arrayI= null;	// sorted
	long[] arrayL= null;
	HashMap histogram= null;
	       
	public Distribution(long[] surface) {
		Arrays.sort(surface);
		this.arrayL= surface;
	}
	
	public Distribution(double[] surface) {
		Arrays.sort(surface);
		this.arrayD= surface;
	}
	
	public Distribution(int[] surface) {
		Arrays.sort(surface);
		this.arrayI= surface;
	}
	
	public double getMin() {
		if (arrayD!= null) {
			double min= java.lang.Double.MAX_VALUE;
			for (int i = 0; i < arrayD.length; i++) 
				if (arrayD[i]< min)
					min= arrayD[i];
			return min;
		} else if (arrayI!= null) {
			int min= java.lang.Integer.MAX_VALUE;
			for (int i = 0; i < arrayI.length; i++) 
				if (arrayI[i]< min)
					min= arrayI[i];
			return min;
		} else if (arrayL!= null) {
			long min= java.lang.Long.MAX_VALUE;
			for (int i = 0; i < arrayL.length; i++) 
				if (arrayL[i]< min)
					min= arrayL[i];
			return min;
		}	
		return 0d;
	}
	
	public double getMax() {
		if (arrayD!= null) {
			double max= java.lang.Double.MIN_VALUE;
			for (int i = 0; i < arrayD.length; i++) 
				if (arrayD[i]> max)
					max= arrayD[i];
			return max;
		} else if (arrayI!= null) {
			int max= java.lang.Integer.MIN_VALUE;
			for (int i = 0; i < arrayI.length; i++) 
				if (arrayI[i]> max)
					max= arrayI[i];
			return max;
		} else if (arrayL!= null) {
			long max= java.lang.Long.MIN_VALUE;
			for (int i = 0; i < arrayL.length; i++) 
				if (arrayL[i]> max)
					max= arrayL[i];
			return (double) max;
		}
		return 0d;
	}
	
	public HashMap getHistogram() {
		if (histogram == null) {
			histogram= new HashMap();
			for (int i = 0; arrayD!= null&& i < arrayD.length; i++) {
				Double key= new Double(arrayD[i]);
				Integer val= (Integer) histogram.get(key);
				if (val== null) 
					val= new Integer(1);
				else
					val= new Integer(val.intValue()+ 1);
				histogram.put(key, val);
			}
			
			for (int i = 0; arrayI!= null&& i < arrayI.length; i++) {
				Integer key= new Integer(arrayI[i]);
				Integer val= (Integer) histogram.get(key);
				if (val== null) 
					val= new Integer(1);
				else
					val= new Integer(val.intValue()+ 1);
				histogram.put(key, val);
			}
		}

		return histogram;
	}
	
	public double getMean() {
		
		if (arrayD!= null) {
			if (arrayD.length== 0)
				return 0d;
			double sum= 0d;
			for (int i = 0; i < arrayD.length; i++) 
				sum+= arrayD[i];
			return (sum/ arrayD.length);
		} else if (arrayI!= null) {
			if (arrayI.length== 0)
				return 0d;
			double sum= 0d;
			for (int i = 0; i < arrayI.length; i++) 
				sum+= arrayI[i];
			return (sum/ arrayI.length);
		}
		
		return 0d;
	}
	
	public double getTotal() {
		if (arrayD!= null) {
			double sum= 0;
			for (int i = 0; i < arrayD.length; i++) 
				sum+= arrayD[i];
			return sum;
		} else if (arrayI!= null) {
			int sum= 0;
			for (int i = 0; i < arrayI.length; i++) 
				sum+= arrayI[i];
			return sum;
		}
		
		return 0d;
	}
	
	double getMedian(int j, int i) {
		if (arrayD!= null) {
			if (arrayD.length== 0)
				return 0d;
			int medPos= j* arrayD.length/ i;
			if ((j*arrayD.length)% i== 0&& (medPos+1)< arrayD.length)
				return ((arrayD[medPos]+ arrayD[medPos+ 1])/ 2d);
			else
				return arrayD[medPos];
		} else if (arrayI!= null) {
			if (arrayI.length== 0)
				return 0d;
			int medPos= j* arrayI.length/ i;
			if ((j* arrayI.length)% i== 0&& (medPos+1)< arrayI.length)
				return ((arrayI[medPos]+ arrayI[medPos+ 1])/ 2d);
			else
				return arrayI[medPos];
		} else if (arrayL!= null) {
			if (arrayL.length== 0)
				return 0d;
			int medPos= j* arrayL.length/ i;
			if ((j* arrayL.length)% i== 0&& (medPos+1)< arrayL.length)
				return ((arrayL[medPos]+ arrayL[medPos+ 1])/ 2d);
			else
				return arrayL[medPos];
		}
		return 0d;
	}
	
	public double getMedian() {
		return getMedian(1,2);
	}
	public double get1stQuart() {
		return getMedian(1,4);
	}
	public double get3rdQuart() {
		return getMedian(3,4);
	}

	public double getSum() {		
		if (arrayD!= null) {
			double sum= 0d;
			for (int i = 0; i < arrayD.length; i++) 
				sum+= arrayD[i];
			return sum;
		} else if (arrayI!= null) {
			int sum= 0;
			for (int i = 0; i < arrayI.length; i++) 
				sum+= arrayI[i];
			return sum;
		}
		return 0d;
	}
	
	public double getStandardDeviation() {
		if (arrayD!= null) {
			if (arrayD.length== 0)
				return 0d;
			double med= getMean();		
			double sum= 0d;
			for (int i = 0; i < arrayD.length; i++) {
				double val= arrayD[i]- med;
				val*= val;
				sum+= val;
			}
			
			sum/= (arrayD.length- 1);
			return java.lang.Math.sqrt(sum);
			
			
		} else if (arrayI!= null) {
			if (arrayI.length== 0)
				return 0d;
			double med= getMean();		
			double sum= 0d;
			for (int i = 0; i < arrayI.length; i++) {
				double val= arrayI[i]- med;
				val*= val;
				sum+= val;
			}
			
			sum/= (arrayI.length- 1);	// ok, mathworld s_N-1 (sample corrected)
			return java.lang.Math.sqrt(sum);
		}
		
		return 0d;
	}
	
	/**
	 * 
	 * @return mean, median, stdDev
	 */
	public String[] toStatString() {
		String[] statStr= new String[3];
		statStr[0]= new Double(getMean()).toString();
		statStr[1]= new Double(getMedian()).toString();
		statStr[2]= new Double(getStandardDeviation()).toString();
		for (int i = 0; i < statStr.length; i++) 
			statStr[i]= statStr[i].substring(0, statStr[i].indexOf('.')+ 2);

		return statStr;
	}
}
