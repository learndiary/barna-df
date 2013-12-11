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

package barna.flux.capacitor.reconstruction;

public class Kernel {

	public static final byte KERNEL_EPANECHNIKOV= 1;
	

	public static void main(String[] args) {
		getEpanechnikovKernel(10);
	}
	
	public static double[] getEpanechnikovKernel(int w) {
		double[] f= new double[w];
		// http://en.wikipedia.org/wiki/Kernel_(statistics)
		// 3/(2* windowsize) * (1- u^2)
		// u= -1+ (idx/ (w-1))* 2
		double csum= 0;
		for (int i = 0; i < f.length; ++i) {
			double u= ((2* i)/ (double) (w- 1d)) -1d;
			f[i]= (1- (u*u))* 3d/ (2d* w);
			csum+= f[i];
		}
		if (csum!= 1) {
			double csum2= 0;
			for (int i = 0; i < f.length; ++i) {
				f[i]/= csum;
				csum2+= f[i];
			}
			csum= csum2;
			//System.err.println(csum2);
		}
		// must be, sometimes its 0.9999999999...
		//assert(csum== 1);
		
		return f;
	}
	/**
	 * 
	 * @param v half the window size (!!!)
	 * @param F
	 * @return
	 */
	public static double[] getTriweightKernel(int v, double F) {
		double[] f= new double[2* v+ 1];
		// http://en.wikipedia.org/wiki/Kernel_(statistics)
		// 3/(2* windowsize) * (1- u^2)
		// u= -1+ (idx/ (w-1))* 2
		double csum= 0;
		for (int i = -v; i <= v; ++i) {
			double uw= -i/(double) v;
			double k= (1d- (uw* uw));
			double k3= k* k* k;
			f[v+ i]= k3;
			csum+= k3;
		}
		if (csum!= F) {
			double csum2= 0;
			double frac= csum/ F;
			for (int i = 0; i < f.length; ++i) {
				f[i]/= frac;
				csum2+= f[i];
			}
			csum= csum2;
			//System.err.println(csum2);
		}
		// must be, sometimes its 0.9999999999...
		//assert(csum== 1);
		
		return f;
	}
	
	public static double smoothen(byte kernel, int w, double[] b) {
		
		double[] f= getKernel(kernel, w);
		if (f== null)
			return -1;
		
		double[] a= new double[b.length];
		for (int i = 0; i < a.length; i++) 
			a[i]= 0;
		
		for (int i = -w/ 2; i < a.length- (w/ 2); ++i) {
			double med= b[i+ (w/ 2)];
			for (int j = 0; j < w; ++j) {
				if (i+ j< 0|| i+ j>= a.length)
					continue;
				a[i+ j]+= f[j]* med;
			}
		}
		
		double sum= 0;
		for (int i = 0; i < a.length; i++) {
			//b[i]= Math.round(a[i]);
			sum+= a[i];
		}
		
		return sum;
	}
	
	public static double[] getKernel(byte kernel, int w) {
		if (kernel== KERNEL_EPANECHNIKOV)
			return getEpanechnikovKernel(w);
		System.err.println("Not implemented kernel "+ kernel);
		return null;
	}
	
}
