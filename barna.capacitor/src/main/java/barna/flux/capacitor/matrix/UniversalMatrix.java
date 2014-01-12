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

package barna.flux.capacitor.matrix;

import barna.commons.log.Log;
import barna.commons.system.OSChecker;
import barna.flux.capacitor.reconstruction.Kernel;
import barna.model.commons.IntVector;
import barna.model.constants.Constants;

import java.util.Arrays;
import java.util.Vector;

public class UniversalMatrix {

	static Vector<IntVector> poolIntVector= null;
	public static void addIntVector(IntVector v) {
		
		if (1== 1)
			return;
		
		if (poolIntVector== null)
			poolIntVector= new Vector<IntVector>();
		if (v.size()< 3000)
			poolIntVector.add(v); 
	}
	static Object lock= new Object();
	public static IntVector getIntVector(int size) {
		
		if (1== 1)
			return null;
		
//		synchronized(lock) {
//			for (int i = 0; poolIntVector!= null&& i < poolIntVector.size(); i++) {
//				if (poolIntVector.elementAt(i).capacity()>= size) {
//					IntVector v= poolIntVector.remove(i);
//					v.reset();
//					return v;
//				}
//			}
//		}
		
		// slow
		IntVector v= null;
		try {
			v= poolIntVector.lastElement();
		} catch (Exception e) {
			return null;
		}
		
		return v;
	}
	
	public double[] sense, asense;
	public double sums;
    public double suma;
    int nrTranscripts;

    /**
     * Constructor initializes values
     * @param length the size of the profile
     */
	public UniversalMatrix(int length) {
		sense= new double[length];
		asense= new double[length];
        // prevent from empty slots
        Arrays.fill(sense, 1);
        Arrays.fill(asense, 1);
		sums= sense.length;
		suma= asense.length;
        nrTranscripts= 0;
	}

    /**
     * Adds a transcript observation to the profile.
     * @param a expression normalized observation, sum(a[i])= 1
     * @param dir directionality
     */
    public void add(double[] a, byte dir) {

        if (a== null)
            throw new IllegalArgumentException("Array to be added cannot be null!");
        if (a.length!= sense.length)
            throw new IllegalArgumentException("Array to be added does not correspond to profile length: "
                    + a.length+ "<>"+ sense.length
            );

        for (int i = 0; i < a.length; i++) {
            if (dir== Constants.DIR_FORWARD) {
                sense[i]+= a[i];
                sums+= a[i];
            } else {
                asense[i]+= a[i];
                suma+= a[i];
            }
        }
    }


	public double get(int p1, int p2, int tlen, byte dir) {
		// 090820 <p2 
		// exclusive regions for back-normalization needed
		// otherwise fracs> transcriptcount for gene (and reads also)
		if (p1== Integer.MIN_VALUE|| p2== Integer.MIN_VALUE|| p1> tlen|| p2> tlen) {
            Log.warn("Read length is too small or exceed transcript length. " + OSChecker.NEW_LINE + "Lengths: Read 1 - " + p1 +
                    " Read 2 - " + p2 + " Transcript - " + tlen);
			return 0;
		}

        // acceptor site= donor site+ 1
        // example: chr9:139685807-139735639W (pos: 139,702,191-2)
        if (p1> p2)
            return 0;

        // p1 and p2 are included
        // transform to:
        // included p1 -> 0..(tlen-1)
        // excluded p2 -> 1..tlen
        double  rPos1= ((p1/ (double) (tlen- 1))* (sense.length- 1)),
                rPos2= (((p2+ 1)/ (double) (tlen- 1))* (sense.length- 1));  // allow last cell shared in "down-casting"

        // prevent overflow, happens during "up-casting" (tx shorter than profile)
        if (rPos2>= sense.length)
            rPos2= sense.length;

        // shared matrix entries
        int iPos1= (int) rPos1, iPos2= (int) rPos2;
        double sum= 0, f= 1;
        if (iPos1== iPos2) {
            f= (rPos2- rPos1);
            if (dir== Constants.DIR_FORWARD|| dir== Constants.DIR_BOTH)
                sum+= f* sense[iPos1];  // dont return here, e.g. rPos1== rPos2== 0
            if (dir== Constants.DIR_BACKWARD|| dir== Constants.DIR_BOTH)
                sum+= f* asense[iPos1];
            if (f!= 0)
                return sum;
        } else {
            if (iPos1< rPos1) {
                f=  (1d- (rPos1- iPos1));
                if (dir== Constants.DIR_FORWARD|| dir== Constants.DIR_BOTH)
                    sum+= f* sense[iPos1];
                if (dir== Constants.DIR_BACKWARD|| dir== Constants.DIR_BOTH)
                    sum+= f* asense[iPos1];
                ++iPos1;    // only increment if there is a float rest
            }
            if (iPos2< rPos2) {
                f=  (rPos2- iPos2);
                if (dir== Constants.DIR_FORWARD|| dir== Constants.DIR_BOTH)
                    sum+= f* sense[iPos2];
                if (dir== Constants.DIR_BACKWARD|| dir== Constants.DIR_BOTH)
                    sum+= f* asense[iPos2];
            }

            --iPos2;    // decrement always, as the position will be included in the next edge
        }

        // intermediate complete matrix entries
        if (dir== Constants.DIR_FORWARD|| dir== Constants.DIR_BOTH) {
            for (int i = iPos1; i <= iPos2; ++i)
                sum+= sense[i];
        }
        if (dir== Constants.DIR_BACKWARD|| dir== Constants.DIR_BOTH) {
            for (int i = iPos1; i <= iPos2; ++i)
                sum+= asense[i];
        }

        return sum;
	}

    public void smooth(int window) {
        sums= Kernel.smoothen(Kernel.KERNEL_EPANECHNIKOV, window, sense);
        suma= Kernel.smoothen(Kernel.KERNEL_EPANECHNIKOV, window, asense);
    }

	public double getFrac(int p1, int p2, int tlen, byte dir) {
		double sum= get(p1, p2, tlen, dir);
		
		long v= 0;
		if (dir== Constants.DIR_FORWARD|| dir== Constants.DIR_BOTH)
			v+= sums;
		if (dir== Constants.DIR_BACKWARD|| dir== Constants.DIR_BOTH)
			v+= suma;
		double frac= sum/ v;
        return frac;
	}

	public int getLength() {
		assert(sense.length== asense.length);
		return sense.length;
	}

	public double getSum(byte dir) {
		if (dir== Constants.DIR_FORWARD)
			return sums;
		else if (dir== Constants.DIR_BACKWARD)
			return suma;
		else if (dir== Constants.DIR_BOTH)
			return sums+ suma;
		System.err.println("[ASSERT] no direction "+ dir);
		return 0;
	}

	@Override
	public String toString() {
		return toStringBuilder().toString();
	}
	
	public StringBuilder toStringBuilder() {
		
		StringBuilder sb= new StringBuilder(sense.length* 3); 
		
		// sense
		for (int i = 0; i < sense.length; i++) {
			sb.append(Double.toString(sense[i]));
			sb.append(",");
		}
		sb.deleteCharAt(sb.length()- 1);
		sb.append(barna.commons.system.OSChecker.NEW_LINE);

		// anti-sense
		for (int i = 0; i < asense.length; i++) {
			sb.append(Double.toString(asense[i]));
			sb.append(",");
		}
		sb.deleteCharAt(sb.length()- 1);
		sb.append(barna.commons.system.OSChecker.NEW_LINE);

		// sum
		for (int i = 0; i < sense.length; i++) {
			sb.append(Double.toString(sense[i]+ asense[i]));
			sb.append(",");
		}
		sb.deleteCharAt(sb.length()- 1);
		sb.append(barna.commons.system.OSChecker.NEW_LINE);

		return sb;
	}
	
	public StringBuilder toStringBuilder(int nrBins) {
		
		if (nrBins> sense.length)
			throw new RuntimeException("Not implemented.");
		
		StringBuilder sb= new StringBuilder(nrBins* 3); 
		float div= sense.length/ (float) nrBins;	// assuming nrBins< 
		
		int[] sBins= new int[nrBins], aBins= new int[nrBins];
		for (int i = 0; i < aBins.length; i++) {
			sBins[i]= 0; 
			aBins[i]= 0;
		}
		for (int i = 0; i < sense.length; i++) {
			int p= (int) (i/ div);
			sBins[p]+= sense[i];
			aBins[p]+= asense[i];
		}
		
		// sense
		for (int i = 0; i < sBins.length; i++) {
			sb.append(Integer.toString(sBins[i]));
			sb.append(",");
		}
		sb.deleteCharAt(sb.length()- 1);
		sb.append(barna.commons.system.OSChecker.NEW_LINE);

		// anti-sense
		for (int i = 0; i < aBins.length; i++) {
			sb.append(Integer.toString(aBins[i]));
			sb.append(",");
		}
		sb.deleteCharAt(sb.length()- 1);
		sb.append(barna.commons.system.OSChecker.NEW_LINE);

		// sum
		for (int i = 0; i < sBins.length; i++) {
			sb.append(Integer.toString(sBins[i]+ aBins[i]));
			sb.append(",");
		}
		sb.deleteCharAt(sb.length()- 1);
		sb.append(barna.commons.system.OSChecker.NEW_LINE);

		return sb;
	}
	
	public byte[] toByteArray() {
		StringBuilder sb= toStringBuilder();
		byte[] b= new byte[sb.length()];
		for (int i = 0; i < b.length; i++) 
			b[i]= (byte) sb.charAt(i);
		
		return b;
	}

	public int fill() {
		int sum= 0;
		for (int i = 0; i < sense.length; i++) {
			++sense[i];
			++sums;
			++asense[i];
			++suma;
			sum+= 2;
		}
		return sum;
	}

	double nfactor= Double.NaN;
	public double getNfactor() {
		return getNfactor(0.2d);
	}
	public double getNfactor(double convFactor) {
		if (Double.isNaN(nfactor)) {
			
			// avg coverage where != 0
/*			long sum= 0;
			int nr= 0;
			for (int i = 0; i < sense.length; i++) {
				int val= sense[i]+ asense[i];
				if (val== 0)
					continue;
				sum+= val;
				++nr;
			}
			double avg= sum/ (double) nr;
*/
			// median
			//IntVector v= new IntVector(sense.length/ 2);
/*			int[] vals= new int[sense.length];
			for (int i = 0; i < sense.length; i++) 
				vals[i]= sense[i]+ asense[i];
			Arrays.sort(vals);
			double med= vals[vals.length/ 2];
			if (med== 0)
				nfactor= 1;
			else {
				double filled= med* sense.length;
				double observed= sums+ suma;
				nfactor= filled/ observed;
			}
*/
			
			// convergence
			int cnt= 0;
			double[] maxis= new double[sense.length/ 2];
			for (int i = 0; i < maxis.length; i++) 
				maxis[i]= 0;
			
			double avg= 0;
			for (int i = 0; i < sense.length; i++) {
				double val= sense[i]+ asense[i];
//				max= Math.max(max, sense[i]);
//				max= Math.max(max, asense[i]);
				avg+= val;
				if (val> maxis[0]) {
					System.arraycopy(maxis, 0, maxis, 1, maxis.length- 1);
					maxis[0]= val;
					++cnt;
				}
				//max= Math.max(max, val);
			}
			avg/= sense.length;
			
			int midIdx= -1;
			cnt= Math.min(cnt, maxis.length);
			for (int i = 1; i < cnt; i++) {
				double delta= Math.abs(maxis[i]- maxis[i-1]);
				double frac= delta/ Math.max(maxis[i], maxis[i- 1]);
				if (frac< convFactor) {
					midIdx= i;
					break;
				}
			}
			if (midIdx< 0)
				midIdx= cnt/ 2;
	
/*			if (maxis[midIdx]== 0)
				nfactor= 1;
			else {
				double filled= maxis[midIdx]* sense.length;
				double observed= sums+ suma;
				nfactor= filled/ observed;
			}
*/						
			// 20101205: choose avg or median according to
			// you compute median and average
			// you compute the normalization factor for the median (call it N_m) and the corresponding one for the average (call it N_a)
			// decision rule:
			// if N_m<1 and N_a<1, then you take the largest one (the one closest to 1)
			// if N_m>1 and N_a>1, then you take the smallest one (the one closest to 1)
			// if N_m>1 and N_a<1, then you take the one such that N_m and 1/N_a is smallest (the one closest to 1)
			// if N_m<1 and N_a>1, then you take the one such that 1/N_m and N_a is smallest (the one closest to 1).
			double obscov= (sums+ suma)/ (2d* sense.length);
			double N_m= (maxis[midIdx]> obscov)? obscov/ maxis[midIdx]:maxis[midIdx]/ obscov, 
					N_a= (avg> obscov)? obscov/ avg:avg/ obscov;
			nfactor= Math.max(N_m, N_a);

			if (Double.isNaN(nfactor)|| Double.isInfinite(nfactor)|| nfactor< 0.1||nfactor> 10) {
				System.err.println("\n\tunusual normalization factor: "+nfactor);
				System.currentTimeMillis();
			}
			
/*			int lastval= -1, minIdx= -1;
			double min= Double.MAX_VALUE, lastDelta= Double.NaN;
			for (int i = 0; i < sense.length; i++) {
				int val= sense[i]+ asense[i];
				if (lastval!= -1&& val+ lastval> 0) {
					double delta= Math.abs(val- lastval)/ (double) (Math.max(val, lastval));
					if (!Double.isNaN(lastDelta)) {
						if (delta+ lastDelta< min) {
							min= delta+ lastDelta;
							minIdx= i;
						}
					}
					lastDelta= delta;
				}
				lastval= val;
			}
			if (min== 0)
				nfactor= 1;
			else {
				double avg= (sense[minIdx]+ asense[minIdx]+ sense[minIdx- 1]+ asense[minIdx- 1])/ 2;
				double filled= avg* sense.length;
				double observed= sums+ suma;
				nfactor= filled/ observed;
			}
*/			
				
		}	
		
		return nfactor;
	}

	public boolean hasEmptyPositions() {
		for (int i = 0; i < sense.length; ++i) 
			if (sense[i]== 0)
				return true;
		for (int i = 0; i < asense.length; ++i) 
			if (asense[i]== 0)
				return true;
		return false;
	}

}
		