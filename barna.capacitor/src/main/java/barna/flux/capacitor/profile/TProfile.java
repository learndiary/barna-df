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

package barna.flux.capacitor.profile;

import barna.flux.capacitor.matrix.Matrix;
import barna.flux.capacitor.matrix.SimpleMatrix;
import barna.flux.capacitor.matrix.UniversalMatrix;
import barna.model.constants.Constants;

import java.util.Comparator;

/**
 * @deprecated
 */
public class TProfile {

	public static class Tuple {
		int x= 0, y= 0, count= 1;
		public Tuple(int x, int y) {
			this.x= x;
			this.y= y; 
		}
	}
	
	static class TupleByXYComparator implements Comparator<Tuple> {
		//@Override
		public int compare(Tuple o1, Tuple o2) {
			if (o1.x< o2.x)
				return -1;
			if (o1.x> o2.x)
				return 1;
			
			if (o1.y< o2.y)
				return -1;
			if (o1.y> o2.y)
				return 1;
			return 0;
		}
	}
	
	static int readLen= 0;
	static TupleByXYComparator defaultTupleByXYComparator= new TupleByXYComparator();
	
	Matrix m;
	boolean strandSpecific= false;
	String ID= Constants.EMPTYSTRING;
	
	public TProfile(int length, boolean strandSpecific, boolean pairedEnd) {
		
		//System.err.println("profile: "+(++profileCtr));
		
		byte dir= strandSpecific? barna.model.constants.Constants.DIR_BOTH: barna.model.constants.Constants.DIR_FORWARD;
		m= (Matrix)	// TODO 
			(pairedEnd?new UniversalMatrix(length):new SimpleMatrix(length, dir));
		
	}
	
	
	public TProfile(String ID, int length, boolean strandSpecific, boolean pairedEnd) {
		this(length, strandSpecific, pairedEnd);
		this.ID= ID;
	}
	
	/*
	 * values [0] readLen, [1] x1 ([2] y1)
	 */
	public void addRead(int exonicStartPos, int readLen, byte dir) {
		
		m.add(exonicStartPos, readLen, dir);
	}
	
	public void addReadPair(int pos1, int pos2, int readLen) {
		m.add(pos1, pos2, readLen);
	}

	public int getReads() {
		return m.getSum();
	}
	
	public int getLength() {
		return m.getLength();
	}

	public int length() {
		return m.getLength();
	}
	
	public void fill(int[] insertSize, int readLen) {
		m.fill(insertSize, readLen);
	}
	
	public long getArea(int[] coords, int readLen, int[] insertMinMax, byte dir) {
		
		assert(coords.length== 2|| coords.length== 4);
		long val= 0;
		if (coords.length== 2) 
			val= m.get(coords[0], coords[1], readLen, insertMinMax, dir); // TODO branch for forward/reverse
		else
			val= m.get(coords[0], coords[1], coords[2], coords[3], readLen);
		
		return val;
	}

	public void addProfile(TProfile profile, int readLen) {
			
		m.merge(profile.m, readLen);
	}
	
	public String printMatrix() {
		return null;
	}
	
	@Override
	public String toString() {
		return m.toString();
	}
	
	public byte[] toByteArray() {
		String s= toString();
		byte[] b= new byte[s.length()];
		for (int i = 0; i < b.length; i++) 
			b[i]= (byte) s.charAt(i);
		return b;
	}

	/**
	 * @deprecated not called currently
	 * @param coords
	 * @param readLen
	 * @param dir
	 * @return
	 */
	public double getAreaFrac(int[] coords, int readLen, int[] insertMinMax, byte dir) {
		
		assert(coords.length== 2|| coords.length== 4);
		double val= 0;
		if (coords.length== 2) 
			val= m.get(coords[0], coords[1], readLen, insertMinMax, dir); // TODO branch for forward/reverse
		else
			val= m.get(coords[0], coords[1], coords[2], coords[3], readLen);
		
		val/= m.getSum();
		return val;
	}


	public String getID() {
		return ID;
	}
}
