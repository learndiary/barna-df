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

import barna.model.constants.Constants;

public class SimpleMatrix implements Matrix {

	int[][] m;
	int s;
	
	public SimpleMatrix(int len, byte dir) {
		if (dir== Constants.DIR_BOTH)
			m= new int[2][];
		else
			m= new int[1][];
		for (int i = 0; i < m.length; i++) {
			m[i]= new int[len];
			for (int j = 0; j < m[i].length; j++) 
				m[i][j]= 0;
		}
	}
	
	public void add(int p, int len, byte dir) {
		
		if (p< 0|| p>= getLength()) {
			System.err.println("Matrix ArrayIndexOutOfBounds "+ p+" in "+ getLength());
			if (p< 0)
				p= 0;
			else
				p= getLength()- 1;
		}
		
		if (dir== Constants.DIR_BACKWARD)
			++m[1][p];
		else
			++m[0][p];
		++s;
	}

	public void add(int p1, int p2, int len) {
		; // :)
	}

	public int getSum() {
		return s;
	}

	public int getLength() {
		return m[0].length;
	}
	
	public void merge(Matrix nn, int readLen) {
		
		SimpleMatrix n= (SimpleMatrix) nn;
		
		double fac= getLength()/ (double) n.getLength();
		for (int x = 0; x < m.length; x++) {
			for (int i = 0; i < n.getLength(); i++) {
//				int p= (int) Math.round(fac* i);
//				if (p+ readLen- 1>= m[x].length)
//					continue; // falls outside
//				p= (p== m[x].length?p-1:p);
				
				// make blurry, dont respect readlength..
				// ..delegate counted ones to TSuperProfile
				double floatPos= fac* i;
				int minPos= (int) Math.floor(floatPos), maxPos= (int) Math.ceil(floatPos);
				try {assert(floatPos> 0&& minPos> 0&& maxPos> 0);} catch (AssertionError err) {
					if (Constants.verboseLevel> Constants.VERBOSE_NORMAL)
						System.err.println("[ASSERTION] "+getClass().getName()+".merge() matrix index negative:\n\t"
							+"length1= "+ getLength()
							+", length2= "+ n.getLength()
							+", floatPos= "+ floatPos
							+", min/max= "+ minPos+"/"+maxPos
						);
					continue;
				}
				for (int p = minPos; p <= maxPos&& p< m[x].length; p++) {
					m[x][p]+= n.m[x][i];
					s+= n.m[x][i];
				}
			}
		}
	}
	
	public void fill(int[] insertSize, int readLen) {
		for (int i = 0; i < m.length; i++) {
			for (int j = 0; j < m[i].length; j++) {
				++m[i][j];
				++s;
			}
		}
	}

	public int get(int p1, int p2, int readLen, int[] insertMinMax, byte dir) {
		int count= 0;
		try {assert(p1>= 0&& p2<= m[0].length);} catch (AssertionError err){
			if (Constants.verboseLevel> Constants.VERBOSE_NORMAL)
				System.err.println("[ASSERTION] "+getClass().getName()+".get():\n\t" +
						"index pair ("+p1+","+p2+") out of matrix bounds (length "+m[0].length+")!");
			return 0;
		};
		// 090820 < p2 
		// exclusive regions for back-normalization needed
		// otherwise fracs> transcriptcount for gene (and reads also)
		for (int i = p1; i < p2; i++) {
			if (dir== Constants.DIR_FORWARD|| dir== Constants.DIR_BOTH)
				count+= m[0][i];
			if (dir== Constants.DIR_BACKWARD|| dir== Constants.DIR_BOTH)
				count+= m[1][i];
		}
		return count;
	}

	public int get(int p1, int p2, int p3, int p4, int readLen) {
		return 0;	// :)
	}

	private StringBuilder toStringBuilder() {
		StringBuilder sb= new StringBuilder();
		for (int i = 0; i < m[0].length; i++) {
			for (int j = 0; j < m.length; j++) {
				sb.append(Integer.toString(m[j][i]));
				sb.append("\t");
			}
			sb.deleteCharAt(sb.length()- 1);
			sb.append("\n");
		}
		return sb;
	}
	
	@Override
	public String toString() {
		return toStringBuilder().toString();
	}
	
	public byte[] toByteArray() {
		StringBuilder sb= toStringBuilder();
		
		byte[] b= new byte[sb.length()];
		for (int i = 0; i < b.length; i++) 
			b[i]= (byte) sb.charAt(i);
		
		return b;
	}
	
	public int project(int[][] b) {
		if (b== null|| b.length!= m.length)
			return -1;
		float len= m[0].length, projLen= b.length;
		int sum= 0;
		for (int i = 0; i < m.length; i++) {
			for (int j = 0; j < m[i].length; j++) {
				int p= (int) Math.round(j* projLen/ len);	// TODO perf
				int v= m[i][j];
				b[i][p]+= v;
				sum+= v;
			}
		}
		assert(s== sum);
		return sum;
	}

	/**
	 * @deprecated stub
	 */
	public void add(int p, byte dir) {
		
	}
}
