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

import barna.model.commons.IntVector;

import java.util.Vector;

public class TupleMatrix implements Matrix {

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
	
	static class Tuple {
		int len, offset;
		public Tuple(int len, int offset) {
			this.len= len;
			this.offset= offset;
		}
		@Override
		protected Object clone() throws CloneNotSupportedException {
			Tuple t= new Tuple(len, offset);		
			return t;
		}
	}
	
	IntVector[] m;
	int s;
	
	public TupleMatrix(int length) {
		m= new IntVector[length];
		s= 0;
	}
	
	public void add(int p, int len, byte dir) {
		; // :)
	}

	public void add(int p, byte dir) {
		; // :)
	}
	
	public void add(int p1, int p2, int readLen) {
		if (p1> p2) {
			int h= p1;
			p1= p2;
			p2= h;
		}
		
		if (m[p1]== null) {
			m[p1]= getIntVector(1);
			if (m[p1]== null)
				m[p1]= new IntVector(1,10);
		}
		//p2+= readLen- 1;
		m[p1].add(p2- p1);
		++s;
	}

	public void fill(int[] insertSize, int readLen) {
		for (int i = 0; i < m.length; i++) {
			m[i]= new IntVector(1,1);
			for (int j = (i+ readLen+ insertSize[0]); 
					j < (i+ readLen+ insertSize[1])
					&& j< (m.length- readLen); j++) {
				int offset= j- i;
				m[i].add(offset);	// absolute offset
				++s;
			}
		}
	}

	/**
	 * for single read fractions, fragment ends are counted twice
	 */
	public int get(int p1, int p2, int readLen, int[] insertMinMax, byte dir) {
		
		int sum= 0;
		for (int i = p1; i < p2; i++) {
			if (m[i]== null)
				continue;
			sum+= m[i].size();
		}
		
		int minUSdist= Math.max(0, p2- (insertMinMax[0]+ readLen));
		int maxUSdist= Math.max(0, p1- (insertMinMax[1]+ readLen));
		for (int i = maxUSdist; i <= minUSdist; ++i) {
			if (m[i]== null)
				continue;
			for (int j = 0; j < m[i].length; j++) {
				int t= i+ m[i].get(j);
				if (t>= p1&& t<= p2)
					++sum;
			}
		}
		
		return sum;
	}

	public int get(int p1, int p2, int p3, int p4, int readLen) {
		int sum= 0;
//		for (int i = p1-1; i>= 0; --i) {
//			if (m[i]== null)
//				continue;
//			int j = 0;
//			for (; j < m[i].size(); j++) {
//				if (i+ m[i].elementAt(j).len< p1)
//					break;
//				for (int k = 0; k < m[i].elementAt(j).len; k++) {
//					int h= i+ k, l= h+ m[i].elementAt(j).offset;
//					if (h>= p1&& h<= p2&& l>= p3&& l<= p4)
//						++s;
//				}
//			}
//			if (j< m[i].size())
//				break;
//		}
		//p3+= readLen- 1;
		//p4+= readLen- 1;
		
		// 090820 <p2 
		// exclusive regions for back-normalization needed
		// otherwise fracs> transcriptcount for gene (and reads also)
		for (int i = p1; i < p2; i++) {
			if (m[i]== null)
				continue;
			for (int j = 0; j < m[i].size(); j++) {
				int t= i+ m[i].get(j);
				if (t>= p3&& t<= p4)
					++sum;
//				for (int h = 0; h < m[i].elementAt(j).len&& i+ h<= p2; h++) {
//					int k= i+ h+ m[i].elementAt(j).offset;
//					if (k>= p3&& k<= p4)
//						++s;
//				}
			}
		}
		
		return sum;
	}

	public int getLength() {	
		return m.length;
	}

	public int getSum() {
		return s;
	}

	public void merge(Matrix nn, int readLen) {
		
		TupleMatrix n= (TupleMatrix) nn; 
		
		double fac= getLength()/ (double) n.getLength();
		for (int i = 0; i < n.getLength(); i++) {
			if (n.m[i]== null)
				continue;
			int p= (int) Math.round(fac* i);
			p=(p== m.length)?p-1:p;
			assert (p< m.length&& p> -1);
			if (m[p]== null) {
				m[p]= TupleMatrix.getIntVector(n.m[i].size());
				if (m[p]== null)
					m[p]= new IntVector(n.m[i].size(), 10);
			}
			for (int j = 0; j < n.m[i].size(); j++) {
				int t= n.m[i].get(j);	// NO: do not correct insert sizes (int) Math.round(fac* n.m[i].get(j));
				if (p+ t+ readLen- 1>= m.length)
					continue;	// but skip the ones that fall outside to achieve 100% possible
				m[p].add(t);
				++s;
			}
			
		}

	}

	@Override
	public String toString() {
		return toStringBuilder().toString();
	}
	
	public StringBuilder toStringBuilder() {
		StringBuilder sb= new StringBuilder(m.length); // not * m.length, can exceed integer bounds
		for (int i = 0; i < m.length; i++) {
			if (m[i]== null)
				continue;
			sb.append(Integer.toString(i));
			sb.append("\t");
			for (int j = 0; j < m[i].size(); j++) {
				sb.append(Integer.toString(m[i].get(j)));
				sb.append("\t");
			}
			//sb.replace(sb.length()- 1, sb.length()- 1, barna.commons.system.OSChecker.NEW_LINE);
			sb.deleteCharAt(sb.length()-1);
			sb.append(barna.commons.system.OSChecker.NEW_LINE);
		}
		return sb;
	}
	
	public byte[] toByteArray() {
		StringBuilder sb= toStringBuilder();
		byte[] b= new byte[sb.length()];
		for (int i = 0; i < b.length; i++) 
			b[i]= (byte) sb.charAt(i);
		
		return b;
	}

	public int project(int[][] b) {
		// TODO Auto-generated method stub
		return 0;
	}
}
