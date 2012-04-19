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

public class DoubleVector {

	static final int DEFAULT_LOAD_FACTOR= 13;
	double[] vector;
	public int length= -1;
	int incrementSize;
	
	public DoubleVector() {
		this(DEFAULT_LOAD_FACTOR);
	}
	
	public DoubleVector(int initialCapacity) {
		this(initialCapacity, DEFAULT_LOAD_FACTOR);
	}
	
	public DoubleVector(int initialCapacity, int loadFactor) {
		vector= new double[initialCapacity];
		incrementSize= loadFactor;
		length= 0;
	}
	
	public double getValue(int pos) {
		if (pos< length)
			return vector[pos];
		return 0;
	}
	
	public void add(double x) {
		if (length== vector.length)
			extendVector();
		vector[length++]= x;
	}
	
	void extendVector() {
		double[] newVector= new double[vector.length+ incrementSize];
		for (int i = 0; i < vector.length; i++) 
			newVector[i]= vector[i];
		vector= newVector;
	}
	
	public int size() {
		return length;
	}
	
	public double remove(int pos) {
		double result= vector[pos];
		--length;
		for (int i = pos; i < length; i++) 
			vector[pos]= vector[pos+1];
		return result;
	}
	
	public double get(int pos) {
		return vector[pos];
	}
	
	public void set(int pos, double val) {
		if (pos< vector.length)
			vector[pos]= val;
	}
	
	public double[] toDoubleArray() {
		double[] result= new double[length];
		for (int i = 0; i < length; i++) 
			result[i]= vector[i];
		return result;
	}
	
	public String toString() {
		StringBuffer sb= new StringBuffer("[");
		for (int i = 0; i < length; i++) {
			sb.append(java.lang.Double.toString(vector[i]));
			if (i< length-1)
				sb.append(",");
		}
		sb.append("]");
		return sb.toString();
	}

	public void insert(double val, int p) {
		
		if (p< 0)
			p= (p+1)* (-1);
		
		double[] newA= new double[vector.length+ 1];
		for (int i = 0; i < p; i++) 
			newA[i]= vector[i];
		newA[p]= val;
		for (int i = p+1; i < newA.length; i++) 
			newA[i]= vector[i-1];
		
		vector= newA;
		length= newA.length;
	}

	public void removeAll() {
		this.length= 0;
	}
}
