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

public class FloatVector {

	static final int DEFAULT_LOAD_FACTOR= 13;
	public float[] vector;
	public int length= -1;
	int incrementSize;
	
	public FloatVector() {
		this(DEFAULT_LOAD_FACTOR);
	}
	
	public FloatVector(int initialCapacity) {
		this(initialCapacity, DEFAULT_LOAD_FACTOR);
	}
	
	public FloatVector(int initialCapacity, int loadFactor) {
		vector= new float[initialCapacity];
		incrementSize= loadFactor;
		length= 0;
	}
	
	public void add(float x) {
		if (length== vector.length)
			extendVector();
		vector[length++]= x;
	}
	
	@Override
	protected Object clone() throws CloneNotSupportedException {
		FloatVector newV= new FloatVector();
		newV.vector= this.vector.clone();
		newV.length= this.length;
		newV.incrementSize= this.incrementSize;
		return newV;
	}
	
	public FloatVector cloneIntVector() {
		try {
			return (FloatVector) clone();
		} catch (CloneNotSupportedException e) {
			e.printStackTrace();
		}
		return null;
	}
	
	void extendVector() {
		float[] newVector= new float[vector.length+ incrementSize];
		for (int i = 0; i < vector.length; i++) 
			newVector[i]= vector[i];
		vector= newVector;
	}
	
	public int size() {
		return length;
	}
	
	public int capacity() {
		return vector.length;
	}
	
	public void reset() {
		incrementSize= DEFAULT_LOAD_FACTOR; 
		length= 0;
		for (int i = 0; i < vector.length; i++) 
			vector[i]= 0;
	}
	
	public float remove(int pos) {
		float result= vector[pos];
		--length;
		// TODO
//		if (pos< length)
//			System.arraycopy(vector, pos+ 1, vector, pos, length- pos);
		for (int i = pos; i < length; i++) 
			vector[pos]= vector[pos+1];
		return result;
	}
	
	public float get(int pos) {
		while(pos>= vector.length)
			extendVector();
		return vector[pos];
	}
	
	public void set(int pos, int val) {
		while (pos>= vector.length)
			extendVector();
		//if (pos< vector.length)
		vector[pos]= val;
		if (length<= pos)
			length= pos;
	}
	public float[] toFloatArray() {
		float[] result= new float[length];
		for (int i = 0; i < length; i++) 
			result[i]= vector[i];
		return result;
	}
	
	public String toString() {
		StringBuffer sb= new StringBuffer("[");
		for (int i = 0; i < length; i++) {
			sb.append(java.lang.Float.toString(vector[i]));
			if (i< length-1)
				sb.append(",");
		}
		sb.append("]");
		return sb.toString();
	}

	public void insert(int val, int p) {
		
		if (p< 0)
			p= (p+1)* (-1);
		
		float[] newA= new float[vector.length+ 1];
		for (int i = 0; i < p; i++) 
			newA[i]= vector[i];
		newA[p]= val;
		for (int i = p+1; i < newA.length; i++) 
			newA[i]= vector[i-1];
		
		vector= newA;
		length= newA.length;
	}
	
	public float getValue(int pos) {
		if (pos< length)
			return vector[pos];
		return 0;
	}
	
	public void putValue(int pos, int val) {
		if (pos>= length) {
			int rounds= (pos+ 1)- size();
			for (int i = 0; i < rounds; i++) 
				add(0);
		}
		vector[pos]= val;
	}

	public void removeAll() {
		length= 0;
	}
}
