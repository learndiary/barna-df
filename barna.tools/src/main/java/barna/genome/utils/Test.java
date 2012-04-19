package barna.genome.utils;/*
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

import barna.model.splicegraph.SplicingGraph;


public class Test {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
//		float f= 0.00001f;
//		String s= String.format("%1$f", f);
//
//		System.out.println(s+"*");
//		
//		System.err.println("3 mod 0"+(3%0));
		
		long[] a= new long[2];
		for (int i = 0; i < a.length; i++) {
			for (int j = 0; j < 64; j++) {
				a[i]|= (int) Math.pow(2, j);
			}
		}
		
		for(int idx= -1; (idx= SplicingGraph.getNextTxIdx(a, idx))!= -1; ){
			System.out.println(idx);
		}
		
		long x= 1;
		for (int i = 0; i < 64; i++, x*= 2) {
			System.err.println(i+"\t"+x);
		}
	}

	void kstest() {
		//StatisticalComparison.compare(null, null);
	}
}
