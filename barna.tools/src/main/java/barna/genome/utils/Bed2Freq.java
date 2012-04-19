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

import barna.commons.ByteArrayCharSequence;
import barna.io.BufferedBACSReader;
import barna.model.Graph;
import barna.model.bed.BEDobject2;

import java.io.File;
import java.io.FileInputStream;



public class Bed2Freq {

	public static void main(String[] args) {
		File genDir= new File("");
		File bedFile= new File("");
		//long[][] freqs= getFreq(bedFile,0,10);
		//printFreq(freqs);
	}

	protected static long[][] getFreq(File bedFile, int fivePrime, int threePrime) throws Exception {
		
		int mLen= (-fivePrime)+ threePrime;
		long[][] f= new long[4][];
		for (int k = 0; k < f.length; k++) {
			f[k]= new long[mLen];
			for (int j = 0; j < f[k].length; j++) 
				f[k][j]= 0;
		}
		
		BEDobject2 bed= null;
		ByteArrayCharSequence cs= new ByteArrayCharSequence(150),
			buf= new ByteArrayCharSequence(mLen);
		BufferedBACSReader reader= new BufferedBACSReader(new FileInputStream(bedFile));
		while (reader.readLine(cs)!= null) {
			bed= new BEDobject2(cs);
			if (bed.getBlockCount()> 1)
				continue;
			Graph.readSequence(bed.getChr(), 
					bed.getStrand()>= 0, 
					bed.getStart()+ 1, 
					bed.getEnd(), 
					buf, 
					0, mLen);
			buf.toUpperCase();
			
		}
		
		return null;
	}
}
