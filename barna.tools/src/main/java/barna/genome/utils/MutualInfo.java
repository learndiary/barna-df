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

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.PrintStream;
import java.util.HashMap;


public class MutualInfo {

	public static final char[] SYMBOLS= new char[] {'A', 'C', 'G', 'T'};
	
	int n= 0;
	String[] tupels= null;
	
	public MutualInfo(int tupelSize) {
		n= tupelSize;
		assert(n>0);
		
		tupels= new String[(int) Math.pow(SYMBOLS.length, n)];
		int p= 0;
		p= generateTupels("", p);
		assert(p== tupels.length);
	}
	
	int generateTupels(String s, int p) {
		if (s.length()== n) {
			tupels[p++]= s;
			return p;
		}
		
		for (int i = 0; i < SYMBOLS.length; i++) 
			p= generateTupels(s+ SYMBOLS[i], p);
		
		return p;
	}
	
	private static HashMap<String, String> stringMap= new HashMap<String, String>();
	private String getKey(String s) {
		if (stringMap.containsKey(s))
			return stringMap.get(s);
		stringMap.put(s, s);
		return s;
	}
	
	double[][] getMatrix(File seqFile) {
		
		try {
			
			// count occurrences
			BufferedReader buffy= new BufferedReader(new FileReader(seqFile));
			HashMap<String, Integer>[][] m= null;
			int cnt= 0;
			for (String s= null; (s= buffy.readLine())!= null;) {
				++cnt;
				String[] ss= s.split("\\s");
				ss[0]= ss[0].toUpperCase();
				//ss[0]= ss[0].substring(10,30);
				if (m== null) {
					m= new HashMap[ss[0].length()- n][];
					for (int i = 0; i < ss[0].length()- n; i++) {
						m[i]= new HashMap[ss[0].length()- n];
						for (int j = i; j < m[i].length; j+= 1) {
							m[i][j]= new HashMap<String, Integer>(tupels.length);
						}
					}
				}
				for (int i = 0; i < ss[0].length()- n; i++) {
					String tup1= ss[0].substring(i, i+ n);
					if (m[i][i].containsKey(tup1)) {
						m[i][i].put(getKey(tup1), m[i][i].get(tup1)+ 1);						
					} else {
						m[i][i].put(getKey(tup1), 1);
					}
					for (int j = i+ n; j < ss[0].length()- n; j+= 1) {
						String tup2= tup1+ ss[0].substring(j, j+ n);
						if (m[i][j].containsKey(tup2)) {
							m[i][j].put(getKey(tup2), m[i][j].get(tup2)+ 1);						
						} else {
							m[i][j].put(getKey(tup2), 1);
						}
					}
				}
			}
			buffy.close();
			
			// build information matrix
			double[][] mm= new double[m.length][];
			for (int i = 0; i < mm.length; i++) {
				mm[i]= new double[m.length];
				for (int j = i+ n; j < mm.length; j+= 1) {
					
					mm[i][j]= 0d;
					
					Object[] oo= m[i][j].keySet().toArray();
					for (int k = 0; k < oo.length; k++) {
						String s= (String) oo[k];	// formula 0 for no joint occurrences
						String tup1= s.substring(0, n), 
								tup2= s.substring(n, s.length());
						int posi= m[i][i].get(tup1),
							posj= m[j][j].get(tup2),
							posij= m[i][j].get(s);
						//int cntij= posi+ posj;
						double pi= posi/ (double) cnt,
								pj= posj/ (double) cnt,
								pij= posij/ (double) cnt;
						double mij= pij* (Math.log(pij/ (pi* pj))/ Math.log(2)); 
						mm[i][j]+= mij;
					}
				}
			}
			
			return mm;
			
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		return null;
	}
	
	public static void main(String[] args) {
		File inFile= new File("/Users/micha/projects/bene/mutual_info/leeF.peak.seq"),
				outFile= new File("/Users/micha/projects/bene/mutual_info/leeF.peak.minf-4");
		int n= 4;
		MutualInfo myMI= new MutualInfo(n);
		double[][] m= myMI.getMatrix(inFile);
		
		try {
			PrintStream p= new PrintStream(outFile);
			for (int i = 0; i < m.length; i++) {
				for (int j = 0; j < i+ n&& j< m.length; j+= 1) {
					p.print("NA\t");
				}
				for (int j = i+ n; j < m.length; j+= 1) {
					p.print(m[i][j]+"\t");
				}
				p.println();
			}
			p.flush();
			p.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
}
