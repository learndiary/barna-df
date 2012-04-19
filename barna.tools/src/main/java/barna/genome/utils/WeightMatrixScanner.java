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

import java.io.*;
import java.util.HashMap;


public class WeightMatrixScanner {

	static HashMap<Character, Integer> mapChar2Pos;
	static {
		mapChar2Pos= new HashMap<Character, Integer>(5);
		mapChar2Pos.put('a', 0);
		mapChar2Pos.put('c', 1);
		mapChar2Pos.put('g', 2);
		mapChar2Pos.put('t', 3);
		mapChar2Pos.put('n', 4);
	}
	static int ctrSeq= 0, ctrFound= 0;
	
	
	/**
	 * 
	 * @param f
	 * @param from 1-based
	 * @param to 1-based, excluded
	 * @return
	 */
	static double[][] loadMatrix(File f, int from, int to) {
		try {
			double[][] m= new double[5][];
			BufferedReader buffy= new BufferedReader(new FileReader(f));
			for (String s; (s= buffy.readLine())!= null; ) {
				if (s.startsWith("#"))
					continue;
				String[] ss= s.split("\\s");
				char c= Character.toLowerCase(ss[0].charAt(0));
				int pos= mapChar2Pos.get(c);
				if (from< 1)
					from= 1;
				if (to< 1)
					to= ss.length;
				int len= to- from;
				m[pos]= new double[len];
				for (int i = from; i < to; i++) 
					m[pos][i-from]= Double.parseDouble(ss[i]);
			}
			buffy.close();
			return m;
		} catch (Exception e) {
			e.printStackTrace();
		}
		return null;
	}
	
	static void scanSequences(String s, int from, int to, String[] tag, double[][] m, double thr, PrintStream p) {
		int mlen= m[0].length;
		for (int i = 0; i <= s.length()- Math.max(mlen,to); i++) {
			double score= 1d;
			for (int j = 0; j < mlen; j++) {
				char c= s.charAt(i+ j);
				score*= m[mapChar2Pos.get(c)][j];
			}
			if (score>= thr&& i+ from> 0&& i+ to< s.length()) {
				++ctrFound;
				String t= s.substring(i+ from, i+ 1+ to);
				p.println(t+ "\t"+ Double.toString(score)+ "\t"+ tag[0]+ "\t"+ (i+1));					
			}
		}
	}
	
	static String getSequenceFromFasta(File f, String[] tag, int from, int to) {
		try {
			FileInputStream fin= new FileInputStream(f);
			fin.skip(bytesRead);
			BufferedReader buffy= new BufferedReader(new InputStreamReader(fin));
			StringBuffer sb= new StringBuffer(500);
			for (String s; (s= buffy.readLine())!= null; ) {
				if (s.startsWith(">")) {	
					if (sb.length()> 0) {
						buffy.close();
						String t= sb.toString();
						if (from> 0)
							t= t.substring(from- 1);
						if (to> 0)
							t= t.substring(0, to);
						++ctrSeq;
						return t;
					}
					tag[0]= s.substring(1);
				} else {
					sb.append(s);
				}
				bytesRead+= s.length()+ 1;
			}
			buffy.close();
			String t= sb.toString();
			if (from> 0)
				t= t.substring(from- 1);
			if (to> 0)
				t= t.substring(0, to);
			++ctrSeq;
			return t;
		} catch (Exception e) {
			e.printStackTrace();
		}
		return null;
	}
	
	static long bytesRead= 0;
	/**
	 * 
	 * @param f
	 * @param tag
	 * @param from 1-based
	 * @param to 1-based, included
	 * @return
	 */
	static String getSequenceFromPeakFile(File f, String[] tag, int from, int to) {
		try {
			FileInputStream in= new FileInputStream(f);
			in.skip(bytesRead);
			BufferedReader buffy= new BufferedReader(new InputStreamReader(in));
			String s= buffy.readLine();
			if (s== null)
				return null;
			bytesRead+= s.length()+ 1;
			buffy.close();
			String[] ss= s.split("\\s");	
			++ctrSeq;
			from= Math.max(1, from);
			if(to<= 0)
				to=ss[0].length();
			return ss[0].substring(from- 1, to);
		} catch (Exception e) {
			e.printStackTrace();
		}
		return null;
	}
	
	public static void scanPeaks(String[] args) {
		File wFile= new File("/Users/micha/projects/bene/weight_matrix_forward_peak.txt");
		File pFile= new File("/Users/micha/projects/bene/new/liverF_2.peak25.seq");
		File oFile= new File("/Users/micha/projects/bene/refMrna.seqScores_gt10-5");
		PrintStream oStream= null;
		try {
			oStream= new PrintStream(oFile);
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		
		double[][] m= loadMatrix(wFile, 11, 21);
		String[] tag= new String[1];
		String s;
		while ((s= getSequenceFromPeakFile(pFile, tag, 11, 21))!= null) {
			scanSequences(s, -10, 100, tag, m, 0, oStream);
		}
		
	}
	
	public static void main(String[] args) {
		
		if (args.length< 3) {
			System.err.println("Usage: WeightedMatrixScanner matrixFile inputFile outputFile [threshold(Default:0) seq-from(DEFAULT:0) seq-to(DEFAULT:length-motif) output-from(DEFAULT:-10) output-to(DEFAULT:+10)]");
			System.exit(-1);
		}
		
//		File wFile= new File("/Users/micha/projects/bene/weight_matrix_forward_peak.txt");
//		File sDir= new File("/Users/micha/projects/bene/refMrna");
//		File oFile= new File("/Users/micha/projects/bene/new/liverF_2.peak25.seqScores-10+100");
		File wFile= new File(args[0]);
		File sDir= new File(args[1]);
		File oFile= new File(args[2]);
		PrintStream oStream= null;
		try {
			oStream= new PrintStream(oFile);
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		int from= -1, to= -1, outFrom= -10, outTo= 10;		
		float threshold= 0f;
		if (args.length> 3) {
			try {
				threshold= Float.parseFloat(args[3]);
				from= Integer.parseInt(args[4]);
				to= Integer.parseInt(args[5]);
				outFrom= Integer.parseInt(args[6]);
				outTo= Integer.parseInt(args[7]);
			} catch (Exception e) {
				e.printStackTrace();
				System.exit(-1);
			}
		}
		
		System.err.println("settings:");
		System.err.println("matrix\t"+wFile);
		System.err.println("input\t"+sDir);
		System.err.println("output\t"+oFile);
		System.err.println("threshold\t"+threshold);
		System.err.println("scan sequences (from,to):\t"+from+","+to);
		System.err.println("output motif flans (from,to):\t"+outFrom+","+outTo);
		
		System.err.print("\nloading Matrix..");
		double[][] m= loadMatrix(wFile, 11, 21);
		if (m== null) {
			System.err.println("error");
			System.exit(-1);
		}
		System.err.println("ok");
		
		if (sDir.isDirectory()) {
			System.err.println("\nscanning files in directory "+sDir.getName()+" ");
			File[] files= sDir.listFiles();
			String[] tag= new String[1];
			for (int i = 0; i < files.length; i++) {
				String name= files[i].getName();
				if (name.contains(".")|| !name.startsWith("NM_"))
					continue;
				tag[0]= name;
				String s= getSequenceFromFasta(files[i], tag, from, to);
				if (s== null)
					break;
				scanSequences(s, -10, 100, tag, m, threshold, oStream);
			}
		} else {
			// check for fasta
			String[] tag= new String[1];
			System.err.println("\nscanning files in file "+sDir.getName()+" ");
			try {
				BufferedReader buffy= new BufferedReader(new FileReader(sDir));
				String s= buffy.readLine();
				buffy.close();
				if (s.startsWith(">")) {
					while ((s= getSequenceFromFasta(sDir, tag, from, to))!= null) {
						scanSequences(s, outFrom, outTo, tag, m, threshold, oStream);
					}
				} else {
					while ((s= getSequenceFromPeakFile(sDir, tag, from, to))!= null) {
						scanSequences(s, outFrom, outTo, tag, m, threshold, oStream);
					}
				}
			} catch (Exception e) {
				e.printStackTrace();
				System.err.println("error during detecting file format");
				System.exit(-1);
			}
		}
		
		System.err.println("found "+ctrFound+" motifs in "+ctrSeq+" sequences.");
		
	}
}
