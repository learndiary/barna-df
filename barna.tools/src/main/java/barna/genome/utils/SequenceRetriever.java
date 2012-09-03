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

public class SequenceRetriever {
	
	static final HashMap<Character, Character> mapRC= new HashMap<Character, Character>(22);
	static{
		mapRC.put('a', 't');
		mapRC.put('c', 'g');
		mapRC.put('g', 'c');
		mapRC.put('t', 'a');
		mapRC.put('A', 'T');
		mapRC.put('C', 'G');
		mapRC.put('G', 'C');
		mapRC.put('T', 'A');
		
		mapRC.put('n', 'n');
		mapRC.put('N', 'N');
		mapRC.put('r', 'y');
		mapRC.put('R', 'Y');
		mapRC.put('y', 'r');
		mapRC.put('Y', 'R');
		mapRC.put('s', 's');
		mapRC.put('S', 'S');
		mapRC.put('w', 'w');
		mapRC.put('W', 'W');
		mapRC.put('k', 'm');
		mapRC.put('K', 'M');
		mapRC.put('m', 'k');
		mapRC.put('M', 'K');
	}
	
	public static void getSequenceFromGTF(File fileGTF, File fastaDir, PrintStream p) {
		try {
			BufferedReader buffy= new BufferedReader(new FileReader(fileGTF)), bufasta= null;
			StringBuffer sb= new StringBuffer();
			FileInputStream istream= null;
			String lchr= null, lstr= null, lline= null;
			int pos= -1, lineLen= -1;
			for (String s= null; (s= buffy.readLine())!= null; ) {
				String[] ss= s.split("\\s");
				int start= Integer.parseInt(ss[3]), end= Integer.parseInt(ss[4]); // incl.
				sb.ensureCapacity(end- start+ 1);
				if ((!ss[0].equals(lchr))|| (!ss[6].equals(lstr))|| start< pos) {
					if (istream!= null) {
						bufasta.close();
						istream.close();
					}
					istream= new FileInputStream(new File(fastaDir.getAbsolutePath()+ File.separator+ ss[0]+ ".fa"));
					bufasta= new BufferedReader(new InputStreamReader(istream));
					bufasta.readLine(); // >
					lline= bufasta.readLine();
					lineLen= lline.length();
					pos= 1;	// 1-based
					lchr= ss[0];
					lstr= ss[6];					
				}

				// skip to start line
				int skipLines= (start- pos)/ lineLen;
				if (skipLines> 0) {
					bufasta.skip((skipLines- 1)* (lineLen+ 1)); // FS len
					pos+= skipLines* lineLen;
					lline= bufasta.readLine();
				}
				// first line
				sb.append(lline.substring(start- pos, Math.min(end- pos+ 1, lineLen)));
				// intermediate
				while(pos+ lineLen<= end) {
					lline= bufasta.readLine();
					pos+= lineLen;
					sb.append(lline.substring(0, Math.min(end- pos+ 1, lineLen)));
				}
				
				if (ss[6].equals("-")) 
					reverseComplement(sb);
				
				p.println(">"+ ss[0]+ ":"+ ss[3]+ "-"+ ss[4]+ ""+ ss[6]);
				p.println(sb.toString());
				sb.setLength(0);
				//p.println(ss[0]+"\t"+ss[3]+"\t"+ss[4]+"\t"+seq.length()+" "+(end- start+ 1)+barna.commons.system.OSChecker.NEW_LINE+seq);
				//System.currentTimeMillis();
			}
			
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	static void reverseComplement(StringBuffer sb) {
		int len= sb.length(), half= len/ 2;		
		for (int i = 0; i < half; ++i) {
			int i2= len- 1- i;
			char c1= sb.charAt(i), c2= sb.charAt(i2);
//			if (Character.isLowerCase(c1)|| Character.isLowerCase(c2)) {
//				char d1= mapRC.get(c1), d2= mapRC.get(c2);
//				System.currentTimeMillis();
//			}
			sb.setCharAt(i, mapRC.get(c2));
			sb.setCharAt(i2, mapRC.get(c1));
		}
		if (len% 2!= 0) 
			sb.setCharAt(half, mapRC.get(sb.charAt(half)));		
	}
	
	public static void main(String[] args) {
		try {
//			getSequenceFromGTF(new File("/Users/micha/annotation/hg19_RefSeq_fromUCSC100615_introns_uniq.gtf"), new File("/genomes/hg19"), 
//					new PrintStream(new File("/Users/micha/annotation/hg19_RefSeq_fromUCSC100615_introns_uniq.mfasta")));
			getSequenceFromGTF(new File("/Users/micha/projects/demassy/download/IP5300109chrall_F.gtf"), new File("/Users/micha/genomes/mm9"), 
					new PrintStream(new File("/Users/micha/projects/demassy/download/IP5300109chrall_F.fasta")));

		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
	}
}
