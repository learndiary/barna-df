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
import barna.io.BufferedByteArrayReader;
import barna.io.FileHelper;
import barna.io.bed.BEDwrapper;
import barna.model.bed.BEDobject2;
import barna.model.commons.IntVector;

import java.io.*;


public class Bed2Wig {

	static void myMain() {
		File dir= new File("/Users/micha/projects/demassy/download_new/regions");
		String[]  ss= dir.list();
		for (int i = 0; i < ss.length; i++) {
			if (!ss[i].endsWith("_reads.sorted.bed"))
				continue;
			File fin= new File(dir.getAbsolutePath()+ File.separator+ ss[i]);
			File fout= FileHelper.replaceSfx(fin, ".bgraph");
			
			System.err.println(fin.getName());
			bed2wig(fin, fout);
		}
	}
	
	public static void main(String[] args) {
		
//		if (1== 1)
//			myMain();
		
		if (args== null|| args.length!= 2) {
			System.err.println("Usage: Bed2Wig [input.bed] [output.wig]");
			System.exit(-1);
		}
		
		
		File bedFile= new File(args[0]);
		File wigFile= new File(args[1]);
		if (!bedFile.exists()) {
			System.err.println(args[0]+ " not found.");
			System.exit(-1);
		}
		
		bed2wig(bedFile, wigFile);
	}
	
	/**
	 * assumes sorted by (1) chr, (2) start position
	 * @param in
	 * @param out
	 */
	static void bed2wig(File in, File out) {
		try {
			
			BEDwrapper bedreader= new BEDwrapper(in);
			if (!bedreader.isApplicable()) {
				File tmp= FileHelper.createTempFile(Bed2Wig.class.getName(), in.getName());
				tmp.deleteOnExit();
				bedreader.sort(tmp);
				in= tmp;
			}
			
			//BufferedReader buffy= new BufferedReader(new FileReader(in));
			if (out.exists())
				out.delete();
			BufferedByteArrayReader buffy= new BufferedByteArrayReader();
			FileInputStream istream= new FileInputStream(in);
			ByteArrayCharSequence cs= new ByteArrayCharSequence(300);
			IntVector v= new IntVector(10000, 10000);
			int start= -1;
			int ctr= 0;
			String chr= null, lastChr= null;			
			while(buffy.readLine(istream, cs)!= null) {
				++ctr;
				if (ctr% 1000== 0) {
					System.err.print("*");
					System.err.flush();
				}
				BEDobject2 obj= new BEDobject2(cs);
				chr= obj.getChr().toString();
				
				if (!chr.equals(lastChr)) {
					writeBedGraph(in.getName(), lastChr, start, v, out);
					v.removeAll();
					start= obj.getStart();
					lastChr= chr;
				}
				if (obj.getBlockCount()< 2) {
					for (int i = obj.getStart(); i < obj.getEnd(); i++) {
						v.set(i- start, v.get(i- start)+ 1);
					}
				} else {
					int bedstart= obj.getStart();
					for (int i = 0; i < obj.getBlockCount(); i++) {
						int bsize= obj.getNextBlockSize();
						int bstart= obj.getNextBlockStart();
						for (int j = 0; j < bsize; j++) {
							int gpos= bedstart+ bstart+ j;
							if (gpos- start< 0)
								System.currentTimeMillis();
							v.set(gpos- start, v.get(gpos- start)+ 1);
						}
					}
				}
			}
			istream.close();

			writeBedGraph(in.getName(), lastChr, start, v, out);

		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	private static void writeFixedStep(String name, String chr, int start, IntVector v, File out) {
		
		// write
		try {
			boolean header= !out.exists();
			BufferedWriter writer= new BufferedWriter(new FileWriter(out, true));
			++start;
			if (header) {
				writer.write("track type=wiggle_0 name="+ name+ barna.commons.system.OSChecker.NEW_LINE);
				writer.write("fixedStep chrom="+chr+" start="+start+" step=1\n");
			}
			for (int i = 0; i < v.size(); i++) {
				writer.write(Integer.toString(v.get(i)));
				writer.write(barna.commons.system.OSChecker.NEW_LINE);
			}
			writer.flush();
			writer.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}
	
	/**
	 * @deprecated not implemented
	 * @param name
	 * @param chr
	 * @param start
	 * @param v
	 * @param out
	 */
	private static void writeVariableStep(String name, String chr, int start, IntVector v, File out) {
		
		// write
		try {
			boolean header= !out.exists();
			BufferedWriter writer= new BufferedWriter(new FileWriter(out, true));
			++start;
			if (header) {
				writer.write("track type=wiggle_0 name="+ name+ barna.commons.system.OSChecker.NEW_LINE);
				writer.write("variableStep chrom="+chr+" start="+start+" step=1\n");
			}
			for (int i = 0; i < v.size(); i++) {
				writer.write(Integer.toString(v.get(i)));
				writer.write(barna.commons.system.OSChecker.NEW_LINE);
			}
			writer.flush();
			writer.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}
	
	private static void writeBedGraph(String name, String chr, int start, IntVector v, File out) {
		
		// write
		try {
			boolean header= !out.exists();
			BufferedWriter writer= new BufferedWriter(new FileWriter(out, true));
			++start;
			if (header) {
				writer.write("track type=bedGraph name=\""+ name+ "\" description=\""+ name+"\" " +
						"visibility=dense color=200,100,0 altColor=0,100,200 priority=20\n");
			}
			for (int i = 0; i < v.size(); i++) {
				int val= v.get(i);
				int j = i+1;
				for (; j < v.size(); j++) 
					if (v.get(j)!= val)
						break;
				
				if (val!= 0) 
					writer.write(chr+ " "+ (start+ i)+ " "+ (start+ j)+" "+ Integer.toString(val)+ barna.commons.system.OSChecker.NEW_LINE);
				
				i= (j-1);
			}
			writer.flush();
			writer.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}


}
