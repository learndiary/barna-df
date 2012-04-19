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

import barna.io.gtf.GTFwrapper;
import barna.model.Gene;
import barna.model.Transcript;
import barna.model.commons.IntVector;

import java.io.*;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;

/**
 * Assigns the closest promotor to a list of transcripts
 * @author micha
 *
 */
public class PromoterAssigner {
	
	public static void main(String[] args) {
		
		args= new String[] {
			"/Volumes/NastyMondays/data/demassy/download_new/analysis/cisgenome/mouse_promoters/regulatory_regionIDs.bed",
			"/Volumes/NastyMondays/annotation/mm9/mm9_refGene_fromUCSC110303.gtf",
			"/Volumes/NastyMondays/annotation/mm9/mm9_refGene_fromUCSC110303_tassign.txt"
		};
		
		File bedFile= new File(args[0]);
		File gtfFile= new File(args[1]);
		File outFile= new File(args[2]);
		
		PromoterAssigner myAssigner= new PromoterAssigner();
		HashMap<String, int[][]> map= myAssigner.getElementMap(bedFile);
		HashMap<String, Integer> tmap= myAssigner.assignElements(map, gtfFile);
		
		PrintWriter p= null;
		try {
			p = new PrintWriter(outFile);
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		myAssigner.output(tmap, map, p);
		p.close();
	}
	
	public void output(HashMap<String, Integer> tmap,
			HashMap<String, int[][]> map, PrintWriter p) {

		Iterator<String> iter= tmap.keySet().iterator();
		while(iter.hasNext()) {
			String key= iter.next();
			String[] ss= key.split("@");
			int[][] a= map.get(ss[0]);
			int q= tmap.get(key);
			p.println(ss[0]+ "\t"+ 
					ss[1]+ "\t"+ 
					Integer.toString(a[0][q])+ "\t"+ 
					Integer.toString(a[1][q]));
		}
		
	}

	/**
	 * Nucleotides downstream of TSS to be considered as promoter
	 */
	protected int tssDstream= 100;
	
	/**
	 * Maximum distance (nt) upstream of TSS to be considered as promoter
	 */
	protected int tssUstream= 10000;
	
	/**
	 * 
	 * @param bedFile
	 * @return
	 */
	public HashMap<String, int[][]> getElementMap(File bedFile) {
		
		HashMap<String, IntVector[]> map= new HashMap<String, IntVector[]>();
		try {
			
			BufferedReader buffy= new BufferedReader(new FileReader(bedFile));
			for(String s; (s= buffy.readLine())!= null;) {
				String[] ss= s.split("\\s");
				IntVector[] a= null;
				if (map.containsKey(ss[0])) {
					a= map.get(ss[0]);
				} else {
					a= new IntVector[2];
					a[0]= new IntVector();
					a[1]= new IntVector();
					map.put(ss[0], a);
				}
				a[0].add(Integer.parseInt(ss[1]));
				a[1].add(Integer.parseInt(ss[2]));
			}
			
		} catch (Exception e) {
			throw new RuntimeException(e);
		}
		
		// convert
		HashMap<String, int[][]> map2= new HashMap<String, int[][]>(map.size(), 1f);
		Iterator<String> iter= map.keySet().iterator();
		while (iter.hasNext()) {
			String key= iter.next();
			int[][] a= new int[2][];
			map2.put(key, a);
			IntVector[] b= map.get(key);
			a[0]= b[0].toIntArray();
			a[1]= b[1].toIntArray();
			Arrays.sort(a[0]);
			Arrays.sort(a[1]);
		}
		return map2;
	}
	
	public HashMap<String, Integer> assignElements(HashMap<String, int[][]> map, File gtfFile) { 
		try {
			GTFwrapper reader= new GTFwrapper(gtfFile);
//			System.err.println("checking");
//			if (!reader.isApplicable()) {
//				System.err.println("not sorted");
//				System.exit(-1);
//			}
			
			Gene[] g= null;
			HashMap<String, Integer> tmap= new HashMap<String, Integer>(map.size());
			for (reader.read(); (g= reader.getGenes())!= null; reader.read()) {
				for (int i = 0; i < g.length; i++) {
					Transcript[] t= g[i].getTranscripts();
					int[][] a= map.get(g[i].getChromosome());
					if (a== null)
						continue;
					for (int j = 0; j < t.length; j++) {
						int p= Math.abs(t[j].get5PrimeEdge());
						if (t[j].getStrand()>= 0) {
							p+= 100;
							int q= Arrays.binarySearch(a[0], p);
							q= (q< 0? -(q+ 1): q);
							if (q== 0)
								continue;
							--q;
							int d= (p- a[1][q]);
							if (d> tssDstream+ tssUstream)
								continue;
							tmap.put(t[j].getChromosome()+ "@"+ t[j].getTranscriptID(), q);
							
						} else {
							p-= 100;
							
							int q= Arrays.binarySearch(a[1], p);
							q= (q< 0? -(q+ 1): q);
							if (q== 0|| q>= a[0].length)
								continue;
							int d= (a[0][q]- p);
							if (d> tssDstream+ tssUstream)
								continue;
							tmap.put(t[j].getChromosome()+ "@"+ t[j].getTranscriptID(), q);
						}
					} /* transcripts */
				} /* genes */
			} /* reader */
			
			return tmap;
			
		} catch (Exception e) {
			throw new RuntimeException(e);
		}
	}
}
