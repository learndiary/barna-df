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

package barna.astalavista.statistics;

import barna.io.gtf.GTFwrapper;
import barna.model.gff.GFFObject;

import java.io.*;
import java.util.*;


/**
 * Annotates events with protein attributes from david database.
 * @author micha
 *
 */
public class ProteinAnnotator {

	public static class Int2Array {
		int[][] a;
		
		public Int2Array(int[][] a) {
			this.a= a;
		}
		
		@Override
		public int hashCode() {
			return toString().hashCode();
		}
		
		@Override
		public String toString() {
			StringBuffer sb= new StringBuffer();
			for (int i = 0; i < a.length; i++) {
				for (int j = 0; a[i]!=null&& j < a[i].length; j++) {
					sb.append(a[i][j]);
					sb.append(",");
				}
				sb.append(";");
			}
			return sb.toString();
		}
		
		
		@Override
		public boolean equals(Object obj) {
			int[][] b= ((Int2Array) obj).a;
			for (int i = 0; i < b.length; i++) {
				if (a[i]!= null&& b[i]!= null) {
					if (a[i].length!= b[i].length)
						return false;
					for (int j = 0; j < b[i].length; j++) {
						if (a[i][j]!= b[i][j])
							return false;
					}
				} else {
					if (a[i]== null^ b[i]== null)
						return false;
				}
			}
			return true;
		}
		
	}

	public static void main(String[] args) {
		//replaceDomainsByNumber();
		//annotateTranscriptDomains();
		annotateWithDavid();
	}
	
	public static void replaceDomainsByNumber() {
		// 109_PFAM_DOMAINS_sorted.tbl
		
		File inFile= new File("I:\\work\\claudia\\sarray2\\sets\\hg16_events_ref_genome_domEvent_clean.gtf"),
		outFile= new File("I:\\work\\claudia\\sarray2\\sets\\hg16_events_ref_genome_domEvent_clean_replaced.gtf"),
		annFile= new File("I:\\work\\claudia\\sarray2\\sets\\109_PFAM_DOMAINS_sorted.tbl");
	
		// read annotation database
		HashMap<String, String> mapAnn= new HashMap<String, String>();	//264222
		String[] colNames= null;
		try {
			System.err.print("reading annotation ");
			System.err.flush();
			BufferedReader buffy= new BufferedReader(new FileReader(annFile));
			int rowCtr= 1, perc= 0;
			long fsize= annFile.length(), bread= 0;
			for (String s= buffy.readLine(); s!= null; s= buffy.readLine(), ++rowCtr) {
				bread+= s.length()+1;
				if (bread*10d/fsize> perc) {
					++perc;
					System.err.print("*");
					System.err.flush();
				}
				String[] tokens= s.split("\t");
				assert(tokens.length==2);
				
				mapAnn.put(tokens[1], tokens[0]);
			}
			System.err.println();
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		
		try {
			System.err.print("annotating ");
			System.err.flush();
			BufferedWriter writer= new BufferedWriter(new FileWriter(outFile));
			GTFwrapper reader= new GTFwrapper(inFile.getAbsolutePath());
			reader.setReadGene(false);
			reader.setReadGTF(true);
			reader.setReadAheadLimit(1);
			reader.setSilent(false);
			reader.setReadFeatures(new String[] {"as_event", "ds_event"});	// GTFObject.FEATURE_ASEVENT
			GFFObject[] obj;
			int objCtr= 0;
			for (reader.read(); (obj= reader.getGtfObj())!= null; reader.read()) {
				for (int i = 0; i < obj.length; i++) {
					String s= (String) obj[i].getAttributes().remove(ATTRIBUTE_DOMAIN_EVENT);
					if (s!= null) {
						String[] sx= s.split("/");
						StringBuffer sb= new StringBuffer();
						for (int j = 0; j < sx.length; j++) {
							StringTokenizer toki= new StringTokenizer(sx[j], ",", true);
							while (toki.hasMoreTokens()) {
								String sxx= toki.nextToken();
								if (sxx.equals(","))
									sb.append(sxx);
								else {
									if (mapAnn.get(sxx)== null)
										System.err.println("NOT FOUND: "+sxx);
									else
										sb.append(mapAnn.get(sxx));
								}
							}
							
							if (j< sx.length-1)
								sb.append("/");
						}
						obj[i].getAttributes().put(ATTRIBUTE_DOMAIN_EVENT, sb.toString());
					}
					writer.write(obj[i].toString());
					writer.write(barna.commons.system.OSChecker.NEW_LINE);
				}
			}
			System.err.println();
			writer.flush();
			writer.close();
			
		} catch (Exception e) {
			e.printStackTrace();
		}
		
	}

	public static void annotateWithDavid() {
		File inFile= new File("M:\\work\\claudia\\sarray3\\source\\hg17_RefSeq-UCSC081117_mRNA-UCSC081117_intronEST-UCSC081117_realign_events_loc.gtf"),
			outFile= new File("M:\\work\\claudia\\sarray3\\source\\hg17_RefSeq-UCSC081117_mRNA-UCSC081117_intronEST-UCSC081117_realign_events_loc_david.gtf"),
			annFile= new File("M:\\work\\claudia\\david_genbank\\save\\annotation_all_108_108.tmp");
		
		// read annotation database
		HashMap<Int2Array,Int2Array> arraysSet= new HashMap<Int2Array,Int2Array>();	// 50000,1f
		HashMap<String, int[][]> mapAnn= new HashMap<String, int[][]>(265000, 1f);	//264222
		String[] colNames= null;
		try {
			System.err.print("reading annotation ");
			System.err.flush();
			BufferedReader buffy= new BufferedReader(new FileReader(annFile));
			int rowCtr= 1, perc= 0;
			long fsize= annFile.length(), bread= 0;
			for (String s= buffy.readLine(); s!= null; s= buffy.readLine(), ++rowCtr) {
				bread+= s.length()+1;
				if (bread*10d/fsize> perc) {
					++perc;
					System.err.print("*");
					System.err.flush();
				}
				String[] tokens= s.split("\t");
				int ep= s.length()-1;
				while (s.substring(ep,ep+1).equals("\t"))
					--ep; 
				if (ep< s.length()-1) {
					String[] tokSave= tokens;
					tokens= new String[tokSave.length+ (s.length()-ep-1)];
					System.arraycopy(tokSave,0,tokens,0,tokSave.length);
					for (int i = tokSave.length; i < tokens.length; i++) {
						tokens[i]= "";
					}
				}
				assert(tokens.length== 108);

				if (rowCtr== 1) {
					colNames= tokens;
					continue;
				}
				int[][] annot= new int[tokens.length-1][];
				for (int i = 0; i < annot.length; i++) {
					if (tokens[i+1].trim().length()== 0)
						continue;
					String[] ints= tokens[i+1].split(",\\s+");
					annot[i]= new int[ints.length];
					for (int j = 0; j < annot[i].length; j++) {
						try {
							annot[i][j]= Integer.parseInt(ints[j]);
						} catch (NumberFormatException e) {
							e.printStackTrace();
						}
					}
				}
				
				Int2Array arr= new Int2Array(annot);
				if (arraysSet.get(arr)!= null) {
					annot= arraysSet.get(arr).a;					
				} else {
					arraysSet.put(arr, arr);
				}
				mapAnn.put(tokens[0], annot);
//				if (rowCtr== 100)
//					break;
			}
			buffy.close();
			System.err.println();
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		System.err.println("Clean up");
		arraysSet= null;
		for (int i = 0; i < colNames.length; i++) {
			colNames[i]= colNames[i].replace(' ','_');
		}
		System.gc();
		try {
			Thread.currentThread().sleep(5000);
		} catch (InterruptedException e1) {
			e1.printStackTrace();
		}
		
		try {
			System.err.print("annotating ");
			System.err.flush();
			BufferedReader buffy= new BufferedReader(new FileReader(inFile));
			BufferedWriter writer= new BufferedWriter(new FileWriter(outFile));
			int objCtr= 0, perc= 0;
			long fsize= inFile.length(), bread= 0;
			for (String s= buffy.readLine(); s!= null; s= buffy.readLine()) {
				
//				System.gc();
//				Thread.currentThread().yield();
				
				bread+= s.length()+1;
				if (bread*10f/fsize> perc) {
					++perc;
					System.err.print("*");
					System.err.flush();
				}
				
				String[] ss= s.split("\\s+");
				String[] ids= ss[9].split("[/,]");
				
				StringBuffer sb= new StringBuffer();
				for (int j = 1; j < colNames.length; j++) {
					HashSet<Integer> littleHash= new HashSet<Integer>(); 
					for (int k = 0; k < ids.length; k++) { 
						int[][] annot= mapAnn.get(ids[k]);
						for (int m = 0; annot!=null&& annot[j-1]!= null&& m < annot[j-1].length; m++) {
							littleHash.add(annot[j-1][m]);
						}
					}
					if (littleHash.size()> 0) {
						sb.append(" ");
						sb.append(colNames[j]);
						sb.append(" \"");
						Integer[] ints= new Integer[littleHash.size()];
						littleHash.toArray(ints);
						Arrays.sort(ints);
						for (int k = 0; k < ints.length; k++) {
							sb.append(Integer.toString(ints[k]));
							if (k< ints.length-1)
								sb.append("/");
						}
						sb.append("\";");
					}
				}
				
				writer.write(s);
				if (sb.length()> 0) {
					writer.write(sb.toString());
					++objCtr;
				}
				writer.write(barna.commons.system.OSChecker.NEW_LINE);
				writer.flush();
				
//				if (objCtr> 10)
//					break;
			}
			writer.flush();
			writer.close();
			System.err.println();
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		if (1==1)
			return;
		
		
		// annotate
		try {
			System.err.print("annotating ");
			System.err.flush();
			BufferedWriter writer= new BufferedWriter(new FileWriter(outFile));
			GTFwrapper reader= new GTFwrapper(inFile.getAbsolutePath());
			reader.setReadGene(false);
			reader.setReadGTF(true);
			reader.setReadAheadLimit(1);
			reader.setSilent(false);
			reader.setReadFeatures(new String[] {"as_event", "ds_event"});	// GTFObject.FEATURE_ASEVENT
			GFFObject[] obj;
			int objCtr= 0;
			for (reader.read(); (obj= reader.getGtfObj())!= null; reader.read()) {
				for (int i = 0; i < obj.length; i++) {
					String s= obj[i].getAttribute(GFFObject.TRANSCRIPT_ID_TAG);
					String[] ids= s.split("[/,]");
					
					for (int j = 1; j < colNames.length; j++) {
						HashSet<Integer> littleHash= new HashSet<Integer>(); 
						for (int k = 0; k < ids.length; k++) { 
							int[][] annot= mapAnn.get(ids[k]);
							for (int m = 0; annot!=null&& annot[j-1]!= null&& m < annot[j-1].length; m++) {
								littleHash.add(annot[j-1][m]);
							}
						}
						if (littleHash.size()> 0) {
							Integer[] ints= new Integer[littleHash.size()];
							littleHash.toArray(ints);
							Arrays.sort(ints);
							StringBuffer sb= new StringBuffer();
							for (int k = 0; k < ints.length; k++) {
								sb.append(Integer.toString(ints[k]));
								if (k< ints.length-1)
									sb.append("/");
							}
							obj[i].addAttribute(colNames[j], sb.toString());
							++objCtr;
						}
					}
					
					writer.write(obj[i].toString());
					writer.write(barna.commons.system.OSChecker.NEW_LINE);
					
//					if (obj[i].getAttributes().size()> 20)
//						++objCtr;
					
				} // for all obj
				
//				if (objCtr> 10)
//					break;
			}
			writer.flush();
			writer.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public static final String ATTRIBUTE_DOMAIN_TRANSCRIPT= "domain_variant";
	public static void annotateTranscriptDomains() {
		
		File inFile= new File("I:\\work\\claudia\\sarray2\\sets\\hg16_events_ref_genome_domEvent_clean_replaced.gtf"),
		outFile= new File("I:\\work\\claudia\\sarray2\\sets\\hg16_events_ref_genome_domEvent_clean_replaced_tdomains.gtf"),
		annFile= new File("I:\\work\\claudia\\sarray2\\sets\\PFAM_DOMAINS.tmp");
	
		// read annotation database
		HashMap<String, int[]> mapAnn= new HashMap<String, int[]>();	//264222
		String[] colNames= null;
		try {
			System.err.print("reading annotation ");
			System.err.flush();
			BufferedReader buffy= new BufferedReader(new FileReader(annFile));
			int rowCtr= 1, perc= 0;
			long fsize= annFile.length(), bread= 0;
			for (String s= buffy.readLine(); s!= null; s= buffy.readLine(), ++rowCtr) {
				bread+= s.length()+1;
				if (bread*10d/fsize> perc) {
					++perc;
					System.err.print("*");
					System.err.flush();
				}
				String[] tokens= s.split("\t");
				assert(tokens.length==2);
				
				String[] ss= tokens[1].split("/");
				int[] ints= new int[ss.length];
				for (int i = 0; i < ints.length; i++) {
					ints[i]= Integer.parseInt(ss[i]);
				}
				
				mapAnn.put(tokens[0], ints);
			}
			System.err.println();
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		
		try {
			System.err.print("annotating ");
			System.err.flush();
			BufferedWriter writer= new BufferedWriter(new FileWriter(outFile));
			GTFwrapper reader= new GTFwrapper(inFile.getAbsolutePath());
			reader.setReadGene(false);
			reader.setReadGTF(true);
			reader.setReadAheadLimit(1);
			reader.setSilent(false);
			reader.setReadFeatures(new String[] {"as_event", "ds_event"});	// GTFObject.FEATURE_ASEVENT
			GFFObject[] obj;
			int objCtr= 0;
			for (reader.read(); (obj= reader.getGtfObj())!= null; reader.read()) {
				for (int i = 0; i < obj.length; i++) {
					String s= obj[i].getTranscriptID();
					
					
					String[] sx= s.split(",");
					StringBuffer sb= new StringBuffer();
					for (int j = 0; j < sx.length; j++) {
						HashSet<Integer> littleMap= new HashSet<Integer>();
						String[] ssx= sx[j].split("/");
						for (int k = 0; k < ssx.length; k++) {
							int[] xxx= mapAnn.get(ssx[k]);
							for (int x = 0; xxx!=null&& x < xxx.length; x++) 
								littleMap.add(xxx[x]);
						}
						
						Iterator<Integer> iter= littleMap.iterator();
						while(iter.hasNext()) {
							sb.append(iter.next());
							if (iter.hasNext())
								sb.append("/");
						}
						if (j< sx.length-1)
							sb.append(",");
					}
										
					if (sb.length()> 1)
						obj[i].addAttribute(ATTRIBUTE_DOMAIN_TRANSCRIPT, sb.toString());
					writer.write(obj[i].toString());
					writer.write(barna.commons.system.OSChecker.NEW_LINE);
				}
			}
			System.err.println();
			writer.flush();
			writer.close();
			
		} catch (Exception e) {
			e.printStackTrace();
		}
		
	}

	public static String ATTRIBUTE_DOMAIN_EVENT= "domain_event";
	
}
