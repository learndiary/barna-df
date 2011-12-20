/*
 * This file is part of the Flux Library.
 *
 * The code of the Flux Library may be freely distributed and modified under the terms of the
 * European Union Public Licence (EUPL) published on the web site <http://www.osor.eu/eupl/european-union-public-licence-eupl-v.1.1>.
 * Copyright for the code is held jointly by the individual authors, who should be listed
 * in @author doc comments. According to Article 5 and Article 11 of the EUPL, publications that
 * include results produced by the Flux Library are liable to reference the Work,
 * see the Flux Library homepage <http://flux.sammeth.net> for more information.
 */

package barna.genome.model.splicegraph;

public class TestPEmapping {
		static boolean output= false;

//		public static void main(String[] args) {
//
//			int[] len= new int[] {128, 178};	//128, 178=  [200;250] - 2* readLen
//			analyzeSplicingGraph(len, false);
//		}
//
//		public static void analyzeSplicingGraph(int[] len, boolean exonFragment) {
//			long t0= System.currentTimeMillis();
//			int lineCtr= 0;
//			try {
//				// /home/msammeth/annotations/hg18_RefSeq_mRNA_ESTs_0710_sorted.gtf
//				// /home/msammeth/annotations/hg18_knownGenes_fromUCSC081019_sorted.gtf
//				GFFReader reader= new GFFReader("/home/msammeth/annotations/hg18_RefSeq_mRNA_ESTs_0710_sorted.gtf");
//				reader.read();
//				Gene[] g= reader.getGenes();
//				BufferedReader buffy= new BufferedReader(new FileReader("/home/msammeth/solexa/pends/moreShit_sorted.txt"));
//				BufferedWriter writer= new BufferedWriter(new FileWriter("/home/msammeth/solexa/pends/moreShit_sorted_out.txt"));
//				String line= buffy.readLine();
//				++lineCtr;
//				while(g!= null&& line!= null) {
//					String[] toki= line.split("\\s");
//					if (toki.length!= 3)
//						System.err.println("line "+lineCtr+" with "+toki.length+
//								" <> 3 fields.");
//					int[] pos= new int[2];
//					pos[0]= Integer.parseInt(toki[1]);
//					pos[1]= Integer.parseInt(toki[2]);
//
//					int last= -1;
//					Graph gr= null;
//					for (int i = 0; i < g.length; i+=(pos[0]>= g[i].getStart()&& pos[0]<= g[i].getEnd())?0:1) {
//
//						if (g[i].getStrand()< 0)
//							continue;
//
//						while (line== null|| pos[0]< g[i].get5PrimeEdge()|| toki[0].compareTo(g[i].getChromosome())< 0) {
//							if (output)
//								System.out.println("skipped line "+line+" w gene "+g[i].getChromosome()+" "+g[i].getStart()+" "+g[i].getEnd());
//							line= buffy.readLine();
//							++lineCtr;
//							toki= line.split("\\s");
//							if (toki.length!= 3)
//								System.err.println("line "+lineCtr+" with "+toki.length+
//										" <> 3 fields.");
//							if (toki[0].compareTo(g[i].getChromosome())> 0)
//								break;
//							pos= new int[2];
//							pos[0]= Integer.parseInt(toki[1]);
//							pos[1]= Integer.parseInt(toki[2]);
//						}
//						if (toki[0].compareTo(g[i].getChromosome())> 0|| pos[0]> g[i].get3PrimeEdge()) {
//							if (output)
//								System.out.println("skipped gene "+g[i].getChromosome()+" "+g[i].getStart()+" "+g[i].getEnd()+ "w line "+line);
//							continue;
//						}
//
//						if (last< i) {
//							gr= new Graph(g[i]);
//							gr.constructGraph();
//							if (exonFragment)
//								gr.transformToFragmentGraph();
//						}
//						last= i;
//
//						if (pos[1]> g[i].getStart()&& pos[1]< g[i].getEnd()) {
//							String s= gr.getPathes(pos[0], pos[1], len[0], len[1]);
//							if (output)
//								System.out.println("MAPPED line "+line+" w gene "+g[i].getChromosome()+" "+g[i].getStart()+" "+g[i].getEnd()+" : "+s);
////							if (s.indexOf("TRIED_ALL")>= 0)
////								System.currentTimeMillis();
//							writer.write(s);
//							writer.write('\n');
//						} else if (output)
//							System.out.println("DROPPED line "+line+" w gene "+g[i].getChromosome()+" "+g[i].getStart()+" "+g[i].getEnd());
//
//						line= buffy.readLine();
//						if (line== null)
//							break;
//						++lineCtr;
//						toki= line.split("\\s");
//						if (toki.length!= 3)
//							System.err.println("line "+lineCtr+" with "+toki.length+
//									" <> 3 fields.");
//						pos= new int[2];
//						pos[0]= Integer.parseInt(toki[1]);
//						pos[1]= Integer.parseInt(toki[2]);
//					}
//
//					reader.read();
//					g= reader.getGenes();
//				}
//				writer.flush();
//				writer.close();
//				buffy.close();
//			} catch (Exception e) {
//				System.err.println(lineCtr+ " lines.");
//				e.printStackTrace();
//			}
//
//			System.err.println(lineCtr+" lines, took "+((System.currentTimeMillis()- t0)/1000)+" sec.");
//
//		}
	
	
}
