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

package barna.model.splicegraph;

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
