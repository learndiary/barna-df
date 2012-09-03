package barna.genome.utils;

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

import barna.astalavista.statistics.David;
import barna.model.commons.DoubleVector;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Vector;


public class Analyzes {

	
	//5islets_woID2__6tissues.tx.inc-exc.txt
	public static void transcriptIncExc_101005() {
		try {
			// 101005: /Users/micha/projects/islets/6x6tissues/5islets_woID2__6tissues.tx.inc-exc.txt
			// 101018: /Users/micha/projects/islets/flux/hg19_RefSeq/5islets-ctl__6islets-cyt.txt
			//5ctl_5cyt_io.txt
			BufferedReader buffy= new BufferedReader(new FileReader("/Users/micha/projects/islets/flux/hg19_RefSeq/5ctl_5cyt_io.txt"));
			String line= null;
			Vector<String> idv= new Vector<String>();
			Vector<double[][]> v= new Vector<double[][]>();
			while((line= buffy.readLine())!= null) {
				String[] ss= line.split("\\s");
				idv.add(ss[0]);
				double[][] m= new double[11][];
				for (int i = 0; i < m.length; i++) {
					m[i]= new double[11];
				}
				for (int i = 1; i <= 20; i+= 2) {
					for (int j = i+ 2; j <= 20; j+=2) {

                        David.FisherExactTest fisher= new David.FisherExactTest();
						fisher.fisherExact(
								Math.round(Float.parseFloat(ss[i])),
								Math.round(Float.parseFloat(ss[i+1])),
								Math.round(Float.parseFloat(ss[j])),
								Math.round(Float.parseFloat(ss[j+1]))
						);
						double p= fisher.getTwotail();
						if (p> 1)
							System.currentTimeMillis();
						m[(i-1)/2][(j-1)/2]= p;
						m[(j-1)/2][(i-1)/2]= p;
					}
				}
				v.add(m);
			}
			buffy.close();
			
			// correct
/*			for (int i = 1; i <= 10; i++) {
				for (int j = i+1; j <= 10; j++) {
					double[] dd= new double[v.size()];
					for (int k = 0; k < dd.length; k++) 
						dd[k]= v.elementAt(k)[i-1][j-1];
					correctBH(dd);
					for (int k = 0; k < dd.length; k++) {
						v.elementAt(k)[i-1][j-1]= dd[k];
						v.elementAt(k)[j-1][i-1]= dd[k];
					}
				}
			}
*/			
			// output
			// 101005: "/Users/micha/projects/islets/6x6tissues/5islets_woID2__6tissues.tx.inc-exc.pval."+i+".txt"
			// 101018: "/Users/micha/projects/islets/flux/hg19_RefSeq/5islets-ctl__6islets-cyt"+i+".txt"
			for (int i = 1; i <= 10; i++) {
				BufferedWriter writer= new BufferedWriter(new FileWriter("/Users/micha/projects/islets/flux/hg19_RefSeq/5ctl_5cyt_io.pval."+i+".txt"));
				for (int x = 0; x < v.size(); x++) {
					StringBuilder sb= new StringBuilder(idv.elementAt(x));
					for (int j = 1; j <= 10; j++) {
						if (i== j)
							continue;
						sb.append("\t");
						sb.append(v.elementAt(x)[i-1][j-1]);
					}
					sb.append(barna.commons.system.OSChecker.NEW_LINE);
					writer.write(sb.toString());
				}
				writer.flush();
				writer.close();
			}
			
			
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	public static double[] correctBH(double[] d) {
		double[] dd= d.clone();
		Arrays.sort(dd);
		int n= dd.length;
		HashMap<Double, Double> map= new HashMap<Double, Double>();
		for (int i = dd.length- 1; i >= 0; --i) {
			if (map.containsKey(dd[i]))
				continue;
			map.put(dd[i],n* dd[i]/(i+ 1d));	// can be >1 (see http://www.silicongenetics.com/Support/GeneSpring/GSnotes/analysis_guides/mtc.pdf
		}
		
		for (int i = 0; i < dd.length; i++) {
			d[i]= map.get(d[i]);
		}
		return d;
	}
	
	public static void main(String[] args) {
		transcriptIncExc_101005();
	}

	public static void specificGeneExpression_100907() {
		try {
			BufferedReader buffy= new BufferedReader(new FileReader("/Users/micha/projects/islets/6x6tissues/brain_liver_contingency.txt"));
			String line= null;
			Vector<String> idv= new Vector<String>();
			DoubleVector v= new DoubleVector();
			while((line= buffy.readLine())!= null) {
				String[] ss= line.split("\\s");
				idv.add(ss[0]);
				double[][] m= new double[12][];
				for (int i = 0; i < m.length; i++) {
					m[i]= new double[12];
				}
				David.FisherExactTest fisher= new David.FisherExactTest();
				fisher.fisherExact(
						Math.round(Float.parseFloat(ss[1])),
						Math.round(Float.parseFloat(ss[2])),
						Math.round(Float.parseFloat(ss[3])),
						Math.round(Float.parseFloat(ss[4]))
				);
				double p= fisher.getTwotail();
				v.add(p);
			}
			buffy.close();
			
			// correct
			correctBH(v.toDoubleArray());
			
			// output
			BufferedWriter writer= new BufferedWriter(new FileWriter("/Users/micha/projects/islets/6x6tissues/brain_liver_contingency.pval.txt"));
			for (int x = 0; x < v.size(); x++) {
				StringBuilder sb= new StringBuilder(idv.elementAt(x));
				sb.append("\t");
				sb.append(v.getValue(x));
				sb.append(barna.commons.system.OSChecker.NEW_LINE);
				writer.write(sb.toString());
			}
			writer.flush();
			writer.close();
			
			
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	//5islets_woID2__6tissues.tx.inc-exc.txt
	public static void specificGeneExpression_100901() {
		try {
			BufferedReader buffy= new BufferedReader(new FileReader("/Users/micha/projects/islets/6x6tissues/6islets_6tissues.locusexpr.contingency.txt"));
			String line= null;
			Vector<String> idv= new Vector<String>();
			Vector<double[][]> v= new Vector<double[][]>();
			while((line= buffy.readLine())!= null) {
				String[] ss= line.split("\\s");
				idv.add(ss[0]);
				double[][] m= new double[12][];
				for (int i = 0; i < m.length; i++) {
					m[i]= new double[12];
				}
				for (int i = 1; i <= 12; i++) {
					for (int j = i+1; j <= 12; j++) {
						David.FisherExactTest fisher= new David.FisherExactTest();
						fisher.fisherExact(
								Math.round(Float.parseFloat(ss[i+i-1])),
								Math.round(Float.parseFloat(ss[i+i])),
								Math.round(Float.parseFloat(ss[j+j-1])),
								Math.round(Float.parseFloat(ss[j+j]))
						);
						double p= fisher.getTwotail();
						m[i-1][j-1]= p;
						m[j-1][i-1]= p;
					}
				}
				v.add(m);
			}
			buffy.close();
			
			// correct
			for (int i = 1; i <= 12; i++) {
				for (int j = i+1; j <= 12; j++) {
					double[] dd= new double[v.size()];
					for (int k = 0; k < dd.length; k++) 
						dd[k]= v.elementAt(k)[i-1][j-1];
					correctBH(dd);
					for (int k = 0; k < dd.length; k++) {
						v.elementAt(k)[i-1][j-1]= dd[k];
						v.elementAt(k)[j-1][i-1]= dd[k];
					}
				}
			}
			
			// output
			for (int i = 1; i <= 12; i++) {
				BufferedWriter writer= new BufferedWriter(new FileWriter("/Users/micha/projects/islets/6x6tissues/6islets_6tissues.locusexpr.contingency.pval."+i+".txt"));
				for (int x = 0; x < v.size(); x++) {
					StringBuilder sb= new StringBuilder(idv.elementAt(x));
					for (int j = 1; j <= 12; j++) {
						if (i== j)
							continue;
						sb.append("\t");
						sb.append(v.elementAt(x)[i-1][j-1]);
					}
					sb.append(barna.commons.system.OSChecker.NEW_LINE);
					writer.write(sb.toString());
				}
				writer.flush();
				writer.close();
			}
			
			
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
