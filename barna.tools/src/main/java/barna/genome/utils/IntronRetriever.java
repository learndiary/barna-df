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

import barna.io.FileHelper;
import barna.io.gtf.GTFwrapper;
import barna.model.Gene;
import barna.model.Transcript;
import barna.model.splicegraph.Node;
import barna.model.splicegraph.SimpleEdge;
import barna.model.splicegraph.SplicingGraph;

import java.io.File;
import java.io.FileOutputStream;
import java.io.PrintStream;

/**
 * builds up splice graph, outputs introns
 * @author msammeth
 *
 */
public class IntronRetriever {

	public static void main(String[] args) {
		if (args== null|| args.length!= 1)
			System.err.println("Usage: IntronRetriever <GTF file>");
		
		// hg19_splicedESTs_UCSC100525.gtf
		File f= new File(args[0]);
		GTFwrapper reader= new GTFwrapper(f.getAbsolutePath());
		if (!reader.isApplicable()) {
			File tmpGTF= reader.sort();
			f= new File(args[0]+ "_sorted.gtf");
			FileHelper.move(tmpGTF, f);
			reader= new GTFwrapper(f.getAbsolutePath());
		}
		Gene[] genes= null;
		try {			
			PrintStream outF= new PrintStream(new FileOutputStream(args[0]+"_introns.gtf"));
			for(reader.read(); (genes= reader.getGenes())!= null; reader.read()) {
				for (int i = 0; i < genes.length; i++) {
					SplicingGraph g= new SplicingGraph(genes[i]);
					g.constructGraph();
					outputIntrons(g, outF);
				}
			}
			outF.flush();
			outF.close();
			System.err.println("output in "+ args[0]+"_introns.gtf");
		} catch (Exception e) {
			e.printStackTrace();
		}
		
	}
	
	public static void outputIntrons(SplicingGraph g, PrintStream p) {
		
		Node[] n= g.getNodesInGenomicOrder();
		for (int i = 0; i < n.length; i++) {
			for (int j = 0; n[i].getOutEdges()!= null&& j < n[i].getOutEdges().size(); j++) {
				SimpleEdge e= n[i].getOutEdges().elementAt(j);
				if (!e.isIntronic())
					continue;
				StringBuilder sb= new StringBuilder(g.trpts[0].getChromosome());
				sb.append("\t");
				sb.append(g.trpts[0].getSource());
				sb.append("\tintron\t");
				if (g.trpts[0].getStrand()>= 0) {
					sb.append(Integer.toString(e.getTail().getSite().getPos()));
					sb.append("\t");
					sb.append(Integer.toString(e.getHead().getSite().getPos()));
				} else {
					sb.append(Integer.toString(Math.abs(e.getHead().getSite().getPos())));
					sb.append("\t");
					sb.append(Integer.toString(Math.abs(e.getTail().getSite().getPos())));
				}
				sb.append("\t");
				sb.append(Integer.toString(g.getTranscriptNb(e.getTranscripts())));
				sb.append("\t");
				if (g.trpts[0].getStrand()> 0) 
					sb.append("+");
				else if (g.trpts[0].getStrand()< 0)
					sb.append("-");
				else
					sb.append(".");
				sb.append("\t.\t");	// TODO get phase from annotation
				sb.append("transcript_id \"");
				Transcript[] tt= g.decodeTset(e.getTranscripts());
				for (int k = 0; k < tt.length; k++) {
					sb.append(tt[k].getTranscriptID());
					sb.append("/");
				}
				sb.deleteCharAt(sb.length()- 1);
				sb.append("\"; gene_id \"");
				sb.append(g.trpts[0].getGene().getGeneID());
				sb.append("\";");
				p.println(sb.toString());
			}
		}
	}
	
}
