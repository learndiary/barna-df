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

package barna.astalavista;

import barna.io.gtf.GTFwrapper;
import barna.model.Gene;
import barna.model.Transcript;
import barna.model.Translation;
import barna.model.gff.GFFObject;
import barna.model.splicegraph.SimpleEdge;
import barna.model.splicegraph.SplicingGraph;

import java.util.HashMap;

public class PrimerDesigner {

	public static void main(String[] args) {
		method001();
	}
	
	static void method001() {
		
		int minAmpl= 150, maxAmpl= 200, 
			minPlen= 18, maxPlen= 18,
			minPovl= -1, seqLen= 35;	//minCommon= primer
		int minIntronLen= 200; // because of artifacts 
		
		try {
            // /users/rg/projects/encode/scaling_up/whole_genome/mouse/Annotation/ensembl65/Long/Element/mm65.long.exon.gtf

			GTFwrapper reader= new GTFwrapper("");
			reader.setReadGTF(true);
			
			Gene[] g;
			for (reader.read(); (g= reader.getGenes())!= null; reader.read()) {
				for (int i = 0; i < g.length; i++) {
					SplicingGraph gr= new SplicingGraph(g[i]);
					gr.constructGraph();
					gr.transformToFragmentGraph();
					//gr.getVariations(minAmpl, maxAmpl, minPlen, maxPlen, minPovl, seqLen, minIntronLen);

                   HashMap<long[],Integer> map= get(gr);
                   int minStart= Integer.MAX_VALUE, maxStart= Integer.MIN_VALUE,
                           minEnd= Integer.MAX_VALUE, maxEnd= Integer.MIN_VALUE;
                    for (int j = 0; j < gr.trpts.length; j++) {
                        int start= gr.trpts[i].get5PrimeEdge(),
                                end= gr.trpts[i].get3PrimeEdge();
                        if(start< minStart)
                            minStart= start;
                        if (start> maxStart)
                            maxStart= start;
                        if (end< minEnd)
                            minEnd= end;
                        if (end> maxEnd)
                            maxEnd= end;
                    }

                    GFFObject obj= GFFObject.createGFFObject(g[i]);
                    StringBuilder sb= new StringBuilder(obj.toString());
                    sb.append(" exon_nr \""+ g[i].getExons().length+ "\";");
                    sb.append(" tx_nr \""+ g[i].getTranscriptCount()+ "\";");

                    sb.append(" pct \"");
                    for (int j = 0; j < gr.trpts.length; j++) {
                        if (gr.trpts[i].isCoding())
                            sb.append(gr.trpts[i].getTranscriptID()+ ",");
                    }
                    sb.deleteCharAt(sb.length()- 1);
                    sb.append("\";");

                    sb.append(" fuzzy \""+ ((maxStart- minStart)+ (maxEnd- minEnd)+ "\";"));

                    for (long[] longs : map.keySet()) {
                        Transcript[] tx= gr.decodeTset(longs);
                        StringBuilder id= new StringBuilder();
                        for (int j = 0; j < tx.length; j++)
                            id.append(tx[j].getTranscriptID()+ "_");
                        id.deleteCharAt(id.length()- 1);
                        sb.append(" "+ id+ " \""+ map.get(longs)+ "\";");
                    }

                    System.err.println(sb);
                }
			}
			
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

    private static HashMap<long[], Integer> get(SplicingGraph gr) {
        // output features of the full graph
        HashMap<long[], Integer> map= new HashMap<long[], Integer>(gr.trpts.length* gr.trpts.length);
        for (SimpleEdge e : gr.getEdgeHash().values()) {
            if (!e.isExonic())
                continue;
            long[] sig= e.getTranscripts();
            Transcript[] tx= gr.decodeTset(sig);
            boolean skip= false;
            if (tx.length== gr.trpts.length)
                continue;
            for (int i = 0; i < tx.length; i++) {
                Transcript transcript = tx[i];
                if (transcript.isCoding()) {
                    if (e.getTail().getSite().getPos()< transcript.getTranslations()[0].get5PrimeEdge()
                            || e.getTail().getSite().getPos()> transcript.getTranslations()[0].get3PrimeEdge())
                    skip= true;
                }
            }
            if (skip)
                continue;

            int len= e.getHead().getSite().getPos()-
                    e.getTail().getSite().getPos();
            if ((!map.containsKey(sig))|| map.get(sig)< len)
                map.put(sig, len);
        }
        return map;  //To change body of created methods use File | Settings | File Templates.
    }


}
