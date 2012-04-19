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
import barna.model.splicegraph.SplicingGraph;

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
			GTFwrapper reader= new GTFwrapper("P:\\annotation\\hg18\\hg18_EnsemblGenes_fromUCSC090615_sorted.gtf");
			reader.setReadGTF(true);
			
			Gene[] g;
			for (reader.read(); (g= reader.getGenes())!= null; reader.read()) {
				for (int i = 0; i < g.length; i++) {
					SplicingGraph gr= new SplicingGraph(g[i]);
					gr.constructGraph();
					gr.transformToFragmentGraph();
					gr.getVariations(minAmpl, maxAmpl, minPlen, maxPlen, minPovl, seqLen, minIntronLen);
				}
			}
			
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	
	
}
