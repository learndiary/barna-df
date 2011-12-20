package barna.genome.astalavista;

import barna.genome.io.gtf.GTFwrapper;
import barna.genome.model.Gene;
import barna.genome.model.splicegraph.SplicingGraph;

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
