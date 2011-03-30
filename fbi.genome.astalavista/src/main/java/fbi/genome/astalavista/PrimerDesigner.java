package fbi.genome.astalavista;

import fbi.genome.io.gff.GFFReader;
import fbi.genome.model.Gene;
import fbi.genome.model.splicegraph.Graph;

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
			GFFReader reader= new GFFReader("P:\\annotation\\hg18\\hg18_EnsemblGenes_fromUCSC090615_sorted.gtf");
			reader.setReadGTF(true);
			
			Gene[] g;
			for (reader.read(); (g= reader.getGenes())!= null; reader.read()) {
				for (int i = 0; i < g.length; i++) {
					Graph gr= new Graph(g[i]);
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
