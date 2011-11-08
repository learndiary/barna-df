import fbi.genome.io.FileHelper;
import fbi.genome.io.gtf.GTFwrapper;
import fbi.genome.model.Gene;
import fbi.genome.model.Transcript;
import fbi.genome.model.splicegraph.SimpleEdge;
import fbi.genome.model.splicegraph.SplicingGraph;
import fbi.genome.model.splicegraph.Node;

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
