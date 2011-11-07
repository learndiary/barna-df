package fbi.genome.sequencing.rnaseq.graph;

import fbi.genome.model.splicegraph.Edge;
import fbi.genome.model.splicegraph.Node;

/**
 * Extends default edges by coverage arrays.
 * 
 * @author Micha Sammeth (gmicha@gmail.com)
 * 
 */
public class EdgeNt extends Edge {

	/**
	 * Coverage of genomic positions.
	 */
	int[] coverage= null;
	
	/**
	 * Creates edge and initializes coverage array.
	 * @param src the source node
	 * @param snk the sink node
	 * @overrides
	 */
	public EdgeNt(Node src, Node snk) {
		super(src, snk);
		coverage= new int[length()];
	}
}
