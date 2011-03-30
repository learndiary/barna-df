package fbi.genome.sequencing.rnaseq.reconstruction;

import fbi.genome.model.splicegraph.Edge;
import fbi.genome.model.splicegraph.Graph;

import java.util.Vector;

public interface ReadStatCalculator {

	double getReads(Vector<Edge> v, byte dir, long[] sigExcl, boolean normalized);
	
	double getReadsAvg(Vector<Edge> v, byte dir, Graph g, long[] sig, boolean exclusive, boolean normalized);
}
