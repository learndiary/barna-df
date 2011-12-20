/*
 * This file is part of the Flux Library.
 *
 * The code of the Flux Library may be freely distributed and modified under the terms of the
 * European Union Public Licence (EUPL) published on the web site <http://www.osor.eu/eupl/european-union-public-licence-eupl-v.1.1>.
 * Copyright for the code is held jointly by the individual authors, who should be listed
 * in @author doc comments. According to Article 5 and Article 11 of the EUPL, publications that
 * include results produced by the Flux Library are liable to reference the Work,
 * see the Flux Library homepage <http://flux.sammeth.net> for more information.
 */

package barna.genome.sequencing.rnaseq.reconstruction;

import barna.genome.model.splicegraph.AbstractEdge;
import barna.genome.model.splicegraph.SplicingGraph;

import java.util.Vector;

public interface ReadStatCalculator {

	double getReads(Vector<AbstractEdge> v, byte dir, long[] sigExcl, boolean normalized);
	
	double getReadsAvg(Vector<AbstractEdge> v, byte dir, SplicingGraph g, long[] sig, boolean exclusive, boolean normalized);
}
