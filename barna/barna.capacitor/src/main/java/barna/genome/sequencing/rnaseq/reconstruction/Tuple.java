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

import java.util.Comparator;

public class Tuple {

	public static class TupleByXComparator implements Comparator<Tuple> {
		public int compare(Tuple o1, Tuple o2) {
			return (o1.x- o2.x);
		}
	}
	
	public TupleByXComparator defaultTupleByXComparator= new TupleByXComparator();
	
	public int x, y;
	public Tuple(int x, int y) {
		this.x= x;
		this.y= y;
	}
}
