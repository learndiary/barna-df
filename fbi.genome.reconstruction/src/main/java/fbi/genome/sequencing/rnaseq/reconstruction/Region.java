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

package fbi.genome.sequencing.rnaseq.reconstruction;

import java.util.Comparator;
import java.util.Vector;

public class Region {
	
	public static class PositionComparator implements Comparator<Region> {
		//@Override
		public int compare(Region o1, Region o2) {
			if (o1.exStart< o2.exStart)
				return -1;
			if (o1.exStart> o2.exStart)
				return 1;
			if (o1.exEnd< o2.exEnd)
				return -1;
			if (o1.exEnd> o2.exEnd)
				return 1;
			return 0;
		}
	}
	
	public static PositionComparator defaultPositionComparator= new PositionComparator(); 
	
	public int exStart, exEnd;
	Vector<Variation> inVarVec= new Vector<Variation>(), outVarVec= new Vector<Variation>();
}
