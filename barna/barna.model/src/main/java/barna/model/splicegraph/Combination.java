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

package barna.model.splicegraph;

import java.util.Arrays;
import java.util.Comparator;

public class Combination {
	public static class CombiSorter implements Comparator<long[]> {
		public int compare(long[] o1, long[] o2) {
			for (int i = 0; i < o1.length; i++) {
				if (o1[i]< o2[i])
					return -1;
				else if (o2[i]< o1[i])
					return 1;
			}
			return 0;
		}
	}
	
	static CombiSorter defaultCombiSorter= new CombiSorter();
	
	long[][] combi;
	
	public Combination(long[][] newCombi) {
		this.combi= newCombi;
		Arrays.sort(this.combi, defaultCombiSorter);
	}
	
	
}
