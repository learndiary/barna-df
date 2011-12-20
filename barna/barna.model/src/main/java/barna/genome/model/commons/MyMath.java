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

package barna.genome.model.commons;

public class MyMath {	
	
	public static long pow(int a, int b) {
		if (b== 0)
			return 1;
		if (b== 1)
			return a;
		long c= a;
		for (int i = 0; i < b-1; i++) 
			c*= a;
		return c;
	}
	
	public static long powSafe(int a, int b) throws IllegalArgumentException {
		if (b== 0)
			return 1;
		if (b== 1)
			return a;
		long c= a;
		for (int i = 0; i < b-1; i++) {
			if (c*a> Long.MAX_VALUE)
				throw new IllegalArgumentException("Exceeding long range!");
			c*= a;
		}
		return c;
	}

}
