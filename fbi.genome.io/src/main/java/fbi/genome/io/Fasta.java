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

package fbi.genome.io;

public class Fasta {

	public static byte getReadDir(String s) {
		//chr16:1500431-1600693C;uc002cma.1;9163;1;49
		int p= 0;
		for (int i = 0; i < 3; i++, p= s.indexOf(';', p+1));
		++p;
		byte dir= Byte.parseByte(s.substring(p, s.indexOf(';', p)));
		return dir;
	}
	
	public static String getReadID(String s) {
		//chr16:1500431-1600693C;uc002cma.1;9163;1;49
		int p= 0;
		for (int i = 0; i < 3; i++, p= s.indexOf(';', p+1));
		
		return s.substring(0,p);
	}

	public static final String QFASTA_PLUS = "+";
}
