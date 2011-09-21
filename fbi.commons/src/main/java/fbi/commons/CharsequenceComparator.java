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

package fbi.commons;

import java.util.Comparator;

public class CharsequenceComparator implements Comparator<CharSequence> {

	public static final CharsequenceComparator DEFAULT_CHARSEQUENCE_COMPARATOR= new CharsequenceComparator();
	
	@Override
	public int compare(CharSequence paramT1, CharSequence paramT2) {
		
		int len= Math.min(paramT1.length(), paramT2.length());
		for (int i = 0; i < len; ++i) {
			char c1= paramT1.charAt(i), 
				c2= paramT2.charAt(i);
			if (c1!= c2)
				return (c1-c2);
		}
				
		return paramT1.length()- paramT2.length();
	}
}
