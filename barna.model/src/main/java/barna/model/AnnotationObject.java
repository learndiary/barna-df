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

package barna.model;

import java.util.Hashtable;

public class AnnotationObject {
	public static String FIELD_DOT= ".";
	public static String FIELD_DASH= "-";
	
	static Hashtable<String,String> map;	// TODO is there not a solution without values?
	static protected String getTerm(String term) {
		String s= getMap().get(term);
		if (s== null)
			getMap().put(term,term);
		else
			return s;
		return term;
	}
	
	static protected Hashtable<String,String> getMap() {
		if (map == null) {
			map = new Hashtable<String,String>();
		}

		return map;
	}
}