package fbi.genome.model;

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
