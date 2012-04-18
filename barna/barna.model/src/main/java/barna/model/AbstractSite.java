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

/*
 * Created on Mar 12, 2006
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
package barna.model;

import java.io.Serializable;
import java.util.Comparator;
import java.util.HashMap;

/**
 * 
 * 
 * @author msammeth
 */
public class AbstractSite implements Serializable {
	
	int pos= -1;
	Transcript[] transcripts = null;
	HashMap attributes = null;
	String id = null;
	float score= Float.NaN;
	Gene gene= null;
	static final long serialVersionUID = 3169139368723074072L;
	public static class PositionToSpliceSiteComparator implements Comparator {

		public int compare(Object arg0, Object arg1) {
			
			AbstractSite s= null;
			Integer i= null;
			try {
				s= (AbstractSite) arg0;
				i= (Integer) arg1;
			} catch (ClassCastException e) {
				try {
					i= (Integer) arg0;
					s= (AbstractSite) arg1;
				} catch (ClassCastException ex) {
					ex.printStackTrace();
					return -1;
				}
			}
			
			if (s.getPos()< i.intValue())
				return -1;
			if (s.getPos()> i.intValue())
				return 1;
			return 0;
		}
		
	}
	
	public static class PositionComparator implements Comparator {
	
		public int compare(Object arg0, Object arg1) {
			
			AbstractSite s1= null, s2= null;
			try {
				s1= (AbstractSite) arg0;
				s2= (AbstractSite) arg1;
			} catch (ClassCastException e) {
				e.printStackTrace();
				return -1;
			}
			
			if (s1.getPos()< s2.getPos())
				return -1;
			if (s1.getPos()> s2.getPos())
				return 1;
			return 0;
		}
		
	}

	public Transcript[] getTranscripts() {
		return transcripts;
	}
	
	public boolean equals(Object obj) {
		if (obj.getClass()!= getClass())	// for bi-symetrical a.equals(b) comparison with descendant classes
			return false;
		
		AbstractSite as= (AbstractSite) obj;	// no class cast
		
			// do not check transcripts here, for node comparison in Graph
		if (pos!= as.getPos())
			return false;
		return true;
	}
	
//	---------------------------------------------- 
//	hashCode 
//	public int hashCode() 
//	Returns a hash code value for the object. This method is supported for the 
//	benefit of hashtables such as those provided by java.util.Hashtable. 
//	The general contract of hashCode is: 
//	* Whenever it is invoked on the same object more than once during an 
//	execution of a Java application, the hashCode method must consistently 
//	return the same integer, provided no information used in equals comparisons 
//	on the object is modified. This integer need not remain consistent from one 
//	execution of an application to another execution of the same application. 
//	* If two objects are equal according to the equals(Object) method, then 
//	calling the hashCode method on each of the two objects must produce the same 
//	integer result. 
//	* It is not required that if two objects are unequal according to the 
//	equals(java.lang.Object) method, then calling the hashCode method on each of 
//	the two objects must produce distinct integer results. However, the 
//	programmer should be aware that producing distinct integer results for 
//	unequal objects may improve the performance of hashtables. 
//
//
//	As much as is reasonably practical, the hashCode method defined by class 
//	Object does return distinct integers for distinct objects. (This is 
//	typically implemented by converting the internal address of the object into 
//	an integer, but this implementation technique is not required by the JavaTM 
//	programming language.) 
//	----------------------------------------------
	public int hashCode() {
		return pos;
	}
	
	public AbstractSite(int newPos) {
		this.pos= newPos;
	}
	public String toString() {
		return Integer.toString(getPos());
	}
	
	public int getPos() {
		return pos;
	}
	public void setPos(int pos) {
		this.pos = pos;
	}
	public boolean isTSS() {
		for (int i = 0; i < getTranscripts().length; i++) {
			if (getPos()== getTranscripts()[i].get5PrimeEdge())
				return true;
		}
		return false;
	}
	
	public void setTranscripts(Transcript[] transcripts) {
		this.transcripts = transcripts;
	}
	public void addAttribute(Object id, Object val) {
		if (attributes== null)
			attributes= new HashMap();
		attributes.put(id, val);
	}

	public Object getAttribute(Object id) {
		if (attributes== null)
			return null;
		return attributes.get(id);
	}

	public HashMap getAttributes() {
		return attributes;
	}

	public String getID() {
		return id;
	}

	public void setID(String id) {
		this.id = id;
	}

	public float getScore() {
		return score;
	}

	public void setScore(float score) {
		this.score = score;
	}

	public Gene getGene() {
		return gene;
	}

	public void setGene(Gene gene) {
		this.gene = gene;
	}
}
