/*
 * Copyright (c) 2010, Micha Sammeth
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * The names of its contributors may be not used to endorse or promote
 *       products derived from this software without specific prior written
 *       permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL MICHA SAMMETH BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

/*
 * Created on Nov 27, 2005
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
package barna.model;

import java.util.Comparator;
import java.util.HashMap;

/**
 * 
 * 
 * @author msammeth
 */
public abstract class AbstractRegion implements Region {

	String id= null;
	HashMap attributes= null;
	
	static final long serialVersionUID=  5443375142823871946L;
	public abstract Species getSpecies();

	public int getLength() {
		return (getEnd()- getStart()+ 1);
	}
	
	public boolean hasInvalidCoordinates() {
		if ((start== 0)|| end== 0|| start== Integer.MAX_VALUE|| end== Integer.MIN_VALUE|| 
				((start< 0)!= (end< 0)))
			return true;
		return false;
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
	
	public abstract String getChromosome();

	/* assumption: exons share same strand
	 * 
	 * @author micha
	 */
	public static class StartComparator implements Comparator {
		
		public int compare(Object o1, Object o2) {

			int start2= ((Region) o2).getStart();
			int start1= ((Region) o1).getStart();
			if (start1== start2)
				return 0;
			if (start1< start2)
				return -1;
			return 1;
		}
	}

	/* assumption: exons share same strand
		 * 
		 * @author micha
		 */
		public static class PositionComparator implements Comparator {
			
			public int compare(Object o1, Object o2) {
	
				AbstractRegion reg1= (DirectedRegion) o1;
				AbstractRegion reg2= (DirectedRegion) o2;
				
				if (reg1.getChromosome()!= null&& reg2.getChromosome()!= null) {
					if (!reg1.getChromosome().equalsIgnoreCase(reg2.getChromosome()))
						return (reg1.getChromosome().compareTo(reg2.getChromosome()));
				}

				int end1= ((Region) o1).getEnd();
				int start2= ((Region) o2).getStart();
				int end2= ((Region) o2).getEnd();
				int start1= ((Region) o1).getStart();
				if (start1== start2&& end1== end2)	// no object identity
					return 0;
				
					// non-overlapping, one before the other
				// cancelled, not working for neg. strand (clustering, sort array asc with start, end pos)
	//			if (end1< start2)
	//				return -1;		// one stops before the other
	//			if (end2< start1)
	//				return 1;
				
					// overlapping: none stops before the other
				if (start1< start2)
					return -1;
				if (start2< start1)
					return 1;
				
					// overlapping and same start position
				if (end1< end2)
					return -1;
				if (end2< end1)
					return 1;
				
				//System.err.println("assertion in abstractregion.positioncomparator failed");
				return 0;	// identical positions --> never reached
				
			}
		}

	/* assumption: exons share same strand
		 * 
		 * @author micha
		 */
		public static class EndComparator implements Comparator {
			
			public int compare(Object o1, Object o2) {
		
				int end2= ((Region) o2).getEnd();
				int end1= ((Region) o1).getEnd();
				if (end1== end2)
					return 0;
				if (end1< end2)
					return -1;
				return 1;
			}
		}

	public boolean overlaps(Region anotherRegion) {
		
		if (getChromosome()!= null&& anotherRegion.getChromosome()!= null
			&& (!anotherRegion.getChromosome().equalsIgnoreCase(getChromosome())))
			return false;
		
		if ((getStart()>= anotherRegion.getStart()&& getStart()< anotherRegion.getEnd())
				|| (anotherRegion.getStart()>= getStart()&& anotherRegion.getStart()< getEnd()))
			return true;
		return false;
	}
	
	/**
	 * Returns <code>true</code> if <code>this</code> region contains <code>anotherRegion</code>.
	 */
	public boolean contains(Region anotherRegion) {

		if (!anotherRegion.getChromosome().equalsIgnoreCase(getChromosome()))
			return false;
		
		if ((getStart()<= anotherRegion.getStart())&& (getEnd()>= anotherRegion.getEnd()))
			return true;
		return false;
	}

	int end = 0;	// Integer.MAX ???? removed
	int start = 0;
	static PositionComparator defaultPositionComparator= new PositionComparator();
	
	
	
	/**
	 * @return Returns the end.
	 */
	public int getEnd() {
		return end;
	}

	/**
	 * @return Returns the start.
	 */
	public int getStart() {
		return start;
	}

	/**
	 * @param end The end to set.
	 */
	public void setEnd(int end) {
		this.end = end;
	}

	/**
	 * @param start The start to set.
	 */
	public void setStart(int start) {
		this.start = start;
	}

	public String toString() {
		return "["+getStart()+";"+getEnd()+"]";
	}
	

	public String getID() {
		return id;
	}

	public void setID(String id) {
		this.id = id;
	}
	
	protected Object clone() throws CloneNotSupportedException {
		return super.clone();
	}

	public HashMap getAttributes() {
		return attributes;
	}

	public static Comparator getDefaultPositionComparator() {
		return defaultPositionComparator;
	}
	
}
