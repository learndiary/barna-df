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

/**
 * 
 * 
 * @author msammeth
 */
public class DefaultRegion extends AbstractRegion implements Region {

	static long serialVersionUID= 6558964245030911377l;
	String chromosome = null;
	Species species = null;
	double score= 0d;
	
	public boolean contains(SpliceSite ss) {
		if (!ss.getGene().getChromosome().equals(chromosome))
			return false;
		
		int apos= Math.abs(ss.getPos());
		if (apos>= getStart()&& apos<= getEnd())
			return true;
		return false;
	}
	
	public Object clone() {
		DefaultRegion clone= new DefaultRegion(getStart(), getEnd());
		clone.setID(getID());
		clone.setChromosome(getChromosome());
		clone.setSpecies(getSpecies());
		return clone;
	}
	
	public DefaultRegion(int newStart, int newEnd) {
		setStart(newStart);
		setEnd(newEnd);
	}
	
	public DefaultRegion() {
	}
	
	/**
	 * @return Returns the chromosome.
	 */
	public String getChromosome() {
		return chromosome;
	}

	/**
	 * @param chromosome The chromosome to set.
	 */
	public void setChromosome(String chromosome) {
		this.chromosome = chromosome;
	}

	/**
	 * @return Returns the species.
	 */
	public Species getSpecies() {
		return species;
	}
	/**
	 * @param species The species to set.
	 */
	public void setSpecies(Species species) {
		this.species = species;
	}
	
	public String toString() {
		String str= super.toString();
		if (getID()!= null)
			str= getID()+ " "+ str;
		if (getChromosome()!= null)
			str= getChromosome()+ ":"+ str;
		if (getSpecies()!= null)
			str= getSpecies().getCommonName()+ "-"+ str;
		if (!Double.toString(getScore()).equals(Double.toString(Double.NaN)))
			str+= ": "+getScore();
			
		return str;
	}

	public double getScore() {
		return score;
	}

	public void setScore(double score) {
		this.score = score;
	}
}
