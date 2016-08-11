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
 */ /*
 * Created on Mar 3, 2005
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
package barna.model;

import java.util.Comparator;
import java.util.Vector;

/**
 * See <code>http://www.genomicglossaries.com/content/gene_def.asp</code>:
 * "What are the rules for deciding whether two partially overlapping mRNAs should be 
 * declared to be  alternative transcripts of the same gene or products of different genes? 
 * We have none."
 * <br><br>
 * Ensembl says in <code>http://whatislife.com/gene-definition.html</code>: 
 * A gene is a set of connected transcripts. A transcript is a set of exons via transcription 
 * followed (optionally) by pre-mRNA splicing. Two transcripts are connected if they share at 
 * least part of one exon in the genomic coordinates. At least one transcript must be expressed 
 * outside of the nucleus and one transcript must encode a protein (see footnotes).  
 * <br><br>
 * Hence, we (Sylvain, Tyler and me) define an exon to belong to strictly one gene.
 * <br><br>
 * 
 * 
 * @author micha
 */
public class Exon extends DirectedRegion {

	public static final String GFF_FEATURE_EXON= "exon";
	public static int MIN_INTRON_LENGTH_HUMAN= 33;
	public static boolean checkAcceptableIntron(SpliceSite donor, SpliceSite acceptor) {
//		   if (($seq1 =~ /^GT/i && $seq2 =~ /AG$/i)
//		   ||($seq1 =~ /^GC/i && $seq2 =~ /AG$/i)
//		   ||($seq1 =~ /^ATATC/i && $seq2 =~ /AG$/i)
//		   ||($seq1 =~ /^ATATC/i && $seq2 =~ /AC$/i)
//		   ||($seq1 =~ /^ATATC/i && $seq2 =~ /AT$/i)
//		   ||($seq1 =~ /^GTATC/i && $seq2 =~ /AT$/i)
//		   ||($seq1 =~ /^ATATC/i && $seq2 =~ /AA$/i) 

		String donSeq= Graph.readSequence(donor, 0, 3).toUpperCase();
		String donSeq2= donSeq.substring(0,2);
		String accSeq= Graph.readSequence(acceptor, 0, 0).toUpperCase();
		if (!((donSeq.startsWith("G")|| donSeq.startsWith("A"))&& accSeq.startsWith("A")))
			return false;
		if (accSeq.equals("AG")) {
			if (donSeq2.equals("GT")|| donSeq2.equals("GC")|| donSeq.equals("GTATC"))
				return true;
			return false;
		} else {
			if (donSeq.equals("ATATC"))
				return true;
			if (donSeq.equals("GTATC")&& accSeq.equals("AT"))
				return true;
			return false;
		}
			
	}
	
	byte frame= -1;
	EqualComparator defaultEqualComparator= new EqualComparator();
	Gene gene= null;
	
	
	
	public static class EqualComparator implements Comparator {
		public int compare(Object o1, Object o2) {
			Exon ex1= (Exon) o1;
			Exon ex2= (Exon) o2;

			int start1= ex1.getStart();
			int end1= ex1.getEnd();
			String chr1= ex1.getChromosome();
			SpliceSite don1= ex1.getDonor();
			SpliceSite acc1= ex1.getAcceptor();
			int start2= ex2.getStart();
			int end2= ex2.getEnd();
			String chr2= ex2.getChromosome();
			SpliceSite don2= ex2.getDonor();
			SpliceSite acc2= ex1.getAcceptor();
			
			if (start1< start2)
				return -1;
			if (start2< start1)
				return 1;
			if (end1< end2)
				return -1;
			if (end2< end1)
				return 1;
			if (acc1.getType()< acc2.getType())
				return -1;
			if (acc2.getType()< acc1.getType())
				return 1;
			if (don1.getType()< don2.getType())
				return -1;
			if (don2.getType()< don1.getType())
				return 1;
			return 0;
		}
	}
	public static class IdentityComparator implements Comparator {

		/**
		 * checks chr, pos (start, end), types (don, acc)
		 */
		public int compare(Object o1, Object o2) {
		
				Exon ex1= (Exon) o1;
				Exon ex2= (Exon) o2;

				int start1= ex1.getStart();
				int end1= ex1.getEnd();
				String chr1= ex1.getChromosome();
				SpliceSite don1= ex1.getDonor();
				SpliceSite acc1= ex1.getAcceptor();
				int start2= ex2.getStart();
				int end2= ex2.getEnd();
				String chr2= ex2.getChromosome();
				SpliceSite don2= ex2.getDonor();
				SpliceSite acc2= ex1.getAcceptor();
				
				if (chr1.equals(chr2)&& start1== start2&& end1== end2&&
						don1== don2&& acc1== acc2)	// no object identity
					return 0;
				
					// non-overlapping, one before the other
				// cancelled, not working for neg. strand (clustering, sort array asc with start, end pos)
	//			if (end1< start2)
	//				return -1;		// one stops before the other
	//			if (end2< start1)
	//				return 1;
				
				if (!chr1.equals(chr2))
					return chr1.compareTo(chr2);
				
					// overlapping: none stops before the other
				if (start1< start2)
					return -1;
				if (start2< start1)
					return 1;
				
					// overlapping and same start position
				if (start1< end2)
					return -1;
				if (end2< start1)
					return 1;
				
				if (don1== null) {
					if (don2!= null)
						return -1; 	// abstract sites before real SSs
				} else {
					if (don2== null)
						return 1;
					else {
						if (don1.getPos()< don2.getPos())
							return -1;
						else if (don2.getPos()< don1.getPos())
							return 1;
					}
				}
				if (acc1== null) {
					if (acc2!= null)
						return 1; 	// abstract sites after real SSs
				} else {
					if (acc2== null)
						return -1;
					else {
						if (acc1.getPos()< acc2.getPos())
							return -1;
						else if (acc2.getPos()< acc1.getPos())
							return 1;
					}
				}
				
				//System.err.println("assertion in abstractregion.positioncomparator failed");
				return 0;	// identical positions --> never reached
					
				}
		
	}
	
	
	public static class PositionSSComparator extends AbstractRegion.PositionComparator {
	
		/**
		 * checks chr, pos (start, end), types (don, acc)
		 */
		public int compare(Object o1, Object o2) {
				
				int res= super.compare(o1, o2);
				if (res!= 0)
					return res;
			
				Exon ex1= (Exon) o1;
				Exon ex2= (Exon) o2;
	
				SpliceSite don1= ex1.getDonor();
				SpliceSite acc1= ex1.getAcceptor();
				SpliceSite don2= ex2.getDonor();
				SpliceSite acc2= ex1.getAcceptor();
				
				if (don1== null&& don2!= null)	// no object identity
					return -1;
				if (don1!= null&& don2== null)
					return 1;
				if (acc1== null&& acc2!= null)
					return -1;
				if (acc1!= null&& acc2== null)
					return 1;
				
				if (don1!= null&& don2!= null) {
					if (don1.getPos()< don2.getPos())
						return -1;
					if (don1.getPos()> don2.getPos())
						return 1;
				}
				if (acc1!= null&& acc2!= null) {
					if (acc1.getPos()< acc2.getPos())
						return -1;
					if (acc1.getPos()> acc2.getPos())
						return 1;
				}
				return 0;
					
			}
		
	}


	public boolean isCodingCompletely() {
		if (get5PrimeCDS()!= 0&& get3PrimeCDS()!= 0)
			return true;
		return false;
	}

	public boolean isCoding() {
		if (get5PrimeCDS()!= 0|| get3PrimeCDS()!= 0)
			return true;
		return false;
	}
	
	@Override
	public Object clone() {
		Exon newExon= new Exon();
		newExon.gene= gene;
		newExon.start= start;
		newExon.end= end;
		
		return newExon;
	}
	
	public boolean setPhase(byte phase) {
		
		byte f= (byte) (((byte) 3)- phase);
		if (frame>= 0&& frame!= f) {
			System.out.println("Exon "+this+" already used in another frame "+frame+"!");
			return false;
		}
		
		frame= f;
		return true;
	}
	
	public int get5PrimeCDS() {
		Transcript[] t= getGene().exonHash.get(this);
		int cds5Prime= 0;
		for (int i = 0; i < getTranscripts().length; i++) {
			if (!getTranscripts()[i].isCoding())
				continue;
			Translation trans= getTranscripts()[i].getTranslations()[0];
			if (this.overlaps(trans)) {
				int p5= Math.max(this.get5PrimeEdge(), trans.get5PrimeEdge());
				if (cds5Prime== 0)
					cds5Prime= p5;
				else
					cds5Prime= Math.min(p5, cds5Prime);
			} 
		}

		return cds5Prime;
	}
	
	public int get3PrimeCDS() {
		
		int cds3Prime= 0;
		for (int i = 0; i < getTranscripts().length; i++) {
			if (!getTranscripts()[i].isCoding())
				continue;
			Translation trans= getTranscripts()[i].getTranslations()[0];
			if (this.overlaps(trans)) {
				int p3= Math.min(this.get3PrimeEdge(), trans.get3PrimeEdge());
				if (cds3Prime== 0)
					cds3Prime= p3;
				else
					cds3Prime= Math.max(p3, cds3Prime);
			} 
		}
		
		return cds3Prime;
	}
	
	public boolean isCoding5Prime() {
		int x= get5PrimeCDS();
		if (x!= 0&& get5PrimeCDS()== get5PrimeEdge())
			return true;
		return false;
	}
	
	/**
	 * 
	 * @return <true> if some transcript shows this exon as internal exon
	 */
	public boolean isInternal() {
		for (int i = 0; i < getTranscripts().length; i++) {
			if (getTranscripts()[i].isInternalExon(this))
				return true;
		}
		return false;
	}
	
	
	public boolean isCodingSomewhere5Prime() {
		for (int i = 0; i < getTranscripts().length; i++) {
			if (!getTranscripts()[i].isCoding())
				continue;
			Translation tln= getTranscripts()[i].getTranslations()[0];
			if (tln.get5PrimeEdge()<= get5PrimeEdge()&& tln.get3PrimeEdge()>= get5PrimeEdge())
				return true;
		}
		return false;	 
	}
	public boolean isUpstreamRegion(DirectedRegion reg) {
		if (reg.get3PrimeEdge()+ 1== get5PrimeEdge())
			return true;
		return false;
	}
	
	public boolean isDownstreamRegion(DirectedRegion reg) {
		if (get3PrimeEdge()+1== reg.get5PrimeEdge())
			return true;
		return false;
	}

	public boolean overlapsCDS() {
		for (int i = 0; i < getTranscripts().length; i++) {
			if (!getTranscripts()[i].isCoding())
				continue;
			Translation tln= getTranscripts()[i].getTranslations()[0];
			if (tln.overlaps(this))
				return true;
		}
		return false;	 
	}
	
	public boolean isCodingSomewhere3Prime() {
		for (int i = 0; i < getTranscripts().length; i++) {
			if (!getTranscripts()[i].isCoding())
				continue;
			Translation tln= getTranscripts()[i].getTranslations()[0];
			if (tln.get5PrimeEdge()<= get3PrimeEdge()&& tln.get3PrimeEdge()>= get3PrimeEdge())
				return true;
		}
		return false;	 
	}
	
	
	public boolean isCoding3Prime() {
		int x= get3PrimeCDS();
		if (x!= 0&& get3PrimeCDS()== get3PrimeEdge())
			return true;
		return false;
	}
	
	public boolean setFrame(byte newFrame) {
		
		if (frame!= 0&& frame!= newFrame) {
			System.out.println("Exon "+this+" already used in another frame "+frame+"!");
			return false;
		}
				
		frame= newFrame;
		return true;
	}
	
	/**
	 * @param b
	 */
	public boolean checkStrand(boolean b) {
		
		return (b== getGene().isForward());
	}
	public boolean checkStrand(String newStrand) {
		
		String nStrand= newStrand.trim();
		if (nStrand.equals("1"))	// || nStrand.equals("forward")
			return checkStrand(true);
		else if (nStrand.equals("-1"))	// || nStrand.equals("reverse")
			return checkStrand(false);
		
		return false; // error			
	}
	
	String exonID= null;	
	protected Exon() {
		// for subclasses
	}
	public Exon(Transcript newTranscript, String stableExonID, int start, int end) {

		this.strand= newTranscript.getStrand();
		this.chromosome= newTranscript.getChromosome();
		this.gene= newTranscript.getGene();
		setStart(start);
		setEnd(end);		
		setID("exon");
		newTranscript.addExon(this);
	}	
	/**
	 * @return the ID of <code>this</code> exon
	 */
	public String getExonID() {
		return exonID;
	}

	/**
	 * @return the gene to which <code>this</code> exon belongs
	 */
	public Gene getGene() {
		return this.gene;
	}
	
	public Transcript[] getTranscripts() {
		return getGene().getTranscripts(this);
	}
	

	public byte getStrand() {
		if (strand!= 0)
			return strand;
		if (getGene()!= null)
			return getGene().getStrand();
		return 0;
	}
	
	public String toString() {
		// return getExonID();
		String res= "";
		if (getAcceptor()!= null)
			res+= getAcceptor();
		else
			res+=0;
		res+="==";
		if (getDonor()!= null)
			res+= getDonor();
		else
			res+=0;
		
		return res;
	}

	public String toPosString() {
		return getGene().getChromosome()+" "+getStart()+" "+getEnd();
	}

	public String getChromosome() {
		if (getGene()== null) {
			if (chromosome!= null)
				return chromosome.toUpperCase();
			return chromosome;
		}
		return getGene().getChromosome();
	}
	
	public Species getSpecies() {
		if (species!= null) 
			return species;
		if (getGene()!= null)
			return getGene().getSpecies();
		
		return null;
	}

	static final long serialVersionUID = 8914674126313232057L;
	public SpliceSite getAcceptor() {
//		if (get5PrimeEdge()== 367659)
//			System.currentTimeMillis();
		
		Vector<SpliceSite> v= getGene().getSpliceSites(new Integer(get5PrimeEdge()));
		if (v== null)
			return null;
		for (int i = 0; i < v.size(); i++) {
			if (v.elementAt(i).isLeftFlank())
				return v.elementAt(i);
		}
		return null;
	}
	public SpliceSite getDonor() {
		Vector<SpliceSite> v= getGene().getSpliceSites(new Integer(get3PrimeEdge()));
		if (v== null)
			System.currentTimeMillis();
		for (int i = 0; i < v.size(); i++) {
			if (v.elementAt(i).isRightFlank())
				return v.elementAt(i);
		}
		return null;
	}
	

	public int get3PrimeFrame() {
		int fr= getFrame();
		int len= get3PrimeCDS()- get5PrimeCDS()+ 1;
		if (fr>= len- 1)	 
			fr-= (len- 1);
		else 
			fr= (fr+ ((len- 1)% 3))% 3;		// simpler?, inverse op of %??
		
		return fr;		
	}
	
	
	public byte getFrame() {
		if (frame== -1) {
			
			for (int i = 0; i < getTranscripts().length; i++) {
				if ((!getTranscripts()[i].isCoding())|| (!getTranscripts()[i].getTranslations()[0].contains(get5PrimeCDS())))
					continue;
				int pos= getTranscripts()[i].getTranslations()[0].getTranslatedPosition(get5PrimeCDS());
				byte f= (byte) (pos% 3);
				if (frame== -1)
					frame= f;
				else
					assert(frame== f);
			}
			
		}
		return frame;
	}

	public int getEndCDS() {
		if (isForward())
			return get3PrimeCDS();
		return get5PrimeCDS();
	}

	public int getStartCDS() {
		if (isForward())
			return get5PrimeCDS();	// here not getter method, avoid init for buildup
		return get3PrimeCDS();
	}

	public EqualComparator getDefaultEqualComparator() {
		return defaultEqualComparator;
	}

	public void setGene(Gene gene) {
		this.gene = gene;
	}

}
