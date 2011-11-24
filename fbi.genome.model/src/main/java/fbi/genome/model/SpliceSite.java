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
 * Created on Feb 22, 2006
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
package fbi.genome.model;

import java.util.Comparator;
import java.util.Vector;


/**
 * 
 * <tt>07/11/14</tt>: soft starts..
 * @author msammeth
 */
public class SpliceSite extends AbstractSite {

	public static final String ATTRIBUTE_ID_SS_TYPE= "ss_type";
	
	public final static byte TYPE_NOT_INITED= -1, TYPE_DONOR= 1, TYPE_ACCEPTOR= 2, TYPE_HARD_START= 3, TYPE_HARD_END= 4, TYPE_SOFT_START=5, TYPE_SOFT_END= 6, TYPE_CODON_START= 7, TYPE_CODON_STOP= 8;	
	
	public final static byte NOTYPE_SS= 0;
	public final static byte CONSTITUTIVE= 1;
	public final static byte ALTERNATIVE= 2;
	
	public final static int DEFAULT_DONOR_5FLANK= 2;
	public final static int DEFAULT_DONOR_3FLANK= 6;
	public final static int DEFAULT_ACCEPTOR_5FLANK= 15;
	public final static int DEFAULT_ACCEPTOR_3FLANK= 2;
	
	static final long serialVersionUID = 422169949942461293L;
	static final int SPLICE_SITE_FLANK_SIZE= 20; // Alternative splicing in the human, mouse and rat genomes is associated with an increased frequency of exon creation and/or loss
												 // Barmak Modrek & Christopher J Lee
	public static final int DELTA_RANGE= 20;
	
	public static final int NOBORU_DON5_EXTENT= 20;
	public static final int NOBORU_DON3_EXTENT= 18;
	public static final int NOBORU_ACC5_EXTENT= 21;
	public static final int NOBORU_ACC3_EXTENT= 17;

	public static final int SHAPIRO_DON5_FIRST= -3;
	public static final int SHAPIRO_DON3_LAST= +6;
	public static final int SHAPIRO_ACC5_FIRST= -13;
	public static final int SHAPIRO_ACC3_LAST= +1;

	public static final String CANONICAL_DONOR = "GT";

	public static final String CANONICAL_ACCEPTOR = "AG";


	static PositionTypeComparator defaultPositionTypeComparator= null;
	static PositionComparator defaultPositionComparator= null;

	public static char SYMBOL_UNDEFINED = '?', SYMBOL_DONOR = '^', SYMBOL_ACCEPTOR = '-', SYMBOL_HARD_START = '[', SYMBOL_HARD_END = ']', SYMBOL_SOFT_START='(', SYMBOL_SOFT_END= ')',
		SYMBOL_CODON_START= '>', SYMBOL_CODON_STOP= '<';

    public final static String FEATURE_ID_ACCEPTOR= "acceptor";

    public final static String FEATURE_ID_DONOR= "donor";

	
	public static byte inverseType(byte ssType) {
		if (ssType== TYPE_HARD_START)
			return TYPE_HARD_END;
		if (ssType== TYPE_HARD_END)
			return TYPE_HARD_START;
		if (ssType== TYPE_DONOR)
			return TYPE_ACCEPTOR;
		if (ssType== TYPE_ACCEPTOR)
			return TYPE_DONOR;
		if (ssType== TYPE_SOFT_START)
			return TYPE_SOFT_END;
		if (ssType== TYPE_SOFT_END)
			return TYPE_SOFT_START;
		if (ssType!= TYPE_NOT_INITED)
			System.err.println("Unknown ss type "+ssType);
		return TYPE_NOT_INITED;
	}
	
	public static final SpliceChainComparator defaultSpliceChainComparator= new SpliceChainComparator();
	
	public static class SpliceChainComparator implements Comparator<SpliceSite[]> {
		public int compare(SpliceSite[] o1, SpliceSite[] o2) {
			if (o1.length< o2.length)
				return -1;
			if (o1.length> o2.length)
				return 1;
			for (int i = 0; i < o2.length; i++) {
				if (o1[i].getPos()< o2[i].getPos())
					return -1;
				if (o1[i].getPos()> o2[i].getPos())
					return 1;
			}
			return 0;
		}
	}
	
	public static PositionTypeComparator getDefaultPositionTypeComparator() {
		if (defaultPositionTypeComparator == null) {
			defaultPositionTypeComparator = new PositionTypeComparator();			
		}
	
		return defaultPositionTypeComparator;
	}

	public static PositionComparator getDefaultPositionComparator() {
		if (defaultPositionComparator == null) {
			defaultPositionComparator = new PositionComparator();			
		}

		return defaultPositionComparator;
	}
	
	
	float scoreU12= Float.NaN;
	

	
	byte type= TYPE_NOT_INITED, sourceType= Transcript.ID_SRC_UNDEFINED;
	byte modality= NOTYPE_SS;
	boolean constitutive= true;
	String dinucleotide= null, stringRep= null;
	public static class PositionComparator implements Comparator {
		public int compare(Object arg0, Object arg1) {
			SpliceSite ss1= (SpliceSite) arg0;
			SpliceSite ss2= (SpliceSite) arg1;
			if (ss1.getPos()< ss2.getPos())
				return -1;
			if (ss2.getPos()< ss1.getPos())
				return 1;
			return 0;
		}
	}
	
	/**
	 * 
	 * @deprecated replaced by PositionEqualSSTypeComparator
	 *
	 */
	
	public static class PositionTypeComparator_old extends PositionComparator {
		public int compare(Object arg0, Object arg1) {
			
			int res= super.compare(arg0, arg1);
			if (res!= 0)
				return res;
			
			// positions equal
			if (arg0 instanceof SpliceSite) {
				if (arg1 instanceof SpliceSite) {
					SpliceSite s0= (SpliceSite) arg0;
					SpliceSite s1= (SpliceSite) arg1;
					if (s0.isDonor()&& !s1.isDonor())
						return 1;	// for 1-lentgth exons, sort acceptor before donor (same pos)
					if (s1.isDonor()&& !s0.isDonor())
						return -1;
					return 0;	// same type and position					
				} else 
					return 1;	// arg1 is no SS					
			} else {
				if (arg1 instanceof SpliceSite) 
					return -1;	// arg0 is no SS
				else {	// 2 no SS
					boolean tss0= ((AbstractSite) arg0).isTSS();
					boolean tss1= ((AbstractSite) arg1).isTSS();
					if (tss0== tss1)
						return 0;
					if (tss0)
						return -1;
					return 1;
				}
			}
		}

		public int compare_old(Object arg0, Object arg1) {
				// positions equal
			if (arg0 instanceof SpliceSite) {
				if (arg1 instanceof SpliceSite) {
					SpliceSite s0= (SpliceSite) arg0;
					SpliceSite s1= (SpliceSite) arg1;
					if (s0.isDonor()&& !s1.isDonor())
						return 1;	// for 1-lentgth exons, sort acceptor before donor (same pos)
					if (s1.isDonor()&& !s0.isDonor())
						return -1;
					return 0;	// same type and position					
				} else 
					return 1;	// arg1 is no SS					
			} else {
				if (arg1 instanceof SpliceSite) 
					return -1;	// arg0 is no SS
				else {	// 2 no SS
					boolean tss0= ((AbstractSite) arg0).isTSS();
					boolean tss1= ((AbstractSite) arg1).isTSS();
					if (tss0== tss1)
						return 0;
					if (tss0)
						return -1;
					return 1;
				}
			}
		}
	}

	public static class PositionEqualSSTypeComparator extends PositionComparator {
		public int compare(Object arg0, Object arg1) {
			
			int res= super.compare(arg0, arg1);
			if (res!= 0)
				return res;
			
			// positions equal: TSS before donor, acc, TES
			SpliceSite s0= (SpliceSite) arg0;
			SpliceSite s1= (SpliceSite) arg1;
			if (s0.isTSS()&& !s1.isTSS())
				return -1;
			if ((!s0.isTSS())&& s1.isTSS())
				return 1;
			if (s0.isAcceptor()&& !s1.isAcceptor())
				return -1;
			if ((!s0.isAcceptor())&& s1.isAcceptor())
				return 1;
			if (s0.isDonor()&& !s1.isDonor())
				return -1;
			if ((!s0.isDonor())&& s1.isDonor())
				return 1;
			if (s0.isTES()&& !s1.isTES())
				return -1;
			if ((!s0.isTES())&& s1.isTES())
				return 1;
			
			return 0; 	// everything equal
		}
	
		public int compare_old(Object arg0, Object arg1) {
				// positions equal
			if (arg0 instanceof SpliceSite) {
				if (arg1 instanceof SpliceSite) {
					SpliceSite s0= (SpliceSite) arg0;
					SpliceSite s1= (SpliceSite) arg1;
					if (s0.isDonor()&& !s1.isDonor())
						return 1;	// for 1-lentgth exons, sort acceptor before donor (same pos)
					if (s1.isDonor()&& !s0.isDonor())
						return -1;
					return 0;	// same type and position					
				} else 
					return 1;	// arg1 is no SS					
			} else {
				if (arg1 instanceof SpliceSite) 
					return -1;	// arg0 is no SS
				else {	// 2 no SS
					boolean tss0= ((AbstractSite) arg0).isTSS();
					boolean tss1= ((AbstractSite) arg1).isTSS();
					if (tss0== tss1)
						return 0;
					if (tss0)
						return -1;
					return 1;
				}
			}
		}
	}

	public static class PositionTypeComparator extends PositionComparator {
		public int compare(Object arg0, Object arg1) {
			
			int res= super.compare(arg0, arg1);
			if (res!= 0)
				return res;
			
			// positions equal: TSS before donor, acc, TES
			SpliceSite s0= (SpliceSite) arg0;
			SpliceSite s1= (SpliceSite) arg1;
			
			// codons after left flank, before right flank
			if (s0.isCodon()|| s1.isCodon()) {
				if ((s0.isLeftFlank()&& s1.isCodon())|| (s0.isCodon()&& s1.isRightFlank()))
					return -1;
				if ((s1.isLeftFlank()&& s0.isCodon())|| (s1.isCodon()&& s0.isRightFlank()))
					return 1;
				if (s0.isCodon()&& s1.isCodon()) {
					if (s0.getType()== s1.getType())
						return 0;
					if (s0.isCodonStart())
						return -1;
				}
			}			
			
			if (s0.isLeftFlank()&& s1.isRightFlank())
				return -1;
			if (s0.isRightFlank()&& s1.isLeftFlank())
				return 1;
			if (s0.isTSS()&& !s1.isTSS())
				return -1;
			if (s1.isTSS()&& !s0.isTSS())
				return 1;
			if ((!s0.isTES())&& s1.isTES())
				return -1;
			if (s0.isTES()&& !s1.isTES())
				return 1;
			
			return 0; 	// everything equal, TODO undefined
		}
	
		public int compare_old(Object arg0, Object arg1) {
				// positions equal
			if (arg0 instanceof SpliceSite) {
				if (arg1 instanceof SpliceSite) {
					SpliceSite s0= (SpliceSite) arg0;
					SpliceSite s1= (SpliceSite) arg1;
					if (s0.isDonor()&& !s1.isDonor())
						return 1;	// for 1-lentgth exons, sort acceptor before donor (same pos)
					if (s1.isDonor()&& !s0.isDonor())
						return -1;
					return 0;	// same type and position					
				} else 
					return 1;	// arg1 is no SS					
			} else {
				if (arg1 instanceof SpliceSite) 
					return -1;	// arg0 is no SS
				else {	// 2 no SS
					boolean tss0= ((AbstractSite) arg0).isTSS();
					boolean tss1= ((AbstractSite) arg1).isTSS();
					if (tss0== tss1)
						return 0;
					if (tss0)
						return -1;
					return 1;
				}
			}
		}
	}

	public static boolean isLeftFlank(byte type) {
		return type== TYPE_HARD_START|| type== TYPE_ACCEPTOR|| type== TYPE_SOFT_START;
	}
	
	public static boolean isRightFlank(byte type) {
		return type== TYPE_HARD_END|| type== TYPE_DONOR|| type== TYPE_SOFT_END;
	}
	
	public static boolean isU2ss(String id) {
		if (id.contains("U2"))
			return true;
		return false;
	}
	
	public static int[] getPositions(SpliceSite[] ss) {
		if (ss== null)
			return null;
		int[] pos= new int[ss.length];
		for (int i = 0; i < pos.length; i++) 
			pos[i]= ss[i].getPos();
		return pos;
	}
	
	public SpliceSite(int pos, byte newType, Gene newGene) {
	
		super(pos);
		setType(newType);		
		setGene(newGene);
	}

	public SpliceSite(int pos, byte newType, Gene newGene, byte newSourceType) {
		
		this(pos, newType, newGene);
		setSourceType(newSourceType);
	}
	
	/**
	 * checks for identity, not for homology
	 */
	public boolean equals(Object obj) {
		
		SpliceSite otherSS= (SpliceSite) obj;
		if (getPos()== otherSS.getPos()&& getType()== otherSS.getType())
			return true;
		return false;
		
//		if (!super.equals(anotherSS))
//			return false;
//		
//		if (!(anotherSS instanceof SpliceSite))
//			return false;
//		
//		SpliceSite s2= (SpliceSite) anotherSS;
//		if (s2.isDonor()!= isDonor())
//			return false;
////		SpliceSite aSS= (SpliceSite) anotherSS;
////		if (gene!= aSS.getGene()|| pos!= aSS.getPos())
////			return false;
//		return true;
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
		// not required, read text
		// bullshit, gotta change otherwise ss overwrite as in hash
		// and the other way
		return getPos(); //super.hashCode();
		
	}
	
	public static byte getTypeBySymbol(char c) {
		if (c== SYMBOL_HARD_START)
			return TYPE_HARD_START;
		if (c== SYMBOL_DONOR)
			return TYPE_DONOR;
		if (c== SYMBOL_ACCEPTOR)
			return TYPE_ACCEPTOR;
		if (c== SYMBOL_HARD_END)
			return TYPE_HARD_END;
		if (c== SYMBOL_SOFT_START)
			return TYPE_SOFT_START;
		if (c== SYMBOL_SOFT_END)
			return TYPE_SOFT_END;
		return TYPE_NOT_INITED;
	}
	
	public boolean isDonor() {
		return (type== TYPE_DONOR);
	}

	public boolean isTSS() {
		return (type== TYPE_HARD_START|| type== TYPE_SOFT_START);
	}
	
	public boolean isTES() {
		return (type== TYPE_HARD_END|| type== TYPE_SOFT_END);
	}
	
	public boolean isHardEdge() {
		return (type== TYPE_HARD_START|| type== TYPE_HARD_END);
	}
	
	public boolean isSoftEdge() {
		return (type== TYPE_SOFT_START|| type== TYPE_SOFT_END);
	}
	
	public boolean isSpliceSite() {
		return (isDonor()|| isAcceptor());
	}

	public boolean isAcceptor() {
		return (type== TYPE_ACCEPTOR);
	}
	
	public boolean isStartOrStopCodon() {
		return (isCodonStart()|| isCodonStop());
	}
	
	public boolean isCodonStart() {
		return (type== TYPE_CODON_START);
	}
	
	public boolean isCodonStop() {
		return (type== TYPE_CODON_STOP);
	}
	
	public boolean isCodon() {
		return (isCodonStart()|| isCodonStop());
	}
	
	public boolean isRealSpliceSite() {
		return (isDonor()|| isAcceptor());
	}
	
	public boolean isLeftFlank() {
		return (type== TYPE_ACCEPTOR|| type== TYPE_HARD_START|| type== TYPE_SOFT_START);
	}
	
	public boolean isRightFlank() {
		return (type== TYPE_DONOR|| type== TYPE_HARD_END|| type== TYPE_SOFT_END);
	}
	
	
	
	public boolean isCanonical() {
		if ((!this.isDonor())&& (!this.isAcceptor()))
			System.err.println("Not canonical() called for non-splicesite");
		
		String seq= getDinucleotide();		
		if ((isDonor()&& seq.equalsIgnoreCase(CANONICAL_DONOR))||
				(isAcceptor()&& seq.equalsIgnoreCase(CANONICAL_ACCEPTOR)))
			return true;
		return false;				
	}
	
	public String getDinucleotide() {
		if (dinucleotide == null) {
			dinucleotide = Graph.readSequence(this, 0, 0);
		}

		return dinucleotide;
	}
	
	public void setType(byte newType) {
		this.type = newType;
		if (isDonor())
			setID(FEATURE_ID_DONOR);
		else if (isAcceptor())
			setID(FEATURE_ID_ACCEPTOR);
	}
	public String toString() {
		if (stringRep == null) {	// || stringRep.charAt(stringRep.length()- 1)== '?'
			stringRep= Integer.toString(getPos())+ getSiteSymbol();
		}

		return stringRep;
	}
	
	public String toOutString() {
		String result= getGene().getChromosome()+ " "+ Math.abs(getPos())+ " "+ getGene().getStrand();
		return result;
	}
	
	public Transcript[] getTranscripts() {
		if (transcripts== null) {
			Vector<Transcript> v= getGene().ssTrptHash.get(this);
			// inconsistent in build up
			Transcript[] transcripts= new Transcript[v.size()];	// TODO 081212 better to always ask gene, use setTranscripts() to fix them
			v.toArray(transcripts);
			return transcripts;
		}
		return transcripts;
	}
	
	public void setTranscripts(Transcript[] tt) {
		this.transcripts= tt;
	}
	
	public boolean isCDSRealTranscript() {
		
				// net gut, twilight-frontier events !
	//		for (int i = 0; i < asVars.length; i++) 
	//			if (asVars[i].isProteinCoding())
	//				return true;
			
			int ctr= 0;
			for (int i = 0; i < transcripts.length; i++) {
				if (transcripts[i].getTranslations()== null|| transcripts[i].getTranslations().length< 1)
					continue;
				if (transcripts[i].getTranslations()[0].contains(this))
					ctr++;	// at least in the CDS of one transcript
				else 
					return false;
			}
			return (ctr> 0);	// if all are non-annotated, ok, its never in CDS
		}

	public boolean isCDSMaxTranscript() {
		
		// net gut, twilight-frontier events !
//		for (int i = 0; i < asVars.length; i++) 
//			if (asVars[i].isProteinCoding())
//				return true;
	
		for (int i = 0; i < transcripts.length; i++) {
			if (transcripts[i].getTranslations()== null|| transcripts[i].getTranslations().length< 1)
				continue;
			if (transcripts[i].getTranslations()[0].contains(getPos()))
				return true;
		}
		return false;	// if all are non-annotated, ok, its never in CDS
	}
	
	public boolean is5UTRRealTranscript() {
		int ctr= 0;
		for (int i = 0; i < transcripts.length; i++) {
			if (transcripts[i].getTranslations()== null|| transcripts[i].getTranslations().length< 1)
				continue;
			if (this.getPos()>= transcripts[i].getTranslations()[0].get5PrimeEdge())
				return false;	// at least outside of 5UTR
			++ctr;
		}
		return (ctr> 0);	// at least one transcript w annotated CDS found
	}
	
	public String getGenicLocation() {
		Transcript[] t= getTranscripts();
		boolean cds= false, utr5= false, utr3= false;
		int pos= getRealSSpos();
		for (int i = 0; i < t.length; i++) {
			if (!t[i].isCoding())
				continue;
			if (t[i].getTranslations()[0].contains(pos))
				cds= true;
			else if (pos < t[i].getTranslations()[0].get5PrimeEdge())
				utr5= true;
			else if (pos > t[i].getTranslations()[0].get3PrimeEdge())
				utr3= true;
		}

		if (cds)
			return "CDS";
		if (utr5&& !utr3)
			return "5UTR";
		if (utr3&& !utr5)
			return "3UTR";
		return "UNDEFINED";
	}
	
	public boolean is5UTRMaxTranscript() {
		for (int i = 0; i < transcripts.length; i++) {
			if (transcripts[i].getTranslations()== null|| transcripts[i].getTranslations().length< 1
					|| transcripts[i].getTranslations()[0].isOpenEnded5())
				continue;
			if (this.getPos()< transcripts[i].getTranslations()[0].get5PrimeEdge())
				return true;	// at least outside of 5UTR
		}
		return false;
	}

	public boolean is3UTRMaxTranscript() {
		for (int i = 0; i < transcripts.length; i++) {
			if (transcripts[i].getTranslations()== null|| transcripts[i].getTranslations().length< 1||
					transcripts[i].getTranslations()[0].isOpenEnded3())
				continue;
			if (this.getPos()> transcripts[i].getTranslations()[0].get3PrimeEdge())
				return true;	// at least outside of 5UTR
		}
		return false;
	}


	/**
	 * 
	 * @param constitutive
	 * @deprecated as to link SSs with ASVars
	 */
	public void setConstitutive(boolean constitutive) {
		this.constitutive = constitutive;
	}
	
	/**
	 * from (-x to + y) number of exonic/intronic positions to
	 * number of positions before and after ss
	 * before and after ss
	 * @param left
	 * @param right
	 * @return
	 */
	public int[] convertToExonIntronPositions(int left, int right) {
		int[] result= new int[2];
		if (isDonor()) {
			result[0]= Math.abs(left);
			result[1]= right- 2;
		} else {
			result[0]= Math.abs(left)- 2;
			result[1]= right;
		}
		return result;
	}
	
	public DirectedRegion getShapiroRegion() {
		
		if (isDonor()) {
			int[] extent= convertToExonIntronPositions(SHAPIRO_DON5_FIRST, SHAPIRO_DON3_LAST);
			return getRegion(extent[0], extent[1]);
		} else {
			int[] extent= convertToExonIntronPositions(SHAPIRO_ACC5_FIRST, SHAPIRO_ACC3_LAST);
			return getRegion(extent[0], extent[1]);
		}
			
	}
	
	public DirectedRegion getRegion(int downstream, int upstream) {
		int start= getPos();
		int end= getPos();
		
		if (isDonor()) {
			start-= downstream- 1;	// -2
			end+= upstream+ 2;	// +4
		} else {
			start-= downstream+ 2;	// +1
			end+= upstream- 1;	// 
		}
			
		DirectedRegion reg= new DirectedRegion(start, end, getTranscripts()[0].getStrand());
		reg.setChromosome(getTranscripts()[0].getChromosome());
		reg.setSpecies(getTranscripts()[0].getSpecies());
		return reg;
	}

	public DirectedRegion getNoboruRegion() {
		int start= getPos();
		int end= getPos();
		if (isDonor()) {
			start-= NOBORU_DON5_EXTENT- 1;	// -2
			end+= NOBORU_DON3_EXTENT+ 2;	// +4
		} else {
			start-= NOBORU_ACC5_EXTENT+ 2;	// +1
			end+= NOBORU_ACC3_EXTENT- 1;	// 
		}
			
		DirectedRegion reg= new DirectedRegion(start, end, getTranscripts()[0].getStrand());
		reg.setChromosome(getTranscripts()[0].getChromosome());
		reg.setSpecies(getTranscripts()[0].getSpecies());
		return reg;
	}

	public float getScoreU12() {
		return scoreU12;
	}

	public void setScoreU12(float scoreU12) {
		this.scoreU12 = scoreU12;
	}

	public byte getType() {
		return type;
	}

	public char getSiteSymbol() {
		if (isDonor())
			return SYMBOL_DONOR;
		if (isAcceptor())
			return SYMBOL_ACCEPTOR;
		if (type== TYPE_HARD_START)
			return SYMBOL_HARD_START;
		if (type== TYPE_HARD_END)
			return SYMBOL_HARD_END;
		if (type== TYPE_SOFT_START)
			return SYMBOL_SOFT_START;
		if (type== TYPE_SOFT_END)
			return SYMBOL_SOFT_END;
		if (type== TYPE_CODON_START)
			return SYMBOL_CODON_START;
		if (type== TYPE_CODON_STOP)
			return SYMBOL_CODON_STOP;
		
		return SYMBOL_UNDEFINED;
	}
	
	/**
	 * +1 for donors, -1 for acceptors
	 * @return
	 */
	public int getRealSSpos() {
		if (isDonor())
			return getPos()+ 1;
		if (isAcceptor())
			return getPos()- 1;
		return getPos();
	}

	public boolean isAlternative() {
		//return (!isConstitutive());
		return (modality== ALTERNATIVE);
	}

	public boolean isConstitutive() {
		//return constitutive;
		//return (asVars== null|| asVars.length< 1);
		return (modality== CONSTITUTIVE);
	}

	public byte getModality() {
		return modality;
	}

	public void setModality(byte modality) {
		this.modality = modality;
	}

	public static SpliceChainComparator getDefaultSpliceChainComparator() {
		return defaultSpliceChainComparator;
	}

	public byte getSourceType() {
		if (sourceType == Transcript.ID_SRC_UNDEFINED&&
				getTranscripts().length> 0) {
			sourceType= Byte.MAX_VALUE;
			for (int j = 0; j < getTranscripts().length; j++) {
				if (getTranscripts()[j].getSourceType()< sourceType)
					sourceType= getTranscripts()[j].getSourceType();
			}
		}

		return sourceType;
	}

	public void setSourceType(byte sourceType) {
		this.sourceType = sourceType;
	}
}
