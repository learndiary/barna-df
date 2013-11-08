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

package barna.model;

import barna.commons.utils.ArrayUtils;
import barna.commons.utils.StringUtils;
import barna.model.commons.MyHashMap;
import barna.model.gff.GFFObject;

import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

//import gphase.graph.Tuple;
// import genome.tools.MyArray;

/*
 * Created on Mar 3, 2005
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */

/**
 * 
 * 
 * @author micha
 */
public class Gene extends DirectedRegion {

	public static final String[] LOCALIZATION_STRING= new String[] {"N/A", "5UTR", "5UTR-CDS", "CDS", "CDS-3UTR", "3UTR", "5UTR-CDS-3UTR"};
	public static final String GFF_FEATURE_LOCUS= "locus";
	
	SpliceSite[] spliceSites= null; // ordered according to position
	HashMap<Integer, Vector<SpliceSite>> spliceHash= new HashMap<Integer, Vector<SpliceSite>>();
	HashMap<SpliceSite, Vector<Transcript>> ssTrptHash= new HashMap<SpliceSite, Vector<Transcript>>();
	barna.model.commons.MyHashMap<DefaultRegion, Transcript[]> exonHash= new barna.model.commons.MyHashMap<DefaultRegion, Transcript[]>(7,1);
	
	Exon[] exons= null;
	boolean construct= true;	// used for hashing during construction
	DirectedRegion[] localizationRef= null;
	
	public boolean isProteinCoding() {
		for (int i = 0; transcripts!= null&& i < transcripts.length; i++) 
			if (!transcripts[i].isNonCoding())
				// transcripts[i].getTranslations()!= null&& transcripts[i].getTranslations().length> 0)
				return true;
		return false;
	}
	
	public int getMinCDSStart() {
		int min= Integer.MAX_VALUE;
		for (int i = 0; i < transcripts.length; i++) {
			Translation[] tl= transcripts[i].getTranslations();
			if (tl== null)
				continue;
			for (int j = 0; j < tl.length; j++) {
				if (Math.abs(tl[j].get5PrimeEdge())< Math.abs(min))
					min= tl[j].get5PrimeEdge();
			}
		}
		if (min== Integer.MAX_VALUE)
			return 0;
		return min;
	}
	
	public int getMinCDSEnd() {
		int min= Integer.MAX_VALUE;
		for (int i = 0; i < transcripts.length; i++) {
			Translation[] tl= transcripts[i].getTranslations();
			if (tl== null)
				continue;
			for (int j = 0; j < tl.length; j++) {
				if (Math.abs(tl[j].get3PrimeEdge())< Math.abs(min))
					min= tl[j].get3PrimeEdge();
			}
		}
		if (min== Integer.MAX_VALUE)
			return 0;
		return min;
	}
	

	
	public Transcript[][] recluster() {
		Arrays.sort(transcripts, new AbstractRegion.PositionComparator());
		Vector v= null;
		int max= Integer.MIN_VALUE;
		Vector clusters= new Vector();
		for (int i = 0; i < transcripts.length; i++) {
			if (Math.abs(transcripts[i].getStart())> Math.abs(max)) {
				if (v!= null)
					clusters.add(v);
				v= new Vector();
			} 
			v.add(transcripts[i]);
			if (Math.abs(transcripts[i].getEnd())> Math.abs(max))
				max= transcripts[i].getEnd();
		}
		if (v!= null)
			clusters.add(v);
		
		return (Transcript[][]) ArrayUtils.toField(clusters);
	}
	
	public DirectedRegion[] getCDSS() {
		
		if (cdssRegs == null) {
			Vector<DirectedRegion> vRegs= new Vector<DirectedRegion>(2,2); 
			for (int i = 0; i < transcripts.length; i++) {
				if ((!transcripts[i].source.contains("refGene"))|| (!transcripts[i].isCoding()))
					continue;
				Translation tln= transcripts[i].getTranslations()[0];
				int j = 0;
				for (; j < vRegs.size(); j++) {
					if (vRegs.elementAt(j).overlaps(tln)) {
						if (tln.get5PrimeEdge()< vRegs.elementAt(j).get5PrimeEdge())
							vRegs.elementAt(j).set5PrimeEdge(tln.get5PrimeEdge());
						if (tln.get3PrimeEdge()> vRegs.elementAt(j).get3PrimeEdge())
							vRegs.elementAt(j).set3PrimeEdge(tln.get3PrimeEdge());
						break;
					}
				}
				if (j== vRegs.size()) 
					vRegs.add((DirectedRegion) tln.clone());
			}
			cdssRegs= new DirectedRegion[vRegs.size()];
			for (int i = 0; i < cdssRegs.length; i++) 
				cdssRegs[i]= vRegs.elementAt(i);
		}

		return cdssRegs;
	}
	
	
	public DirectedRegion getReal5UTR() {
		DirectedRegion reg;
		if (isForward()) {
			int x= getMinCDSStart();
			if (x== 0)
				return null;
			reg= new DirectedRegion(getStart(), x- 1, getStrand());	// utr starts outside CDS			
		} else {
			int x= getMaxCDSStart();
			if (x== 0)
				return null;
			reg= new DirectedRegion(x- 1, getEnd(), getStrand());
		}
		reg.setChromosome(getChromosome());
		return reg;
		
	}
	
	public DirectedRegion getRealCDS() {
		
		DirectedRegion reg;
		if (isForward())
			reg= new DirectedRegion(getMaxCDSStart(), getMinCDSEnd(), getStrand());
		else 
			reg= new DirectedRegion(getMaxCDSEnd(), getMinCDSStart(), getStrand());
		reg.setChromosome(getChromosome());
		return reg;

	}

	
	public DirectedRegion[] getRegionAllCDS() {
		
		HashMap<Integer, Integer> hashStart= new HashMap<Integer, Integer>(), hashEnd= new HashMap<Integer, Integer>();
		
		for (int i = 0; i < getTranscripts().length; i++) {
			Transcript trpt= getTranscripts()[i];
			if (!trpt.isCoding())
				continue;
			Translation tln= getTranscripts()[i].getTranslations()[0];
			hashStart.put(new Integer(tln.get5PrimeEdge()), null);
			hashEnd.put(new Integer(tln.get3PrimeEdge()), null);
		}
		
		int[] starts= new int[hashStart.size()], ends= new int[hashEnd.size()];
		Iterator<Integer> iter= hashStart.keySet().iterator();
		int cnt= 0;
		while (iter.hasNext())
			starts[cnt++]= iter.next();
		iter= hashEnd.keySet().iterator();
		cnt= 0;
		while (iter.hasNext())
			ends[cnt++]= iter.next();
		Arrays.sort(starts);
		Arrays.sort(ends);
		
		// make overlap
		int ptrS= 0, ptrE= 0;
		Vector<DirectedRegion> v= new Vector<DirectedRegion>();
		while (ptrS< starts.length|| ptrE< ends.length) {
				
			int start= starts[ptrS];// min start
			int advS= 1, advE= 0;
			while (advS!= advE&& (ptrS+advS)< starts.length&& (ptrE+advE)< ends.length) {
				if (starts[ptrS+advS]< ends[ptrE+ advE])
					++advS;
				else
					++advE;
			}
			// max end
			int end= ends[ends.length-1];
			if (((ptrS+advS))< starts.length&& ((ptrE+ advE)< ends.length))			// ((ptrE+ advE)< ends.length)
				end= ends[ptrE+advE];
			else
				advE= ends.length- ptrE- 1;
			
			DirectedRegion reg= new DirectedRegion(start, end, getStrand());
			reg.setChromosome(getChromosome());
			reg.setID("allCDS");
			v.add(reg);
			
			ptrS+= advS;
			ptrE+= advE+1;
		}
		
		DirectedRegion[] regs= new DirectedRegion[v.size()];
		for (int i = 0; i < regs.length; i++) 
			regs[i]= v.elementAt(i);
		
		return regs;
	} 

	
	public DirectedRegion getMaxCDS() {
		
		DirectedRegion reg;
		if (isForward())
			reg= new DirectedRegion(getMinCDSStart(), getMaxCDSEnd(), getStrand());
		else 
			reg= new DirectedRegion(getMinCDSEnd(), getMaxCDSStart(), getStrand());
		reg.setChromosome(getChromosome());
		return reg;

	} 
	
	public boolean isRealCDS(DirectedRegion reg) {
		if (getRealCDS().contains(reg))
			return true;
		return false;
	}
	
	public boolean isReal5UTR(DirectedRegion reg) {
		if (getReal5UTR().contains(reg))
			return true;
		return false;
	}
	
	public DirectedRegion getReal3UTR() {
		
		DirectedRegion reg;
		if (isForward()) {
			int x= getMaxCDSEnd();
			if (x== 0)
				return null;
			reg= new DirectedRegion(x+1, getEnd(), getStrand());	// utr starts outside CDS			
		} else {
			int x= getMinCDSEnd();
			if (x== 0)
				return null;
			reg= new DirectedRegion(getStart(), x+ 1, getStrand());	// neg strand -(-1)
		}
		
		reg.setChromosome(getChromosome());
		return reg;
		
	}
	
	public int getMaxCDSEnd() {
		int max= 0;
		for (int i = 0; i < transcripts.length; i++) {
			Translation[] tl= transcripts[i].getTranslations();
			if (tl== null)
				continue;
			for (int j = 0; j < tl.length; j++) {
				if (Math.abs(tl[j].get3PrimeEdge())> Math.abs(max))
					max= tl[j].get3PrimeEdge();
			}
		}
		if (max== 0)
			return 0;
		return max;
	}

	public Exon[] getConstitutiveInternalExons() {
		
		HashMap<Integer, Integer> yet= new HashMap<Integer, Integer>();	// start to end of each investigated
		for (int i = 0; i < getTranscripts().length; i++) {
			Exon[] e= getTranscripts()[i].getExons();
			for (int j = 0; j < e.length; j++) {
				// first check if this exon is valid
				
			} 
		}
		// TODO
		return null;
	}
	
	
	public int getMaxCDSStart() {
		int max= 0;
		for (int i = 0; i < transcripts.length; i++) {
			Translation[] tl= transcripts[i].getTranslations();
			if (tl== null)
				continue;
			for (int j = 0; j < tl.length; j++) {
				if (Math.abs(tl[j].get5PrimeEdge())> Math.abs(max))
					max= tl[j].get5PrimeEdge();
			}
		}
		if (max== 0)
			return 0;
		return max;
	}
	
	public AbstractRegion getCDSRegion() {
		return new DefaultRegion(getMinCDSStart(), getMaxCDSEnd());
	}
	
	@Override
	public Object clone() {
		Gene newGene= new Gene(Gene.getUniqueID());
		newGene.start= start;
		newGene.end= end;
		newGene.chromosome= chromosome;
		newGene.construct= construct;
		
		newGene.transcripts= new Transcript[transcripts.length];
		for (int i = 0; i < transcripts.length; i++) 
			newGene.transcripts[i]= transcripts[i];
		
		if (attributes!= null)
			newGene.attributes= (HashMap) attributes.clone();
		if (exonHash!= null)
			newGene.exonHash= (MyHashMap<DefaultRegion, Transcript[]>) exonHash.clone();
		if (exons!= null) {
			newGene.exons= new Exon[exons.length];
			for (int i = 0; i < exons.length; i++) 
				newGene.exons[i]= exons[i];
		}
		if (spliceSites!= null) {
			newGene.spliceSites= new SpliceSite[spliceSites.length];
			for (int i = 0; i < spliceSites.length; i++) 
				newGene.spliceSites[i]= spliceSites[i];
		}
		if (spliceHash!= null)
			newGene.spliceHash= (HashMap) spliceHash.clone();
		if (ssTrptHash!= null)
			newGene.ssTrptHash= (HashMap) ssTrptHash.clone();

		if (localizationRef!= null) {
			newGene.localizationRef= new DirectedRegion[localizationRef.length];
			for (int i = 0; i < localizationRef.length; i++) 
				newGene.localizationRef[i]= localizationRef[i];
		}
		
		return newGene;
	}


    /**
     * Comparator to compare gene boundaries, either amongst genes or between genes and integer values.
     */
    public static class BoundaryComparator implements Comparator {

        /**
         * Flag to indicate whether the start or the end of genes are compared.
         */
        boolean start= true;

        /**
         * Constructs a comparator that compares the indicated boundary
         * of genes with an integer value.
         * @param start
         */
        public BoundaryComparator(boolean start) {
            this.start= start;
        }

        /**
         * Compares genes with genes respectively integer values describing the cooridnate of one of the boundaries.
         * @param o1 a <code>gene</code> or an integer value
         * @param o2 another <code>gene</code> or an integer value
         * @return <code>0</code> if both arguments are equal, otherwise the distance between the first and the
         * second argument
         */
        @Override
        public int compare(Object o1, Object o2) {
            int v1, v2;
            if (o1 instanceof Gene) {
                v1= (start? Math.abs(((Gene) o1).getStart()): Math.abs(((Gene) o1).getEnd()));
            } else if (o1 instanceof Integer)
                v1= ((Integer) o1).intValue();
            else throw new RuntimeException("Cannot compare "+ o1.getClass().getSimpleName());
            if (o2 instanceof Gene) {
                v2= (start? Math.abs(((Gene) o2).getStart()): Math.abs(((Gene) o2).getEnd()));
            } else if (o2 instanceof Integer)
                v2= ((Integer) o2).intValue();
            else throw new RuntimeException("Cannot compare "+ o2.getClass().getSimpleName());

            if (v1== v2)
                return 0;
            return (v1- v2);
        }
    }


        /**
     * Comparator to compare two genes by their position in the reference genome,
     * disregarding the transcription directionality.
     */
    public static class OverlapComparator implements Comparator {

        /**
         * Compares the reference sequence sequence name, and the start/end position
         * of two loci. The strand (i.e., the transcription directionality) is not
         * considered in the comparison.
         * @param arg0 a gene locus
         * @param arg1 another gene locus
         * @return in case the two loci are from different reference sequences (e.g.,
         * chromosomes or scaffolds), the return value corresponds to the value of
         * <code>String.compareTo()</code> between the chromosome names. Otherwise
         * the method returns the distance between start respectively end positions
         * (latter in the case of equal start positions of both of the compared genes).
         */
        public int compare(Object arg0, Object arg1) {
            Gene g1= (Gene) arg0, g2= (Gene) arg1;
            int val= g1.getChromosome().compareTo(g2.getChromosome());
            if (val!= 0)
                return val;

            int b1= Math.abs(g1.getStart()), b2= Math.abs(g2.getStart()), e1= Math.abs(g1.getEnd()), e2= Math.abs(g2.getEnd());
            if (b1!= b2)
                return (b1- b2);
            return (e1- e2);
        }

        /**
         * Determines genomic overlap of two gene loci, disregarding their transcription
         * directionality.
         * @param g1 a gene locus
         * @param g2 another gene locus
         * @return <code>true</code> if both genes intersect in their genomic coordinates,
         * <code>false</code> otherwise.
         */
        public boolean overlaps(Gene g1, Gene g2) {
            int b1= Math.abs(g1.getStart()), b2= Math.abs(g2.getStart()), e1= Math.abs(g1.getEnd()), e2= Math.abs(g2.getEnd());
            if ((b1>= b2&& b1<= e2)|| (b2>= b1&& b2<= e1))
                return true;   // overlap
            return false;
        }
    }

    /**
     * Comparator for comparing arbitrary two instances of String or Gene
     * @author msammeth
     */
	public static class StableIDComparator implements Comparator {

		public int compare(Object arg0, Object arg1) {
			
			String sID1= null;
			try {
				sID1= ((Gene) arg0).getStableID();
			} catch (ClassCastException e) {
				sID1= (String) arg0;
			}
			
			String sID2= null;
			try {
				sID2= ((Gene) arg1).getStableID();
			} catch (ClassCastException e) {
				sID2= (String) arg1;
			}
			
			return sID1.compareToIgnoreCase(sID2);
		}
	}
	
	final static int getArrayPosition(String query, String[] array) {
 
		if (query== null)
			return -1;
		
		int i;
		for (i = 0; i< array.length; i++) 
			if (query.equalsIgnoreCase(array[i]))
				break;
		
		if (i>= array.length) 
			return -1;
		
		return i;
	}
	
	public Exon getExon(String exonID) {
		
		if (transcripts== null)
			return null;
		
		Exon e= null;
		for (int i = 0; i < transcripts.length; i++) {
			e= transcripts[i].getExon(exonID);
			if (e!= null)
				return e;
		}
		return e;
	}

	public DirectedRegion[] getRegionMaxCDS() {
		int min= Integer.MAX_VALUE;
		int max= Integer.MIN_VALUE;
		for (int i = 0; i < transcripts.length; i++) {
			if (!transcripts[i].isCoding())
				continue;
			if (transcripts[i].getTranslations()[0].get5PrimeEdge()< min)
				min= transcripts[i].getTranslations()[0].get5PrimeEdge();
			if (transcripts[i].getTranslations()[0].get3PrimeEdge()> max)
				max= transcripts[i].getTranslations()[0].get3PrimeEdge();
		}
		
		if (min== Integer.MAX_VALUE|| max== Integer.MIN_VALUE)
			return null;
		
		DirectedRegion reg= new DirectedRegion(min, max, getStrand());
		reg.setChromosome(getChromosome());
		reg.setID("maxCDS");
		return new DirectedRegion[] {reg};
	}
	
	public DirectedRegion[] getRegionMinCDS() {
		int min= Integer.MIN_VALUE;
		int max= Integer.MAX_VALUE;
		for (int i = 0; i < transcripts.length; i++) {
			if (!transcripts[i].isCoding())
				continue;
			if (transcripts[i].getTranslations()[0].get5PrimeEdge()> min)
				min= transcripts[i].getTranslations()[0].get5PrimeEdge();
			if (transcripts[i].getTranslations()[0].get3PrimeEdge()< max)
				max= transcripts[i].getTranslations()[0].get3PrimeEdge();
		}
		
		if (min== Integer.MIN_VALUE|| max== Integer.MAX_VALUE)
			return null;
		
		DirectedRegion reg= new DirectedRegion(min, max, getStrand());
		reg.setChromosome(getChromosome());
		reg.setID("maxCDS");
		return new DirectedRegion[] {reg};
	}
	

	public DirectedRegion[] getLocalizationRef() {
		if (localizationRef == null) {
			// TODO take into account loc reference
			localizationRef= null; 	// getRegionAllCDS();		// getRegionMaxCDS(); 
		}

		return localizationRef;
	}
	
	/**
	 * handles multiple ref regions
	 * 1= 5UTR, 2= 5UTR+CDS, 3= CDS, 4= CDS+3UTR, 5= 3UTR, 6= 5UTR+CDS+3UTR
	 * @param reg
	 * @return
	 */
	public byte[] getLocalization(DirectedRegion reg, DirectedRegion[] refRegs) {
		if (refRegs.length== 0)
			return new byte[] {0};
		
		byte[] regCodes= new byte[refRegs.length];
		for (int i = 0; i < refRegs.length; i++) {
			if (refRegs[i].contains(reg))
				regCodes[i]= 3;
			else {
				if (reg.get5PrimeEdge()< refRegs[i].get5PrimeEdge()) {
					if (reg.get3PrimeEdge()< refRegs[i].get5PrimeEdge())
						regCodes[i]= 1;
					else {
						if (reg.get3PrimeEdge()> refRegs[i].get3PrimeEdge())
							regCodes[i]= 6;
						else
							regCodes[i]= 2;
					}
				} else {
					if (reg.get5PrimeEdge()> refRegs[i].get3PrimeEdge())
						regCodes[i]= 5;
					else
						regCodes[i]= 4;
				}
			}
		}
		return regCodes;
	}
	
	
	public Exon getExon(int start, int end) {
		Exon[] exons= getExons();
		for (int i = 0; i < exons.length; i++) 
			if (exons[i].getStart()== start&& exons[i].getEnd()== end)
				return exons[i];
		return null;
	}
	
	public SpliceSite checkSpliceSite(SpliceSite ss) {
		if (spliceSites== null)
			return null;
		
		int p= Arrays.binarySearch(spliceSites, ss, new SpliceSite.PositionComparator());	// no abstract site, check don/acc !!
		if (p>= 0) {
			if (spliceSites[p].isDonor()!= ss.isDonor()) 
				System.err.println("Mismatching splice site type (donor/acceptor): "+ getStableID()+" pos "+ss.getPos());
			return spliceSites[p]; 	// bullshit, add first and remove afterwards when filtering!
		}
		return null;
	}

	public Transcript[] getTranscripts(int tlnInit) {
		Vector v= new Vector();
		for (int i = 0; i < getTranscripts().length; i++) {
			if (getTranscripts()[i].isNonCoding())
				continue;
			Translation tln= getTranscripts()[i].getTranslations()[0];
			if (tln.isOpenEnded5())	// exclude truncated 
				continue;
			if (tln.get5PrimeEdge()== tlnInit)
				v.add(getTranscripts()[i]);
		}
		return ((Transcript[]) ArrayUtils.toField(v));
	}
	
	public int countCodingSpliceForms() {
		int cnt= 0;
		for (int i = 0; i < getTranscripts().length; i++) 
			if (getTranscripts()[i].isCoding())
				++cnt;
		return cnt;
	}
	
	public Transcript[] getNonCodingTranscripts() {
		Vector v= new Vector();
		for (int i = 0; i < getTranscripts().length; i++) 
			if (getTranscripts()[i].isNonCoding())
				v.add(getTranscripts()[i]);
		return ((Transcript[]) ArrayUtils.toField(v));
	}
	
	/**
	 * @param ss
	 * @return
	 */
	public boolean addSpliceSite(SpliceSite ss, Vector<Transcript> trptV) {
		
		ss.setGene(this);
		Vector<SpliceSite> v= spliceHash.get(ss.pos);
		if (v== null) {
			v= new Vector<SpliceSite>(1,1);
//			if (ss.getPos()== 368595|| ss.getPos()== 367659)
//				System.currentTimeMillis();
			spliceHash.put(new Integer(ss.getPos()), v);
		}
		int i = 0;
		for (; i < v.size(); i++) {
			if (v.elementAt(i).isLeftFlank()== ss.isLeftFlank()|| 
					v.elementAt(i).isRightFlank()== ss.isRightFlank()) {
				if (ss.isRealSpliceSite()) {
					if (v.elementAt(i).isRealSpliceSite()) { // add to old real ss
						Vector<Transcript> vv= ssTrptHash.get(v.elementAt(i));
						if (vv== null) {
							vv= new Vector<Transcript>();
							ssTrptHash.put(v.elementAt(i), vv);
						}
						ArrayUtils.addAllUniqueSorted(vv, trptV, Transcript.getDefaultIDComparator());	// 080822: uniq add
						v.elementAt(i).setSourceType((byte) Math.min(v.elementAt(i).getSourceType(), ss.getSourceType())); 
						return false;
					} else { // replace old not real ss
						// byte oldSrc= v.elementAt(i).getSourceType();	// get here, afterwards trpts no longer in gene hash
						Vector<Transcript> vv= ssTrptHash.remove(v.elementAt(i));	
						if (vv== null) 
							vv= new Vector<Transcript>();
						ArrayUtils.addAllUniqueSorted(vv, trptV, Transcript.getDefaultIDComparator());
						ssTrptHash.put(ss,vv);	// better not the other way, V reused in Transcript.addExon(), do not have to update spliceHash: pos-based
						//ss.setSourceType((byte) Math.min(oldSrc, ss.getSourceType()));
						//ss.setSourceType(ss.getSourceType());	// TODO eliminated 081212
						v.set(i, ss);
						return true;
					}
					
				} else { // no real splice site -> edge or codon
					if (ss.isCodon()|| v.elementAt(i).isCodon()) {
						if (v.elementAt(i).isCodon()&& !ss.isCodon()) {
							Vector<Transcript> vv= ssTrptHash.remove(v.elementAt(i));
							if (vv== null) 
								vv= new Vector<Transcript>();
							ArrayUtils.addAllUniqueSorted(vv, trptV, Transcript.getDefaultIDComparator());
							ssTrptHash.put(ss,vv);	
							v.set(i, ss);
						} else
							ssTrptHash.put(ss, trptV);
						return true;
					} else if (ss.isHardEdge()&& v.elementAt(i).isSoftEdge()
							&& ((ss.isLeftFlank()&& v.elementAt(i).isLeftFlank())
									|| (ss.isRightFlank()&& v.elementAt(i).isRightFlank()))) {	// override soft edges by hard edges
						Vector<Transcript> vv= ssTrptHash.remove(v.elementAt(i));
						if (vv== null) 
							vv= new Vector<Transcript>();
						ArrayUtils.addAllUniqueSorted(vv, trptV, Transcript.getDefaultIDComparator());
						ssTrptHash.put(ss,vv);	// better not the other way, V reused in Transcript.addExon(), do not have to update spliceHash: pos-based
						// quatsch: ss source immer < v.elementAt(i), weil ersteres hard und letzteres softedge
						//ss.setSourceType((byte) Math.min(v.elementAt(i).getSourceType(), ss.getSourceType()));
						v.set(i, ss);
						return true;
					} else {
						Vector<Transcript> vv= ssTrptHash.get(v.elementAt(i));
						if (vv== null) {
							vv= new Vector<Transcript>();
							ssTrptHash.put(v.elementAt(i), vv);
						}
						ArrayUtils.addAllUniqueSorted(vv, trptV, Transcript.getDefaultIDComparator());
						//v.elementAt(i).setSourceType((byte) Math.min(v.elementAt(i).getSourceType(), ss.getSourceType()));
						//v.elementAt(i).setSourceType(v.elementAt(i).getSourceType());	// TODO eliminated 081212
						return false;
					}
				}
			}
		}
		Vector<Transcript> vv= new Vector<Transcript>(2,2);
		ArrayUtils.addAllUniqueSorted(vv, trptV, Transcript.getDefaultIDComparator());	//vv.addAll(trptV);
		ssTrptHash.put(ss, vv);
		v.add(ss);
		return true;
	}
	
	public final static Gene[] toGeneArray(Vector v) {

		if (v== null)
			return null;
		return toGeneArray(v.toArray());
	}
	
	public final static Gene[] toGeneArray(Object[] o) {

		Gene[] result= new Gene[o.length];
		for (int i = 0; i < result.length; i++) 
			result[i]= (Gene) o[i];
		return result;
	}
		public boolean isExonicPosition(int absPos) {
			
			Exon[] exons= getExons();
			for (int i = 0; i < exons.length; i++) 
				if (exons[i].contains(absPos))
					return true;
			
			return false;
		}
		
		// all exons identical, one/some missing		
		Exon[][] filterAlternativeExons(Transcript t1, Transcript t2) {
			
				// get& check
			Exon[] e1= t1.getExons();
			Exon[] e2= t2.getExons();
			if (e1== null|| e2== null|| e1.length< 1|| e2.length< 1)
				return null;
			
				// compare
			Vector e1Vec= new Vector();
			for (int i = 0; i < e1.length; i++) 
				e1Vec.add(e1[i]);
			Vector e2Vec= new Vector();
			for (int i = 0; i < e2.length; i++) 
				e2Vec.add(e2[i]);
			for (int i = 0; i < e1Vec.size(); i++)	// remove identical exon pairs 
				for (int j = 0; i < e2Vec.size(); j++) 
					if (e1Vec.elementAt(i).equals(e2Vec.elementAt(j))) {
						e1Vec.remove(i);
						e2Vec.remove(j);
						--i; --j;
					}
			
			Exon[][] concat= new Exon[2][];
			concat[0]= new Exon[e1Vec.size()];
			for (int i = 0; i < e1Vec.size(); i++) 
				concat[0][i]= (Exon) e1Vec.elementAt(i);
			for (int i = 0; i < e2Vec.size(); i++) 
				concat[1][i]= (Exon) e2Vec.elementAt(i);
            return concat;
		}
		
	String geneID= null;
		
	
	Transcript[] transcripts= null;
	DirectedRegion[] cdssRegs= null;

	public Gene(String newGeneID) {
		geneID= newGeneID;
		setID("gene");
		//setStrand(getStrand());	// hae?
	}
	
	public Gene(Species spec, String newGeneID) {
		this (newGeneID);
		setSpecies(spec);
	}
	
	public String toString() {
		return getStableID();
	}

	public String toStringSSPattern() {
		String s= "";
		for (int i = 0; getSpliceSites()!= null&& i < getSpliceSites().length; i++) {
			if (getSpliceSites()[i].isDonor())
				s+= "^";
			else
				s+= "-";
		}
		return s;
	}

	public void repairAlignmentErrors() {
		SpliceSite[] ss= getSpliceSites();
		for (int i = 0; i < ss.length; i++) {
			
		}
	}

    /**
     * Retrieves a list of exonic regions where overlapping
     * exons are merged to form one Super-Exon.
     * @return non-redundant exons of the gene
     */
    public Exon[] getSuperExons() {



        return null;
    }

        /**
       * Retrieves a non-redundant set of exons based on their flanking sites.
       * @return non-redundant exons of the gene
       */
	public Exon[] getExons() {

        if (exons== null) {

            // take map, do not change equals() method of exon
            HashMap<String, Exon> set= new HashMap<String, Exon>();

            // create non-redundant list of exons
            for (int i = 0; i < getTranscriptCount(); i++) {
                Exon[] ex= getTranscripts()[i].getExons();
                for (int j = 0; j < ex.length; j++) {
                    set.put(ex[j].toString(), ex[j]);
                }
            }

            exons= (Exon[]) ArrayUtils.toField(set.values());

        }

		return exons;
	}
	
	public void merge(Gene anotherGene) {

		// replace/add splice sites; before exons, to keep determine the good ones
		Iterator<Vector<SpliceSite>> iter= anotherGene.spliceHash.values().iterator();
		while (iter.hasNext()) {
			Vector<SpliceSite> v= iter.next();
			for (int i = 0; i < v.size(); i++) 
				addSpliceSite(v.elementAt(i), anotherGene.ssTrptHash.get(v.elementAt(i)));
		}

		//spliceSites= (SpliceSite[]) ArrayUtils.addAll(spliceSites, anotherGene.getSpliceSites());
		for (int i = 0; i < anotherGene.getTranscripts().length; i++) 
			addTranscript(anotherGene.getTranscripts()[i]);	// adds transcripts and exons
		
	}
	
	public DirectedRegion[] getExonicRegions() {
	
		DirectedRegion[] superExons= DirectedRegion.unite((DirectedRegion[]) getExons());
		// if (superExons.length== getExons().length)
		//	System.err.println("No AS");
		return superExons;
	}
	
	public String getExonicRegionsSplicedSequence() {
		
		DirectedRegion[] superExons= getExonicRegions();
		
		String result= "";
		for (int i = 0; i < superExons.length; i++) {
			String tmp= Graph.readSequence(getSpecies(), getChromosome(), isForward(),
					Math.abs(superExons[i].getStart()), Math.abs(superExons[i].getEnd()));
			if (!isForward())
				tmp= StringUtils.reverseComplement(tmp);
			result+= tmp;
		}
		return result;
	}
	
	public int[] getExonicRegionsSSCoords() {
		SpliceSite[] ss= getSpliceSites();
		Comparator compi= new SpliceSite.PositionComparator();
		Arrays.sort(ss, compi);
		
		int[] result= new int[ss.length];
		DirectedRegion[] sExons= getExonicRegions();
		int j= 0;
		int pos= 0;
		for (int i = 0; j< ss.length&& i < sExons.length; i++) {
			while (j< ss.length&& sExons[i].contains(ss[j].getPos())) {	// check for complained SSs
				int aPos= pos+ ss[j].getPos()- sExons[i].get5PrimeEdge();
					// forward/rev does not matter in this case since sequences are already reversed
				result[j++]= aPos;
			}
			pos+= sExons[i].getLength();
		}
	
		return result;
	}
	
	/**
	 * 
	 * @param seqs	return by parameter !!!
	 * @return
	 * @deprecated too mem-intensive
	 */
	int align_qalign(String[] seqs) {
//		
//		QAlign qalign= new QAlign();
//		try {
//			qalign.setAll(
//			    0,	// weighting tree
//			    1,	// output console
//			    1,	// simultaneous alignment
//			    seqs,
//			    CostTable.DNARNA,
//			    5,	// minimal epsilon
//			    null,
//			    null,
//			    null);
//			qalign.run();
//			seqs[0]= qalign.getSimultaneousLayout()[0];
//			seqs[1]= qalign.getSimultaneousLayout()[1];
//		} catch (CancelException e) {
//			e.printStackTrace();
//		}
//		
//		return qalign.getCost();
		return -1;
	}
	
	/**
	 * @return
	 */
	public String getLocusID() {
		//return geneID;
		StringBuffer sb= new StringBuffer(getChromosome());
		sb.append(":");
		int start= Math.abs(getStart());
		int end= Math.abs(getEnd());
		if (start> end) {
			int h= start;
			start= end;
			end= h;
		}
		sb.append(Integer.toString(start));
		sb.append("-");
		sb.append(Integer.toString(end));
		sb.append(getWatsonCrickStrandSymbol());
		return sb.toString();
	}

    public String getGeneID() {
        HashSet<String> geneIDs = new HashSet<String>();
        for (Transcript t : transcripts) {
            String id = t.getAttribute(GFFObject.GENE_ID_TAG).toString();
            if (!geneIDs.contains(id))
                geneIDs.add(id);
        }
        StringBuilder out = new StringBuilder();
        for (String id : geneIDs) {
            if (!out.toString().isEmpty())
                out.append(",");
            out.append(id);
        }
        return out.toString();
    }

    public String getGeneID(Transcript[] transcripts) {
        HashSet<String> geneIDs = new HashSet<String>();
        for (Transcript t : transcripts) {
            String id = t.getAttribute(GFFObject.GENE_ID_TAG).toString();
            if (!geneIDs.contains(id))
                geneIDs.add(id);
        }
        StringBuilder out = new StringBuilder();
        for (String id : geneIDs) {
            if (!out.toString().isEmpty())
                out.append(",");
            out.append(id);
        }
        return out.toString();
    }

    public String getGeneID(String tid) {
        for (Transcript t : transcripts) {
            if (t.getTranscriptID().toString().equals(tid)) {
                return t.getAttribute(GFFObject.GENE_ID_TAG).toString();
            }
        }
        return "";
    }
	
	/**
	 * @return
	 */
	public Transcript[] getTranscripts() {
		return transcripts;
	}
	
	public Transcript[] getTranscripts(Exon e) {
		if(exonHash== null)
			return null;
		return exonHash.get(e);
	}
	
	public int getTranscriptCount() {
		if (transcripts== null)
			return 0;
		return transcripts.length;
	}
	
	/**
	 * gets non-redundant set of coding exons
	 * @param completelyCoding
	 * @return
	 */
	public Exon[] getCodingExons(boolean completelyCoding) {
		Exon[] ex= getExons();
		Vector<Exon> resEx= new Vector<Exon>();
		for (int i = 0; i < ex.length; i++) {
			if (completelyCoding&& ex[i].isCodingSomewhere5Prime()&& 
					ex[i].isCodingSomewhere3Prime())
				resEx.add(ex[i]);
			if ((!completelyCoding)&& ex[i].overlapsCDS())
				resEx.add(ex[i]);
		}
		return (Exon[]) ArrayUtils.toField(resEx);
	}
	
	public int getTranscriptNbFromSource(String source) {
		int nb= 0;
		for (int i = 0; i < transcripts.length; i++) {
			if (transcripts[i].getTranscriptID().contains(source))
				++nb;
		}
		return nb;
	}
	
	/**
	 * @deprecated does not work that easy for constitutive exons, 
	 * bette take non-involvement in ASVariation as criterion
	 * @param constitutive
	 * @return
	 */
	public Exon[] getExons(boolean constitutive) {
		
		Vector v= new Vector();
		for (int i = 0; i < transcripts.length; i++) {
			for (int j = 0; j < transcripts.length; j++) {
				Exon[] ex1= transcripts[i].getExons();
				for (int k = 1; k < ex1.length- 1; k++) {	// no terminal exons for the comparison
					if (!transcripts[j].contains(ex1[k]))
						continue;
					Exon[] ex2= transcripts[j].getExons();
					int m;
					for (m = 1; m < ex2.length- 2; m++) {
						if (ex1[k].overlaps(ex2[m])) {
							//if (ex1[k].get)
						}
					}
					if (m== ex2.length&& !constitutive)
						v.add(ex1);
				}
			}
		}
		return null;
	}
	
	HashMap addTrptHash(HashMap trptHash, SpliceSite[] schain, Transcript trpt) {
		Object[] o= trptHash.keySet().toArray();
		int i;
		for (i = 0; i < o.length; i++) {
			SpliceSite[] keyChain= (SpliceSite[]) o[i];
			if (keyChain.length!= schain.length)
				continue;
			int j;
			for (j = 0; j < keyChain.length; j++) {
				if (keyChain[j].getPos()!= schain[j].getPos())
					break;
			}
			if (j>= keyChain.length) 
				break;
		}
		Vector v= new Vector();
		if (i< o.length) {
			v= (Vector) trptHash.remove(o[i]);
		}
		v.add(trpt);
		trptHash.put(schain, v);
		return trptHash;
	}
	
	/**
	 * (marked for deletion)
	 * @param spliceChains
	 * @param tt
	 * @param omitStart
	 * @param omitEnd
	 * @return
	 * @deprecated uses deprecated API
	 */
	static Vector tokenizeASEvents(SpliceSite[][] spliceChains, Transcript[] tt, boolean omitStart, boolean omitEnd) {
/*		
			// find splice sites common in ALL transcripts
		SpliceSite.PositionComparator c= new SpliceSite.PositionEqualSSTypeComparator();	// here! SpliceSite.PositionComparator
		Vector tokens= new Vector();
		for (int i = 0; i < spliceChains[0].length; i++) {
			AbstractSite s= spliceChains[0][i];
			int[] pos= new int[spliceChains.length];
			pos[0]= i;
			int j;
			for (j = 1; j < spliceChains.length; j++) {
				int p= Arrays.binarySearch(spliceChains[j], s, c);
				if (p< 0)
					break;	// as soon as not found in one schain, give up..
				else
					pos[j]= p;
			}
			if (j>= spliceChains.length)	// found
				tokens.add(pos);
		}
		
			// tokenize by these conserved splice sites
		int[] posOld= new int[spliceChains.length];
		for (int i = 0; i < posOld.length; i++) 
			posOld[i]= (-1);	// anchor before sequence start	
		int[] pos;
		Vector result= new Vector();
		for (int i = 0; i < tokens.size(); i++) {
			pos= (int[]) tokens.elementAt(i);
			if (omitStart&& i== 0) {	// omit splice chains containing start
				posOld= pos;
				continue;
			}
			SpliceSite[][] ass= new SpliceSite[pos.length][];
			for (int j = 0; j < ass.length; j++) {
				int len= pos[j]- posOld[j]- 1;
				if (len< 0)
					break;
				ass[j]= new SpliceSite[len];
				for (int k = 0; k < ass[j].length; k++) { 
					ass[j][k]= spliceChains[j][posOld[j]+ k+ 1];
				}
			}
			int j;
			for (j = 0; j < ass.length; j++) 	// at least one transcript needs to provide alternative SSs
				if (ass[j]!= null&& ass[j].length> 0) 
					break;
			
			if (j< ass.length) {
				ASVariation var= new ASVariation(tt[0], tt[1], ass[0], ass[1]);
				AbstractSite flank5= null;
				if (i== 0)
					flank5= null;
				else
					flank5= spliceChains[j][posOld[j]];
				AbstractSite flank3= spliceChains[j][pos[j]];
				var.setAnchors(flank5, flank3);
				result.add(var);
				
			}
			
			posOld= pos;
		}
		
		if (!omitEnd) {
			pos= new int[spliceChains.length];		// anchor after sequence end
			for (int i = 0; i < pos.length; i++) 
				pos[i]= spliceChains[i].length;
			SpliceSite[][] ass= new SpliceSite[pos.length][];
			for (int j = 0; j < ass.length; j++) {
				ass[j]= new SpliceSite[pos[j]- posOld[j]- 1];
				for (int k = 0; k < (pos[j]- posOld[j]- 1); k++) 
					ass[j][k]= spliceChains[j][posOld[j]+ k+ 1];				
			}
			int j;
			for (j = 0; j < ass.length; j++) 	// at least one transcript needs to provide alternative SSs 
				if (ass[j].length> 0) 
					break;
			if (j< ass.length) {
				ASVariation var= new ASVariation(tt[0], tt[1], ass[0], ass[1]);
				AbstractSite flank5= null;
				if (tokens.size()== 0)
					flank5= null;
				else
					flank5= spliceChains[j][posOld[j]];
				AbstractSite flank3= null;
				var.setAnchors(flank5, flank3);
				result.add(var);
				
			}

		}

		return result;
*/		
		return null;
	}

	static Vector getASVariations(SpliceSite[][] spliceChains) {
		
			// find splice sites common in ALL transcripts
		SpliceSite.PositionComparator c= new SpliceSite.PositionEqualSSTypeComparator();	// here! SpliceSite.PositionComparator
		for (int i = 0; i < spliceChains.length; i++) 	// just to make sure
			Arrays.sort(spliceChains[i], c);
		
		Vector tokens= new Vector();
		for (int i = 0; i < spliceChains[0].length; i++) {
			SpliceSite s= spliceChains[0][i];
			int[] pos= new int[spliceChains.length];
			pos[0]= i;
			int j;
			for (j = 1; j < spliceChains.length; j++) {
				int p= Arrays.binarySearch(spliceChains[j], s, c);
				if (p< 0)
					break;	// as soon as not found in one schain, give up..
				else
					pos[j]= p;
			}
			if (j>= spliceChains.length)	// found
				tokens.add(pos);
		}
		
			// tokenize by these conserved splice sites
		SpliceSite[][] ss= new SpliceSite[spliceChains.length][];
		int[] posOld= new int[spliceChains.length];
		for (int i = 0; i < posOld.length; i++) 
			posOld[i]= (-1);	// anchor before sequence start	
		int[] pos;
		Vector result= new Vector();
		for (int i = 0; i < tokens.size(); i++) {
			pos= (int[]) tokens.elementAt(i);
						
			SpliceSite[][] ass= new SpliceSite[pos.length][];
			for (int j = 0; j < ass.length; j++) {
				int len= pos[j]- posOld[j]- 1;
				if (len< 0)
					break;
				ass[j]= new SpliceSite[len];
				for (int k = 0; k < ass[j].length; k++) { 
					ass[j][k]= spliceChains[j][posOld[j]+ k+ 1];
				}
			}
			for (int j = 0; j < ass.length; j++) 	// at least one transcript needs to provide alternative SSs
				if (ass[j]!= null&& ass[j].length> 0) {
					result.add(ass);
					break;
				}
			posOld= pos;
		}
		
			// last
		pos= new int[spliceChains.length];		// anchor after sequence end
		for (int i = 0; i < pos.length; i++) 
			pos[i]= spliceChains[i].length;
		SpliceSite[][] ass= new SpliceSite[pos.length][];
		for (int j = 0; j < ass.length; j++) {
			ass[j]= new SpliceSite[pos[j]- posOld[j]- 1];
			for (int k = 0; k < (pos[j]- posOld[j]- 1); k++) 
				ass[j][k]= spliceChains[j][posOld[j]+ k+ 1];				
		}
		for (int j = 0; j < ass.length; j++) 	// at least one transcript needs to provide alternative SSs 
			if (ass[j].length> 0) {
				result.add(ass);
				break;
			}
		
		
		return result;
	}

	static Vector tokenizeASClusters(SpliceSite[][] spliceChains, boolean omitStart, boolean omitEnd) {
		
			// find splice sites common in ALL transcripts
		SpliceSite.PositionComparator c= new SpliceSite.PositionEqualSSTypeComparator();	// here! SpliceSite.PositionComparator
		for (int i = 0; i < spliceChains.length; i++) 	// just to make sure
			Arrays.sort(spliceChains[i], c);
		
		Vector tokens= new Vector();
		for (int i = 0; i < spliceChains[0].length; i++) {
			SpliceSite s= spliceChains[0][i];
			int[] pos= new int[spliceChains.length];
			pos[0]= i;
			int j;
			for (j = 1; j < spliceChains.length; j++) {
				int p= Arrays.binarySearch(spliceChains[j], s, c);
				if (p< 0)
					break;	// as soon as not found in one schain, give up..
				else
					pos[j]= p;
			}
			if (j>= spliceChains.length)	// found
				tokens.add(pos);
		}
		
			// tokenize by these conserved splice sites
		SpliceSite[][] ss= new SpliceSite[spliceChains.length][];
		int[] posOld= new int[spliceChains.length];
		for (int i = 0; i < posOld.length; i++) 
			posOld[i]= (-1);	// anchor before sequence start	
		int[] pos;
		Vector result= new Vector();
		for (int i = 0; i < tokens.size(); i++) {
			pos= (int[]) tokens.elementAt(i);
			if (omitStart&& i== 0) {	// omit splice chains containing start
				posOld= pos;
				continue;
			}
			SpliceSite[][] ass= new SpliceSite[pos.length][];
			for (int j = 0; j < ass.length; j++) {
				int len= pos[j]- posOld[j]- 1;
				if (len< 0)
					break;
				ass[j]= new SpliceSite[len];
				for (int k = 0; k < ass[j].length; k++) { 
					ass[j][k]= spliceChains[j][posOld[j]+ k+ 1];
				}
			}
			for (int j = 0; j < ass.length; j++) 	// at least one transcript needs to provide alternative SSs
				if (ass[j]!= null&& ass[j].length> 0) {
					result.add(ass);
					break;
				}
			posOld= pos;
		}
		
		if (!omitEnd) {
			pos= new int[spliceChains.length];		// anchor after sequence end
			for (int i = 0; i < pos.length; i++) 
				pos[i]= spliceChains[i].length;
			SpliceSite[][] ass= new SpliceSite[pos.length][];
			for (int j = 0; j < ass.length; j++) {
				ass[j]= new SpliceSite[pos[j]- posOld[j]- 1];
				for (int k = 0; k < (pos[j]- posOld[j]- 1); k++) 
					ass[j][k]= spliceChains[j][posOld[j]+ k+ 1];				
			}
			for (int j = 0; j < ass.length; j++) 	// at least one transcript needs to provide alternative SSs 
				if (ass[j].length> 0) {
					result.add(ass);
					break;
				}
		}
		
		return result;
	}

	/**
	 * @param i
	 */
	public void setGeneID(String i) {
		geneID= i; 
	}

	/**
	 * called during clustering process
	 * @param newExon
	 * @return
	 */
	public boolean addExon(Exon newExon, Transcript[] trpts) {

		newExon.setGene(this);	// can now be before, exon can get its transcripts via old gene 
		Transcript[] t= exonHash.get(newExon);
		if (t== null) {
			exonHash.put(newExon, trpts);
			return true;
		}
		
		Transcript[] nT= new Transcript[t.length+ newExon.getTranscripts().length];
		for (int i = 0; i < trpts.length; i++) {
			if (i< t.length)	// TODO check redundancy?? actually it cannot happen in the current gene merging
				nT[i]= t[i];
			else
				nT[i]= newExon.getTranscripts()[i-t.length];
		}
		
			// force updating with the good exons (the ones with the overriding SSs)
		exonHash.put(newExon, nT);	// replacing not necessary, redundant objects, but consistent transcript mapping
		return false;
		
	}
	
	// adjusts the first/last exon of ESTs
	// if there is a refseq exon overlapping with it, the "TSS/TTS" boundary will be adjusted to the corresponding flank there (extension or shortening)
	// if more than one ref exon is overlapping with the EST exon, take the one with the shortest distance for flank extension (if not, then shorten)
	// if there are no refseq exons overlapping, try mRNA (as above)
	// if there are no mRNA exons overlapping, take other EST exons into account
	// if overlapping with internal EST exon with a good 5' intron: trust this flank
	// if overlapping with terminal EST exon(s) -> extend all of them to the first TSS 
	public void trimESTextremities() {
		DefaultRegion[] regs= new DefaultRegion[exonHash.size()];
		Iterator<DefaultRegion> iter= exonHash.keySet().iterator();
		int cnt= 0;
		while(iter.hasNext()) 
			regs[cnt++]= iter.next();		
		Arrays.sort(regs, DefaultRegion.getDefaultPositionComparator());
		
		
		for (int i = 0; i < getTranscripts().length; i++) {
			if (!getTranscripts()[i].getSource().contains("EST"))
				continue;
			Vector<DefaultRegion> v= new Vector<DefaultRegion>();
			int p= Arrays.binarySearch(regs, getTranscripts()[i].getExons()[0], DefaultRegion.getDefaultPositionComparator());
		}
	}

    /**
     * Marks the splice sites of a locus by <code>ALTERNATIVE</code> when
     * there is at least one transcript that contains that site but does not
     * employ it for splicing. Otherwise the splice site is marked
     * <code>CONSTITUTIVE</code>.<br>
     *
     * <b>Warning:</b> complexity O(T x S) for a locus with T transcripts and
     * S splice sites.
     */
	public void markAlternativeSpliceSites() {
		SpliceSite[] sites= getSpliceSites();
		Comparator compi= new SpliceSite.PositionTypeComparator_old();
        // for each splice site
		for (int i = 0; i < sites.length; i++) {

            // efficiency, handle the trivial case first
			if (sites[i].getTranscripts().length== getTranscripts().length) {
				sites[i].setModality(SpliceSite.CONSTITUTIVE);
				continue;
            }

            // look for a transcript that shows alternative ussage
			for (int j = 0; j < getTranscripts().length; j++) {
				SpliceSite[] ss= getTranscripts()[j].getSpliceSitesAll();	// not schain
				int p= Arrays.binarySearch(ss, sites[i], compi);
				if (p< 0) {
					p= -(p+ 1);
					if (p!= 0&& p!= ss.length) {
						sites[i].setModality(SpliceSite.ALTERNATIVE);
						break;
					}
				}
			}


            if (sites[i].getModality()== SpliceSite.NOTYPE_SS)
                sites[i].setModality(SpliceSite.CONSTITUTIVE);

		}
	}
	
	/**
     * Adds a transcript to a locus.
	 * @param newTranscript the transcript to be added to the locus
	 */
	public boolean addTranscript(Transcript newTranscript) {

			// search transcript
//		for (int i = 0; transcripts!= null&& i < transcripts.length; i++) 
//			if (transcripts[i].getStableID().equalsIgnoreCase(newTranscript.getStableID()))
//				return false;

		newTranscript.gene= this;
		updateBoundaries(newTranscript);
		
		if(transcripts== null) 
			transcripts= new Transcript[] {newTranscript};
		else {
			Transcript[] nTranscripts= new Transcript[transcripts.length+ 1];
			for (int i= 0; i < transcripts.length; i++) 
				nTranscripts[i]= transcripts[i];
			nTranscripts[nTranscripts.length- 1]= newTranscript;
			newTranscript.setGene(this);
			transcripts= nTranscripts;
		}
		
		for (int i = 0;newTranscript.getExons()!= null && 
					i < newTranscript.getExons().length; i++) {
			addExon(newTranscript.getExons()[i], newTranscript.getExons()[i].getTranscripts());
			if (exonHash.getKey(newTranscript.getExons()[i])!= newTranscript.getExons()[i])
				newTranscript.getExons()[i]= (Exon) exonHash.getKey(newTranscript.getExons()[i]);	// force update 
		}
		//sites= null;	// are to be re-inited
		
		return true;
	}
	
	/**
	 * takes the transcript with the most 5' start to name the locus.
	 * dangerous when they omitted from future releases..
	 * @deprecated
	 */
	public Transcript getNameTranscript() {
		int minStart= Integer.MAX_VALUE;
		for (int i = 0; i < transcripts.length; i++) {
			if (transcripts[i].get5PrimeEdge()< minStart)
				minStart= transcripts[i].get5PrimeEdge();
		}
		Transcript t= null;
		for (int i = 0; i < transcripts.length; i++) {
			if ((transcripts[i].get5PrimeEdge()== minStart) &&(t== null ||
				(transcripts[i].getTranscriptID().compareTo(t.getTranscriptID())< 0)))
					t= transcripts[i];
		}
		return t;
	}
	
	/**
	 * takes the transcripts with the highest confidence level, and of those
	 * the longest one (yeah, length matters !!) or the one with the most exons (yeah, 
	 * quantity also matters !!!).
	 * @return
	 */
	public Transcript getReferenceTranscript() {
		HashMap<Integer, Vector<Transcript>> map= new HashMap<Integer, Vector<Transcript>>();
		for (int i = 0; i < transcripts.length; i++) {
			int x= transcripts[i].getConfidenceLevel();
			Integer key= new Integer(x);
			Vector<Transcript> v= map.get(key);
			if (v== null) {
				v= new Vector<Transcript>();
				map.put(key, v);
			}
			v.add(transcripts[i]);
		}
		
		Integer[] keys= new Integer[map.size()];
		Iterator<Integer> iter= map.keySet().iterator();
		int cnt= 0;
		while(iter.hasNext())
			keys[cnt++]= iter.next();
		
		Arrays.sort(keys);
		Vector<Transcript> v= map.get(keys[0]);
		Transcript ref= v.elementAt(0);
		for (int i = 1; i < v.size(); i++) {
			if (v.elementAt(i).getLength()> ref.getLength()||
					((v.elementAt(i).getLength()== ref.getLength())&& v.elementAt(i).getExons().length> ref.getExons().length))
				ref= v.elementAt(i);
		}
		
		return ref;
	}
	
	public void updateBoundaries(DirectedRegion reg) {
		if (strand== 0)
			setStrand(reg.getStrand());
		if (chromosome== null)
			setChromosome(reg.getChromosome());
		super.updateBoundaries(reg);
	}

	/**
	 * @param transcripts
	 */
	public void setTranscripts(Transcript[] transcripts) {
		this.transcripts= transcripts;
	}
	
	public static final String getSpeciesPfx(String ensemblString) {
		
		Pattern patty= Pattern.compile("^(\\D+)\\d+$");
		Matcher matty= patty.matcher(ensemblString);
		
		if (matty.matches())
			return matty.group(1).substring(0, matty.group(1).length()- 1);	// chop marker ("G", "T", "E",..)
		return null;	// error
	}
	
	public static final String getStableID(String pfx, char marker, int nb) {
		
		StringBuffer result= new StringBuffer(""+ nb); // the int nb
		
			// add leading '0's
		char[] c= new char[11- result.length()];
		Arrays.fill(c, '0');
		result.insert(0,c);		
		
		result.insert(0, marker);	// add marker (gene, transcript or exon)
		result.insert(0, pfx);	// add species pfx
		
		return result.toString();
	}
	
	/**
	 * @return
	 */
	public String getStableID() {
		
		return geneID;
	}	

	static final long serialVersionUID = 8933737248601221991L;


	public String printGeneStructure() {
		
		String result= getStart()+ "----";
		for (int i = 0; i < spliceSites.length; i++) {
			if (!spliceSites[i].isDonor())
				result+= "<";
			result+= spliceSites[i].getPos();
			if (spliceSites[i].isDonor())
				result+= ">";
			result+= "---";
		}
		result+= getEnd();
		
		return result;
	}

	public Exon[] getExons(int region) {
		Vector v= new Vector();
		Comparator compi= new DirectedRegion.PositionComparator();
		for (int i = 0; i < transcripts.length; i++) 
			ArrayUtils.addAllUniqueSorted(v, transcripts[i].getExons(), compi);
			
		DirectedRegion reg= getRegion(region);
		if (reg== null)
			return null;
		for (int i = 0; i < v.size(); i++) 
			if (!reg.contains((Exon) v.elementAt(i))) {
				v.remove(i--);
				break;
			}
		
		return  (Exon[]) ArrayUtils.toField(v);
	}
	
	/**
	 * just returning global regions, eg real/max/transcript utr/cds
	 * intronic/exonic arrays delegated to submethods..
	 * @param regionID
	 * @return
	 */ 
	public DirectedRegion getRegion(int regionID) {
		
		Object o= null;
		if (regionID== REGION_REAL_5UTR)
			o= getReal5UTR();
		else if (regionID== REGION_REAL_CDS)
			o= getRealCDS();
		else if (regionID== REGION_REAL_3UTR)
			o= getReal3UTR();
		else if (regionID== REGION_MAX_CDS)
			o= getMaxCDS();
		else if (regionID== REGION_COMPLETE_GENE)
			o= this;	// new DirectedRegion(getStart(), getEnd(), getStrand());
		
		return (DirectedRegion) o;
	}
	

	String chromosome = null;
	public static final String[] REGION_ID= 
	{"REGION_COMPLETE_GENE", "REGION_REAL_5UTR", "REGION_REAL_CDS", "REGION_REAL_3UTR",
		"REGION_MAX_5UTR", "REGION_MAX_CDS", "REGION_MAX_3UTR", "REGION_TRANSCRIPT_5UTR",
		"REGION_TRANSCRIPT_CDS", "REGION_TRANSCRIPT_3UTR"};
	public static final byte REGION_COMPLETE_GENE= 0, REGION_REAL_5UTR= 1, REGION_REAL_CDS= 2, REGION_REAL_3UTR= 3, REGION_MAX_5UTR= 4, REGION_MAX_CDS= 5, REGION_MAX_3UTR= 6,
		REGION_TRANSCRIPT_5UTR= 7, REGION_TRANSCRIPT_CDS= 8, REGION_TRANSCRIPT_3UTR= 9;

	/**
	 * @return
	 */
	public String getChromosome() {
		return chromosome;
	}

	/**
	 * @param string
	 */
	public void setChromosome(String string) {
//		String stringU= string.toUpperCase();
//		if (stringU.startsWith("SCAFFOLD")|| stringU.startsWith("REFTIG")
//				|| stringU.startsWith("CONTIG")|| stringU.startsWith("CHR")
//				|| stringU.startsWith("GROUP"))
//			chromosome= string;
//		else {
//			if (stringU.equals("MT"))	// correct Ensembl to GRIB jargon
//				chromosome= "M";
//			else {						// add cher
//				if (string.startsWith("0"))
//					string= string.substring(1, string.length());
//				chromosome= "chr"+ string;
//			}
//		}
		chromosome= string;
	}

    /**
     * Creates an array of sites in this gene, sorted by their position
     * and type. If the list has already been generated, it is lazily
     * returned.
     * @return an array of sites
     */
	public SpliceSite[] getSpliceSites() {

        if (spliceSites== null|| spliceSites.length< spliceHash.size()) {

            int n= 0;
            for (Vector<SpliceSite> siteV : spliceHash.values())
                n+= siteV.size();
            spliceSites= new SpliceSite[n];
            n= 0;
            for (Vector<SpliceSite> siteV : spliceHash.values())
                for (SpliceSite spliceSite : siteV)
                    spliceSites[n++]= spliceSite;

            Arrays.sort(spliceSites, new SpliceSite.PositionTypeComparator());
        }

		return spliceSites;	// evtl lazy extractions from exons
	}
	
	public Vector<SpliceSite> getSpliceSites(Integer pos) {
		return spliceHash.get(pos);
	}

	private static int uniqueID= 0;
	
	public static String getUniqueID() {
		
		String nb= Integer.toString(uniqueID++);
		for (int i = 0; i < 13- nb.length(); i++) 
			nb= "0"+ nb;
		
		return nb;
	}
	
	public static void resetUniqueID() {
		uniqueID= 0;
	}

	public boolean isConstruct() {
		return construct;
	}

	public void setConstruct(boolean construct) {
		this.construct = construct;
	}

	public void setLocalizationRef(DirectedRegion[] localizationRef) {
		this.localizationRef = localizationRef;
	}

	public void setLocalizationRef(HashMap<String, String> trptMap) {
		for (int i = 0; i < getTranscripts().length; i++) 
			if (trptMap.get(getTranscripts()[i].getTranscriptID())!= null
					&& getTranscripts()[i].isCoding())
				addLocalizationRef(getTranscripts()[i].getTranslations()[0]);
	}

	public void addLocalizationRef(DirectedRegion locRef) {
		if (localizationRef== null)
			localizationRef= new DirectedRegion[] {locRef};
		else {
			DirectedRegion[] h= localizationRef;
			localizationRef= new DirectedRegion[h.length+1];
			for (int i = 0; i < h.length; i++) 
				localizationRef[i]= h[i];
			localizationRef[localizationRef.length-1]= locRef;
		}
	
			
	}

}

