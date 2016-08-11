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
 * Created on Mar 3, 2005
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
package barna.model;

//import genome.NMDSimulator;
//import genome.SpliceSiteConservationComparator;

import barna.commons.log.Log;
import barna.commons.utils.ArrayUtils;
import barna.model.tools.NMDSimulator;

import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Vector;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

//import genome.tools.MyArray;
//import genome.tools.ENCODE;

/**
 * 
 * 
 * @author micha
 */
public class Transcript extends DirectedRegion {

	public static final String GFF_FEATURE_TRANSCRIPT= "transcript";

	private static Comparator exonOrderCompi= null; 
	
	public static boolean removeGaps= true;
	
	public static byte maxLengthIntronIsGap = 4;
	public static final String ID_TAG= "id_tag";
	public static final String[] POLY_A_SITES= new String[] {
		"AATAAA","ATTAAA", 	// canonical 
		"AGTAAA", "TATAAA", "CATAAA", "GATAAA", "AATATA", "AATACA", "AATAGA", "AATGAA", "ACTAAA", // Beaudoing et al. 2000  
		"AACAAA", "TTTAAA"	// F. Lopez and D. Gautheret, unpubl.
	};
	public static final String GTF_ATTRIBUTE_NMD_TAG= "nmd";
	public static final String GTF_ATTRIBUTE_3COMPLETE_TAG= "complete3",
		GTF_ATTRIBUTE_TAG_ENDDIST= "end_dist", GTF_ATTRIBUTE_TAG_LAST_EXON= "last_exon",
		GTF_ATTRIBUTE_SOURCE= "source";
	
	public static final byte ID_SRC_NEG= -1, ID_SRC_REFSEQ= 0, ID_SRC_UCSC= 1, ID_SRC_MRNA= 2, ID_SRC_EST= 3, ID_SRC_UNDEFINED= 4, ID_SRC_MOST_INCONFIDENT= Byte.MAX_VALUE;
	
	public static class SpliceChainComparator implements Comparator {
		public int compare(Object o1, Object o2) {
//			SpliceSite[] sc1= ((Transcript) o1).getSpliceChain();
//			SpliceSite[] sc2= ((Transcript) o2).getSpliceChain();
//			
//			if (sc1== null|| sc1.length< sc2.length)
//				return -1;
//			if (sc2== null|| sc2.length< sc1.length)
//				return 1;
//			
//			for (int i = 0; i < sc2.length; i++) {
//				if (sc1[i].getPos()!= sc2[i].getPos())
//					return -1;
//			}
			return 0;
		}
	}
	
	
	static final long serialVersionUID = 2863324463934791891L;
	static byte edgeConfidenceLevel= ID_SRC_MOST_INCONFIDENT; // now default trust all, ID_SRC_NEG;	// the global max confidence level for trusting transcript edges
	byte nmd= 0;
	byte srcType= ID_SRC_UNDEFINED;	
	String source= null; 
	DirectedRegion[] introns= null;
	
	public static int REGION_5UTR= 1;
	public static int REGION_3UTR= 2;
	public static int REGION_CDS= 3;
	public static int REGION_COMPLETE_TRANSCRIPT= 4;
	public static String[] TYPE_TO_ID= new String[] {"RefSeq", "knownGene", "mRNA", "EST", "Undefined"};
	
	public final String[][] CONFIDENCE_LEVELS= new String[][] {{"refGene"}, {"knownGene"}, {"mrna"},{"EST"},{"Undefined"}};

	public static String getSourceStr(byte level) {
		if (level< 0)
			return "negative";
		if (level >= TYPE_TO_ID.length)
			return "most_inconfident";
		
		return TYPE_TO_ID[level];
	}
	
	static Comparator<Transcript> defaultIDComparator= null;
	
	public static Comparator<Transcript> getDefaultIDComparator() {
		if (defaultIDComparator == null) {
			defaultIDComparator = new Comparator<Transcript>() {
				//@Override
				public int compare(Transcript o1, Transcript o2) {
					return o1.getTranscriptID().compareTo(o2.getTranscriptID());
				}
			};
		}

		return defaultIDComparator;
	}
	/**
	 * Gives 0-based position relative to transcript start.
	 * Negative values indicate position before the transcription start site,
     * counted 1-based, i.e., (-3) is 3 positions upstream of the transcript
     * start site.
	 * Values >= (transcript length) are positions after the cleavage site,
     * also counted 1based, i.e., the distance to the cleavage site can
     * by reconstructed by (exon pos) - (transcript length).
     * Consequently, the value of the transcript length is never returned by
     * the method.
	 * Integer.MIN_VALUE indicates that the position falls in an intron.
	 * 
	 * @param pos genomic position, positive or negative
	 * @return transcript-relative position of the genomic location,
     *          or <code>Integer.MIN_VALUE</code> if the position is in an
     *          intron
	 */
	public int getExonicPosition(int pos) {
		
		// safety correction
		if (pos> 0&& !isForward())
			pos= -pos;
		
		// before transcription start, 1-based distance
		if (pos< exons[0].get5PrimeEdge())
			return pos- exons[0].get5PrimeEdge();
		
		// find containing exon, 0-based coordinates
		int dist= 0;
		for (int x = 0; x < exons.length; x++) {
			if (exons[x].get5PrimeEdge()> pos) 
				return Integer.MIN_VALUE;	// pos was in preceeding intron
			if (exons[x].contains(pos))
                return dist+ (pos- exons[x].get5PrimeEdge());
			else
				dist+= exons[x].getLength();
		}
	
		// after transcript end, 1-based distance
		return dist+ (pos- exons[exons.length- 1].get3PrimeEdge());
	}

    /**
     * Gives 0-based position relative to transcript start of the next
     * transcriptomic coordinate after the given genomic position.
     *
     * @param pos genomic position, positive or negative
     * @return transcript-relative position of the genomic location,
     *          or <code>(-1)</code> if the next position is after
     *          the transcript's 3'-end
     */
    public int getNextExonicPosition(int pos) {

        // safety correction
        if (pos> 0&& !isForward())
            pos= -pos;

        // before transcription start, next position is 0
        if (pos< exons[0].get5PrimeEdge())
            return 0;

        // find containing exon, 0-based coordinates
        int dist= 0;
        for (int x = 0; x < exons.length; x++) {
            if (exons[x].get5PrimeEdge()> pos)
                return dist+ 1;	// pos was in preceeding intron
            if (exons[x].contains(pos)) {
                if (x== exons.length&& pos== exons[x].get3PrimeEdge())
                    return (-1);    // position after pos is outside transcript
                return dist+ (pos- exons[x].get5PrimeEdge());
            }
            // else
            dist+= exons[x].getLength();
        }

        return (-1);    // pos after transcript end
    }

    /**
     * Gives 0-based position relative to transcript start of the next
     * transcriptomic coordinate after the given genomic position.
     *
     * @param pos genomic position, positive or negative
     * @return transcript-relative position of the genomic location,
     *          or <code>(-1)</code> if the next position is before
     *          the transcript's 5'-end
     */
    public int getPrevExonicPosition(int pos) {

        // safety correction
        if (pos> 0&& !isForward())
            pos= -pos;

        // before transcription start, next position is 0
        if (pos< exons[0].get5PrimeEdge())
            return (-1);    // pos is before 5'-start

        // find containing exon, 0-based coordinates
        int dist= 0;
        for (int x = 0; x < exons.length; x++) {
            if (exons[x].get5PrimeEdge()> pos)
                return dist;	// pos was in preceeding intron
            if (exons[x].contains(pos)) {
                return dist+ (pos- exons[x].get5PrimeEdge())- 1;
            }
            // else
            dist+= exons[x].getLength();
        }

        return (dist- 1);    // pos after transcript end
    }


    /**
     * Gives 0-based position relative to transcript cds start of the next
     * transcriptomic coordinate after the given genomic position.
     *
     * @param pos genomic position, positive or negative
     * @return transcript cds-relative position of the genomic location,
     *          or <code>(-1)</code> if the next position is before
     *          the transcript's 5'-end
     */
    public int getPrevCDSPosition(int pos) {

        // safety correction
//        if (pos> 0&& !isForward())
//            pos= -pos;

        // before transcription start, next position is 0
        if (pos< getCDSRegions()[0].get5PrimeEdge())
            return (-1);    // pos is before 5'-start

        // find containing exon, 0-based coordinates
        int dist= 0;
        for (int x = 0; x < getCDSRegions().length; x++) {
            if (getCDSRegions()[x].get5PrimeEdge()> pos)
                return dist;	// pos was in preceeding intron
            if (getCDSRegions()[x].contains(pos)) {
                return dist+ (pos- getCDSRegions()[x].get5PrimeEdge())- 1;
            }
            // else
            dist+= getCDSRegions()[x].getLength();
        }

        return (dist- 1);    // pos after transcript end
    }


    /**
	 * @param exonPos 0-based coordinate
	 * @return genomic pos 1-based
	 */
	public int getGenomicPosition(int exonPos) {
		
			// find containing exon
		int x;
		int dist= 0;
		for (x = 0; dist<= exonPos&& x < exons.length; x++) 
			dist+= exons[x].getLength();
		if (x> 0) {
			--x;
			dist-= exons[x].getLength();
		}
		
		int genPos= exons[x].get5PrimeEdge()+ (exonPos- dist);		
		return genPos;
	}
	
	public String toString() {
		return getTranscriptID();
	}
	
	public String toStringNMD() {
//		String s= transcriptID+ "\t";
//		if ((translations!= null)&& (!translations[0].isOpenEnded()))
//			s+= translations[0].getStart()+ "\t"+ translations[0].getEnd()+ "\t";
//		else
//			s+= ".\t.\t";
//		if (predORF!= null)
//			s+= predORF.getStart()+ "\t"+ predORF.getEnd()+ "\t";
//		else
//			s+= ".\t.\t";
//		s+= "NMD"+ nmd;
//		return s;
		return null;
	}
	

	public void updateBoundaries(DirectedRegion reg) {
		super.updateBoundaries(reg);
		getGene().updateBoundaries(reg);
	}
	
	public static String toStringSChain(SpliceSite[] spliceChain) {
		String s= "{";
		for (int i = 0; spliceChain!= null&& i < spliceChain.length; i++) 
			s+= spliceChain[i].getPos()+ ",";
		if (spliceChain!= null&& spliceChain.length> 0)
			s= s.substring(0, s.length()- 1);
		s+= "}";
		return s;
	}

	Translation[] translations= null;
	SpliceSite[] spliceSites= null;	// sorted!!
	public final static Transcript[] toTranscriptArray(Vector v) {
		if (v== null)
			return null;
		Transcript[] result= new Transcript[v.size()];
		for (int i = 0; i < result.length; i++) 
			result[i]= (Transcript) v.elementAt(i);
		return result;
	}	
	
	/**
	 * 
	 * @return the first translation registered for <code>this</code> transcript.
	 */
	public Translation getDefaultTranslation() {
		if (translations== null) {
			translations= new Translation[]{new Translation(this)};
		}
		return translations[0];
	}
	
	public SpliceSite[] getSpliceSitesAll() {
		return getSpliceSitesAll(false);
	}
	public SpliceSite[] getSpliceSitesAll(boolean includeCDS) {
		if (spliceSites== null|| includeCDS) {
			spliceSites= new SpliceSite[exons.length* 2];
			int pos= 0;
			for (int i = 0; i < exons.length; i++) {
				spliceSites[pos++]= exons[i].getAcceptor();
				spliceSites[pos++]= exons[i].getDonor();
			}
			
			// TODO: never tested
			if (includeCDS) {
				int pStart= Integer.MIN_VALUE, pStop= Integer.MIN_VALUE, n= 0;
				if (getTranslations()!= null) {
					SpliceSite codon= getTranslations()[0].getCodonStart();
					if (codon!= null&& getGene().ssTrptHash.get(codon)!= null) {
						SpliceSite ss= getTranslations()[0].getCodonStart();
						pStart= Arrays.binarySearch(spliceSites, ss, SpliceSite.getDefaultPositionTypeComparator());
						if (pStart< 0)
							pStart= -(pStart+ 1);
						++n;
					}
					codon= getTranslations()[0].getCodonStop();
					if (codon!= null&& getGene().ssTrptHash.get(codon)!= null) {
						SpliceSite ss= getTranslations()[0].getCodonStop();
						pStop= Arrays.binarySearch(spliceSites, ss, SpliceSite.getDefaultPositionTypeComparator());
						if (pStop< 0)
							pStop= -(pStop+ 1);
						if (pStart>= 0)
							++pStop;
						++n;
					}
					if (n> 0) {
						SpliceSite[] allSitesNew= Arrays.copyOf(spliceSites, spliceSites.length+ n);
						codon= getTranslations()[0].getCodonStart();
						if (codon!= null&& getGene().ssTrptHash.get(codon)!= null) {
							System.arraycopy(allSitesNew, pStart, allSitesNew, pStart+ 1, spliceSites.length- pStart);
							allSitesNew[pStart]= getTranslations()[0].getCodonStart();
						}
						codon= getTranslations()[0].getCodonStop();
						if (codon!= null&& getGene().ssTrptHash.get(codon)!= null) {
							System.arraycopy(allSitesNew, pStop, allSitesNew, pStop+ 1, allSitesNew.length- pStop- 1);
							allSitesNew[pStop]= getTranslations()[0].getCodonStop();
						}
						spliceSites= allSitesNew;
					}
				}			
				
			}
		}
		return spliceSites;
	}
	
	public SpliceSite getSpliceSite(SpliceSite ss) {
		int low = 0;
		int high = exons.length- 1;

		while (low <= high) {
		    int mid = (low + high) >>> 1;
		    Exon midVal = exons[mid];

		    if (midVal.get3PrimeEdge() < ss.getPos())
		    	low = mid + 1;
		    else if (midVal.get5PrimeEdge() > ss.getPos())
		    	high = mid - 1;
		    else {
		    	if (ss== midVal.getDonor()|| ss== midVal.getAcceptor())
		    		return ss; // key found
		    	else
		    		break;
		    }
		}
		return null;
	}
	
	public SpliceSite[] getSpliceSitesBetween(SpliceSite ss1, SpliceSite ss2) {
		return getSitesBetween(ss1, ss2, false);
	}
	public SpliceSite[] getSitesBetween(SpliceSite ss1, SpliceSite ss2, boolean includeCDS) {
		SpliceSite[] ssa= getSpliceSitesAll(includeCDS);
		int p1= Arrays.binarySearch(ssa, ss1, SpliceSite.getDefaultPositionTypeComparator());
		int p2= Arrays.binarySearch(ssa, ss2, SpliceSite.getDefaultPositionTypeComparator());
		if (p1<0) // root
			p1= -(p1+2);	// limit outside
		if (p2<0)
			p2= getSpliceSitesAll().length;
		
		++p1; --p2;
		if (p1> p2)
			return new SpliceSite[0];
		
		SpliceSite[] ss= new SpliceSite[p2-p1+1];
		for (int i = 0; i < ss.length; i++) {
			ss[i]= getSpliceSitesAll()[i+p1];
		}
		return ss;
	}
	
	public DirectedRegion[] get5UTRRegion(boolean beforeSplicing) {
		Translation[] trans= getTranslations();
		if (trans== null|| trans.length< 1)
			return null;
		if (beforeSplicing) {
			return new DirectedRegion[] {
					new DirectedRegion(
						Math.min(Math.abs(trans[0].get5PrimeEdge()), Math.abs(get5PrimeEdge())),
						Math.max(Math.abs(trans[0].get5PrimeEdge()), Math.abs(get5PrimeEdge())),
						getStrand()
					)};
		} else {
			int atg= trans[0].get5PrimeEdge();
			Vector v= new Vector();
			for (int i = 0; i < exons.length; i++) {
				if (exons[i].get3PrimeEdge()>= atg) {	// atg containing exon
					int min= exons[i].get5PrimeEdge();
					int max= atg- 1;
					if (!isForward()) {
						int h= min;
						min= max;
						max= h;
					}
					DirectedRegion dir= null;
					if (Math.abs(max)- Math.abs(min)>= 0) {	// non-empty region
						dir= new DirectedRegion(min, max, getStrand());
						dir.setChromosome(getChromosome());
						v.add(dir);
					}
					break;
				}
				DirectedRegion dir= new DirectedRegion(exons[i].getStart(), exons[i].getEnd(), exons[i].getStrand());
				dir.setChromosome(getChromosome());
				v.add(dir);
			}
			return (DirectedRegion[]) ArrayUtils.toField(v);
		}
	}

	public int get5UTR(boolean beforeSplicing) {
		Translation[] trans= getTranslations();
		if (trans== null|| trans.length< 1)
			return -1;
		if (beforeSplicing) {
			return (trans[0].get5PrimeEdge()- get5PrimeEdge()+ 1);
		} else {
			int atg= trans[0].get5PrimeEdge();
			int sum= 0;
			for (int i = 0; i < exons.length; i++) {
				if (exons[i].get3PrimeEdge()> atg) {	// atg containing exon
					sum+= atg- exons[i].get5PrimeEdge()+ 1;	// add rest
					break;
				}
				sum+= exons[i].getLength();
			}
			return sum;
		}
	}
	
	/**
	 * what is the difference to exons?
	 * @return
	 */
	public DirectedRegion[] getExonicRegions() {
		if (exons== null)
			return null;
		DirectedRegion[] regs= new DirectedRegion[exons.length];
		for (int i = 0; i < regs.length; i++) {
			regs[i]= new DirectedRegion(exons[i].getStart(), exons[i].getEnd(), exons[i].getStrand());
			regs[i].setChromosome(getChromosome());
			regs[i].setSpecies(getSpecies());
		}
		return regs;
	}
	
	public String getSplicedSequence() {
		DirectedRegion[] regs= getExonicRegions();	// not sorted
		if (regs== null)
			return "";
		java.util.Arrays.sort(regs, new DirectedRegion.DirectedPositionComparator());
		StringBuffer sb= new StringBuffer();
		for (int i = 0; i < regs.length; i++) { 
			String s= Graph.readSequence(regs[i]);
			if (s!= null)
				sb.append(s);
		}
		return sb.toString();
	}


    public String getSplicedSequence(int leftFlank, int rightFlank, String leftFlankChar, String rightFlankChars) {
        DirectedRegion[] regs= getExonicRegions();	// not sorted
        if (regs== null)
            return "";
        java.util.Arrays.sort(regs, new DirectedRegion.DirectedPositionComparator());
        StringBuffer sb= new StringBuffer();
        for (int i = 0; i < regs.length; i++) {
            String s= Graph.readSequence(regs[i]);
            if (s!= null)
                sb.append(s);
        }

        // append flanks
        if(leftFlank > 0){
            // append start
            int start = get5PrimeEdge() - leftFlank;
            String ss = "";
            ss= Graph.readSequence(
                        getSpecies(),
                        getChromosome(),
                        isForward(),
                        start,
                        get5PrimeEdge()-1
                );
            while(ss.length() < leftFlank){
                ss = leftFlankChar+ss;
            }
            sb.insert(0, ss);
        }

        if(rightFlank > 0){
            int end = get3PrimeEdge() + rightFlank ;
            String se= Graph.readSequence(
                    getSpecies(),
                    getChromosome(),
                    isForward(),
                    get3PrimeEdge(),
                    end
            );
            while(se.length() < rightFlank){
                se = se+rightFlankChars;
            }
            sb.append(se);
        }

        return sb.toString();
    }

	
	private int elength= -1;
	public void setExonicLength(int newLen) {
		elength= newLen;
	}

    /**
     * Get the length of exonic region
     * @return Sum of all exons lenght of transcript
     */
	public int getExonicLength() {
		if (elength == -1) {
			elength= 0;
			for (int i = 0; i < exons.length; i++) 
				elength+= exons[i].getLength();
		}

		return elength;
	}
	
	public int getIntronicLength() {
		int sum= 0;
		DirectedRegion[] regs= getIntrons();
		for (int i = 0; regs!= null&& i < regs.length; i++) 
			sum+= regs[i].getLength();
		return sum;
	}
		
	public int getCDSLength(boolean spliced) {
		if (translations== null|| translations.length< 1)
			return -1;
		if (spliced) {
			DirectedRegion[] reg= getCDSRegions();
			int sum= 0;
			for (int i = 0; i < reg.length; i++) 
				sum+= reg[i].getLength();
			return sum;
		} else {
			return translations[0].get3PrimeEdge()- translations[0].get5PrimeEdge()+ 1;
		}
	}
	
	public int get5UTRLength(boolean spliced) {
		if (translations== null|| translations.length< 1)
			return -1;
		if (spliced) {
			DirectedRegion[] reg= get5UTRRegion(false);
			int sum= 0;
			for (int i = 0; reg!= null&& i < reg.length; i++) 
				sum+= reg[i].getLength();
			return sum;
		} else {
			return translations[0].get5PrimeEdge()- getStart();
		}
	}
	
	/**
	 * @deprecated repair
	 * @return
	 */
	public String get5UTRSequence() {
//		DirectedRegion[] reg= get5UTRRegion(false);
//		
//		String region= "";
//		int min= 0, max= 0;
//		for (int j = 0; reg!= null&& j < reg.length; j++) {
//			if (j== 0)
//				min= Math.abs(reg[j].get5PrimeEdge());
//			if (j== reg.length- 1)
//				max= Math.abs(reg[j].get3PrimeEdge());
//			DirectedRegion regSav= (DirectedRegion) reg[j].clone();
//			reg[j]= ENCODE.convertToEncodeCoord(reg[j]);
//			if (reg[j]== null) {
//				System.err.println("outside encode, skipping transcript "+getTranscriptID());
//				return null;
//			}
//			String s= SpliceSiteConservationComparator.getSubstring( 
//					reg[j].getID(),
//					"human", 
//					Math.abs(reg[j].getStart()), 
//					Math.abs(reg[j].getEnd())+ 1);	// end excl
//			if (!reg[j].isForward())
//				s= gphase.tools.Arrays.reverseComplement(s);
//			region+= s;
//		}
//		
//		return region;
		return null;
	}

	/**
	 * @deprecated repair
	 * @param regCode
	 * @return
	 */
	public String getSequence(int regCode) {
//		DirectedRegion[] reg= null;
//		if (regCode== REGION_5UTR)
//			reg= get5UTRRegion(false);
//		else if (regCode== REGION_COMPLETE_TRANSCRIPT)
//			reg= getExonicRegions();
//		else if (regCode== REGION_CDS)
//			reg= getCDSRegions();
//		
//		String region= "";
//		int min= 0, max= 0;
//		for (int j = 0; reg!= null&& j < reg.length; j++) {
//			if (j== 0)
//				min= Math.abs(reg[j].get5PrimeEdge());
//			if (j== reg.length- 1)
//				max= Math.abs(reg[j].get3PrimeEdge());
//			DirectedRegion regSav= (DirectedRegion) reg[j].clone();
//			reg[j]= ENCODE.convertToEncodeCoord(reg[j]);
//			if (reg[j]== null) {
//				System.err.println("outside encode, skipping transcript "+getTranscriptID());	// skip fragments
//				return null;	// continue?
//			}
//			String s= SpliceSiteConservationComparator.getSubstring( 
//					reg[j].getID(),
//					"human", 
//					Math.abs(reg[j].getStart()), 
//					Math.abs(reg[j].getEnd())+ 1);	// end excl
//			if (!reg[j].isForward())
//				s= gphase.tools.Arrays.reverseComplement(s);
//			region+= s;
//		}
//		
//		return region;
		return null;
	}
	

	
	public static SpliceSite getSuccSpliceSite(SpliceSite[] ss, int pos) {
		int p= Arrays.binarySearch(ss, new Integer(pos+1), new AbstractSite.PositionToSpliceSiteComparator());	// abstract site for only comparing pos
		if(p>= 0) {	// found ?!
			return ss[p];
		}
		p= -(p+ 1);	// insertionpoint= point of successor
		if (p< ss.length&& ss[p].getPos()> pos)
			return ss[p];
		else
			return null;
	}
	
	public static int getSuccPos(int[] ss, int pos) {
		int p= Arrays.binarySearch(ss, pos+1);	// abstract site for only comparing pos
		if(p>= 0) {
			return ss[p];
		}
		p= -(p+ 1);	// insertionpoint= point of successor
		if (p< ss.length&& ss[p]> pos)
			return ss[p];
		else
			return 0;
	}
	
	public static SpliceSite getPredSpliceSite(SpliceSite[] ss, int pos) {
		int p= Arrays.binarySearch(ss, new Integer(pos-1), new AbstractSite.PositionToSpliceSiteComparator());	// abstract site for only comparing pos
		if(p>= 0) {	// found ?!
			return ss[p];
		}
		p= -(p+ 1)- 1;	// insertion point- 1= point of predecessor
		if (p>= 0&& p< ss.length&& ss[p].getPos()< pos)
			return ss[p];
		else
			return null;
	}
	
	public SpliceSite[] getLastUTRIntron() {
		if (translations== null|| exons.length< 3|| exons[1].get3PrimeEdge()>= translations[0].get5PrimeEdge())
			return null;
		int i;
		for (i = 0; i < exons.length; i++) 
			if (exons[i].get3PrimeEdge()>= translations[0].get5PrimeEdge())
				break;
		SpliceSite[] ss= new SpliceSite[2];
		ss[0]= exons[i-1].getDonor();
		ss[1]= exons[i].getAcceptor();
		return ss;
	}
	
	/**
	 * @deprecated repair
	 * @param intron
	 * @return
	 */
	public Exon getUpstreamExon(DirectedRegion intron) {
//		for (int i = 0; i < exons.length; i++) {
//			if (intron.isUpstreamRegion(exons[i]))
//				return exons[i];
//		}
		return null;
	}
	
	/**
	 * @deprecated repair
	 * @param intron
	 * @return
	 */
	public Exon getDownstreamExon(DirectedRegion intron) {
//		for (int i = 0; i < exons.length; i++) {
//			if (intron.isDownstreamRegion(exons[i]))
//				return exons[i];
//		}
		return null;
	}
	
	public static int getPredPos(int[] ss, int pos) {
		int p= Arrays.binarySearch(ss, pos-1);	// abstract site for only comparing pos
		if(p>= 0) {
			return ss[p];
		}
		p= -(p+ 1)- 1;	// insertion point- 1= point of predecessor
		if (p>= 0&& p< ss.length&& ss[p]< pos)
			return ss[p];
		else
			return 0;
	}
	
	public static SpliceSite getSpliceSiteByPos(SpliceSite[] ss, int pos) {
		int p= Arrays.binarySearch(ss, new Integer(pos), new AbstractSite.PositionToSpliceSiteComparator());	// abstract site for only comparing pos
		if(p>= 0) 
			return ss[p];
		return null;
	}
	
	public static int getPos(int[] ss, int pos) {
		int p= Arrays.binarySearch(ss, pos);	// abstract site for only comparing pos
		if(p>= 0) 
			return ss[p];
		return -1;
	}
	
	public boolean addTranslation(Translation newTrans) {

		// search transcript for same translation, not necessary
		
			// new transcipt array
		if (translations== null) {
			translations= new Translation[] {newTrans};
		} else {	// add
			Translation[] nTranslations= new Translation[translations.length+ 1];
			for (int i= 0; i < translations.length; i++) 
				nTranslations[i]= translations[i];
			nTranslations[nTranslations.length- 1]= newTrans;
			translations= nTranslations;
		}
		
		return true;
	}
	
	public void setTranslations(Translation[] tln) {
		translations= tln;
	}
	
	public Translation[] getTranslation(int start, int end) {
		if (translations== null)
			return null;
		Vector result= new Vector();
		for (int i = 0; i < translations.length; i++) {
			if ((translations[i].get5PrimeEdge()<= start)&&
				(translations[i].get3PrimeEdge()>= end))
				result.add(translations[i]);
		}
		
		Object o= ArrayUtils.toField(result);
		if (o== null|| result.size()< 1)
			return null;
		return (Translation[]) o;
	}
	
	public Translation[] getTranslations() {
		return translations;
	}
	
	public boolean isNonCoding() {
		return (translations== null|| translations.length== 0);
	}
	/**
	 * see  Lopez1 et al., RNA 2006
	 */
	public boolean isInternallyPrimed() {
		final int FLANK_REGION= 50;
		int end3= getExons()[getExons().length- 1].get3PrimeEdge();
		DirectedRegion flankReg= new DirectedRegion(
			end3,end3+ FLANK_REGION,getStrand());	// 50 nt downstream
		flankReg.setChromosome(getChromosome());
		flankReg.setSpecies(getSpecies());
		
		return isArichRegion(flankReg);
	}
	
	public boolean isInternalExon(Exon e) {
		int p= Arrays.binarySearch(exons, e, getExonOrderComparator());
		if (p> 0&& p< exons.length- 1)
			return true;
		return false;
	}
	
	private static Comparator getExonOrderComparator() {
		// dont use DirectedRegion.PositionComparator, iterate - strand exons from 5' to 3'
		if (exonOrderCompi == null) {
			exonOrderCompi = new AbstractRegion.PositionComparator();
		}

		return exonOrderCompi;
	}
	
	public static boolean isArichRegion(DirectedRegion reg) {
		final int WIN_SIZE= 10;
		String seq= Graph.readSequence(reg);
		seq= seq.toUpperCase();
		int cntA= 0;
		boolean aRich= false;
		for (int i = 0; i < reg.getLength()- WIN_SIZE+ 1; i++) {
			if (i== 0) {
				for (int j = i; j < i+ WIN_SIZE; j++) 
					if (seq.charAt(j)== 'A')
						++cntA;
			} else {
				if (seq.charAt(i-1)== 'A')
					--cntA;
				if (seq.charAt(i+ WIN_SIZE- 1)== 'A')
					++cntA;
			}
			if (cntA> 8) {
				break;
			}
		}
		
		return aRich;
	}

	/**
	 * see  Lopez1 et al., RNA 2006
	 * check complete3PSandro()
	 * @return
	 */
	public boolean is3Pcomplete() {
		
		final int windowSize= 30;	// window to look for pA
		
			// get seq for winSize
		int nb= getExons().length- 1;		
		String seq= "";
		while (seq.length()< windowSize&& nb> 0) {
			seq= Graph.readSequence(getExons()[nb])+ seq;
			--nb;
		}
		
		seq= seq.toUpperCase();
		int x;
		for (x = 0; x < POLY_A_SITES.length; x++) 
			if (seq.contains(POLY_A_SITES[x]))
				break;
		if (x< POLY_A_SITES.length)
			return true;
		return false;
	}
	
	public int[][] countOccurrences(DirectedRegion[] regs, String motif) {
		if (regs== null)
			return new int[0][];		
		motif= motif.toUpperCase();
		int cnt= 0, pos= 0;
		Vector v= new Vector();
		for (int i = 0; i < regs.length; i++) {
			String seq= Graph.readSequence(getGene().getSpecies(), getGene().getChromosome(), isForward(), regs[i].getStart(), regs[i].getEnd());
			seq= seq.toUpperCase();
			while (true) {
				pos= seq.indexOf(motif, pos);
				if (pos< 0)
					break;
				v.add(new int[] {pos,seq.length()-pos});
				// else
				++pos;
			}
		}
		int[][] vv= new int[v.size()][];
		for (int i = 0; i < vv.length; i++) 
			vv[i]= (int[]) v.elementAt(i);
		return vv;
	}


	public boolean hasNonGTAGIntron() {
		DirectedRegion[] regs= getIntrons();
		if (regs== null)
			return false;
		for (int i = 0; i < regs.length; i++) {
			String seq= Graph.readSequence(regs[i]);
			if (seq.length()< 4)// GTAG
				return true;
			String don= seq.substring(0, 2);
			String acc= seq.substring(seq.length()- 2, seq.length());
			if ((!don.equalsIgnoreCase("GT"))|| (!acc.equalsIgnoreCase("AG")))
				return true;
		}
		return false;
	}
	
	public boolean getNonGTAGIntron(HashMap chrMap) {
		DirectedRegion[] regs= getIntrons();
		boolean hasNonGTAG= false;
		if (regs== null)
			return hasNonGTAG;
		
		for (int i = 0; i < regs.length; i++) {
			int[] ratio= (int[]) chrMap.get(getChromosome());
			++ratio[1];
			String seq= Graph.readSequence(regs[i]);
			if (seq== null)
				return true;	// throw away
			if (seq.length()< 4) {// GTAG
				++ratio[0];
				hasNonGTAG|= true;
				continue;
			}
			String don= seq.substring(0, 2);
			String acc= seq.substring(seq.length()- 2, seq.length());
			if ((!don.equalsIgnoreCase("GT"))|| (!acc.equalsIgnoreCase("AG"))) { 
				++ratio[0];
				hasNonGTAG|= true;	
			}
		}
		return hasNonGTAG;
	}
	
	public boolean is5UTR(int pos) {
		
		if (translations== null)
			return false;
		
		int i;
		for (i = 0; translations!= null&& i < translations.length; i++) {
			if (translations[i].get5PrimeEdge()> pos)
				break;
		}
		if (translations!= null&& i< translations.length)
			return true;
		return false;
	}
	
	public boolean is3UTR(int pos) {
		
		if (translations== null)
			return false;
		
		int i;
		for (i = 0; translations!= null&& i < translations.length; i++) {
			if (translations[i].get3PrimeEdge()< pos)
				break;
		}
		if (translations!= null&& i< translations.length)
			return true;
		return false;
	}

	/**
	 * when 
	 * @param pos
	 * @return
	 */
	public boolean isCDS(int pos) {
		
		if (translations== null)
			return false;
		
		int i;
		for (i = 0; translations!= null&& i < translations.length; i++) 
			if (pos>= translations[i].get5PrimeEdge()&& pos<= translations[i].get3PrimeEdge())
				break;
		
		if (translations!= null&& i< translations.length)
			return true;
		return false;
	}
	
	public boolean isExonic(int pos) {
		for (int i = 0; exons!= null&& i < exons.length; i++) 
			if (pos>= exons[i].get5PrimeEdge()&& pos<= exons[i].get3PrimeEdge())
				return true;
		return false;
	}
	
	public boolean isCoding() {
		return (translations!= null&& translations.length> 0);
	}
	
	public boolean isForward() {
		if (gene!= null)	// strand== 0&& 
			return gene.isForward();
		return super.isForward();
	}
	
	Gene gene= null;

    public void setTranscriptID(String transcriptID) {
        this.transcriptID = transcriptID;
    }

    String transcriptID= null;
	Exon[] exons= new Exon[0];	// sorted !!
	public Transcript(Gene newGene, String stableTranscriptID) {

		this.strand= newGene.getStrand();
		this.gene= newGene;
		this.transcriptID= stableTranscriptID;
		setID("transcript");
	}
	
	public Transcript(String newID) {
		this.transcriptID= newID;
		setID("transcript");
	}
	
	/** Get array with exons of transcript
	 * @return Array of exons
	 */
	public Exon[] getExons() {
		return exons;
	}
	
	
	public boolean isATGStart() {
		if (!isCoding())
			return false;
		int tis= getExonicPosition(getTranslations()[0].get5PrimeEdge());
		String startCodon= getSplicedSequence().substring(tis, tis+3);
		if (startCodon.equalsIgnoreCase(Translation.START_CODON))
			return true;
		return false;
	}
	

	
	public DirectedRegion[] getIntrons() {
		if (exons== null|| exons.length< 2)
			return null;
		if (introns== null) {
			Vector<DirectedRegion> intronV= new Vector<DirectedRegion>();
			Exon[] ex= getExons();
			for (int j = 1; j < ex.length; j++) {
				DirectedRegion intron= new DirectedRegion(
						ex[j-1].get3PrimeEdge()+1, 
						ex[j].get5PrimeEdge()- 1,
						getStrand());
				intron.setChromosome(getChromosome());
				intronV.add(intron);
			}
			introns= (DirectedRegion[]) ArrayUtils.toField(intronV);
		}
		return introns;
	}
	
	public byte predictDirectionality() {
		SpliceSite[] sites= getSpliceSitesAll();
		int cntDir= 0, cntCtrDir= 0;
		Gene revGene= new Gene(Gene.getUniqueID());
		revGene.setSpecies(getGene().getSpecies());
		revGene.setChromosome(getGene().getChromosome());
		revGene.setStrand((byte)-getGene().getStrand());
		revGene.setStart(-getGene().getEnd());
		revGene.setEnd(-getGene().getStart());
		for (int i = 0; i < sites.length-1; i++) {
			if (sites[i].isDonor()&& sites[i+1].isAcceptor()) {
				
				if (sites[i+1].getPos()- sites[i].getPos()+ 1< Exon.MIN_INTRON_LENGTH_HUMAN)
					continue;
				
				if (Exon.checkAcceptableIntron(sites[i], sites[i+1]))
					++cntDir;
				
				SpliceSite revDon= new SpliceSite(
						-sites[i+1].getPos(), 
						SpliceSite.TYPE_DONOR,
						revGene);
				SpliceSite revAcc= new SpliceSite(
						-sites[i].getPos(), 
						SpliceSite.TYPE_ACCEPTOR,
						revGene);
				if (Exon.checkAcceptableIntron(revDon, revAcc))
					++cntCtrDir;
			}
		}
		if (cntDir> cntCtrDir)
			return (byte) getStrand();
		if (cntDir== cntCtrDir)
			return 0;
		return ((byte) -getStrand());
		
	}
	
	/**
	 * @return
	 */
	public Gene getGene() {
		return gene;
	}

	/** Get the identifier of transcript in annotation
	 * @return Transcript ID in annotation
	 */
	public String getTranscriptID() {
		return transcriptID;
	}
	
	public String getPureTranscriptID() {
		String id= getTranscriptID();
		id= stripDup(id);
		id= stripGBversion(id);
		return id;
	}
	
	/**
	 * @return
	 */
	public String getStableID() {
		
		return transcriptID;
	}	

	/**
	 * @param exons
	 */
	public void setExons(Exon[] exons) {
		this.exons= exons;
	}

	/**
	 * @param gene
	 */
	public void setGene(Gene gene) {
		this.gene= gene;
		setStrand(getGene().getStrand());
	}

	public int getConfidenceLevel() {
		for (int i = 0; i < CONFIDENCE_LEVELS.length; i++) {
			for (int j = 0; j < CONFIDENCE_LEVELS[i].length; j++) 
				if (getSource().contains(CONFIDENCE_LEVELS[i][j]))
					return i;
		}
		return CONFIDENCE_LEVELS.length+ 1;
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
		if (getGene()== null)
			return species;
		return getGene().getSpecies();
	}
	
	/**
	 * @param b
	 */
	public boolean checkStrand(boolean b) {
		
		return (b== getGene().isForward());
	}
	
	public boolean checkStrand(String newStrand) {
		
		String nStrand= newStrand.trim();
		if (nStrand.equals("1")) 	// || nStrand.equals("forward")
			return checkStrand(true); 
		else if (nStrand.equals("-1"))	// || nStrand.equals("reverse")
			return checkStrand(false);
		
		return false; // error			
	}


	public void addCDS(int start, int end) {
		if (!isForward()) {
			start= -start;
			end= -end;
		}
			
		if (translations== null) {
			translations= new Translation[] {new Translation(this)};
			translations[0].setStart(start);
			translations[0].setEnd(end);
			translations[0].setChromosome(getChromosome());
			translations[0].setSpecies(getSpecies());
			updateBoundaries(translations[0]);
			return;
		}
		
			// else
		if (Math.abs(start)< Math.abs(translations[0].getStart()))
			translations[0].setStart(start);
		if (Math.abs(end)> Math.abs(translations[0].getEnd()))
			translations[0].setEnd(end);
		
		updateBoundaries(translations[0]);
	}

    private static boolean outputExonMerged = false;

	public boolean removeExon(Exon e, boolean removeAcc, boolean removeDon) {
		
		int p= Arrays.binarySearch(exons, e, AbstractRegion.getDefaultPositionComparator());	// search for identical exon (HERE necessary)
		if (p< 0) 
			return false;	// not contained
		
			// remove SSs from transcription hash
		if (removeAcc) {
			Vector<Transcript> vt= getGene().ssTrptHash.remove(e.getAcceptor());
			if (vt!= null&& vt.size()> 1) {
				vt.remove(this);
				getGene().ssTrptHash.put(e.getAcceptor(), vt);
			} else { // no longer supported, remove SSs from position hash
				Integer pos= new Integer(e.get5PrimeEdge());
				SpliceSite s= e.getAcceptor();	// save splice site
				Vector<SpliceSite> vss= getGene().spliceHash.remove(pos);
				if (vss!= null) {				
					vss.remove(s);
					if (vss.size()> 0)
						getGene().spliceHash.put(pos, vss);
				}
			}
		}
		
		if (removeDon) {
			Vector<Transcript> vt= getGene().ssTrptHash.remove(e.getDonor());
			if (vt!= null&& vt.size()> 1) {
				vt.remove(this);
				getGene().ssTrptHash.put(e.getDonor(), vt);
			} else {
				Integer pos= new Integer(e.get3PrimeEdge());
				SpliceSite s= e.getDonor();
				Vector<SpliceSite> vss= getGene().spliceHash.remove(pos);
				if (vss!= null) {
					vss.remove(s);
					if (vss.size()> 0)
						getGene().spliceHash.put(pos, vss);
				}
	
			}
		}
		

			// remove from gene exon hash
		Transcript[] t= getGene().exonHash.remove(e);
		if (t.length> 1) {
			Transcript[] tt= new Transcript[t.length-1];
			for (int i = 0, j= 0; i < tt.length; ++i, ++j) {
				if (t[i]== this)
					--j;
				else
					tt[j]= t[i];
			}
			getGene().exonHash.put(e, tt);
		}
		
			// kick out of this position array
		Exon[] newE= new Exon[exons.length-1];
		for (int i = 0; i < p; i++) 
			newE[i]= exons[i];
		for (int i = p+1; i < exons.length; i++) 
			newE[i-1]= exons[i];
		exons= newE;
		
		return true;
	}
	
			/**
			 * Inserts the exons in an array sorted according to ascending order
			 * of their start/stop position. <b>IMPORTANT</b>: add exons AFTER adding 
			 * transcripts to ensure the correct init of AS types.
			 * 
			 * @param newExon
			 * @return the exon already contained or <code>newExon</code> case of the exon was added successfully
			 */
			public boolean addExon(Exon newExon) {
		
					// generate splice sites, but do not yet insert into transcript
				int p= Arrays.binarySearch(exons, newExon, AbstractRegion.getDefaultPositionComparator());	// search for identical exon (HERE necessary)
				if (p>= 0) 
					return false;	// already contained, not added - also no need to check gene
				p= -(p+1);	// insertion point
			
					// too short introns
				if (removeGaps) {
					if ((p-1>= 0)&& (exons[p-1].get3PrimeEdge()+ maxLengthIntronIsGap + 1>= newExon.get5PrimeEdge()))  {
						String s= "Merging exon ("+newExon.start+","+newExon.end+") with exon ("+ exons[p-1].start+","+exons[p-1].end+")"+
								" in transcript "+ getTranscriptID()+ " because intervening intron has "+ maxLengthIntronIsGap +" or less nt.";
                        if (outputExonMerged)
                            Log.debug(s);
                        else {
                            Log.warn(s);
                            Log.warn("Further merged exons are sent at debug verbosity level.");
                            outputExonMerged = true;
                        }
						newExon.set5PrimeEdge(exons[p-1].get5PrimeEdge());
						removeExon(exons[p-1], false, true);
						--p;
						//return false;
					}
					if ((p< exons.length)&& (exons[p].get5PrimeEdge()<= newExon.get3PrimeEdge()+ maxLengthIntronIsGap + 1))  {
						String s= "Merging exon (" + newExon.start + "," + newExon.end + ") with exon (" + exons[p].start + "," + exons[p].end + ")" +
                                " in transcript " + getTranscriptID() + " because intervening intron has " + maxLengthIntronIsGap + " or less nt.";
                        if (outputExonMerged)
                            Log.debug(s);
                        else {
                            Log.warn(s);
                            Log.warn(" Further merged exons are sent at debug verbosity level.");
                            outputExonMerged = true;
                        }
						newExon.set3PrimeEdge(exons[p].get3PrimeEdge());
						removeExon(exons[p], true, false);
						; // nothing, insertion point is the same
						//return false;
					}
				}

					// insert
				Vector<Transcript> trptV= new Vector<Transcript>(1,1);
				trptV.add(this);
				if (p== 0) {	// new first exon: ex-first exon now has an acceptor, newExon has tss and donor
					byte type= SpliceSite.TYPE_SOFT_START;
					if (getConfidenceLevel()<= edgeConfidenceLevel)
						type= SpliceSite.TYPE_HARD_START;
					getGene().addSpliceSite(new SpliceSite(newExon.get5PrimeEdge(), type, getGene(), getSourceType()), trptV);
					if (exons.length> 0)
						getGene().addSpliceSite(new SpliceSite(exons[0].get5PrimeEdge(), SpliceSite.TYPE_ACCEPTOR, getGene(), getSourceType()), trptV);
				} 
				 
				if (p> 0) 	// has acceptor
					getGene().addSpliceSite(new SpliceSite(newExon.get5PrimeEdge(), SpliceSite.TYPE_ACCEPTOR, getGene(), getSourceType()), trptV);
				
				if (p< this.exons.length) 			// has donor
					getGene().addSpliceSite(new SpliceSite(newExon.get3PrimeEdge(), SpliceSite.TYPE_DONOR, getGene(), getSourceType()), trptV);
		
				if (p== exons.length) {	// ex-last exon now has an donor
					byte type= SpliceSite.TYPE_SOFT_END;
					if (getConfidenceLevel()<= edgeConfidenceLevel)
						type= SpliceSite.TYPE_HARD_END;
					getGene().addSpliceSite(new SpliceSite(newExon.get3PrimeEdge(), type, getGene(), getSourceType()), trptV);
					if (exons.length> 0)
						getGene().addSpliceSite(new SpliceSite(exons[p-1].get3PrimeEdge(), SpliceSite.TYPE_DONOR, getGene(), getSourceType()), trptV);
				}

				// NOW insert
				getGene().addExon(newExon, new Transcript[] {this});
				if (exons== null) 
					exons= new Exon[] {newExon};
				else
					exons= (Exon[]) ArrayUtils.insert(this.exons, newExon, p);
				
				updateBoundaries(newExon);
				return true;
			}
			
			public void setStart(int v) {
			super.setStart(v);
		}
			
		/**
		 * warning, destroys gene cluster
		 */
			@Override
		public void setStrand(byte strand) {
			byte oldStrand= getStrand();
			super.setStrand(strand);
			if (oldStrand!= strand&& oldStrand!= STRAND_NA) {
				
				// too much mem, cheat
//				Gene newGene= (Gene) getGene().clone();
//				newGene.setStrand(strand);
//				Exon[] newExons= new Exon[exons.length];
//				for (int i = 0; i < exons.length; i++) {
//					Exon e= (Exon) exons[i].clone();
//					e.setStrand(strand);
//					e.setGene(newGene);
//					newExons[exons.length-1-i]= e;
//				}
//				exons= newExons;
				
//				SpliceSite[] newSites= new SpliceSite[spliceSites.length];
//				for (int i = 0; i < spliceSites.length; i++) {
//					SpliceSite ss= new SpliceSite(-spliceSites[i].getPos(), SpliceSite.inverseType(spliceSites[i].getType()), newGene);
//					newSites[spliceSites.length-1-i]= ss;
//				}
//				spliceSites= newSites;
				
				if (getTranslations()!= null&& getTranslations().length> 0) {
					getTranslations()[0].setStrand(strand);
				}
			}
		}
		
		public void setEnd(int v) {
			super.setEnd(v);
		}

	/**
	 * Finds the first exon containing the corresponding position
	 * 
	 * deprecated inconsistent for overlapping exons // why?!		 
	 * @param absPos
	 * @return
	 */
	public int getExonIdx(int absPos) {
		
		for (int i = 0; i < exons.length; i++) 
			if (exons[i].contains(absPos))
				return i;
		
		return -1;
	}
	public Exon getExon(int absPos) {
		int p= getExonIdx(absPos);
		if (p>= 0)
			return exons[p];
		return null;
	}
	
	/**
	 * gets exons in between e1 and e2
	 * @param e1
	 * @param e2
	 * @return
	 */
	public Exon[] getExons(Exon e1, Exon e2) {
		Comparator compi= new AbstractRegion.PositionComparator();
		int p1= Arrays.binarySearch(exons, e1, compi);
		int p2= Arrays.binarySearch(exons, e2, compi);
		
		if (p1== p2)
			return new Exon[0];
		Exon[] res= new Exon[p2-p1-1];
		for (int i = p1+1; i < p2; i++) 
			res[i-p1-1]= exons[i];
		return res;
	}
	
	public Exon getExon(DirectedRegion reg) {
		
		for (int i = 0; i < exons.length; i++) 
			if (exons[i].overlaps(reg))
				return exons[i];
		
		return null;
	}

	public int getOtherSideOfExon(int pos) {
		for (int i = 0; i < exons.length; i++) {
			if (exons[i].getStart()== pos)
				return exons[i].getEnd();
			else if (exons[i].getEnd()== pos)
				return exons[i].getStart();
		}
		
		return -1;
	}
	
	public Translation[][] findORFs(boolean openEnded) {
		String seq= getSplicedSequence();
		Translation[][] orfs= new Translation[3][];
		for (int i = 0; i < 3; i++) {
			orfs[i]= findORFs(seq, i, openEnded);
		}
		// skipped for the moment
//		int cntORFs= 0;
//		for (int i = 0; i < orfs.length; i++) 
//			if (orfs[i]!= null)
//				cntORFs+= orfs[i].length;	// TODO min threshold
//		if (cntORFs== 0) {
//			for (int i = 0; i < orfs.length; i++) {
//				orfs[i]= forceORFs(seq, i);
//			}
//		}
		return orfs;
	}
	
	/**
	 * Use the most upstream ATG as start codon. If there is no stop codon 
	 * upstream of the first ATG, check genomic sequence upstream beyond the 
	 * gene for in-frame stop. If one is found and the start of the gene is
	 * in a CpG island, use the ATG, otherwise start the open-ended CDS at the
	 * start of the gene.
	 * @return
	 */
	public Translation findHavanaORF() {
		
		Translation[][] tlns= findORFs(false);
		if (tlns== null)
			return null;
		
			// find annotated start atgs in the gene
		Transcript[] trpts= getGene().getTranscripts();
		Vector atgPosV= new Vector();
		for (int i = 0; i < trpts.length; i++) {
			if (trpts[i]== this)
				continue;
			if (trpts[i].isCoding()) {
				Translation tmpTln= trpts[i].getTranslations()[0];
				// now, predict open ended, if already reported
				if (tmpTln.isOpenEnded())	
					continue;
				Integer pos= new Integer(tmpTln.get5PrimeEdge());
				int j; 
				for (j = 0; j < atgPosV.size(); j++) 
					if (atgPosV.elementAt(j).equals(pos))
						break;
				if (j== atgPosV.size())
					atgPosV.add(pos);
			}
		}
		
			// find ORFs at the corresponding starts
		Vector tlnV= new Vector();
		for (int i = 0; i < tlns.length; i++) {
			for (int j = 0; tlns[i]!= null&& j < tlns[i].length; j++) {
				int k;
				for (k = 0; k < atgPosV.size(); k++) 
					if (tlns[i][j].get5PrimeEdge()== ((Integer) atgPosV.elementAt(k)).intValue())
						break;
				if (k< atgPosV.size())
					tlnV.add(tlns[i][j]);
			}
		}

			// if coincidence with annotated ORF, get longest of this
		if (tlnV.size()< 1)
			return findLongestORF(tlns);
		else
			return findLongestORF(new Translation[][] {(Translation[]) ArrayUtils.toField(tlnV)});
	}
	
	public Translation findLongestORF() {
		return findLongestORF(findORFs(false));
	}
	
	Translation findLongestORF(Translation[][] predORFs) {
		Translation longestORF= null;
		if (predORFs== null)
			return null;
		for (int i = 0; i < predORFs.length; i++) 
			for (int j = 0; predORFs[i]!= null&& j < predORFs[i].length; j++) 
				if (longestORF== null|| predORFs[i][j].getSplicedLength()> longestORF.getSplicedLength())
					longestORF= predORFs[i][j];
		
		if (longestORF!= null&& longestORF.getSplicedLength()< (NMDSimulator.MIN_ORF_LENGTH_AA* 3))
			return null;
		return longestORF;
	}
	
	public Translation[] getAllORFs() {
		Translation[][] predORFs= findORFs(false);
		if (predORFs== null)
			return null;
		Vector allOrfV= new Vector();
		for (int i = 0; i < predORFs.length; i++) 
			for (int j = 0; predORFs[i]!= null&& j < predORFs[i].length; j++) 
				allOrfV.add(predORFs[i][j]);
		
		Translation[] allORFs= (Translation[]) ArrayUtils.toField(allOrfV);
		return allORFs;
	}
	
	public Translation[] findORFs(String seq, int frame, boolean openEnded) {
//			if (seq== null|| seq.length()< 1)
//				return null;
//			final int minORFlen= NMDSimulator.MIN_ORF_LENGTH_AA* 3;
//			seq= seq.substring(frame).toUpperCase();
//			Vector startV= new Vector(), stopV= new Vector();
//			int pos= 0;
//			for (int i = 0; (i+3) < seq.length(); i+= 3) {
//				String codon= seq.substring(i, i+3);
//				if (codon.equals(Translation.START_CODON))
//					startV.add(new Integer(i));
//				else if (codon.equals(Translation.STOP_CODONS[0])||
//						codon.equals(Translation.STOP_CODONS[1])||
//						codon.equals(Translation.STOP_CODONS[2]))
//					stopV.add(new Integer(i));
//			}
//			int[] startPos= gphase.tools.Arrays.toPrimitive((Integer[]) gphase.tools.Arrays.toField(startV));
//			int[] stopPos= gphase.tools.Arrays.toPrimitive((Integer[]) gphase.tools.Arrays.toField(stopV));
//	
//				// determine ORFs
//			Vector v= new Vector();
//			Translation reg;
//			for (int i = 0; startPos!= null&& stopPos!= null&& i < startPos.length; i++) {
//				int j= Arrays.binarySearch(stopPos, startPos[i]);
//				if (j>= 0)
//					System.err.println("Assertion failed: start/stop at same pos!");
//				j= -j- 1;
//				if (j== stopPos.length) {
//					reg= new Translation(this, getGenomicPosition(frame+ startPos[i]), get3PrimeEdge(), getStrand());
//					reg.setSplicedLength(seq.length()- startPos[i]+ 1);
//				} else {
//					reg= new Translation(this, getGenomicPosition(frame+ startPos[i]), getGenomicPosition(frame+ stopPos[j]+ 2), getStrand());
//					reg.setSplicedLength(stopPos[j]+ 2- startPos[i]+ 1);
//				}
//				reg.setChromosome(getChromosome());
//				reg.setSpecies(getSpecies());
//				//if (reg.getLength()> minORFlen)
//					v.add(reg);
//			}
//			
//			if (!openEnded)
//				return (Translation[]) gphase.tools.Arrays.toField(v);
//
//				// 3'OpenEnded
//			if (stopPos!= null&& stopPos.length> 0) {	// ORF starting outside transcript
//				reg= new Translation(this, get5PrimeEdge(), getGenomicPosition(frame+ stopPos[0]+ 2), getStrand());	// +2 to include stop_codon
//				reg.setSplicedLength(stopPos[0]+ 2+ frame+ 1);	// here, q si frame ... goes outside transcript..
//				reg.setChromosome(getChromosome());
//				reg.setSpecies(getSpecies());
//				v.add(reg);
//			} else {	// bi-OpenEnded
//				reg= new Translation(this, get5PrimeEdge(), get3PrimeEdge(), getStrand());	// +2 to include stop_codon
//				reg.setSplicedLength(seq.length());	// here, q si frame ... goes outside transcript..
//				reg.setChromosome(getChromosome());
//				reg.setSpecies(getSpecies());
//				v.add(reg);
//			}
//				// 5'OpenEnded
//			if (startPos!= null&& startPos.length> 0){	
//				reg= new Translation(this, getGenomicPosition(frame+ startPos[0]), get3PrimeEdge(), getStrand());
//				reg.setSplicedLength(seq.length()- startPos[0]+ 1);
//				reg.setChromosome(getChromosome());
//				reg.setSpecies(getSpecies());
//				v.add(reg);
//			} 
//			
//			return (Translation[]) gphase.tools.Arrays.toField(v);
		return null;
		}
	
	/**
	 * 
	 * @param seq
	 * @param frame
	 * @deprecated not in use, check inframe stop
	 * @return
	 */
	public Translation[] forceORFs(String seq, int frame) {
//		seq= seq.substring(frame).toUpperCase();
//		Vector startV= new Vector(), stopV= new Vector();
//		for (int i = 0; (i+3) < seq.length(); i+= 3) {
//			String codon= seq.substring(i, i+3);
//			if (codon.equals(Translation.START_CODON))
//				startV.add(new Integer(i));
//			else if (codon.equals(Translation.STOP_CODONS[0])||
//					codon.equals(Translation.STOP_CODONS[1])||
//					codon.equals(Translation.STOP_CODONS[2]))
//				stopV.add(new Integer(i));
//		}
//		int[] startPos= gphase.tools.Arrays.toPrimitive((Integer[]) gphase.tools.Arrays.toField(startV));
//		int[] stopPos= gphase.tools.Arrays.toPrimitive((Integer[]) gphase.tools.Arrays.toField(stopV));
//
//			// determine ORFs
//		Vector v= new Vector();
//		Translation reg;
//		if (stopPos!= null&& stopPos.length> 0) {	// ORF starting outside transcript
//			reg= new Translation(this, get5PrimeEdge(), getGenomicPosition(frame+ stopPos[0]+ 2), getStrand());	// +2 to include stop_codon
//			reg.setSplicedLength(stopPos[0]+ 2+ frame+ 1);	// here, q si frame ... goes outside transcript..
//			reg.setChromosome(getChromosome());
//			reg.setSpecies(getSpecies());
//			
//			v.add(reg);
//		} 
//		
//		if (startPos!= null&& startPos.length> 0){	
//			reg= new Translation(this, getGenomicPosition(frame+ startPos[0]), get3PrimeEdge(), getStrand());
//				// TODO check for no inframe stop
//			reg.setSplicedLength(seq.length()- startPos[0]+ 1);
//			reg.setChromosome(getChromosome());
//			reg.setSpecies(getSpecies());
//			v.add(reg);
//		}
//		
//		return (Translation[]) gphase.tools.Arrays.toField(v);
		return null;
	}
	

	
	public Exon getExon(String stableID) {
		
		if (exons== null)
			return null;
		
		for (int i = 0; i < exons.length; i++) 
			if (exons[i].getExonID().equals(stableID))
				return exons[i];
	
		return null;
	}

	public Exon getLastExon() {
		
		if (exons== null|| exons.length< 1)
			return null;
		return exons[exons.length- 1];
	}

	public int getDistFromATG(int gPos) {

		if (getTranslations()== null|| getTranslations().length!= 1)
			return -1;
		
		int atg= getTranslations()[0].get5PrimeEdge();
		int min= Math.min(gPos, atg);
		int max= Math.max(gPos, atg);
		
		int i;
		for (i = 0; i < exons.length; i++) 
			if(exons[i].get3PrimeEdge()>= min)
				break;
		if (i== exons.length)	// outa range
			return -1;

		int j;
		for (j = i; j < exons.length; j++) 
			if(exons[j].get3PrimeEdge()>= max)
				break;
		if (j== exons.length)	// outa range
			return -1;
		
		int dist= 0;
		for (int k = i+1; k < j; k++) 
			dist+= exons[k].getLength();
		
		if (i== j)
			dist+= max- min;
		else {
			dist+= exons[i].get3PrimeEdge()- min;	// 1st exon, not +-1
			dist+= max- exons[j].get5PrimeEdge();	// last exon, not +-1
			++dist;	// and in both cases +1
		}

		return dist;
	}
	/**
	 * working with cds in exons
	 * @return
	 */
	public DirectedRegion[] getCDSRegions() {
		Exon[] exons= getExons();
		Vector v= new Vector(exons.length);
		DirectedRegion reg;
		for (int i = 0; i < exons.length; i++) {
			if (!exons[i].isCoding())
				continue;
			reg= new DirectedRegion(exons[i].get5PrimeCDS(), exons[i].get3PrimeCDS(), exons[i].getStrand());
			reg.setChromosome(getChromosome());
			reg.setSpecies(getSpecies());
			v.add(reg);
		}
		
		DirectedRegion[] regs= (DirectedRegion[]) ArrayUtils.toField(v);
		// no longer necessary, stop in VEGA st included, st not
		//regs[regs.length- 1].set3PrimeEdge(regs[regs.length- 1].get3PrimeEdge()+ 3);	// to include stop codon
		return regs;
	}

    /**
     * Retrieve the nucleotide sequence of CDS
     * @return CDS sequence
     */
	
	public String getCDSSequenceNt() {
		DirectedRegion[] regs= getCDSRegions();	// not sorted
		java.util.Arrays.sort(regs, new DirectedRegion.DirectedPositionComparator());
		StringBuffer sb= new StringBuffer();
		for (int i = 0; i < regs.length; i++) 
			sb.append(Graph.readSequence(regs[i]));
		return sb.toString();
	}
	public byte getNmd() {
		return nmd;
	}
	public void setNmd(byte nmd) {
		this.nmd = nmd;
	}
	public String getSource() {
		return source;
	}
	
	public byte getSourceType() {
		if (srcType== ID_SRC_UNDEFINED) 
			srcType= Transcript.getSourceType(getSource());
		//else
		return srcType;
	}
	
	public static byte getSourceType(String source) {
		if (source.contains("Est"))
			return ID_SRC_EST;
		else if (source.contains("mrna"))
			return ID_SRC_MRNA;
		else if (source.contains("knownGene"))
			return ID_SRC_UCSC;
		else if (source.contains("refGene"))
			return ID_SRC_REFSEQ;
		return ID_SRC_UNDEFINED;
	}
	public void setSource(String source) {
		this.source = source;
	}

	public static byte getEdgeConfidenceLevel() {
		return edgeConfidenceLevel;
	}

	static final Pattern versionSfx= Pattern.compile("(.+)\\.\\d$");
	public static String stripGBversion(String id) {
		return strip (id, versionSfx);
	}
	
	private static String strip(String id, Pattern pattern) {
		Matcher m= pattern.matcher(id);
		if (m.matches()) 
			return m.group(1);
		return id;
	}

	static final Pattern dupSfx= Pattern.compile("(.+)_dup\\d$");
	public static String stripDup(String id) {
		return strip(id, dupSfx);
	}
	
	public static void setEdgeConfidenceLevel(byte edgeConfidenceLevel) {
		Transcript.edgeConfidenceLevel = edgeConfidenceLevel;
	}
	
	public static String getSource(Transcript[] tt) {
		String src= null;
		byte source= Transcript.ID_SRC_UNDEFINED;
		for (int i = 0; i < tt.length; i++) {
			if (src== null) 
				src= tt[i].getSource();
			else if (!src.equals(tt[i].getSource())) {
				src= "Mixed";// todo: was GFFObject.SOURCE_MIXED make it static again
				break;
			}
			source= (byte) Math.max(source, tt[i].getSourceType());
		}
		
		if (source== Transcript.ID_SRC_UNDEFINED)
			return src;
		else
			return Transcript.TYPE_TO_ID[source];
	}

    public boolean containsSS(SpliceSite ss) {
        for (int i = 0; i < spliceSites.length; i++) {
            if (ss.getPos()== spliceSites[i].getPos()&&
                    ss.isDonor()== spliceSites[i].isDonor())
                return true;
        }
        return false;
    }

}