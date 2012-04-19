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
 * Created on May 4, 2005
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
package barna.model;

import barna.commons.utils.ArrayUtils;
import barna.model.commons.IntVector;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Vector;

//import gphase.StopCodons;

/**
 * 
 * 
 * @author micha
 */
public class Translation extends DirectedRegion {

	static final long serialVersionUID= 8996021902187779155L;
	public static final int UNKNOWN_ID= 0;
	public static final int REFSEQ_ID= 1;
	public static final int ENSEMBL_ID= 2;
	public static final int SWISSPROT_ID= 3;
	public static String ORDER_AA= "ILVFMCAGPTSYWQNHEDKR";
	public static String AA_STERIC= "FWHYP";
	public static String AA_POLAR= "DECNQTYSGHKR";	// hydrophob are the others
	public static String AA_POS_CHARGE= "HKR";
	public static String AA_NEG_CHARGE= "DE";
		
	public static String[][] CODONS_AA= new String[][] {
		new String[] {"ATT", "ATC", "ATA"},
		new String[] {"CTT", "CTC", "CTA", "CTG", "TTA", "TTG"},
		new String[] {"GTT", "GTC", "GTA", "GTG"},
		new String[] {"TTT", "TTC"},
		new String[] {"ATG"},
		new String[] {"TGT", "TGC"},
		new String[] {"GCT", "GCC", "GCA", "GCG"},
		new String[] {"GGT", "GGC", "GGA", "GGG"},
		new String[] {"CCT", "CCC", "CCA", "CCG"},
		new String[] {"ACT", "ACC", "ACA", "ACG"},
		new String[] {"TCT", "TCC", "TCA", "TCG", "AGT", "AGC"},
		new String[] {"TAT", "TAC"},
		new String[] {"TGG"},
		new String[] {"CAA", "CAG"},
		new String[] {"AAT", "AAC"},
		new String[] {"CAT", "CAC"},
		new String[] {"GAA", "GAG"},
		new String[] {"GAT", "GAC"},
		new String[] {"AAA", "AAG"},
		new String[] {"CGT", "CGC", "CGA", "CGG", "AGA", "AGG"}
	};
	
	public static HashMap<String, String> CODON_HASH= new HashMap<String, String>();
	{
		CODON_HASH.put("TTT", "F");
		CODON_HASH.put("TTC", "F");
		CODON_HASH.put("TTA", "L");
		CODON_HASH.put("TTG", "L");
		CODON_HASH.put("CTT", "L");
		CODON_HASH.put("CTC", "L");
		CODON_HASH.put("CTA", "L");
		CODON_HASH.put("CTG", "L");
		CODON_HASH.put("ATT", "I");
		CODON_HASH.put("ATC", "I");
		CODON_HASH.put("ATA", "I");
		CODON_HASH.put("ATG", "M");
		CODON_HASH.put("GTT", "V");
		CODON_HASH.put("GTC", "V");
		CODON_HASH.put("GTA", "V");
		CODON_HASH.put("GTG", "V");
		
		CODON_HASH.put("TCT", "S");
		CODON_HASH.put("TCC", "S");
		CODON_HASH.put("TCA", "S");
		CODON_HASH.put("TCG", "S");
		CODON_HASH.put("CCT", "P");
		CODON_HASH.put("CCC", "P");
		CODON_HASH.put("CCA", "P");
		CODON_HASH.put("CCG", "P");
		CODON_HASH.put("ACT", "T");
		CODON_HASH.put("ACC", "T");
		CODON_HASH.put("ACA", "T");
		CODON_HASH.put("ACG", "T");
		CODON_HASH.put("GCT", "A");
		CODON_HASH.put("GCC", "A");
		CODON_HASH.put("GCA", "A");
		CODON_HASH.put("GCG", "A");
		
		CODON_HASH.put("TAT", "Y");
		CODON_HASH.put("TAC", "Y");
		CODON_HASH.put("TAA", "");	// stop
		CODON_HASH.put("TAG", "");	// stop
		CODON_HASH.put("CAT", "H");
		CODON_HASH.put("CAC", "H");
		CODON_HASH.put("CAA", "Q");
		CODON_HASH.put("CAG", "Q");
		CODON_HASH.put("AAT", "N");
		CODON_HASH.put("AAC", "N");
		CODON_HASH.put("AAA", "K");
		CODON_HASH.put("AAG", "K");
		CODON_HASH.put("GAT", "D");
		CODON_HASH.put("GAC", "D");
		CODON_HASH.put("GAA", "E");
		CODON_HASH.put("GAG", "E");
		
		CODON_HASH.put("TGT", "C");
		CODON_HASH.put("TGC", "C");
		CODON_HASH.put("TGA", "");	// stop
		CODON_HASH.put("TGG", "W");
		CODON_HASH.put("CGT", "R");
		CODON_HASH.put("CGC", "R");
		CODON_HASH.put("CGA", "R");
		CODON_HASH.put("CGG", "R");
		CODON_HASH.put("AGT", "S");
		CODON_HASH.put("AGC", "S");
		CODON_HASH.put("AGA", "R");
		CODON_HASH.put("AGG", "R");
		CODON_HASH.put("GGT", "G");
		CODON_HASH.put("GGC", "G");
		CODON_HASH.put("GGA", "G");
		CODON_HASH.put("GGG", "G");

	}
	public static HashMap codonHash= null;
	
	public static final int FRAME_BITCDS0= 0, FRAME_BITCDS1= 1, FRAME_BITCDS2= 2, FRAME_BIT3UTR= 3, FRAME_BITSTART= 4, FRAME_BIT5UTR= 5, FRAME_BITSTOP= 6, FRAME_BITNC= 7, FRAME_BITNI= 8;
	// byte[i]=1<<i, i=[0..8]
	public static final int[] FRAME_BYTEVAL= new int[] {1,2,4,8,16,32,64,128,0};
	public static final String[] FRAME_VERBOSE= new String[] {"CDS0", "CDS1", "CDS2", "3UTR", "START", "5UTR", "STOP", "NC", "NI"};
	public static final HashMap<Integer, String> mapFrameVerbose= new HashMap<Integer, String>();
	static {
		mapFrameVerbose.put((int) FRAME_BYTEVAL[FRAME_BIT5UTR], FRAME_VERBOSE[FRAME_BIT5UTR]);
		mapFrameVerbose.put((int) FRAME_BYTEVAL[FRAME_BITCDS0], FRAME_VERBOSE[FRAME_BITCDS0]);
		mapFrameVerbose.put((int) FRAME_BYTEVAL[FRAME_BITCDS1], FRAME_VERBOSE[FRAME_BITCDS1]);
		mapFrameVerbose.put((int) FRAME_BYTEVAL[FRAME_BITCDS2], FRAME_VERBOSE[FRAME_BITCDS2]);
		mapFrameVerbose.put((int) FRAME_BYTEVAL[FRAME_BIT3UTR], FRAME_VERBOSE[FRAME_BIT3UTR]);
		mapFrameVerbose.put((int) FRAME_BYTEVAL[FRAME_BITSTART], FRAME_VERBOSE[FRAME_BITSTART]);
		mapFrameVerbose.put((int) FRAME_BYTEVAL[FRAME_BITSTOP], FRAME_VERBOSE[FRAME_BITSTOP]);
		mapFrameVerbose.put((int) FRAME_BYTEVAL[FRAME_BITNC], FRAME_VERBOSE[FRAME_BITNC]);
		mapFrameVerbose.put((int) FRAME_BYTEVAL[FRAME_BITNI], FRAME_VERBOSE[FRAME_BITNI]);
		
		mapFrameVerbose.put(getCombinedFrame(FRAME_BYTEVAL[FRAME_BIT5UTR], FRAME_BYTEVAL[FRAME_BIT5UTR]), "5UTR");
		mapFrameVerbose.put(getCombinedFrame(FRAME_BYTEVAL[FRAME_BIT5UTR], FRAME_BYTEVAL[FRAME_BITSTART]), "5UTR");
		mapFrameVerbose.put(getCombinedFrame(FRAME_BYTEVAL[FRAME_BIT5UTR], FRAME_BYTEVAL[FRAME_BITCDS0]), "5UTR-CDS");
		mapFrameVerbose.put(getCombinedFrame(FRAME_BYTEVAL[FRAME_BIT5UTR], FRAME_BYTEVAL[FRAME_BITCDS1]), "5UTR-CDS");
		mapFrameVerbose.put(getCombinedFrame(FRAME_BYTEVAL[FRAME_BIT5UTR], FRAME_BYTEVAL[FRAME_BITCDS2]), "5UTR-CDS");
		mapFrameVerbose.put(getCombinedFrame(FRAME_BYTEVAL[FRAME_BIT5UTR], FRAME_BYTEVAL[FRAME_BITSTOP]), "5UTR-CDS");
		mapFrameVerbose.put(getCombinedFrame(FRAME_BYTEVAL[FRAME_BIT5UTR], FRAME_BYTEVAL[FRAME_BIT3UTR]), "5UTR-CDS-3UTR");
		
		mapFrameVerbose.put(getCombinedFrame(FRAME_BYTEVAL[FRAME_BITSTART], FRAME_BYTEVAL[FRAME_BITSTART]), "CDS");
		mapFrameVerbose.put(getCombinedFrame(FRAME_BYTEVAL[FRAME_BITSTART], FRAME_BYTEVAL[FRAME_BITCDS0]), "CDS");
		mapFrameVerbose.put(getCombinedFrame(FRAME_BYTEVAL[FRAME_BITSTART], FRAME_BYTEVAL[FRAME_BITCDS1]), "CDS");
		mapFrameVerbose.put(getCombinedFrame(FRAME_BYTEVAL[FRAME_BITSTART], FRAME_BYTEVAL[FRAME_BITCDS2]), "CDS");
		mapFrameVerbose.put(getCombinedFrame(FRAME_BYTEVAL[FRAME_BITSTART], FRAME_BYTEVAL[FRAME_BITSTOP]), "CDS");
		mapFrameVerbose.put(getCombinedFrame(FRAME_BYTEVAL[FRAME_BITSTART], FRAME_BYTEVAL[FRAME_BIT3UTR]), "CDS-3UTR");
		
		mapFrameVerbose.put(getCombinedFrame(FRAME_BYTEVAL[FRAME_BITCDS0], FRAME_BYTEVAL[FRAME_BITSTART]), "CDS");
		mapFrameVerbose.put(getCombinedFrame(FRAME_BYTEVAL[FRAME_BITCDS0], FRAME_BYTEVAL[FRAME_BITCDS0]), "CDS");
		mapFrameVerbose.put(getCombinedFrame(FRAME_BYTEVAL[FRAME_BITCDS0], FRAME_BYTEVAL[FRAME_BITCDS1]), "CDS");
		mapFrameVerbose.put(getCombinedFrame(FRAME_BYTEVAL[FRAME_BITCDS0], FRAME_BYTEVAL[FRAME_BITCDS2]), "CDS");
		mapFrameVerbose.put(getCombinedFrame(FRAME_BYTEVAL[FRAME_BITCDS0], FRAME_BYTEVAL[FRAME_BITSTOP]), "CDS");
		mapFrameVerbose.put(getCombinedFrame(FRAME_BYTEVAL[FRAME_BITCDS0], FRAME_BYTEVAL[FRAME_BIT3UTR]), "CDS-3UTR");
		
		mapFrameVerbose.put(getCombinedFrame(FRAME_BYTEVAL[FRAME_BITCDS1], FRAME_BYTEVAL[FRAME_BITSTART]), "CDS");
		mapFrameVerbose.put(getCombinedFrame(FRAME_BYTEVAL[FRAME_BITCDS1], FRAME_BYTEVAL[FRAME_BITCDS0]), "CDS");
		mapFrameVerbose.put(getCombinedFrame(FRAME_BYTEVAL[FRAME_BITCDS1], FRAME_BYTEVAL[FRAME_BITCDS1]), "CDS");
		mapFrameVerbose.put(getCombinedFrame(FRAME_BYTEVAL[FRAME_BITCDS1], FRAME_BYTEVAL[FRAME_BITCDS2]), "CDS");
		mapFrameVerbose.put(getCombinedFrame(FRAME_BYTEVAL[FRAME_BITCDS1], FRAME_BYTEVAL[FRAME_BITSTOP]), "CDS");
		mapFrameVerbose.put(getCombinedFrame(FRAME_BYTEVAL[FRAME_BITCDS1], FRAME_BYTEVAL[FRAME_BIT3UTR]), "CDS-3UTR");

		mapFrameVerbose.put(getCombinedFrame(FRAME_BYTEVAL[FRAME_BITCDS2], FRAME_BYTEVAL[FRAME_BITSTART]), "CDS");
		mapFrameVerbose.put(getCombinedFrame(FRAME_BYTEVAL[FRAME_BITCDS2], FRAME_BYTEVAL[FRAME_BITCDS0]), "CDS");
		mapFrameVerbose.put(getCombinedFrame(FRAME_BYTEVAL[FRAME_BITCDS2], FRAME_BYTEVAL[FRAME_BITCDS1]), "CDS");
		mapFrameVerbose.put(getCombinedFrame(FRAME_BYTEVAL[FRAME_BITCDS2], FRAME_BYTEVAL[FRAME_BITCDS2]), "CDS");
		mapFrameVerbose.put(getCombinedFrame(FRAME_BYTEVAL[FRAME_BITCDS2], FRAME_BYTEVAL[FRAME_BITSTOP]), "CDS");
		mapFrameVerbose.put(getCombinedFrame(FRAME_BYTEVAL[FRAME_BITCDS2], FRAME_BYTEVAL[FRAME_BIT3UTR]), "CDS-3UTR");

		mapFrameVerbose.put(getCombinedFrame(FRAME_BYTEVAL[FRAME_BITSTOP], FRAME_BYTEVAL[FRAME_BITSTART]), "CDS");
		mapFrameVerbose.put(getCombinedFrame(FRAME_BYTEVAL[FRAME_BITSTOP], FRAME_BYTEVAL[FRAME_BITCDS0]), "CDS");
		mapFrameVerbose.put(getCombinedFrame(FRAME_BYTEVAL[FRAME_BITSTOP], FRAME_BYTEVAL[FRAME_BITCDS1]), "CDS");
		mapFrameVerbose.put(getCombinedFrame(FRAME_BYTEVAL[FRAME_BITSTOP], FRAME_BYTEVAL[FRAME_BITCDS2]), "CDS");
		mapFrameVerbose.put(getCombinedFrame(FRAME_BYTEVAL[FRAME_BITSTOP], FRAME_BYTEVAL[FRAME_BITSTOP]), "CDS");
		mapFrameVerbose.put(getCombinedFrame(FRAME_BYTEVAL[FRAME_BITSTOP], FRAME_BYTEVAL[FRAME_BIT3UTR]), "CDS-3UTR");
		
		mapFrameVerbose.put(getCombinedFrame(FRAME_BYTEVAL[FRAME_BIT3UTR], FRAME_BYTEVAL[FRAME_BIT3UTR]), "3UTR");

		mapFrameVerbose.put(getCombinedFrame(FRAME_BYTEVAL[FRAME_BITNC], FRAME_BYTEVAL[FRAME_BITNC]), "NC");

		mapFrameVerbose.put(getCombinedFrame(FRAME_BYTEVAL[FRAME_BITNI], FRAME_BYTEVAL[FRAME_BITNI]), "NA");
		

	}

	
	public static void main(String[] args) {
		int f5= FRAME_BITNI, f3= FRAME_BITNI;
		System.out.println(FRAME_BYTEVAL[f5]+" "+FRAME_BYTEVAL[f3]);
		int x= getCombinedFrame(FRAME_BYTEVAL[f5], FRAME_BYTEVAL[f3]);
		int new5= get5Frame(x), new3= get3Frame(x);
		System.out.println(new5+" "+new3);
	}
	
	/**
	 * @deprecated	
	 * @param prime5
	 * @param prime3
	 * @return
	 */
	public static int cdsCombine(byte prime5, byte prime3) {
		int combined= prime5 << FRAME_BYTEVAL.length;
		combined+= prime3;
		return combined;
	}
	
	public static void cdsAppend(StringBuilder sb, int combined) {
		byte prime3= (byte) (combined& (FRAME_BYTEVAL.length* 2- 1));
		byte prime5= (byte) (combined- prime3);
		int p5= Arrays.binarySearch(FRAME_BYTEVAL, prime5),
			p3= Arrays.binarySearch(FRAME_BYTEVAL, prime3);
		sb.append(FRAME_VERBOSE[p5]);
		sb.append("-");
		sb.append(FRAME_VERBOSE[p3]);		
	}
	
	
	public static boolean isAAposCharge(char aa) {
		if (AA_POS_CHARGE.indexOf(aa)>= 0)
			return true;
		return false;
	}
	public static boolean isAAnegCharge(char aa) {
		if (AA_NEG_CHARGE.indexOf(aa)>= 0)
			return true;
		return false;
	}
	public static boolean isAAnoCharge(char aa) {
		if (AA_POS_CHARGE.indexOf(aa)< 0&& AA_NEG_CHARGE.indexOf(aa)< 0)
			return true;
		return false;
	}
	public static boolean isAApolar(char aa) {
		if (AA_POLAR.indexOf(aa)>= 0)
			return true;
		return false;
	}
	public static boolean isAAunpolar(char aa) {
		return !isAApolar(aa);
	}
	public static boolean isAAsteric(char aa) {
		if (AA_STERIC.indexOf(aa)>= 0)
			return true;
		return false;
	}
	
	public static boolean isStop(char aa) {
		if (aa== '.')
			return true;
		return false;
	}

	public static char getTranslatedAA(String codon) {
		HashMap map= getCodonHash();
		Character c= ((Character) map.get(codon));
		if (c== null) {
			for (int i = 0; i < STOP_CODONS.length; i++) 
				if (STOP_CODONS[i].equals(codon))
					return '.';
			return '?';		// whatever
		}
		
			// else
		return c.charValue();
	}
	
	public static HashMap getCodonHash() {
		if (codonHash == null) {
			codonHash = new HashMap(ORDER_AA.length());
			for (int i = 0; i < CODONS_AA.length; i++) 
				for (int j = 0; j < CODONS_AA[i].length; j++) 
					codonHash.put(CODONS_AA[i][j], new Character(ORDER_AA.charAt(i)));
		}

		return codonHash;
	}
	
	public static int getProteinID(String someID) {
		someID= someID.toUpperCase();
		
			// RefSeq
		if (someID.startsWith("NP_")|| someID.startsWith("XP_")) {
			try {
				Integer.parseInt(someID.substring(3, someID.length()));
				return REFSEQ_ID;
			} catch (NumberFormatException e) {
				; //:)
			}
		}
		
			// Swissprot
		if (someID.startsWith("P")|| someID.startsWith("Q")|| someID.startsWith("O")) {
			try {
				Integer.parseInt(someID.substring(1, someID.length()));
				return SWISSPROT_ID;
			} catch (NumberFormatException e) {
				; //:)
			}
		}
		
		return UNKNOWN_ID;
	}
	
	public static final String START_CODON= "ATG";
	public static final String[] STOP_CODONS= new String[] {"TAA", "TGA", "TAG"};
	
	Vector proteinIDs= null;
	Transcript transcript= null;
	int splicedLength= -1;
	byte frame= -1;
	SpliceSite codonStart= null, codonStop= null;
	
	public static String[] extractCodons(String seq, int frame) {
		Vector codons= new Vector(seq.length()/ 3);
		int start= frame;
		if (frame> 0)
			--start;	// frame is 1-based
		for (int i = start; i <= seq.length()- 3; i+= 3) 
			codons.add(seq.substring(i, i+3));
		
		return (String[]) ArrayUtils.toField(codons);
	}
	
	public static int[] getCodonCount(String[] codons, String seq) {
		int[] res= new int[3];
		for (int i = 0; i < res.length; i++) 
			res[i]= getCodonCount(codons, seq, i);
		return res;
	}
	
	public static int getCodonCount(String[] codons, String seq, int frame) {
		seq= seq.toUpperCase();
		int cnt= 0;
		for (int i = frame; i < seq.length()- 2; i+= 3) {
			String cod= seq.substring(i, i+3);
			int j;
			for (j = 0; j < codons.length; j++) 
				if (codons[j].equals(cod))
					break;
			if (j< codons.length)
				++cnt;
		}
		return cnt;
	}

	/**
	 * condenses all 3 frames to one position array
	 * @param codons
	 * @param seq
	 * @return
	 */
	public static int[] getCodonPositions(String[] codons, String seq) {
		IntVector res= new IntVector();
		for (int i = 0; i < 3; i++) { 
			int[] pos= getCodonPositions(codons, seq, i);
			for (int j = 0; j < pos.length; j++) {
				int ins= java.util.Arrays.binarySearch(res.toIntArray(), pos[j]);
				if (ins< 0) 
					res.insert(pos[j], ins);
			}
		}
		return res.toIntArray();
	}
	
	public static int[] getCodonPositions(String[] codons, String seq, int frame) {
		seq= seq.toUpperCase();
		IntVector pos= new IntVector();
		for (int i = frame; i < seq.length()- 2; i+= 3) {
			String cod= seq.substring(i, i+3);
			int j;
			for (j = 0; j < codons.length; j++) 
				if (codons[j].equals(cod))
					break;
			if (j< codons.length)
				pos.add(i);
		}
		return pos.toIntArray();
	}
	
	public static int[] getStartCount(String seq) {
		return getCodonCount(new String[] {START_CODON}, seq);
	}
	
	public static int[] getStopCount(String seq) {
		return getCodonCount(STOP_CODONS, seq);
	}
	
	public Translation(Transcript newTranscript, int newStart, int newEnd, int newStrand) {
		super(newStart, newEnd, newStrand);
		this.transcript= newTranscript;
	}
	
	public Translation(Transcript newTranscript) {
		this.transcript= newTranscript;
		this.strand= getTranscript().getStrand();
		setID("CDS");
		setStrand(getTranscript().getStrand());
	}
	
	public int getSplicedLength() {
		if (splicedLength< 0) {
			DirectedRegion[] regs= transcript.getCDSRegions();
			splicedLength= 0;
			for (int i = 0; i < regs.length; i++) 
				splicedLength+= regs[i].getLength();
		}
		return splicedLength;
	}
	
	public DirectedRegion[] getExonicRegions() {
		Exon[] ex= getTranscript().getExons();
		Vector regV= new Vector(ex.length);
		for (int i = 0; i < ex.length; i++) {
			if (!ex[i].isCodingCompletely())
				continue;
			if (ex[i].isCoding5Prime()&& ex[i].isCoding3Prime()) {
				regV.add(ex[i]);
				continue;
			}
				
			DirectedRegion reg= new DirectedRegion(
					ex[i].get5PrimeCDS(), ex[i].get3PrimeCDS(), ex[i].getStrand());
			reg.setChromosome(getChromosome());
			reg.setSpecies(getSpecies());
			regV.add(reg);
		}
		return (DirectedRegion[]) ArrayUtils.toField(regV);
	}
	

	public Translation(Transcript newTranscript, String stableTranslationID) {

		this(newTranscript);
		addProteinID(stableTranslationID);
	}

	public String getChromosome() {
		return getTranscript().getChromosome();
	}
	
	/**
	 * @return Returns the transcript.
	 */
	public Transcript getTranscript() {
		return transcript;
	}
	
	public Species getSpecies() {
		return getTranscript().getGene().getSpecies();
	}
	/**
	 * @return Returns the translationID.
	 */
	public String[] getProteinID(int idCode) {
		Vector v= new Vector();
		for (int i = 0; proteinIDs!= null&& i < proteinIDs.size(); i++) {
			if (getProteinID((String) proteinIDs.elementAt(i))== idCode)
				v.add(proteinIDs.elementAt(i));
		}
		
		if (v.size()== 0)
			return null;
		return (String[]) ArrayUtils.toField(v);
	}
	
	public String[] getProteinIDsAll() {
		if (proteinIDs== null|| proteinIDs.size()== 0)
			return null;
		return (String[]) ArrayUtils.toField(proteinIDs);
	}
	/**
	 * @param translationID The translationID to set.
	 */
	public void addProteinID(String newTranslationID) {
		
		if (proteinIDs== null)
			proteinIDs= new Vector();
		proteinIDs.add(newTranslationID);
	}

	public String translate() {
		String s= getSplicedSequence();
		if (s.length()%3!= 0)
			return null;
		StringBuffer b= new StringBuffer(s.length()/3);
		for (int i = 0; i < s.length()- 3; i+= 3) {
			String c= CODON_HASH.get(s.substring(i, i+3).toUpperCase());
			if (c== null|| c.equals(""))	// stop
				return null;
			b.append(c);
		}
		if (s.length()< 3)
			return null;
		String c= CODON_HASH.get(s.substring(s.length()- 3, s.length()).toUpperCase());
		if (c!= null&& !c.equals(""))	// stop
			b.append(c);
		return b.toString();
	}
	
	public void setSplicedLength(int splicedLength) {
		this.splicedLength = splicedLength;
	}

	/**
	 * Never annotate an ATG starting internal of another CDS > 35 aa upstream
	 * of the ATG as is subject to NMD. [HAVANA]
	 * 
	 * @param trans
	 * @param maxDistAA
	 * @return
	 */
	public Translation[] getUsORF() {
		
		Translation[] trns= getTranscript().getAllORFs();
		Vector uOrfV= new Vector();
		for (int i = 0; trns!= null&& i < trns.length; i++) {
			if (trns[i].get3PrimeEdge()< get5PrimeEdge())
				uOrfV.add(trns[i]);
		}
		
		Translation[] uOrf= (Translation[]) ArrayUtils.toField(uOrfV);
		return uOrf;
	}

	/**
	 * @return genomic pos of premature stops, sorted
	 */
	public int[] getPrematureStops() {
		Exon startEx= getFirstCodingExon();
		String seq= getSplicedSequence();
		int[] pos= getCodonPositions(STOP_CODONS, seq, getFrame(startEx));
		
			// correct to genomic coordinates
		for (int i = 0; i < pos.length; i++) 
			pos[i]= getGenomicPosition(pos[i]+ 2);	// genomic pos of end of stop-codon 
		return pos;
	}

	public byte getFrame(Exon ex) {
		int cdsStart= getTranscript().getExonicPosition(get5PrimeEdge());
		int exStart= getTranscript().getExonicPosition(ex.get5PrimeEdge());
		byte frame= getFrame();
		if (exStart> cdsStart) {	// cds starts before exon 5'boundary
			frame+= (exStart- cdsStart)% 3;	// frameshift
			if (frame> 2)
				frame%= 3;
		} //else
			//assert(frame>= 0);	// frame detection must have worked!
		return frame;
	}
	
	/**
	 * @return starting frame of translation
	 */
	public byte getFrame() {
		
		if (frame< 0) {
			if (getTranscript().getGene().getSpecies()== null&& Graph.overrideSequenceDirPath== null)
				return -1;
			frame= getFirstCodingExon().getFrame();	// not reliable, also not in Havana
			if (frame< 0) 	// if initialized..
				frame= 0;	// guess 0, good for predicted reading frames
			
			int tlnExStart= getTranscript().getExonicPosition(get5PrimeEdge())+ frame;
			String exSeq= null;
			try {
				exSeq= getTranscript().getSplicedSequence();
			} catch (Exception e) {
				; // :)
			}
			if (exSeq== null)
				return -1;
			String startCodon= null;
			try {
				startCodon= exSeq.substring(tlnExStart, tlnExStart+3);
			} catch (Exception e) {
				return frame; // no sequence available, have to trust the annotated frame
			}
			
			if (startCodon.equalsIgnoreCase(START_CODON))
				return frame;	// if there is a start codon here, trust the annotated frame
			
				// otherwise have to guess, look for frames without stop codons
			tlnExStart-= frame;
			int tlnExEnd= getTranscript().getExonicPosition(get3PrimeEdge());
			
			exSeq= exSeq.substring(tlnExStart, tlnExEnd+ 1);	// max tln seq, evtlly includes term stop
			int[][] stops= new int[3][];
			for (int i = 0; i < stops.length; i++) 
				stops[i]= getCodonPositions(STOP_CODONS, exSeq, i);	// stops in all 3 frames
			
				// and now? take the one with least stops? the one with longest ORF? ???
			int min= Integer.MAX_VALUE;
			IntVector v= new IntVector();
			for (int i = 0; i < stops.length; i++) {
				if (stops[i]== null) {
					if (min> 0)
						v= new IntVector();
					min= 0;
					v.add(i);
					continue;
				}
				if (stops[i].length< min) {
					min= stops[i].length;
					v= new IntVector();
					v.add(i);
				} else if (stops[i].length== min) 
					v.add(i);
			}
			
				// unique one
			if (v.size()== 1) {
				frame= (byte) v.get(0);
				return frame;
			}
			
			for (int i = 0; i < v.size(); i++) {
				if (v.get(i)== frame)
					return frame;	// trust frame annotation, if no inframe stop
			}
			
			assert(true);
				// what if (the 2 not as "frame" annotated) frames have equal count of stop codons?		
//			int maxORFFrame= -1;
//			int maxORFEnd= Integer.MIN_VALUE;
//			for (int i = 0; i < v.size(); i++) 
//				if (stops[v.get(i)][0]> maxORFEnd) {	// actually, this should be the longest reading frame.. 
//					maxORFEnd= stops[v.get(i)][0];			// maybe not, we keep for the moment this tln, so assume start here/us
//					maxORFFrame= v.get(i);
//				}
//			return maxORFFrame;
		}
		
		return frame;
	}
	
	public String getSplicedSequence() {
		DirectedRegion[] regs= getExonicRegions();	// not sorted?
		if (regs== null)
			return "";
		java.util.Arrays.sort(regs, new DirectedRegion.DirectedPositionComparator());
		StringBuffer sb= new StringBuffer();
		for (int i = 0; i < regs.length; i++) 
			sb.append(Graph.readSequence(regs[i]));
		return sb.toString();
	}
	
	public int getGenomicPosition(int exonPos) {

		Exon[] exons= getTranscript().getExons();
		
			// find containing exon
		int x;
		int dist= 0;
		for (x = 0; dist<= exonPos&& x < exons.length; x++) 
			dist+= getCDSLength(exons[x]);
		if (x> 0) {
			--x;
			dist-= getCDSLength(exons[x]);
		}
		
		int genPos= Math.max(exons[x].get5PrimeEdge(), get5PrimeEdge())+ (exonPos- dist);		
		return genPos;
	}
	
	/**
	 * to also deal with predicted reading frames
	 */
	public int getCDSLength(Exon ex) {
		if (!this.overlaps(ex))
			return 0;
		DirectedRegion reg= this.intersect(ex);
		return reg.getLength();
	}
	

	public Exon getFirstCodingExon() {
		Exon[] ex= getTranscript().getExons();
		int i;
		for (i = 0; i < ex.length; i++) 
			if (ex[i].get3PrimeEdge()> get5PrimeEdge())		// better than exon.isCodin()
				break;
		if (i< ex.length)
			return ex[i];
		return null;
	}
	
	/**
	 * 
	 * @return 0-based position in the reading frame
	 */
	public int getTranslatedPosition(int genomicPos) {
		int a= getTranscript().getExonicPosition(genomicPos),
			b= getTranscript().getExonicPosition(get5PrimeEdge());
		int x= (a- b);
		return x;
	}
	
	/**
	 * 
	 * @param genomicPos
	 * @return 0,1, or 2
	 */
	public int getFrameAtPosition(int genomicPos) {
		return ((getTranslatedPosition(genomicPos))% 3);	// translated pos is 0-based, 20101028, killed +1 for trans.pos	
	}
	
	public int getFrameOrRegion(int genomicPos) {
		//int diff5= genomicPos- get5PrimeEdge() TODO efficiency 
		if (genomicPos< get5PrimeEdge())
			return FRAME_BYTEVAL[FRAME_BIT5UTR];
		if (genomicPos< get5PrimeEdge())
			return FRAME_BYTEVAL[FRAME_BITSTART];
		if (genomicPos> get3PrimeEdge())
			return FRAME_BYTEVAL[FRAME_BIT3UTR];
		if (genomicPos< get5PrimeEdge())
			return FRAME_BYTEVAL[FRAME_BITSTOP];
		int cdsFrame= getFrameAtPosition(genomicPos);
		return FRAME_BYTEVAL[cdsFrame];
	}
	
	public boolean isOpenEnded5() {
		return (get5PrimeEdge()== transcript.get5PrimeEdge());
	}

	public boolean isOpenEnded3() {
		return (get3PrimeEdge()== transcript.get3PrimeEdge());
	}
	
	public boolean isOpenEnded() {
		// TODO: check seq for ATG/stop
		if (isOpenEnded5()|| isOpenEnded3())
			return true;
		return false;
	}
	
	public int get3PrimeCDS(Exon ex) {
		
		int x= 0;
		for (x = 0; x < getTranscript().getExons().length; x++) 
			if (getTranscript().getExons()[x]== ex)
				break;
		if (x== getTranscript().getExons().length)
			return -1;	// exon not found
		
		
		if (ex.get3PrimeEdge()>= get5PrimeEdge()&& ex.get3PrimeEdge()<= get3PrimeEdge())
			return ex.get3PrimeEdge();
		if (get3PrimeEdge()>= ex.get5PrimeEdge()&& get3PrimeEdge()<= ex.get3PrimeEdge())
			return get3PrimeEdge();
		return -1;
	}
	
	public int get5PrimeCDS(Exon ex) {
		
		int x= 0;
		for (x = 0; x < getTranscript().getExons().length; x++) 
			if (getTranscript().getExons()[x]== ex)
				break;
		if (x== getTranscript().getExons().length)
			return -1;	// exon not found
		
		
		if (ex.get5PrimeEdge()>= get5PrimeEdge()&& ex.get5PrimeEdge()<= get3PrimeEdge())
			return ex.get5PrimeEdge();
		if (get5PrimeEdge()>= ex.get5PrimeEdge()&& get5PrimeEdge()<= ex.get3PrimeEdge())
			return get5PrimeEdge();
		return -1;
	}
	
	public int getStartCDS(Exon ex) {
		if (getTranscript().isForward())
			return get5PrimeCDS(ex);	// here not getter method, avoid init for buildup
		return get3PrimeCDS(ex);
	}

	public int getEndCDS(Exon ex) {
		if (getTranscript().isForward())
			return get3PrimeCDS(ex);	// here not getter method, avoid init for buildup
		return get5PrimeCDS(ex);
	}

	public SpliceSite getCodonStart() {
		return codonStart;
	}

	public void setCodonStart(SpliceSite codonStart) {
		this.codonStart = codonStart;
	}

	public SpliceSite getCodonStop() {
		return codonStop;
	}

	public void setCodonStop(SpliceSite codonStop) {
		this.codonStop = codonStop;
	}

	public static int getCombinedFrame(int frame5, int frame3) {
		int result= frame5;
		result= (frame5 << FRAME_BYTEVAL.length) | frame3;
		return result;
	}
	
	public static int get5Frame(int combined) {
		if (combined== 131328)
			return 256;	// overflow, 2x 256
		int result= combined>> FRAME_BYTEVAL.length;		
		return result & 0xFF;
	}

	public static int get3Frame(int combined) {
		if (combined== 131328)
			return 256;	// overflow, 2x 256
		int result= combined & 0xFF;
		return result;
	}
	
	public static String getFrameVerbose(byte frameByte) {
		int ii= (int) frameByte;	// dont rely on automatic widening!
		String s= mapFrameVerbose.get(ii);
		return s;
	}
	
	public static String getFrameVerbose(int frame) {
		String s=  mapFrameVerbose.get(frame);
		return s;
	}
	
	public static String getFrameVerbose(int currFrame5, int frame3) {
		int combined53= getCombinedFrame(currFrame5, frame3);
		assert(mapFrameVerbose.containsKey(combined53));
		String s=  mapFrameVerbose.get(combined53);
		return s;
	}

	public static boolean isCDS(int frame) {
		if (frame== FRAME_BYTEVAL[FRAME_BITCDS0]
			|| frame== FRAME_BYTEVAL[FRAME_BITCDS1]
			|| frame== FRAME_BYTEVAL[FRAME_BITCDS2]
			|| frame== FRAME_BYTEVAL[FRAME_BITSTART])
			return true;
		return false;
	}

	public static int findStop(String seq) {
		for (int i = 0; i < seq.length()- 3; i+=3) {
			String cc= seq.substring(i, i+ 3);
			for (int j = 0; j < STOP_CODONS.length; j++) {
				if (cc.equals(STOP_CODONS[j]))
					return i;
			}
		}
		return -1;
	}

	public static int findStart(String seq) {
		int lastStart= -2;
		for (int i = seq.length()- 3; i>= 0; i-=3) {
			String cc= seq.substring(i, i+ 3);
			if (cc.equals(START_CODON)) {
				lastStart= i;
				continue;
			}
			for (int j = 0; j < STOP_CODONS.length; j++) {
				if (cc.equals(STOP_CODONS[j]))
					return lastStart;
			}
		}
		return -1;
	}



	
}
