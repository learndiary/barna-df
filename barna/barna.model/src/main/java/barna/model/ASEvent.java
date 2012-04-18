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
 * Created on Feb 23, 2006
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
package barna.model;

import barna.commons.utils.ArrayUtils;
import barna.commons.utils.StringUtils;
import barna.model.commons.IntVector;
import barna.model.tools.NMDSimulator;

import java.util.*;

/**
 * 
 * <tt>07/11/15</tt>: 	<code>getScore()</code> implemented.
 * 						<code>getEventRegion()</code> adapted to correct for infinity flanks.
 * <tt>07/11/22</tt>:	implemented <code>getType()</code>, 
 * 						also bugfix in <code>isASevent()</code>, compared vs region of trpt w var SS
 * <tt>07/11/26</tt>:	new getRegions(byte ID)
 * 
 * 
 * @author msammeth
 */
public class ASEvent {

	public static final String GTF_FEATURE_ASEVENT= "as_event",
		GTF_FEATURE_DSEVENT= "ds_event",
		GTF_FEATURE_VSEVENT= "vs_event",
		GTF_FEATURE_NOEVENT= "no_event";
	public static final String GTF_ATTRIBUTE_TAG_STRUCTURE= "structure";
	public static final String GTF_ATTRIBUTE_TAG_SPLICE_CHAIN= "splice_chain";
	public static final String GTF_ATTRIBUTE_TAG_SOURCES= "sources";
	public static final String GTF_ATTRIBUTE_TAG_DIMENSION= "dimension";
	public static final String GTF_ATTRIBUTE_TAG_DEGREE= "degree";
	public static final String GTF_ATTRIBUTE_TAG_FLANKS= "flanks",
	GTF_ATTRIBUTE_SITE_SEQ= "seqsite",
	GTF_ATTRIBUTE_FLANK_MODE= "flank_mode",
	GTF_ATTRIBUTE_CDS_LOC= "cdsLoc",	// TODO not optimal, maybe rename
	GTF_ATTRIBUTE_CDS_VARIANTS= "cdsVariants",	// TODO not optimal, maybe rename
	GTF_ATTRIBUTE_CDS_IMPACT= "cdsImpact";
    public final static String ATTRIBUTE_TAG_LOCALIZATION= "localization";
    public final static char STRAND_SYMBOL_POS= '+', STRAND_SYMBOL_NEG= '-', SYMBOL_NOT_INITED= '.';
    public final static String TRANSCRIPT_ID_TAG= "transcript_id";
    public final static String GENE_ID_TAG= "gene_id", LOCUS_ID_TAG= "locus_id";
	
	public static final byte TYPE_NOT_INITED= -1, TYPE_UNDEFINED= 0, 
		TYPE_VS_EVENT= 1, TYPE_DS_EVENT= 2, TYPE_AS_EVENT= 3;

	public static final byte TYPE_NOTINITED= 0,  
		TYPE_SPLICING= 1, TYPE_INITIATION= 1<<1, TYPE_CLEAVAGE= 1<<2, TYPE_ALTERNATIVE= 1<<3, TYPE_PARTIAL= 1<<4, TYPE_FUSION= 1<<5, TYPE_NOTDEFINED= 1<<6;
	
	public static final HashMap<String, String> EVENT_DESCRIPTION_MAP= new HashMap<String, String>(5,1f);
	static {
		EVENT_DESCRIPTION_MAP.put("0,1-2^","skip 1 exon");
		EVENT_DESCRIPTION_MAP.put("1-,2-","alt acceptor");
		EVENT_DESCRIPTION_MAP.put("1^,2^","alt donor");
		EVENT_DESCRIPTION_MAP.put("1-2^,3-4^","mutually exclusive");
		EVENT_DESCRIPTION_MAP.put("0,1^2-","intron retention");
		EVENT_DESCRIPTION_MAP.put("0,1-2^3-4^","skip 2 exons");
	};

	public static boolean checkDimension= true, checkNMD= false, checkLocalization= false, check3Pcomplete= false, outputFlankMode= false;
	
	
	byte bubbleDimension= 0;
	HashMap<String, String> attributeMap= null; 
	float score= -1f;
	
	int[] frame3, frame5;
	boolean[] cdsValid;

    public void setType(byte type) {
        this.type = type;
    }

    public static class SpliceChainComparator implements Comparator<SpliceSite[]> {
		public int compare(SpliceSite[] o1, SpliceSite[] o2) {
			if (o1== null)	// empty chains at the end
				return 1;
			if (o2== null)
				return -1;

			int min= Math.min(o1.length, o2.length);
			for (int i = 0; i < min; i++) {
				if (o1[i].getPos()< o2[i].getPos())
					return -1;
				else if (o2[i].getPos()< o1[i].getPos())
					return 1;
			}
			
			if (o1.length> min)
				return 1;
			if (o2.length> min)
				return -1;
				// equal length and equal ss
			return 0;
		}
	}
	
	public static SpliceChainComparator defaultSpliceChainComparator= new SpliceChainComparator();

	
	Transcript[][] trpts;
	SpliceSite[][] spliceChains;	// sorted!
	String stringRep= null;
	SpliceSite src, snk;
	byte type= TYPE_NOT_INITED;
	
	public ASEvent() {
	}
	

	public void collapse(SpliceSite[][] reference) {
		Vector<SpliceSite[]> v= new Vector<SpliceSite[]>(reference.length);
		for (int i = 0; i < reference.length; i++) {
			boolean add= true;
			for (int j = i-1; j> 0; --j) {
				if (defaultSpliceChainComparator.compare(reference[i], reference[j])== 0) {
					add= false;
					break;
				}
			}
			if (add) 
				v.add(reference[i]);			
		}
				
	}
	
	
	public void addAttribute(String key, String val) {
		if (attributeMap== null)
			attributeMap= new HashMap<String, String>();
		attributeMap.put(key, val);
	}
	
	public String getAttribute(String key) {
		if (attributeMap== null)
			return null;
		return attributeMap.get(key);
	}
	

	public ASEvent(Transcript[][] newTrpts, SpliceSite[][] ss) {
		SpliceSite[] firsts= new SpliceSite[ss.length];
		HashMap<SpliceSite[], Transcript[]> map= new HashMap<SpliceSite[], Transcript[]>(ss.length, 1f);
		for (int i = 0; i < ss.length; i++) {
			map.put(ss[i], newTrpts[i]);	// null still can only be one, we hope
		}
		Arrays.sort(ss, defaultSpliceChainComparator);
		
		trpts= new Transcript[newTrpts.length][];
		for (int i = 0; i < ss.length; i++) {
			trpts[i]= map.get(ss[i]);
		}
		spliceChains= ss;
	}
	
	public boolean isComplete() {
		for (int i = 0; i < spliceChains.length; i++) {
			if (spliceChains[i].length== 0)
				continue;
			if ((!spliceChains[i][0].isSpliceSite())|| (!spliceChains[i][spliceChains[i].length- 1].isSpliceSite()))
				return false;
		}
		return true;
	}
	
	public void init() {
		
		if (checkLocalization) {
			addAttribute(ATTRIBUTE_TAG_LOCALIZATION, toStringLocalization());
		}
		
		
		if (checkNMD) {
			StringBuffer buf= new StringBuffer();
			for (int i = 0; i < trpts.length; i++) {
				for (int j = 0; j < trpts[i].length; j++) {
					if (!trpts[i][j].isCoding())
						continue;
					NMDSimulator nmdSim= new NMDSimulator(trpts[i][j]);
					if (nmdSim.isNMD(trpts[i][j].getTranslations()[0])) {
						buf.append(trpts[i][j].getTranscriptID());
						buf.append('/');
					}
				}
				if (buf.length()> 0&& buf.charAt(buf.length()-1)== '/')
					buf.deleteCharAt(buf.length()- 1);
				buf.append(',');
			}
			buf.deleteCharAt(buf.length()-1);
			addAttribute(Transcript.GTF_ATTRIBUTE_NMD_TAG, buf.toString());
		}
		
		if (check3Pcomplete) {
			StringBuffer buf= new StringBuffer();
			for (int i = 0; i < trpts.length; i++) {
				for (int j = 0; j < trpts[i].length; j++) {
					if (trpts[i][j].is3Pcomplete()) {
						buf.append(trpts[i][j].getTranscriptID());
						buf.append('/');
					}
				}
				if (buf.length()> 0&& buf.charAt(buf.length()-1)== '/')
					buf.deleteCharAt(buf.length()- 1);
				buf.append(',');
			}
			buf.deleteCharAt(buf.length()-1);
			addAttribute(Transcript.GTF_ATTRIBUTE_3COMPLETE_TAG, buf.toString());
		}
	}
	
	/**
	 * 0= not included
	 * 1= overlapping
	 * 2= included
	 * @deprecated get rid of that
	 */
	public int inCDS() {
		int min= Integer.MAX_VALUE;
		int max= Integer.MIN_VALUE;
		for (int i = 0; i < spliceChains.length; i++) {
			if (spliceChains[i].length== 0)
				continue;
			if (spliceChains[i][0].getPos()< min)
				min= spliceChains[i][0].getPos();
			if (spliceChains[i][spliceChains[i].length-1].getPos()> max)
				max= spliceChains[i][spliceChains[i].length-1].getPos();
		}
		
		DirectedRegion[] cdss= trpts[0][0].getGene().getCDSS();	// assume non-overlapping
		for (int j = 0; j < cdss.length; j++) {
			boolean minContained= cdss[j].contains(min);
			boolean maxContained= cdss[j].contains(max);
			if (minContained&& maxContained)
				return 2;
			if (minContained|| maxContained || (min< cdss[j].get5PrimeEdge()&& max> cdss[j].get3PrimeEdge()))
				return 1;
			
		}
		return 0;
	}
	
	public void setAnchors(SpliceSite src, SpliceSite snk) {
		this.src= src;
		this.snk= snk;
	}
	
	public String getTypeDescription() {
		if (getType()== TYPE_AS_EVENT)
			return GTF_FEATURE_ASEVENT;
		if (getType()== TYPE_DS_EVENT)
			return GTF_FEATURE_DSEVENT;
		if (getType()== TYPE_VS_EVENT)
			return GTF_FEATURE_VSEVENT;

		return GTF_FEATURE_NOEVENT;
	}
	
	public String toStringStructure_old() {
		int[] p= new int[spliceChains.length];
		StringBuffer[] sb= new StringBuffer[spliceChains.length];
		for (int i = 0; i < p.length; i++) { 
			p[i]= 0;
			sb[i]= new StringBuffer();
		}
		int cnt= 1;
		
		while (true) {
			IntVector nextI= new IntVector();
			int nextVal= Integer.MAX_VALUE;
			for (int i = 0; i < spliceChains.length; i++) {
				if (p[i]== spliceChains[i].length)
					continue;
				if (spliceChains[i][p[i]].getPos()< nextVal) {						
					nextI= new IntVector();
					nextI.add(i);
					nextVal= spliceChains[i][p[i]].getPos();
				} else if (spliceChains[i][p[i]].getPos()== nextVal)
					nextI.add(i);
			}
			
			for (int i = 0; i < nextI.size(); i++) {
				sb[nextI.get(i)].append(cnt);
				sb[nextI.get(i)].append(spliceChains[nextI.get(i)][p[nextI.get(i)]++].getSiteSymbol());
			}
			++cnt;	// now increment, for identical positions having same index
			
			int x= 0;
			for (; x < p.length; x++) 
				if (p[x]< spliceChains[x].length)
					break;
			if (x== p.length)
				break;
		}
	
		StringBuffer stringBuf= new StringBuffer();
		for (int i = 0; i < sb.length; i++) {
			if (sb[i].length()> 0)
				stringBuf.append(sb[i]);
			else 
				stringBuf.append('0');
			stringBuf.append(",");
		}
		stringBuf.deleteCharAt(stringBuf.length()- 1);
		
		return stringBuf.toString();
	}

	public String toStringStructure() {
		
		HashMap<SpliceSite, Integer> littleMap= new HashMap<SpliceSite, Integer>(getDegree(),1f);
		SpliceSite[] su= getSpliceUniverse(false);
		Arrays.sort(su, SpliceSite.getDefaultPositionTypeComparator());
		for (int i = 0; i < su.length; i++) 
			littleMap.put(su[i], (i+1));
		
		StringBuffer[] sb= new StringBuffer[spliceChains.length];		
		for (int i = 0; i < spliceChains.length; i++) {
			sb[i]= new StringBuffer();
			for (int j = 0; j < spliceChains[i].length; j++) {
				sb[i].append(littleMap.get(spliceChains[i][j]));
				sb[i].append(spliceChains[i][j].getSiteSymbol());
			}
		}

		StringBuffer stringBuf= new StringBuffer();
		for (int i = 0; i < sb.length; i++) {
			if (sb[i].length()> 0)
				stringBuf.append(sb[i]);
			else 
				stringBuf.append('0');
			stringBuf.append(",");
		}
		stringBuf.deleteCharAt(stringBuf.length()- 1);
		
		return stringBuf.toString();
	}
	
	/**
	 * One line schematical representation of the splicing variation.
	 */
	public String toString() {
		
		if (stringRep== null) {
			stringRep= toStringStructure();
		}
		
		return stringRep;
	}

	public String[] getGTFSources() {
		byte[] stypes= getSources();
		String[] sources= new String[stypes.length];
		for (int i = 0; i < sources.length; i++) 
			sources[i]= Transcript.TYPE_TO_ID[stypes[i]];		//Transcript.TYPE_TO_ID[type];		
		return sources;
	}
	
	public byte[] getSources() {
		byte[] types= new byte[trpts.length];
		for (int i = 0; i < types.length; i++) {
			types[i]= trpts[i][0].getSourceType();
			for (int j = 1; types[i]!= 0&& j < trpts[i].length; j++) 
				if (trpts[i][j].getSourceType()< types[i]) {
					types[i]= trpts[i][j].getSourceType();
				}
		}
		return types;
	}
	
	public int getDegree() {
		int deg= 0;
		for (int i = 0; i < spliceChains.length; i++) 
			deg+= spliceChains.length;
		return deg;
	}
	
	/**
	 * One line schematical representation of the splicing variation.
	 */
	public String toStringASTA() {
	
			// build final string
		StringBuffer result= new StringBuffer(Integer.toString(trpts.length));
		String evCode= toString();
		result.append("\t");
		result.append(evCode);
		result.append("\t");
		result.append("CDS=");
		result.append(inCDS());
		result.append("\t");
		result.append(trpts[0][0].getChromosome());
		result.append(":");
		int srcPos= src.getPos();
		if (srcPos== Integer.MIN_VALUE) {
			srcPos= Integer.MAX_VALUE;
			for (int i = 0; i < spliceChains.length; i++) 
				if (spliceChains[i].length> 0&& spliceChains[i][0].getPos()< srcPos)
					srcPos= spliceChains[i][0].getPos();
		}
		srcPos-= 1;
		
		int snkPos= snk.getPos();
		if (snkPos== Integer.MAX_VALUE) {
			snkPos= Integer.MIN_VALUE;
			for (int i = 0; i < spliceChains.length; i++) 
				if (spliceChains[i].length> 0&& spliceChains[i][spliceChains[i].length-1].getPos()> snkPos)
					snkPos= spliceChains[i][spliceChains[i].length-1].getPos();
		}
		snkPos+= 1;
			
		if (srcPos< 0) {
			result.append(Math.abs(snkPos));
			result.append("-");
			result.append(Math.abs(srcPos));
			result.append("\t");
			result.append("-");
		} else {
			result.append(srcPos);
			result.append("-");
			result.append(snkPos);
			result.append("\t");
			result.append("+");
		}
		
		
		for (int i = 0; i < spliceChains.length; i++) {
			result.append("\t");
			if (spliceChains[i].length> 0) {
				for (int x = 0; x < spliceChains[i].length; x++) {	// do not reverse order of SSs <-> breaks syntheny with rel. positions 
					if (spliceChains[i][x].getPos()>= 0)
						result.append(spliceChains[i][x].getPos());
					else
						result.append(Math.abs(spliceChains[i][x].getPos()));
					result.append(",");
				}
				result.deleteCharAt(result.length()- 1);
			}
		}
	
		for (int i = 0; i < trpts.length; i++) {
			result.append("\t");
			for (int j = 0; j < trpts[i].length; j++) { 
				result.append(trpts[i][j].getTranscriptID());
				result.append(",");
			}
			result.deleteCharAt(result.length()- 1);
		}
		
		String[] sources= getGTFSources();
		for (int i = 0; i < sources.length; i++) {
			result.append("\t");
			result.append(sources[i]);
		}
		
		return result.toString();
	}

	public Gene getGene() {
		return trpts[0][0].getGene();
	}
	
	public int getTranscriptSupport() {
		int sum= 0;
		for (int i = 0; i < trpts.length; i++) 
			sum+= trpts[i].length;
		return sum;
	}
	
	public SpliceSite getFirstVarSS() {
		for (int i = 0; i < getSpliceChains().length; i++) {
			if (getSpliceChains()[i].length!= 0)
				return getSpliceChains()[i][0];
		}
		return null;
	}

	public SpliceSite getLastVarSS() {
		int max= Integer.MIN_VALUE;
		SpliceSite last= null;
		for (int i = 0; i < getSpliceChains().length; i++) {
			if (getSpliceChains()[i].length== 0)
				continue;
			if (getSpliceChains()[i][getSpliceChains()[i].length-1].getPos()> max) {
				max= getSpliceChains()[i][getSpliceChains()[i].length-1].getPos();
				last= getSpliceChains()[i][getSpliceChains()[i].length-1];
			}
		}
		return last;
	}
	public String getGTFSource() {
		String[] sources= getGTFSources();
		int  max= Integer.MIN_VALUE, maxIdx= 0;
		for (int i = 0; i < sources.length; i++) {
			int type= Transcript.getSourceType(sources[i]);
			if (type> max) {
				max= type;
				maxIdx= i;
			}
		}
		return sources[maxIdx];
	}
	
	public String toStringGTF() {
		return toStringGTF(false);
	}
	/**
	 * new GTF compliance
	 */
	public String toStringGTF(boolean outSeq) {
	
			// build final string
		StringBuffer result= new StringBuffer(getGene().getChromosome());
		result.append('\t');

		result.append(getGTFSource());
		result.append('\t');

		byte type= getType();
		if (type== TYPE_AS_EVENT)
			result.append(GTF_FEATURE_ASEVENT);
		else if (type== TYPE_DS_EVENT)
			result.append(GTF_FEATURE_DSEVENT);
		else if (type== TYPE_VS_EVENT)
			result.append(GTF_FEATURE_VSEVENT);
		else 
			result.append(GTF_FEATURE_NOEVENT);
		result.append('\t');


		int strand= getGene().getStrand();
		DirectedRegion reg= getRegionEvent();
		int srcPos= reg.get5PrimeEdge();
		int snkPos= reg.get3PrimeEdge();
		if (strand> 0) {
			result.append(srcPos);
			result.append('\t');
			result.append(snkPos);
			result.append('\t');
		} else {
			result.append(Math.abs(snkPos));
			result.append('\t');
			result.append(Math.abs(srcPos));
			result.append('\t');
		}
		
			// score
		float score= -1f; 	// getScore(); 091208 SLOW !!! check overlap()
		if (score== -1f)
			result.append(".");
		else
			result.append(StringUtils.fprint(getScore(), 5));
		result.append('\t');
		
		if (strand> 0)
			result.append(STRAND_SYMBOL_POS);
		else
			result.append(STRAND_SYMBOL_NEG);
		result.append('\t');
		
		result.append('.');
		result.append('\t');
			
			// attributes
		result.append(TRANSCRIPT_ID_TAG);
		result.append(' ');
		result.append('\"');
		for (int i = 0; i < trpts.length; i++) {
			for (int j = 0; j < trpts[i].length; j++) {
				result.append(trpts[i][j].getTranscriptID());
				result.append('/');
			}
			result.replace(result.length()- 1, result.length(), ",");
		}
		result.deleteCharAt(result.length()- 1);
		result.append('\"');
		result.append(';');
		result.append(' ');

		result.append(GENE_ID_TAG);
		result.append(' ');
		result.append('\"');
		result.append(getGene().getGeneID());	// getGene().getReferenceTranscript()
		result.append('\"');
		result.append(';');
		result.append(' ');

		result.append(GTF_ATTRIBUTE_TAG_FLANKS);
		result.append(' ');
		result.append('\"');
		if (getSrc().getPos()== Integer.MIN_VALUE)
			result.append("null");
		else {
			result.append(Math.abs(getSrc().getPos()));
			result.append(getSrc().getSiteSymbol());
		}
		result.append(",");
		if (getSnk().getPos()== Integer.MAX_VALUE)
			result.append("null");
		else {
			result.append(Math.abs(getSnk().getPos()));
			result.append(getSnk().getSiteSymbol());
		}
		result.append('\"');
		result.append(';');
		result.append(' ');

		
		result.append(GTF_ATTRIBUTE_TAG_STRUCTURE);
		result.append(' ');
		result.append('\"');
		String evCode= toString();
		result.append(evCode);
		result.append('\"');
		result.append(';');
		result.append(' ');

		result.append(GTF_ATTRIBUTE_TAG_SPLICE_CHAIN);	// 5' -> 3'
		result.append(' ');
		result.append('\"');
		for (int i = 0; i < spliceChains.length; i++) {
			for (int j = 0; j < spliceChains[i].length; j++) {
				result.append(Math.abs(spliceChains[i][j].getPos()));
				result.append(spliceChains[i][j].getSiteSymbol()); //'/'
			}
//			if (spliceChains[i].length> 0)
//				result.deleteCharAt(result.length()- 1);
			result.append(',');
		}
		result.deleteCharAt(result.length()- 1);
		result.append('\"');
		result.append(';');
		result.append(' ');

		result.append(GTF_ATTRIBUTE_TAG_SOURCES);
		result.append(' ');
		result.append('\"');
		String[] sources= getGTFSources();
		for (int i = 0; i < sources.length; i++) {
			result.append(sources[i]);
			result.append(',');
		}
		result.deleteCharAt(result.length()-1);
		result.append('\"');
		result.append(';');
		result.append(' ');

		addAttribute(GTF_ATTRIBUTE_TAG_DEGREE, new Integer(getDegree()).toString()); 
		if (getAttribute(GTF_ATTRIBUTE_TAG_DIMENSION)== null)
			addAttribute(GTF_ATTRIBUTE_TAG_DIMENSION, new Integer(trpts.length).toString());
		
		if (outSeq) {
			SpliceSite[] su= getSpliceUniverse(true);
			StringBuffer sb= new StringBuffer(su.length* 25);
			for (int i = 0; i < su.length; i++) {
				sb.append(i);
				sb.append(su[i].getSiteSymbol());
				if (su[i].isTSS()|| su[i].isTES()) {
					sb.append(",");
					continue;
				}
				String s= null;
				if (su[i].isDonor())
				 s= Graph.readSequence(su[i], 15, 4); // us: 3exonic ss+ 4 codons
				else
				 s= Graph.readSequence(su[i], 15, 13);	// ds: 1exonic ss+ 4 codons
				sb.append(s);
				sb.append(",");
			}
			sb.deleteCharAt(sb.length()-1);
			addAttribute(GTF_ATTRIBUTE_SITE_SEQ, sb.toString());
		}
		
		if (outputFlankMode) {
			result.append(GTF_ATTRIBUTE_FLANK_MODE);
			result.append(" \"");
			if (getSrc().getPos()== Integer.MIN_VALUE)
				result.append("null");
			else {
				String mode= "constitutive";
				for (int i = 0; i < getGene().getTranscripts().length; i++) {
					Transcript t= getGene().getTranscripts()[i];
					if (t.get5PrimeEdge()<= getSrc().getPos()&&
							t.get3PrimeEdge()>= getSrc().getPos()) {
						int x;
						for (x = 0; x < getSrc().getTranscripts().length; x++) 
							if (getSrc().getTranscripts()[x]== t)
								break;
						if (x== getSrc().getTranscripts().length) {
							mode= "alternative";
							break;
						}
					}
				}
				result.append(mode);
			}
			result.append(",");
			if (getSnk().getPos()== Integer.MAX_VALUE)
				result.append("null");
			else {
				String mode= "constitutive";
				for (int i = 0; i < getGene().getTranscripts().length; i++) {
					Transcript t= getGene().getTranscripts()[i];
					if (t.get5PrimeEdge()<= getSnk().getPos()&&
							t.get3PrimeEdge()>= getSnk().getPos()) {
						int x;
						for (x = 0; x < getSnk().getTranscripts().length; x++) 
							if (getSnk().getTranscripts()[x]== t)
								break;
						if (x== getSnk().getTranscripts().length) {
							mode= "alternative";
							break;
						}
					}
				}
				result.append(mode);
			}
			result.append("\"; ");
		}
		
		if (cdsImpact!= null) {
			result.append(GTF_ATTRIBUTE_CDS_IMPACT);
			result.append(" \"");
			result.append(cdsImpact);
			result.append("\"; ");
			
			result.append(GTF_ATTRIBUTE_CDS_VARIANTS);
			result.append(" \"");
			result.append(cdsVariants);
			result.append("\"; ");
		}
		
		if (cdsLoc!= null) {
			result.append(GTF_ATTRIBUTE_CDS_LOC);
			result.append(" \"");
			for (int j = 0; j < cdsLoc.length; j++) {
				if (cdsLoc[j]==null)
					result.append("-,");
				else
					result.append(cdsLoc[j]+",");
			}
			result.deleteCharAt(result.length()- 1);
			result.append("\"; ");
			
/*			result.append(GTF_ATTRIBUTE_CDS_IMPACT);
			result.append(" \"");
			HashMap<String,String> littleMap= new HashMap<String, String>(8);
			for (int i = 0; i < cdsLoc.length; i++) 
				littleMap.put(cdsLoc[i], null);
			Object[] oo= littleMap.keySet().toArray();
			Arrays.sort(oo);
			for (int i = 0; i < oo.length; i++) 
				result.append(oo[i]+",");
			result.deleteCharAt(result.length()- 1);
			result.append("\"; ");
*/			
		}
		
		
		if (frame5!= null) {	
			// boundary frames
			result.append("FRAME5");
			result.append(" \"");
			for (int j = 0; j < frame5.length; j++) {
				String s= Translation.getFrameVerbose(frame5[j]);
				if (!cdsValid[j])
					s= s.toLowerCase();
				if (s== null|| s.contains("null"))
					System.currentTimeMillis();
				result.append(s+ ",");
			}
			result.deleteCharAt(result.length()- 1);
			result.append("\"; ");
		}
		
		if (frame3!= null) {
			
			result.append("FRAME3");
			result.append(" \"");
			for (int j = 0; j < frame3.length; j++) {
				String s= Translation.getFrameVerbose(frame3[j]);
				if (!cdsValid[j])
					s= s.toLowerCase();
				if (s== null|| s.contains("null"))
					System.currentTimeMillis();
				result.append(s+ ",");
			}
			result.deleteCharAt(result.length()- 1);
			result.append("\"; ");
		}
		
		result.append("NMorMRNA \""+getClaudiaNMorMRNAref()+"\"; ");
		
		
		Iterator<String> iter= getAttributeMap().keySet().iterator();
		while (iter.hasNext()) {
			String key= iter.next();
			result.append(key);
			result.append(' ');
			result.append('\"');
			result.append(getAttribute(key));
			result.append('\"');
			result.append(';');
			result.append(' ');
		}
		
		
		
		return result.toString();
	}
	
	public SpliceSite getSnk() {
		return snk;
	}

	public void setSnk(SpliceSite snk) {
		this.snk = snk;
	}

	public SpliceSite getSrc() {
		return src;
	}
	
	/**
	 * Adapted from Jim Kent:
	 * Confidence note from the function that returns confidence:
	 * Return the score for this cassette exon. Want to have cassette exons
	 * that are present in multiple transcripts and that are not present in multiple
	 * transcripts. We want to see both forms of the cassette exon, we don't want to have
	 * one outlier be chosen. Thus we count the times that the exon is seen, we
	 * count the times that the exon isn't seen and we calculate a final score by:
	 * (seen + notseen + prior)/(abs(seen - notSeen+prior) + 1) . Thus larger scores are
	 * better, but tend to range between 1-2:
	 * 
	 * conf = (float)(seen + notSeen + prior)/(float)(abs(seen - notSeen) + prior);
	 * 
	 * Here:
	 * 
	 * conf= (float) (seen+ notseen)+ ovl_trpts / max_t1,t2(seen- notseen)+ ovl_trpts 
	 * 
	 * @return
	 */
	public float getScore() {
		if (score == -1f&& getGene()!= null) {
			DirectedRegion reg= getRegionEvent();
			int cnt= 0;
			for (int i = 0; i < getGene().getTranscripts().length; i++) 
				if (getGene().getTranscripts()[i].overlaps(reg))	// overlaps ok, not includes -> edge events!
					++cnt;
			
			int max= 0;
			for (int i = 0; i < trpts.length; i++) 
				for (int j = i+1; j < trpts.length; j++) 
					max= Math.max(max, Math.abs(trpts[i].length- trpts[j].length));				
			
			score= ((float) (getTranscriptSupport()+ cnt))/ ((float) max+ cnt);			
		}

		return score;
	}

	public void setSrc(SpliceSite src) {
		this.src = src;
	}

	public int getDimension() {
		
		if (spliceChains!= null)
			return spliceChains.length;
		
		if (getAttribute(GTF_ATTRIBUTE_TAG_DIMENSION)== null)
			return 0;
		
		String dim= getAttribute(GTF_ATTRIBUTE_TAG_DIMENSION);
		dim= dim.substring(dim.indexOf('_')+ 1);
		return Integer.parseInt(dim);
	}

	public void setDimension(int bubbleDimension) {
		if (checkDimension)
			addAttribute(GTF_ATTRIBUTE_TAG_DIMENSION, trpts.length+"_"+Integer.toString(bubbleDimension));
	}
	
	/**
	 * Gets the genomic region the event covers, extended up to the flanks.
	 * Corrects for infinity boundaries.
	 * @return
	 */
	public DirectedRegion getRegionEvent() {
		int srcPos= src.getPos();	// infinity has to be eliminated
		if (srcPos== Integer.MIN_VALUE)
			srcPos= getFirstVarSS().getPos();
		int snkPos= snk.getPos();		
		if (snkPos== Integer.MAX_VALUE)
			snkPos= getLastVarSS().getPos();

		DirectedRegion reg= new DirectedRegion(srcPos, snkPos, getGene().getStrand());
		reg.setChromosome(getGene().getChromosome());
		reg.setSpecies(getGene().getSpecies());
		return reg;
	}
	
	public static byte REGION_EVENT= 1, REGION_VAR_REGIONS= 2, REGION_EXONIC_REGIONS= 3, REGION_INTRONIC_REGIONS= 4, REGION_VAR_REGIONS_NOT_MISSING= 5;
	
	public Transcript getClaudiaNMorMRNAref() {
		Transcript txNM= null, txMRNA= null;
		for (int i = 0; txNM== null&& i < trpts.length; i++) {
			for (int j = 0; txNM== null&& j < trpts[i].length; j++) {
				byte type= trpts[i][j].getSourceType();
				if (txNM== null&& type== Transcript.ID_SRC_REFSEQ
						//&& trpts[i][j].getTranslations()!= null 
						&& trpts[i][j].getTranscriptID().startsWith("NM")) {
					txNM= trpts[i][j];
					break;
				}
				if (txMRNA== null&& type== Transcript.ID_SRC_MRNA) 
					txMRNA= trpts[i][j];
			}
		}
		if (txNM!= null)
			return txNM;
		//else
		return txMRNA;
	}
	
	
	/**
	 * returns every region between src and snk, if not +/- infinity
	 * @param id
	 * @return
	 */
	// TODO make more efficient ??!
	public DirectedRegion[] getRegion(byte id) {
		
		if (id== REGION_EVENT)
			return new DirectedRegion[] {getRegionEvent()};
		
		byte[] exonic= new byte[getSpliceChains().length];	// 0 out of space, -1 intronic, +1 exonic
		Vector<DirectedRegion> v= new Vector<DirectedRegion>();
		
			// 5' of first
		if (getSrc().getType()== SpliceSite.TYPE_NOT_INITED) {
			for (int i = 0; i < exonic.length; i++) 
				exonic[i]= 0;			
		} else {
			SpliceSite firstSS= getFirstVarSS();
			byte init= 1;
			if (firstSS.isLeftFlank())
				init= -1;
			for (int i = 0; i < exonic.length; i++) 
				exonic[i]= init;			
			if (((init== 1)&& id== REGION_EXONIC_REGIONS)|| ((init== -1)&& id== REGION_INTRONIC_REGIONS)) {				
				DirectedRegion reg= new DirectedRegion(getSrc().getPos()+1, firstSS.getPos()-1, getGene().getStrand());
				reg.setChromosome(getGene().getChromosome());
				reg.setSpecies(getGene().getSpecies());
				v.add(reg);
			}
		}
		
			// all var splice sites
		int[] p= new int[spliceChains.length];
		for (int i = 0; i < p.length; i++)  
			p[i]= 0;
		SpliceSite lastSite= null;
		Comparator compi= SpliceSite.getDefaultPositionTypeComparator();
		while (true) {					
			IntVector nextI= new IntVector();
			SpliceSite nextSite= null;
			for (int i = 0; i < spliceChains.length; i++) {				
				if (p[i]== spliceChains[i].length)
					continue;
				
				int res= -1;
				if (nextSite!= null)
					res= compi.compare(spliceChains[i][p[i]], nextSite);
				if (res< 0) {						
					nextI= new IntVector();
					nextI.add(i);
					nextSite= spliceChains[i][p[i]];
				} else if (res== 0)
					nextI.add(i);
			}
			
			byte init= 1;
			if (nextSite.isRightFlank()) {
				if (nextSite.isTES())
					init= 0;
				else
					init= -1;
			}
			
			if (lastSite!= null) {
				boolean collect= true;
				if (id== REGION_VAR_REGIONS|| id== REGION_VAR_REGIONS_NOT_MISSING) {
					boolean exonicOnce= false, intronicOnce= false, missingOnce= false;
					for (int i = 0; i < exonic.length; i++) {
						if (exonic[i]== 1)
							exonicOnce= true;
						else if (exonic[i]== -1)
							intronicOnce= true;
						else if (exonic[i]== 0)
							missingOnce= true;
					}
					if (id== REGION_VAR_REGIONS)
						collect= exonicOnce& (intronicOnce|| missingOnce);
					else if (id== REGION_VAR_REGIONS_NOT_MISSING)
						collect= exonicOnce& intronicOnce&& (!missingOnce);
				} else {
					for (int i = 0; i < exonic.length; i++) {
						if (id== REGION_EXONIC_REGIONS)
							collect&= exonic[i]== 1;
						else if (id== REGION_INTRONIC_REGIONS)
							collect&= exonic[i]== -1;
					}
				}
				
				if (collect) {
					int start= lastSite.getPos();
					int end= nextSite.getPos();
					if (lastSite.isRightFlank()) 
						++start;
					if (nextSite.isLeftFlank())
						--end;
					DirectedRegion reg= new DirectedRegion(start, end, getGene().getStrand());
					reg.setChromosome(getGene().getChromosome());
					reg.setSpecies(getGene().getSpecies());
					v.add(reg);
				}

			}
			
			for (int i = 0; i < nextI.size(); i++) { 
				exonic[nextI.get(i)]= init;
				++p[nextI.get(i)];
			}			
			lastSite= nextSite;
			
			int x= 0;
			for (; x < p.length; x++) 
				if (p[x]< spliceChains[x].length)
					break;
			if (x== p.length)
				break;
		}

		// 5' of last
		if (getSnk().getType()!= SpliceSite.TYPE_NOT_INITED) {
			SpliceSite nextSite= getSnk();
			boolean collect= true;
			if (id== REGION_VAR_REGIONS) {
				boolean exonicOnce= false, intronicOnce= false;
				for (int i = 0; i < exonic.length; i++) {
					if (exonic[i]== 1)
						exonicOnce= true;
					else if (exonic[i]== -1)
						intronicOnce= true;
				}
				collect= exonicOnce& intronicOnce;
			} else {
				for (int i = 0; i < exonic.length; i++) {
					if (id== REGION_EXONIC_REGIONS)
						collect&= exonic[i]== 1;
					else if (id== REGION_INTRONIC_REGIONS)
						collect&= exonic[i]== -1;
				}
			}
			
			if (collect) {
				int start= lastSite.getPos();
				int end= nextSite.getPos();
				if (lastSite.isRightFlank()) 
					++start;
				if (nextSite.isLeftFlank())
					--end;
				DirectedRegion reg= new DirectedRegion(start, end, getGene().getStrand());
				reg.setChromosome(getGene().getChromosome());
				reg.setSpecies(getGene().getSpecies());
				v.add(reg);
			}

		}
		
		DirectedRegion[] regs= new DirectedRegion[v.size()];
		for (int i = 0; i < regs.length; i++) 
			regs[i]= v.elementAt(i);
		return regs;
	}

	public HashMap<String, String> getAttributeMap() {
		return attributeMap;
	}

	public SpliceSite[][] getSpliceChains() {
		return spliceChains;
	}

	public void setSpliceChains(SpliceSite[][] spliceChains) {
		for (int i = 0; i < spliceChains.length; i++) {
			Arrays.sort(spliceChains[i], SpliceSite.getDefaultPositionTypeComparator());
		}
		this.spliceChains = spliceChains;
	}

	public Transcript[][] getTranscripts() {
		return trpts;
	}
	
	public Transcript[] getTranscripts(SpliceSite[] schain) {
		for (int i = 0; i < getSpliceChains().length; i++) {
			if (getSpliceChains()[i]== schain)
				return getTranscripts()[i];
		}
		return null;
	}

	public void setTranscripts(Transcript[][] trpts) {
		this.trpts = trpts;
	}
	
	public SpliceSite getSpliceSite(int pos) {
		for (int i = 0; i < getSpliceChains().length; i++) {
			SpliceSite dummy= new SpliceSite(pos, SpliceSite.TYPE_NOT_INITED, null);
			int p= Arrays.binarySearch(getSpliceChains()[i], dummy, SpliceSite.getDefaultPositionComparator());
			if (p> 0)
				return getSpliceChains()[i][p];
		}
		return null;
	}
	
	/**
	 * 
	 * @deprecated replaced by new method..
	 */
	public SpliceSite[] getSpliceUniverseWithFlanks() {
		HashMap<SpliceSite, Object> map= new HashMap<SpliceSite, Object>();
		map.put(getSrc(), null);
		map.put(getSnk(), null);
		for (int i = 0; i < spliceChains.length; i++) 
			for (int j = 0; j < spliceChains[i].length; j++) 
				map.put(spliceChains[i][j], null);
		
		Iterator<SpliceSite> iter= map.keySet().iterator();
		int x= 0;
		SpliceSite[] su= new SpliceSite[map.size()];
		while (iter.hasNext())
			su[x++]= iter.next();
		Arrays.sort(su, SpliceSite.getDefaultPositionTypeComparator());
		
		return su;
	}

	public SpliceSite[] getSpliceUniverse(boolean includeFlanks) {
		HashMap<SpliceSite, Object> map= new HashMap<SpliceSite, Object>();
		if (includeFlanks) {
			map.put(getSrc(), null);
			map.put(getSnk(), null);
		}
		for (int i = 0; i < spliceChains.length; i++) 
			for (int j = 0; j < spliceChains[i].length; j++) 
				map.put(spliceChains[i][j], null);
		
		Iterator<SpliceSite> iter= map.keySet().iterator();
		int x= 0;
		SpliceSite[] su= new SpliceSite[map.size()];
		while (iter.hasNext())
			su[x++]= iter.next();
		Arrays.sort(su, SpliceSite.getDefaultPositionTypeComparator());
		
		return su;
	}
	
	/**
	 * the regions where there is at least one exon
	 * @return
	 */
	public DirectedRegion[] getExonicRegions() {
		int[] p= new int[spliceChains.length];
		boolean[] exonic= new boolean[spliceChains.length];
		SpliceSite firstSS= getFirstVarSS();
		for (int i = 0; i < p.length; i++) { 
			p[i]= 0;
			if (firstSS.isLeftFlank())
				exonic[i]= false;
			else
				exonic[i]= true;
		}
		int cnt= 1;
		
		boolean lastState= false;
		for (int i = 0; i < exonic.length; i++) 
			lastState|= exonic[i];
		
		DirectedRegion reg= null;
		if (lastState) {
			reg= new DirectedRegion(getSrc().getPos(), getSrc().getPos(), getGene().getStrand());
			reg.setChromosome(getGene().getChromosome());
		}
		Vector<DirectedRegion> exRegV= new Vector<DirectedRegion>();
		
		while (true) {
				// get next(s)
			IntVector nextI= new IntVector();
			int nextVal= Integer.MAX_VALUE;
			for (int i = 0; i < spliceChains.length; i++) {
				if (p[i]== spliceChains[i].length)
					continue;
				if (spliceChains[i][p[i]].getPos()< nextVal) {						
					nextI= new IntVector();
					nextI.add(i);
					nextVal= spliceChains[i][p[i]].getPos();
				} else if (spliceChains[i][p[i]].getPos()== nextVal)
					nextI.add(i);
			}

				// new exonic state
			SpliceSite site= spliceChains[nextI.get(0)][p[nextI.get(0)]];
			for (int i = 0; i < nextI.size(); i++) {
				exonic[nextI.get(i)]= !exonic[nextI.get(i)];
				p[nextI.get(i)]++;
			}
							
			boolean state= false;	// new state in all
			for (int i = 0; i < exonic.length; i++) 
				state|= exonic[i];

				// state change ?
			if (state!= lastState) {
				if (state) { 	// new start of exonic region
					reg= new DirectedRegion(site.getPos(), site.getPos(), getGene().getStrand());
					reg.setChromosome(getGene().getChromosome());
				} else { // intron
					reg.set3PrimeEdge(site.getPos());
					exRegV.add(reg);
				}
			}
			lastState= state;
			
			int x= 0;
			for (; x < p.length; x++) 
				if (p[x]< spliceChains[x].length)
					break;
			if (x== p.length)
				break;
		}
		
			// close last exonic region
		if (getLastVarSS().isLeftFlank()) {
			reg.set3PrimeEdge(getSnk().getPos());
			exRegV.add(reg);
		}

		DirectedRegion[] exRegs= new DirectedRegion[exRegV.size()];
		for (int i = 0; i < exRegs.length; i++) 
			exRegs[i]= exRegV.elementAt(i);
		return exRegs;
	}
	
	public DirectedRegion[] getVariableExonicRegions() {
		return getRegion(REGION_VAR_REGIONS);
	}
	
	public DirectedRegion[] getVariableRegions() {
		return getRegion(REGION_VAR_REGIONS);
	}
	
	public DirectedRegion[] getVariableExonicRegions_last_how_did_this_shit_ever_work() {
		DirectedRegion[] reg1= getVariableRegions();
		DirectedRegion[] reg2= getExons();
		DirectedRegion[] regs= DirectedRegion.intersect(reg1, reg2);
		return regs;
	}
	
	
	public DirectedRegion[] getExons(boolean withFlanks) {
		
		Vector<DirectedRegion> v= new Vector<DirectedRegion>(2,2);
		SpliceSite[] su= getSpliceUniverse(withFlanks);
		for (int i = 0; i < spliceChains.length; i++) {
			int leftFlankPos= su[0].getPos();
			if (su[0].isRightFlank())
				++leftFlankPos;
			for (int j = 0; j < spliceChains[i].length; j++) {
				if (spliceChains[i][j].isLeftFlank())
					leftFlankPos= spliceChains[i][j].getPos();
				else if (spliceChains[i][j].getPos()>= leftFlankPos) {
					DirectedRegion reg= new DirectedRegion(leftFlankPos, spliceChains[i][j].getPos(), getGene().getStrand());
					reg.setChromosome(getGene().getChromosome());
					v.add(reg);	// not unique/sorted - see getExons()
				}
			}
			if ((spliceChains[i].length== 0&& su[0].isRightFlank())|| 
					(spliceChains[i].length> 0&& spliceChains[i][spliceChains[i].length-1].isLeftFlank()&&
					spliceChains[i][spliceChains[i].length-1].getPos()< su[su.length-1].getPos())) {
				int rightFlankPos= su[su.length-1].getPos();
				if (su[su.length-1].isLeftFlank())
					--rightFlankPos;
				DirectedRegion reg= new DirectedRegion(leftFlankPos, rightFlankPos, getGene().getStrand());
				reg.setChromosome(getGene().getChromosome());
				v.add(reg);	// not unique/sorted - see getExons()
			}				
		}
		
		DirectedRegion[] regs= new DirectedRegion[v.size()];
		for (int i = 0; i < regs.length; i++) 
			regs[i]= v.elementAt(i);
		
		return regs;
	}

	public DirectedRegion[] getExons() {
		
		Vector<DirectedRegion> v= new Vector<DirectedRegion>(2,2);
		SpliceSite[] su= getSpliceUniverse(false);
		for (int i = 0; i < spliceChains.length; i++) {
			int leftFlankPos= su[0].getPos();
			if (su[0].isRightFlank())
				++leftFlankPos;
			for (int j = 0; j < spliceChains[i].length; j++) {
				if (spliceChains[i][j].isLeftFlank())
					leftFlankPos= spliceChains[i][j].getPos();
				else if (spliceChains[i][j].getPos()>= leftFlankPos) {
					DirectedRegion reg= new DirectedRegion(leftFlankPos, spliceChains[i][j].getPos(), getGene().getStrand());
					reg.setChromosome(getGene().getChromosome());
					ArrayUtils.addUniqueSorted(v, reg, DirectedRegion.getDefaultPositionComparator());	// for uniqueness, sorting according to pos
				}
			}
			if ((spliceChains[i].length== 0&& su[0].isRightFlank())|| 
					(spliceChains[i].length> 0&& spliceChains[i][spliceChains[i].length-1].isLeftFlank()&&
					spliceChains[i][spliceChains[i].length-1].getPos()< su[su.length-1].getPos())) {
				int rightFlankPos= su[su.length-1].getPos();
				if (su[su.length-1].isLeftFlank())
					--rightFlankPos;
				DirectedRegion reg= new DirectedRegion(leftFlankPos, rightFlankPos, getGene().getStrand());
				reg.setChromosome(getGene().getChromosome());
				ArrayUtils.addUniqueSorted(v, reg, DirectedRegion.getDefaultPositionComparator());		// for uniqueness, sorting according to pos
			}				
		}
		
		DirectedRegion[] regs= new DirectedRegion[v.size()];
		for (int i = 0; i < regs.length; i++) 
			regs[i]= v.elementAt(i);
		
		return regs;
	}

	public DirectedRegion[] getIntrons() {
		
		Vector<DirectedRegion> v= new Vector<DirectedRegion>(2,2);
		
		for (int i = 0; i < spliceChains.length; i++) {
			int leftFlankPos= getSrc().getPos();
			for (int j = 0; j < spliceChains[i].length; j++) {
				if (spliceChains[i][j].isRightFlank()) 	// left side of an intron
					leftFlankPos= spliceChains[i][j].getPos();
				else if (spliceChains[i][j].getPos()>= leftFlankPos) {
					int rightFlankPos= spliceChains[i][j].getPos();
					++leftFlankPos;
					--rightFlankPos;
//					if (getGene().getStrand()> 0) {	// narrow intron 
//						++leftFlankPos;
//						--rightFlankPos;
//					} else {
//						--leftFlankPos;
//						++rightFlankPos;
//					}
					DirectedRegion reg= new DirectedRegion(leftFlankPos, rightFlankPos, getGene().getStrand());
					reg.setChromosome(getGene().getChromosome());
					// dont take that, we want the introns sorted to their occurence in the event
					//v= gphase.tools.Arrays.addAllUniqueSorted(v, reg, DirectedRegion.getDefaultPositionComparator());
					v.add(reg);
				}
			}
			if ((spliceChains[i].length== 0&& getSrc().isRightFlank())|| 
					(spliceChains[i].length> 0&& spliceChains[i][spliceChains[i].length-1].isRightFlank()&&
					spliceChains[i][spliceChains[i].length-1].getPos()< getSnk().getPos())) {
				int rightFlankPos= getSnk().getPos();
				++leftFlankPos;
				--rightFlankPos;
//				if (getGene().getStrand()>= 0) {	// narrow intron 
//					++leftFlankPos;
//					--rightFlankPos;
//				} else {
//					--leftFlankPos;
//					++rightFlankPos;
//				}
				DirectedRegion reg= new DirectedRegion(leftFlankPos, rightFlankPos, getGene().getStrand());
				reg.setChromosome(getGene().getChromosome());
				// dont take that, we want the introns sorted to their occurence in the event
				//v= gphase.tools.Arrays.addAllUniqueSorted(v, reg, DirectedRegion.getDefaultPositionComparator());
				v.add(reg);
			}				
		}
		
		DirectedRegion[] regs= new DirectedRegion[v.size()];
		for (int i = 0; i < regs.length; i++) 
			regs[i]= v.elementAt(i);
		
		return regs;
	}

	public DirectedRegion[] getVariableExonicRegions_not_working() {
		
		SpliceSite[] su= getSpliceUniverseWithFlanks();
		Vector<DirectedRegion> regV= new Vector<DirectedRegion>();
		for (int i = 1; i < su.length; i++) {	// TODO: actually from the 2nd onward
			byte ssType= SpliceSite.TYPE_NOT_INITED;
			for (int j = 0; j < spliceChains.length; j++) {
				int p= Arrays.binarySearch(spliceChains[j], su[i], SpliceSite.getDefaultPositionComparator());
				if (p< 0)
					p= -(p+1);
				SpliceSite lastSite= getSrc();
				if (p> 0)
					lastSite= spliceChains[j][--p]; // the last site before insertion point
				if (ssType== SpliceSite.TYPE_NOT_INITED)
					ssType= lastSite.getType();
				else if (ssType!= lastSite.getType()&& su[i].isRightFlank()) {	// diff to getVarRegs()
					int start= su[i-1].getPos();
					if (su[i-1].isRightFlank())
						--start;
					int end= su[i].getPos();
					if (su[i].isLeftFlank())
						--end;
					
					DirectedRegion reg= new DirectedRegion(start, end, getGene().getStrand());
					reg.setChromosome(getGene().getChromosome());
					regV.add(reg);
				}
			}
		}
		
		DirectedRegion[] regs= new DirectedRegion[regV.size()];
		for (int i = 0; i < regs.length; i++) 
			regs[i]= regV.elementAt(i);
		
		return regs;
	}

	public DirectedRegion[] getVariableRegions_last_how_did_this_shit_ever_work() {
		
		SpliceSite[] su= getSpliceUniverseWithFlanks();
		Vector<DirectedRegion> regV= new Vector<DirectedRegion>();
		for (int i = 1; i < su.length; i++) {	// TODO: actually from the 2nd onward
			byte ssType= SpliceSite.TYPE_NOT_INITED;
			for (int j = 0; j < spliceChains.length; j++) {
				int p= Arrays.binarySearch(spliceChains[j], su[i], SpliceSite.getDefaultPositionTypeComparator());
				if (p< 0)
					p= -(p+1);
				SpliceSite lastSite= getSrc();
				if (p> 0)
					lastSite= spliceChains[j][--p]; // the last site before insertion point
				if (ssType== SpliceSite.TYPE_NOT_INITED)
					ssType= lastSite.getType();
				else if (ssType!= lastSite.getType()) {
					int start= su[i-1].getPos();
					if (su[i-1].isRightFlank())
						--start;
					int end= su[i].getPos();
					if (su[i].isLeftFlank())
						--end;
					
					DirectedRegion reg= new DirectedRegion(start, end, getGene().getStrand());
					reg.setChromosome(getGene().getChromosome());
					regV.add(reg);
				}
			}
		}
		
		DirectedRegion[] regs= new DirectedRegion[regV.size()];
		for (int i = 0; i < regs.length; i++) 
			regs[i]= regV.elementAt(i);
		
		return regs;
	}

	public DirectedRegion[] getRegionEventCDS() {
		
		HashMap<Integer, Integer> hashStart= new HashMap<Integer, Integer>(), hashEnd= new HashMap<Integer, Integer>();
		
		for (int i = 0; i < getTranscripts().length; i++) {
			for (int j = 0; j < getTranscripts()[i].length; j++) {
				Transcript trpt= getTranscripts()[i][j];
				if (!trpt.isCoding())
					continue;
				Translation tln= trpt.getTranslations()[0];
				hashStart.put(new Integer(tln.get5PrimeEdge()), null);
				hashEnd.put(new Integer(tln.get3PrimeEdge()), null);
			}
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
			
			DirectedRegion reg= new DirectedRegion(start, end, getGene().getStrand());
			reg.setChromosome(getGene().getChromosome());
			reg.setID("eventCDS");
			v.add(reg);
			
			ptrS+= advS;
			ptrE+= advE+1;
		}
		
		DirectedRegion[] regs= new DirectedRegion[v.size()];
		for (int i = 0; i < regs.length; i++) 
			regs[i]= v.elementAt(i);
		
		return regs;
	} 

	
	public String toStringLocalization() {
		String stringLoc= null;
		if (getAttribute(ATTRIBUTE_TAG_LOCALIZATION) == null) {
			
				byte[] regCodes= null; // 1= 5UTR, 2= 5UTR+CDS, 3= CDS, 4= CDS+3UTR, 5= 3UTR, 6= 5UTR+CDS+3UTR
				DirectedRegion[] regs= null;
				if (getGene().getLocalizationRef()!= null&& getGene().getLocalizationRef().length> 0) {
					regs= getGene().getLocalizationRef();
					regCodes= getGene().getLocalization(getRegionEvent(), regs);
				} else {
					regs= getRegionEventCDS();
					regCodes= getGene().getLocalization(getRegionEvent(), regs);
					byte code= 0;
					if (regCodes!= null&& regCodes.length> 0) {
						if (regCodes.length== 1) 
							code= regCodes[0]; 
						else {
							boolean utr5= false, cds= false, utr3= false; 
							if (regCodes!= null|| regCodes.length> 0) { 
								for (int i = 0; i < regCodes.length; i++) {
									if (regCodes[i]== 1|| regCodes[i]== 2|| regCodes[i]== 6)
										utr5= true;
									if (regCodes[i]== 2|| regCodes[i]== 3|| regCodes[i]== 4|| regCodes[i]== 6)
										cds= true;
									if (regCodes[i]>= 4)
										utr3= true;
								}
								code= utr5?(cds?(utr3?(byte)6:(byte)2):(byte)1):(cds?(utr3?(byte)4:(byte)3):(utr3?(byte)5:(byte)0));
//								if (code!= 0&& regCodes.length> 1)
//									System.currentTimeMillis();
								if (code== 0)
									regCodes= new byte[0];
								else
									regCodes= new byte[] {code};
							}
						}
					}
				}
				if (regs== null|| regs.length== 0|| regCodes== null|| regCodes.length== 0)
					stringLoc= "N/A";
				else {
					StringBuffer sb= new StringBuffer(Gene.LOCALIZATION_STRING[regCodes[0]]);
					sb.append("_");
					sb.append(regs[0].getID());
					for (int i = 1; i < regCodes.length; i++) {
						sb.append(",");
						sb.append(Gene.LOCALIZATION_STRING[regCodes[i]]);
						sb.append("_");
						sb.append(getGene().getLocalizationRef()[i].getID());
					}
					stringLoc= sb.toString();
					addAttribute(ATTRIBUTE_TAG_LOCALIZATION, stringLoc);
				}
				
		}
	
		return stringLoc;
	}

	public byte getType_100113() {
	
			if (type== TYPE_NOT_INITED) {
				
				type= TYPE_UNDEFINED;
				byte[] stypes= getSources();
				int[] min= new int[trpts.length], max= new int[trpts.length];
				
				// min/max transcript extension per variant
				for (int h = 0; h < trpts.length; h++) {
					min[h]= Integer.MAX_VALUE;
					max[h]= Integer.MIN_VALUE;
					for (int m = 0; m < trpts[h].length; m++) {
						min[h]= Math.min(min[h], trpts[h][m].get5PrimeEdge());
						max[h]= Math.max(max[h], trpts[h][m].get3PrimeEdge());
					}
				}
	
	//			 another loop over all transcripts needed! for soft starts..			
				for (int i = 0; type!= TYPE_AS_EVENT&& i < spliceChains.length; i++) {
					for (int j = 0; type!= TYPE_AS_EVENT&& j < spliceChains[i].length; j++) {
	//					if (!spliceChains[i][j].isSpliceSite()
	//							|| spliceChains[i][j].getPos()== min[i]			// truncated start/end splice sites cannot decide
	//							|| spliceChains[i][j].getPos()== max[i]) {
	////						if (type< TYPE_VS_EVENT)
	////							type= TYPE_VS_EVENT;
	//						continue;
	//					}
						
						for (int h = 0; h < trpts.length; h++) {
							if (h== i)
								continue;
							
							// not >=, <= ... truncated trpts dont contain the respective donor/acceptor
							if (spliceChains[i][j].getPos()> min[h]&& 
									spliceChains[i][j].getPos()< max[h]) {	 
								if (spliceChains[i][j].isSpliceSite()) {
									type= TYPE_AS_EVENT;	// as soon as AS w.r.t. ONE other spliceform, stop
									break; // highest in hierachy
								} else 
									type= TYPE_DS_EVENT;
								// [i] has a SS outsite other trpt [h], check unsure level
							} else if(type!= TYPE_DS_EVENT&&
									((src.getPos()>= min[h]&& src.getPos()<= max[h])
									|| (snk.getPos()>= min[h]&& snk.getPos()<= max[h]))) {
								type= TYPE_VS_EVENT;	// lowest in hierachy
							}
	//						else 
								// old nomenclature
	//							if (stypes[i]> stypes[h]|| // EST cannot extend other EST
	//								(stypes[i]== stypes[h]&& stypes[i]< Transcript.ID_SRC_EST)) { // ..nor each other .. 	
	//							if (src.isSpliceSite()|| snk.isSpliceSite())
	//								type= TYPE_DS_EVENT;
								// more variable events..
								// e.g., 1-2^3-4^5-7],6[8^9-10]
	//							if ((min[h]>= min[i]&& min[h]<= max[i])
	//									|| (min[i]>= min[h]&& min[i]<= max[h]))	// have to overlap								
	//						}
	
						}
					}
				}
			}
			
			return type;
		}

	public static boolean isOutputFlankMode() {
		return outputFlankMode;
	}

	public static void setOutputFlankMode(boolean outputFlankMode) {
		ASEvent.outputFlankMode = outputFlankMode;
	}

	String[] cdsLoc= null;
	public String cdsImpact, cdsVariants;
	
	public void setCDSloc(int ttpos, String loc) {
		if (cdsLoc== null) {
			cdsLoc= new String[trpts.length];
		}
		cdsLoc[ttpos]= loc;
	}
	
	public void set53Valid(int ttpos, boolean valid) {
		if (cdsValid== null) {
			cdsValid= new boolean[trpts.length];
		}
		cdsValid[ttpos]= valid;
	}
	public void setFrame3(int ttpos, int frame) {
		if (frame3== null) {
			frame3= new int[trpts.length];
		}
		frame3[ttpos]= frame;
	}
	
	public void setFrame5(int ttpos, int frame) {
		if (frame5== null) {
			frame5= new int[trpts.length];
		}
		frame5[ttpos]= frame;
	}
	
	/**
	 * @deprecated does not consider cds frames
	 * @param ttpos
	 */
	public void setCDSloc(int ttpos) {
		
		Translation refCDS= null;
		for (int i = 0; i < trpts[ttpos].length; i++) {
			if (trpts[ttpos][i].getTranslations()!= null) {
				refCDS= trpts[ttpos][i].getTranslations()[0];
				break;
			}
		}
		
		if (refCDS== null) {
			setCDSloc(ttpos, "NA");
			return;
		}
		
		int cdsStart= refCDS.get5PrimeEdge(), cdsEnd= refCDS.get3PrimeEdge(),
			evStart= src.getPos(), evEnd= snk.getPos();
		if (evStart< cdsStart) {
			if (evEnd> cdsEnd)
				setCDSloc(ttpos, "5UTR-CDS-3UTR");
			else if (evEnd>= cdsStart)
				setCDSloc(ttpos, "5UTR-CDS");
			else
				setCDSloc(ttpos, "5UTR");
		} else if (evStart<= cdsEnd) {
			if (evEnd> cdsEnd)
				setCDSloc(ttpos, "CDS-3UTR");
			else
				setCDSloc(ttpos, "CDS");
		} else // evStart > cdsEnd
			setCDSloc(ttpos, "3UTR");
	}
	
	/**
	 * @deprecated incomplete
	 * @param tpos
	 */
	public void setCDSloc(int[] tpos) {
		
		HashMap[] mapLocVar= new HashMap[trpts.length];
		int evStart= src.getPos(), evEnd= snk.getPos();
		Translation refCDS= null;

		// get localizations
		for (int x = 0; x < tpos.length; x++) {
			if (tpos[x]< 0)
				break;
			int p= tpos[x];
			HashMap map= mapLocVar[p];
			if (map== null) {
				map= new HashMap<String, Object>(trpts[p].length);
				mapLocVar[p]= map;
			}
			// check all CDSs	
			for (int i = 0; i < trpts[p].length; i++) {
				if (trpts[p][i].getTranslations()!= null) {
					refCDS= trpts[p][i].getTranslations()[0];
					int cdsStart= refCDS.get5PrimeEdge(), cdsEnd= refCDS.get3PrimeEdge();
					if (evStart< cdsStart) {
						if (evEnd> cdsEnd)
							map.put("5UTR-CDS-3UTR", null);
						else if (evEnd>= cdsStart)
							map.put("5UTR-CDS", null);
						else
							map.put("5UTR", null);
					} else if (evStart<= cdsEnd) { // event starts inside CDS
						int frame= refCDS.getFrameAtPosition(evStart);
						if (evEnd> cdsEnd)
							map.put("CDS"+frame+"-3UTR", null);
						else
							map.put("CDS"+frame, null);
					} else // evStart > cdsEnd
						map.put("3UTR", null);
				}
			}
		}
		
		// global marriage
		HashMap<String, Integer> globalMap= new HashMap<String, Integer>();
		Iterator<String> iter= null;
		String id= null;
		for (int i = 0; i < mapLocVar.length; i++) {
			if (mapLocVar[i]== null)
				continue;
			iter= mapLocVar[i].keySet().iterator();
			while (iter.hasNext()) {
				id= iter.next();
				if (id.startsWith("CDS"))
					id= id.substring(0, 4);	// CDS0.. CDS1.. CDS2..
				if(globalMap.containsKey(id))
					globalMap.put(id, globalMap.get(id)+ 1);
				else
					globalMap.put(id, 1);
			}
		}
		iter= globalMap.keySet().iterator();
		while(iter.hasNext()) {
			id= iter.next();
			//if (id.startsWith("CDS")&& globalMap.get(id))
			// continue here...
		}
		
	}

	public String getCDSloc(int i) {
		if (cdsLoc== null)
			return null;
		return cdsLoc[i];
	}

	public byte getType() {
	
			if (type== TYPE_NOT_INITED) {
				
				type= TYPE_UNDEFINED;
				int[] min= new int[trpts.length], max= new int[trpts.length];
				
				// min/max transcript extension per variant
				for (int h = 0; h < trpts.length; h++) {
					min[h]= Integer.MAX_VALUE;
					max[h]= Integer.MIN_VALUE;
					for (int m = 0; m < trpts[h].length; m++) {
						min[h]= Math.min(min[h], trpts[h][m].get5PrimeEdge());
						max[h]= Math.max(max[h], trpts[h][m].get3PrimeEdge());
					}
				}
	
	//			 another loop over all transcripts needed! for soft starts..			
				for (int i = 0; type!= TYPE_AS_EVENT&& i < spliceChains.length; i++) {
					for (int j = 0; type!= TYPE_AS_EVENT&& j < spliceChains[i].length; j++) {
	//					if (!spliceChains[i][j].isSpliceSite()
	//							|| spliceChains[i][j].getPos()== min[i]			// truncated start/end splice sites cannot decide
	//							|| spliceChains[i][j].getPos()== max[i]) {
	////						if (type< TYPE_VS_EVENT)
	////							type= TYPE_VS_EVENT;
	//						continue;
	//					}
						
						for (int h = 0; h < trpts.length; h++) {
							if (h== i)
								continue;
							
							// not >=, <= ... truncated trpts dont contain the respective donor/acceptor
							if (spliceChains[i][j].isSpliceSite()) {
								if (spliceChains[i][j].getPos()> min[h]&& 
										spliceChains[i][j].getPos()< max[h]) {	 
									type= TYPE_AS_EVENT;	// as soon as AS w.r.t. ONE other spliceform, stop
									break; // highest in hierachy
								} else if ((min[h]>= min[i]&& min[h]<= max[i])
										|| (min[i]>= min[h]&& min[i]<= max[h]))	// tx have to overlap
									type= TYPE_DS_EVENT;
								// [i] has a SS outsite other trpt [h], check unsure level
							} else if(type!= TYPE_DS_EVENT&&
									((min[h]>= min[i]&& min[h]<= max[i])
									|| (min[i]>= min[h]&& min[i]<= max[h]))) {	// tx have to overlap
	//								((src.getPos()>= min[h]&& src.getPos()<= max[h])
	//								|| (snk.getPos()>= min[h]&& snk.getPos()<= max[h]))) {
								type= TYPE_VS_EVENT;	// lowest in hierachy
							}
	//						else 
								// old nomenclature
	//							if (stypes[i]> stypes[h]|| // EST cannot extend other EST
	//								(stypes[i]== stypes[h]&& stypes[i]< Transcript.ID_SRC_EST)) { // ..nor each other .. 	
	//							if (src.isSpliceSite()|| snk.isSpliceSite())
	//								type= TYPE_DS_EVENT;
								// more variable events..
								// e.g., 1-2^3-4^5-7],6[8^9-10]
	//							if ((min[h]>= min[i]&& min[h]<= max[i])
	//									|| (min[i]>= min[h]&& min[i]<= max[h]))	// have to overlap								
	//						}
	
						}
					}
				}
			}
			
			return type;
		}
}
