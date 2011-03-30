/*
 * Created on Mar 2, 2005
 *
 * TODO To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package fbi.genome.model.gff;

import fbi.genome.model.bed.BEDobject;
import fbi.genome.model.AbstractSite;
import fbi.genome.model.DirectedRegion;
import fbi.genome.model.Exon;
import fbi.genome.model.Gene;
import fbi.genome.model.SpliceSite;
import fbi.genome.model.Transcript;
import fbi.genome.model.commons.MyArrays;
import fbi.genome.model.constants.Constants;

import java.util.Comparator;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.Vector;

/**
 * see http://genes.cs.wustl.edu/GTF2.html for description
 * @author msammeth
 *
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
public class GFFObject {
	
	public static String GTF_ATTRIBUTE_EXON_ALTERNATIVE= "alternative";
	public static String[] GTF_ATTRIBUTE_EXON_ALTERNATIVE_VALUES= 
		new String[] {"none", "left", "right", "left,right", "skipped", "skipped,left", "skipped,right", "skipped,left,right"};
	public static byte GTF_ATTRIBUTE_EXON_ALTERNATIVE_NOT= 0,
		GTF_ATTRIBUTE_EXON_ALTERNATIVE_LEFT= 1,
		GTF_ATTRIBUTE_EXON_ALTERNATIVE_RIGHT= 2,
		GTF_ATTRIBUTE_EXON_ALTERNATIVE_BOTH= 3,
		GTF_ATTRIBUTE_EXON_ALTERNATIVE_SKIPPED= 4,
		GTF_ATTRIBUTE_EXON_ALTERNATIVE_SKIPLEFT= 5,
		GTF_ATTRIBUTE_EXON_ALTERNATIVE_SKIPRIGHT= 6,
		GTF_ATTRIBUTE_EXON_ALTERNATIVE_SKIPBOTH= 7;
	
	public final static char STRAND_SYMBOL_POS= '+', STRAND_SYMBOL_NEG= '-', SYMBOL_NOT_INITED= '.';
	public final static String ATTRIBUTE_TAG_LOCALIZATION= "localization";
	public final static String ATTRIBUTE_TAG_LOCAL_REF= "local_ref";
	public final static String FEATURE_BP= "branch_point", FEATURE_BP_U12= "U12_branch_point", FEATURE_PP= "pp_tract", FEATURE_PP_U12= "U12_pp_tract",FEATURE_SS= "splice_site", FEATURE_SS_U12GTAG= "U12gtag_splice_site", FEATURE_SS_U12CTAC= "U12ctac_splice_site";
	public final static String FEATURE_ASEVENT= "asEvent";
	public final static String FEATURE_TAG= "tag";
	public final static byte STRAND_UNDEFINED= 0, STRAND_POS=1, STRAND_NEG= -1;
	public final static int POSITION_UNDEFINED= 0;
	
	public final static String ATTRIBUTE_EVENT_STRUCTURE= "structure";

	/**
	 * For hashing equal string objects to non-redundant objects.
	 */
	public static Hashtable<String,Object> redundancyHash= null;
	
	public static String getTranscriptID(String s) {
		int p= s.indexOf(TRANSCRIPT_ID_TAG)+ TRANSCRIPT_ID_TAG.length();
		if (p< 0)
			return null;
		while (p< s.length()&& Character.isWhitespace(s.charAt(p++)));
		int p2= p;
		while (p2< s.length()&& s.charAt(p2++)!= '\"');
		return s.substring(p, p2- 1);
	}
	
	public static String getField(int x, String s) {
		int p= 0; --x;
		for (int cnt = 0; p< s.length()&& cnt < x; cnt++) {
			//while (p< s.length()&& s.charAt(p)== Constants.TAB.charAt(0));
			while (p< s.length()&& s.charAt(p++)!= Constants.TAB.charAt(0));
			--p;
			while (p< s.length()&& s.charAt(p++)== Constants.TAB.charAt(0));
			--p;
		}
		int p2= p;
		while (p2< s.length()&& s.charAt(p2++)!= Constants.TAB.charAt(0));
		--p2;
		return s.substring(p, p2);
	}
	
	
	public static GFFObject[] fromBed(BEDobject[] beds) {
		Vector<GFFObject> v= new Vector<GFFObject>();
		for (int i = 0; i < beds.length; i++) {
			GFFObject[] obs= fromBed(beds[i]);
			for (int j = 0; j < obs.length; j++) 
				v.add(obs[j]);
		}
		
		GFFObject[] obs= new GFFObject[v.size()];
		for (int i = 0; i < obs.length; i++) 
			obs[i]= v.elementAt(i);
		
		return obs;
	}
	
	public static GFFObject[] fromBed(BEDobject bed) {
		Vector<GFFObject> v= new Vector<GFFObject>();

		if (bed.getBlockCount()== 0) {	// for clowns.. make only one block
			System.out.println("WARNING: improper bed object converted single region:");
			System.out.println(bed.toString());
			GFFObject o= new GFFObject();
			o.setSeqname(bed.getChrom().toString());
			o.setFeature(bed.getName().toString());
			o.setStrand((byte) bed.getStrand());
			o.setStart(bed.getStart()+1);	// 0-based bed
			o.setEnd(bed.getEnd()+1);
			if (bed.getScore()!= Integer.MIN_VALUE)
				o.setScore(bed.getScore());
//			if (bed.getAttributes()!= null)
//				o.attributes= (HashMap<String, String>) bed.getAttributes().clone();
			v.add(o);
		} else {
			for (int i = 0; i < bed.getBlockCount(); i++) {
				GFFObject o= new GFFObject();
				o.setSeqname(bed.getChrom().toString());
				o.setFeature(bed.getName().toString());
				o.setStrand((byte) bed.getStrand());
				int start= bed.getStart()+1+bed.getBlockStart(i);	//+1 for bed is 0-based
//				System.err.println("Check neg. strand bed->gtf conversion.");	// was not good for psl
				if (bed.getStrand()< 0)
					start= -(Math.abs(bed.getEnd())+ bed.getBlockStart(i)+ 1); 
				o.setStart(start);
				int end= bed.getStart()+bed.getBlockStart(i)+bed.getBlockSize(i); 	// +1-1
//				System.err.println("Check neg. strand bed->gtf conversion.");	// was not good for psl
				if (bed.getStrand()< 0)
					end= -(Math.abs(start)+ bed.getBlockSize(i)- 1);	// already 1-based
				o.setEnd(end);
				if (bed.getScore()!= -1)
					o.setScore(bed.getScore());
//				if (bed.getAttributes()!= null&& bed.getAttributes().size()> 0)
//					o.attributes= (HashMap<String, String>) bed.getAttributes().clone();
				v.add(o);
			}
		}
		
		GFFObject[] obs= new GFFObject[v.size()];
		for (int i = 0; i < obs.length; i++) 
			obs[i]= v.elementAt(i);
		
		return obs;
	}
	
	public static class PositionComparator implements Comparator {
		public int compare(Object o1, Object o2) {
			GFFObject gtf1= (GFFObject) o1;
			GFFObject gtf2= (GFFObject) o2;
			
			int res= gtf1.getSeqname().compareTo(gtf2.getSeqname());
			if (res!= 0)
				return res;
			
			if (gtf1.getStrand()< gtf2.getStrand())
				return -1;
			else if (gtf1.getStrand()> gtf2.getStrand())
				return 1;
			
			if (gtf1.getStart()< gtf2.getStart())
				return -1;
			else if (gtf1.getStart()> gtf2.getStart())
				return 1;
			
			if (gtf1.getEnd()< gtf2.getEnd())
				return -1;
			else if (gtf1.getEnd()> gtf2.getEnd())
				return 1;
			
			
			return 0;
		}
	}
	
	static PositionComparator defaultPositionComparator= new PositionComparator();
	
	public static String ID_ATTRIBUTE_SEQUENCE= "seq";
	
	boolean gff= false;
	
	@Override
	public String toString() {
		
		String x= 
			seqname+ "\t"+
			source+ "\t"+ 
			feature+ "\t"+ 
			start+ "\t"+ 
			end+ "\t";
		x+= (new Float(score).isNaN())?".":Float.toString(score); // == Float.NaN fails..
		x+= "\t";
		if (strand== 1)
			x+="+";
		else if (strand== -1)
			x+="-";
		else 
			x+=".";
		x+= "\t";
		x+= (frame> 0)?Integer.toString(frame):".";
		
		if (attributes!= null) {
			x+= "\t";
			if (attributes.get(GFFObject.TRANSCRIPT_ID_TAG)!= null)
				x+= GFFObject.TRANSCRIPT_ID_TAG+" \""+attributes.get(GFFObject.TRANSCRIPT_ID_TAG)+"\";";
			
			Iterator iter= attributes.keySet().iterator();
			while (iter.hasNext()) {
				Object k= iter.next();
				if (k.equals(GFFObject.TRANSCRIPT_ID_TAG))
					continue;
				x+= " "+ k+ " \""+ attributes.get(k)+"\";";
			}
		}
		if (comments!= null)
			x+= "\t"+ comments;
		return x;	
	}

	public static final String SOURCE_RefSeq= "refGene";
	public static final String SOURCE_MRNA= "mrna";
	public static final String SOURCE_EST= "Est";
	public static final String SOURCE_MIXED= "Mixed";
	public static final String SOURCE_UNDEFINED= "Undefined";
	

	public final static String[] FEATURE_VALID= {"ATG", "CDS", "start_codon", "stop_codon", "exon", "intron", "splice_site", "gene", "mRNA", "5UTR", "3UTR", "CDS"};
	public final static String GENE_ID_TAG= "gene_id", LOCUS_ID_TAG= "locus_id";	
	public final static String TRANSCRIPT_ID_TAG= "transcript_id";
	public final static String EXON_ID_TAG= "exon_id";
	public final static String GENE_ALIAS_TAG= "gene_alias";
	public final static String INTRON_ID_TAG= "intron_id";
	public final static String EXON_FEATURE_TAG= "exon";
	public final static String CDS_FEATURE_TAG= "CDS";
	/**
	 * Parses a string that represents a strand assignment, either in the form 
	 * <code>+/-</code> or \"1\"/\"-1\".
	 * @param strandStr the string representation of the strand
	 * @return the byte representation of the string
	 */
	public static byte parseStrand(CharSequence strandStr) {
		
		byte strand= STRAND_UNDEFINED;
		//strandStr= strandStr.trim();
		if (Character.isDigit(strandStr.charAt(strandStr.length()- 1))) {	// "-1"
			if (strandStr.length()== 1&& strandStr.charAt(0)== '1')
				return 1;
			else if (strandStr.length()== 2&& strandStr.charAt(0)== '-'&& 
					strandStr.charAt(1)== '1')
				return -1;
//			try {
//				strand= Byte.parseByte(strandStr);
//			} catch (Exception e) {
//				e.printStackTrace();
//			}
		} else {
			char c= strandStr.charAt(0);
			if (c== STRAND_SYMBOL_POS)
				strand= 1;
			else if (c== STRAND_SYMBOL_NEG)
				strand= -1;
		}
		
		return strand;
	}
	
	public static char getStrandSymbol(int strand) {
		if (strand==1)
			return '+';
		if (strand==-1)
			return '-';
		return '.';
	}
	
	public static void parseGTF(String line, GFFObject obj) {
		String[] tokens= line.split("\t");
		if (tokens.length< 8) 
			System.out.println("WARNING: invalid GTF line "+line);

		obj.setStart( Integer.parseInt(tokens[3]));
		obj.setEnd(Integer.parseInt(tokens[4]));
		try {
			obj.setStrand(tokens[6]);
		} catch (Exception e) {
			System.out.println(e.getMessage());
		}
		obj.setFeature(tokens[2]);
		obj.setSource(tokens[1]);
		obj.setSeqname(tokens[0]);
		obj.setScore(tokens[5]);
		try {
			obj.setFrame(tokens[7]);
		} catch (Exception e) {
			System.out.println(e.getMessage());
		}
		String h= line.substring(line.indexOf(tokens[8]), line.length()).trim();	// attributes, comments
		String[] attrTokens= h.split(";");
		for (int i = 0; i < attrTokens.length; i++) {
			h= attrTokens[i].trim();
			h= h.replaceAll("\\s+", " ");
			int sep= Math.max(h.indexOf(' '), h.indexOf('='));	// = in ASD gff
			if (sep < 0) 						// comments
				break;
			if (sep>= 0) {
				String id= h.substring(0, sep);
				String val= h.substring(sep+ 1, h.length());
				if (val.charAt(0)== '\"')
					val= val.substring(1, val.length()- 1);
				obj.addAttribute(id, val);
			}
		}
		
	}
	
	public static GFFObject createGTFObject(AbstractSite site) {
	//		<seqname> <source> <feature> <start> <end> <score> <strand> <frame> [attributes] [comments] 
			GFFObject gtf= new GFFObject();
			gtf.setSeqname(site.getGene().getChromosome());
			gtf.setStart(Math.abs(site.getPos()));
			gtf.setEnd(Math.abs(site.getPos()));
			try {
				gtf.setFeature("site");
				gtf.setStrand(site.getGene().getStrand());
			} catch (Exception e) {
				e.printStackTrace();
			}
			return gtf;
		}

	public static GFFObject[] createGTFObjects(AbstractSite site, Transcript trpt) {
	//		<seqname> <source> <feature> <start> <end> <score> <strand> <frame> [attributes] [comments] 
			GFFObject gtf= new GFFObject();
			gtf.setSeqname(trpt.getChromosome());
			gtf.setStart(Math.abs(site.getPos()));
			gtf.setEnd(Math.abs(site.getPos()));
			gtf.addAttribute(TRANSCRIPT_ID_TAG, trpt.getTranscriptID());
			gtf.setSource(trpt.getSource());
			gtf.setScore((float) site.getScore());
			try {
				gtf.setFeature(site.getID());
				gtf.setStrand(trpt.getStrand());
			} catch (Exception e) {
				e.printStackTrace();
			}
			
			Vector gtfV= new Vector();
			gtfV.add(gtf);
			createGTFObjectsAddAttributes(site.getAttributes(), gtf, gtfV, trpt);
			
			return (GFFObject[]) MyArrays.toField(gtfV);
		}

	public static GFFObject[] createGTFObjects(SpliceSite site, Transcript trpt) {
//		<seqname> <source> <feature> <start> <end> <score> <strand> <frame> [attributes] [comments]
		GFFObject gtf= new GFFObject();
		gtf.setSeqname(trpt.getChromosome());
		int start= Math.abs((site.isDonor())?site.getPos():site.getPos()-1); 
		int end= Math.abs((site.isDonor())?site.getPos()+1:site.getPos()); 
		gtf.setStart(start);
		gtf.setEnd(end);
		gtf.setSource(trpt.getSource());
		gtf.setScore((float) site.getScore());
		try {
			gtf.setFeature(site.getID());
			gtf.setStrand(trpt.getStrand());
			gtf.addAttribute(TRANSCRIPT_ID_TAG, trpt.getTranscriptID());
		} catch (Exception e) {
			e.printStackTrace();
		}
		Vector gtfV= new Vector();
		
		createGTFObjectsAddAttributes(site.getAttributes(), gtf, gtfV, trpt);
		gtfV.add(gtf);
		
		return (GFFObject[]) MyArrays.toField(gtfV);
	}
	
	public static void createGTFObjectsAddAttributes(HashMap map, GFFObject baseObj, Vector v, Transcript trpt) {
		if (map== null) 
			return;
		
		Object[] keys= map.keySet().toArray();
		java.util.Arrays.sort(keys);
		for (int i = 0; i < keys.length; i++) {
			if (map.get(keys[i]) instanceof AbstractSite)
				MyArrays.addAll(v, createGTFObjects((AbstractSite) map.get(keys[i]), trpt));
			if (map.get(keys[i]) instanceof DirectedRegion)
				MyArrays.addAll(v, createGTFObjects((DirectedRegion) map.get(keys[i]), trpt));
			else if (map.get(keys[i]) instanceof String&& keys[i] instanceof String)
				baseObj.addAttribute((String) keys[i], (String) map.get(keys[i]));
		}

	}
	
	
	
	public final static String INTRON_STATUS_TAG = "intron_status";

	public static GFFObject createGFFObject(DirectedRegion reg) {
		GFFObject gtf= new GFFObject();
		gtf.setGff(true);
		gtf.setSeqname(reg.getChromosome());
		gtf.setStart(Math.abs(reg.getStart()));
		gtf.setEnd(Math.abs(reg.getEnd()));
		try {
			gtf.setFeature(reg.getID());
			gtf.setStrand(reg.getStrand());
		} catch (Exception e) {
			e.printStackTrace();
		}
		return gtf;
	}
	
	public static GFFObject createGFFObject(DirectedRegion reg, String[] attributes) {
		GFFObject gtf= createGFFObject(reg);
		for (int i = 0; attributes!= null&& i < attributes.length; i++) 
			gtf.addAttribute(Integer.toString(i), attributes[i]);
		return gtf;
	}
	
	public static GFFObject createGFFObject(AbstractSite site) {
//		<seqname> <source> <feature> <start> <end> <score> <strand> <frame> [attributes] [comments] 
		GFFObject gtf= new GFFObject();
		gtf.setGff(true);
		gtf.setSeqname(site.getGene().getChromosome());
		gtf.setStart(Math.abs(site.getPos()));
		gtf.setEnd(Math.abs(site.getPos()));
		try {
			gtf.setFeature("site");
			gtf.setStrand(site.getGene().getStrand());
		} catch (Exception e) {
			e.printStackTrace();
		}
		return gtf;
	}
	
	/**
	 * Uses <code>redundancyHash</code> to filter equal string objects for redundancy. 
	 * @param inString a string to be checked for redundancy
	 * @return consistently the same string instance that equals the input string
	 */
	public static String filterRedundantStrings(String inString) {
		String outString = (String) getRedundancyHash().get(inString);	// memory efficiency: kill redundant string objects 
		if (outString == null) {
			outString = inString;
			GFFObject.getRedundancyHash().put(outString, outString);
		}

		return outString;
	}
	
	public static GFFObject createGFFObject(SpliceSite site, Transcript t) {
		GFFObject gtf= createGFFObject(site);
		
		if (t.getTranscriptID().equals("AK095183"))
			System.currentTimeMillis();
		
		if (site.isCodon()) {			
			int b= t.getExonicPosition(site.getPos()),
				x= t.getGenomicPosition(b),
				y= t.getGenomicPosition(b+ t.getStrand()* (site.isCodonStart()?2:-2));
			int len= t.getExonicLength();
			gtf.setStart(x);
			gtf.setEnd(y);
		}
		gtf.addAttribute(TRANSCRIPT_ID_TAG, t.getTranscriptID());
		gtf.addAttribute(GENE_ID_TAG, t.getGene().getGeneID());
		return gtf;
	}
	
	
	public static GFFObject createGFFObject(SpliceSite site) {
		GFFObject gtf= createGFFObject((AbstractSite) site);
		try {
			if (site.isDonor())
				gtf.setFeature("donor");
			else if (site.isAcceptor())
				gtf.setFeature("acceptor");
			else if (site.isTSS())
				gtf.setFeature("tss");
			else if (site.isTES())
				gtf.setFeature("tes");
			else if (site.isCodonStart())
				gtf.setFeature("start_codon");
			else if (site.isCodonStop())
				gtf.setFeature("stop_codon");

		} catch (Exception e) {
			e.printStackTrace();
		}
		
		if (site.isSpliceSite()) {
		
			if (site.isConstitutive()) 
				gtf.addAttribute("modality", "constitutive");
			else
				gtf.addAttribute("modality", "alternative");
			
			Transcript[] t= site.getTranscripts();
			StringBuffer parents= new StringBuffer();
			String source= null;
			for (int i = 0; i < t.length; i++) { 
				parents.append(t[i].getTranscriptID());
				parents.append(',');
				
				if (source== null)
					source= t[i].getSource();
				else if (!source.equals(t[i].getSource())) {
					String srcUp= source.toUpperCase();
					String nSrcUp= t[i].getSource().toUpperCase();
					if ((srcUp.contains("MRNA")&& nSrcUp.contains("GENE"))|| 
							(srcUp.contains("EST")&& (!nSrcUp.contains("EST"))))
						source= t[i].getSource();
				}
			}
			gtf.setSource(source);
			parents.deleteCharAt(parents.length()-1);
			gtf.addAttribute("parent", parents.toString());
			gtf.addAttribute("location", site.getGenicLocation());
		
		}
		
		return gtf;
	}

	public static GFFObject[] createGTFObjects(DirectedRegion reg) {
		GFFObject obj= new GFFObject();
		obj.setSeqname(reg.getChromosome());
		obj.setStart(reg.getStart());
		obj.setEnd(reg.getEnd());
		try {
			obj.setStrand(reg.getStrand());
			obj.setFeature(reg.getID());
		} catch (Exception e) {
			e.printStackTrace();
		}
		if (reg.getAttributes()!= null) {
			Iterator<String> iter= reg.getAttributes().keySet().iterator();
			while (iter.hasNext()) {
				String k= iter.next();
				if (k.equals(GFFObject.TRANSCRIPT_ID_TAG)||
						k.equals(GFFObject.GENE_ID_TAG))
					obj.addAttribute(k, (String) reg.getAttribute(k));
				if (k.equals(Transcript.GTF_ATTRIBUTE_SOURCE))
					obj.setSource((String) reg.getAttribute(k)); 
			}
		}
	
		Vector addObjV= new Vector();
		if (reg instanceof Gene) {
			obj.addAttribute(GENE_ID_TAG, ((Gene) reg).getGeneID());
		} else if (reg instanceof Transcript) {
			obj.addAttribute(GENE_ID_TAG, ((Transcript) reg).getGene().getGeneID());
			obj.addAttribute(TRANSCRIPT_ID_TAG, ((Transcript) reg).getTranscriptID());
			obj.setSource(((Transcript) reg).getSource());
		} else if (reg instanceof Exon) {
			Exon exon= ((Exon) reg);
			
			for (int i = 0; i < exon.getTranscripts().length; i++) {
				GFFObject ex= (GFFObject) obj.clone();
				ex.addAttribute(GENE_ID_TAG, exon.getGene().getGeneID());
				ex.addAttribute(EXON_ID_TAG, exon.getExonID());
				ex.setSource(exon.getTranscripts()[i].getSource());
				
				ex.addAttribute(TRANSCRIPT_ID_TAG, exon.getTranscripts()[i].getTranscriptID());
				addObjV.add(ex);
				
				GFFObject cds= (GFFObject) ex.clone();
				cds.setFeature("CDS");
				cds.setStart(exon.getStartCDS());
				cds.setEnd(exon.getEndCDS());
				if (exon.getStartCDS()!= 0&& exon.getEndCDS()!= 0)
					addObjV.add(cds);
			}
	
		} 
		
		if (addObjV.size()== 0)
			return new GFFObject[] {obj};
		else
			return (GFFObject[]) MyArrays.toField(addObjV);
	}

	public static GFFObject[] createGTFObjects(DirectedRegion reg, Transcript trpt) {
		GFFObject obj= new GFFObject();
		obj.setSeqname(trpt.getChromosome());
		obj.setStart(reg.getStart());
		obj.setEnd(reg.getEnd());
		obj.addAttribute(TRANSCRIPT_ID_TAG, trpt.getTranscriptID());
		obj.setSource(trpt.getSource());
		obj.setScore((float) reg.getScore());
		try {
			obj.setStrand(reg.getStrand());
			obj.setFeature(reg.getID());
		} catch (Exception e) {
			e.printStackTrace();
		}
		Vector gtfV= new Vector();
		gtfV.add(obj);
		
		createGTFObjectsAddAttributes(reg.getAttributes(), obj, gtfV, trpt);
		return (GFFObject[]) MyArrays.toField(gtfV);
	}

	public static GFFObject createGTFObject(Exon exon, Transcript trpt) {
		return createGTFObjects(exon, trpt)[0];
	}
	
	public static GFFObject[] createGTFObjects(Exon exon, Transcript trpt) {
		return createGTFObjects(exon, trpt, false);
	}
	
	public static GFFObject createGFFObject(Gene g) {
		GFFObject obj= new GFFObject();
		
		obj.setSeqname(g.getChromosome());
		obj.setSource(Transcript.getSource(g.getTranscripts()));
		obj.setStart(g.getStart());
		obj.setEnd(g.getEnd());
		try {
			obj.setStrand((byte) g.getStrand());
			obj.setFeature(Gene.GFF_FEATURE_LOCUS);
		} catch (Exception e) {
			e.printStackTrace();
		}
		StringBuffer sb= new StringBuffer();
		for (int i = 0; i < g.getTranscripts().length; i++) {
			sb.append(g.getTranscripts()[i].getTranscriptID());
			if (i< g.getTranscripts().length-1)
				sb.append(",");
		}
		obj.addAttribute(TRANSCRIPT_ID_TAG, sb.toString());
		
		obj.addAttribute(GENE_ID_TAG, g.getGeneID());

		return obj;
	}
	
	public static GFFObject createGFFObject(Transcript t) {
		GFFObject obj= new GFFObject();
		obj.setSeqname(t.getChromosome());
		obj.setSource(t.getSource());
		obj.setStart(t.getStart());
		obj.setEnd(t.getEnd());
		try {
			obj.setStrand((byte) t.getStrand());
			obj.setFeature(Transcript.GFF_FEATURE_TRANSCRIPT);
		} catch (Exception e) {
			e.printStackTrace();
		}
		obj.addAttribute(TRANSCRIPT_ID_TAG, t.getTranscriptID());
		
		obj.addAttribute(GENE_ID_TAG, t.getGene().getGeneID());
		return obj;
	}
	
	public static GFFObject[] createGTFObjects(Exon exon, Transcript trpt, boolean attrAlternative) {
		GFFObject obj= new GFFObject();
		obj.setSeqname(exon.getChromosome());
		obj.setStart(exon.getStart());
		obj.setEnd(exon.getEnd());
		try {
			obj.setStrand((byte) trpt.getStrand());
			obj.setFeature(exon.getID());
		} catch (Exception e) {
			e.printStackTrace();
		}

		Vector addObjV= new Vector();
		GFFObject ex= (GFFObject) obj.clone();
		ex.addAttribute(TRANSCRIPT_ID_TAG, trpt.getTranscriptID());
		if (trpt.getGene()!= null)
			ex.addAttribute(GENE_ID_TAG, trpt.getGene().getGeneID());
		if (exon.getExonID()!= null)
			ex.addAttribute(EXON_ID_TAG, exon.getExonID());
		ex.setSource(trpt.getSource());
			// revise this
		if (attrAlternative) {
			boolean skipped= false, varLeft= false, varRight= false;
			for (int i = 0; (skipped&& varLeft&& varRight)!= true&& i < trpt.getGene().getTranscripts().length; i++) {
				Transcript t= trpt.getGene().getTranscripts()[i];
				if (t== trpt|| !t.overlaps(exon))
					continue;
				for (int j = 0; j < t.getExons().length; j++) {	// TODO binSearch
					if (t.getExons()[j]== exon)
						break;
					if (t.getExons()[j].get5PrimeEdge()> exon.get3PrimeEdge()) {
						skipped= true;
						break;
					}
					if (t.getExons()[j].overlaps(exon)) {
						if (t.getExons()[j].get5PrimeEdge()!= exon.get5PrimeEdge()) {
							varLeft= true;
								// no more..
//							if (t.getExons()[j].get3PrimeEdge()== exon.get3PrimeEdge())	// DEBUG
//								System.out.println("exon ident error");
						} 
						if (t.getExons()[j].get3PrimeEdge()!= exon.get3PrimeEdge()) {
							varRight= true;
//							if (t.getExons()[j].get5PrimeEdge()== exon.get5PrimeEdge())	// DEBUG
//								System.out.println("exon ident error");
						} 
						break;
					} 
					System.currentTimeMillis();
				}
			}
			byte exonStatus= GTF_ATTRIBUTE_EXON_ALTERNATIVE_NOT;
			if (varLeft) 
				exonStatus= GTF_ATTRIBUTE_EXON_ALTERNATIVE_LEFT;
			else if (varRight)
				exonStatus= GTF_ATTRIBUTE_EXON_ALTERNATIVE_RIGHT;
			if (varLeft&& varRight)
				exonStatus= GTF_ATTRIBUTE_EXON_ALTERNATIVE_BOTH;
			if (skipped)
				exonStatus+= 4;
			ex.addAttribute(GTF_ATTRIBUTE_EXON_ALTERNATIVE, GTF_ATTRIBUTE_EXON_ALTERNATIVE_VALUES[exonStatus]);
		}
		addObjV.add(ex);
		
		if (trpt.getTranslations()!= null) {
			GFFObject cds= (GFFObject) ex.clone();
			cds.setFeature("CDS");
			DirectedRegion exonReg= (DirectedRegion) exon.clone();
			exonReg.setStrand(trpt.getStrand());
			DirectedRegion[] cdsRegs= 
				DirectedRegion.intersect(new DirectedRegion[] {exonReg}, new DirectedRegion[] {trpt.getTranslations()[0]});
			if (cdsRegs!= null) {
				cds.setStart(cdsRegs[0].getStart());
				cds.setEnd(cdsRegs[0].getEnd());
				byte frame= (byte) trpt.getTranslations()[0].getFrame(exon);
				if (frame>= 0)
					cds.setFrame(frame);	// needs genomic sequence ??!
				addObjV.add(cds);
			}
		}

		return (GFFObject[]) MyArrays.toField(addObjV);
	}
	
	/**
	 * all attributes set on pattern (ie, not <code>null</code> or 0) are compared against
	 * objects in array. Objects that match all criteria are added to the solution. Attribute
	 * is a minimum score threshold.
	 * @param objs
	 * @param pattern
	 * @return
	 */
	public static GFFObject[] filterGTFObjects(GFFObject[] objs, GFFObject pattern) {
		Vector v= new Vector();
		for (int i = 0; i < objs.length; i++) {
			if (pattern.getSeqname()!= null&& !pattern.getSeqname().equals(objs[i].getSeqname()))
				continue;
			if (pattern.getSource()!= null&& !pattern.getSource().equals(objs[i].getSource()))
				continue;
			if (pattern.getFeature()!= null&& !pattern.getFeature().equals(objs[i].getFeature()))
				continue;
			if (pattern.getStrand()!= 0&& pattern.getStrand()!= objs[i].getStrand())
				continue;
			if (pattern.getStart()!= 0&& pattern.getStart()!= objs[i].getStart())
				continue;
			if (pattern.getEnd()!= 0&& pattern.getEnd()!= objs[i].getEnd())
				continue;
			if (pattern.getStrand()!= 0&& pattern.getStrand()!= objs[i].getStrand())
				continue;
			if (pattern.getScore()!= 0&& pattern.getScore()< objs[i].getScore())
				continue;
			
			HashMap map= pattern.getAttributes();
			Object[] keys= map.keySet().toArray();
			HashMap map2= objs[i].getAttributes();
			int j;
			for (j = 0; keys!= null&& j < keys.length; j++) {
				if (map2== null|| !map.get(keys[j]).equals(map2.get(keys[j])))
					break;
			}
			if (keys!= null&& j< keys.length)
				continue;
			
			v.add(objs[i]);
		}
		
		return (GFFObject[]) MyArrays.toField(v);
	}
	
	public static GFFObject createGFFObject(SpliceSite site, String source) {
		GFFObject gtf= createGFFObject(site);
		gtf.setSource(source);
		return gtf;
	}
	public static GFFObject createGFFObject(DirectedRegion reg, String src, String[] attributes) {
		GFFObject gtf= createGFFObject(reg, attributes);
		gtf.setSource(src);
		return gtf;
	}
	public static GFFObject createGFFObject(SpliceSite site, String[] attributes) {
		GFFObject gtf= createGFFObject(site);
		for (int i = 0; attributes!= null&& i < attributes.length; i++) 
			gtf.addAttribute(Integer.toString(i), attributes[i]);
		return gtf;
	}
	
	public static GFFObject createGFFObject(SpliceSite site, String src, String[] attributes) {
		GFFObject gtf= createGFFObject(site, src);
		for (int i = 0; attributes!= null&& i < attributes.length; i++) 
			gtf.addAttribute(Integer.toString(i), attributes[i]);
		return gtf;
	}
	
	
	//	<seqname> <source> <feature> <start> <end> <score> <strand> <frame> [attributes] [comments] 
	String seqname= null;	// The FPC contig ID from the Golden Path. 
	// The source column should be a unique label indicating where the annotations came from --- typically the name of either a prediction program or a public database.
	String source= null;
//	The following feature types are required: "CDS", "start_codon", "stop_codon". The feature "exon" is optional, since this project will not evaluate predicted splice sites outside of protein coding regions. All other features will be ignored.
//	CDS represents the coding sequence starting with the first translated codon and proceeding to the last translated codon. Unlike Genbank annotation, the stop codon is not included in the CDS for the terminal exon.
	String feature= "";

	//	Integer start and end coordinates of the feature relative to the beginning of the sequence named in <seqname>.  
	// <start> must be less than or equal to <end>. Sequence numbering starts at 1.
	// Values of <start> and <end> that extend outside the reference sequence are technically acceptable, 
	// but they are discouraged for purposes of this project.
	// If the strand is '-', then the first base of the region is value of <end>, 
	// because the corresponding coding region will run from <end> to <start> on the reverse strand.
	int start= -1, end= -1;
	
	//	[attributes]
//	 All four features have the same two mandatory attributes at the end of the record:
//	     * gene_id value;     A globally unique identifier for the genomic source of the transcript
//	     * transcript_id value;     A globally unique identifier for the predicted transcript.
//	 These attributes are designed for handling multiple transcripts from the same genomic region. Any other attributes or comments must appear after these two and will be ignored.
//	 Attributes must end in a semicolon which must then be separated from the start of any subsequent attribute by exactly one space character (NOT a tab character).
//	 Textual attributes should be surrounded by doublequotes.
	HashMap<String,String> attributes= null;	// GTF2, generic attributes (eg, geneAlias, exonID, intronID, intronStatus, ..)

	
	
	
	// The score field will not be used for this project, so you can either provide a meaningful float or replace it by a dot. 
	float score= Float.NaN;
	
//	0 indicates that the first whole codon of the reading frame is located at 5'-most base. 1 means that there is one extra base before the first codon and 2 means that there are two extra bases before the first codon. Note that the frame is not the length of the CDS mod 3.
//	Here are the details excised from the GFF spec. Important: Note comment on reverse strand.
//	'0' indicates that the specified region is in frame, i.e. that its first base corresponds to the first base of a codon. '1' indicates that there is one extra base, i.e. that the second base of the region corresponds to the first base of a codon, and '2' means that the third base of the region is the first base of a codon. If the strand is '-', then the first base of the region is value of <end>, because the corresponding coding region will run from <end> to <start> on the reverse strand.
	byte strand= 0;
	byte frame= -1;
	
	
	String comments= null; // GTF2, not described

	public boolean isExon() {
		return getFeature().equals("exon"); 
	}
	public boolean isIntron() {
		return getFeature().equals("intron"); 
	}
	public boolean isCoding() {
		return (isCDS()|| isStartCodon()|| isStopCodon());
	}
	public boolean isCDS() {
		return getFeature().equals("CDS"); 
	}
	public boolean isStartCodon() {
		return getFeature().equals("start_codon"); 
	}
	public boolean isStopCodon() {
		return getFeature().equals("stop_codon"); 
	}
	public GFFObject(String line) {
		this();
		parseGTF(line, this);
	}
	
	public GFFObject() {
	}
	
	public void addAttribute(String name, String value) {
		
		if (value== null) 
			return;
		
		value= value.trim();
		if (value.startsWith("\"")) 
			value= value.substring(1, value.length()- 1);	// remove quota
		
		getAttributes().put(name, value);
	}
	
	public Object clone() {
		GFFObject obj= new GFFObject();
		if (attributes!= null)
			obj.attributes= (HashMap) attributes.clone(); 
		obj.comments= comments;
		obj.end= end;
		obj.feature= feature;
		obj.frame= frame;
		obj.gff= gff;
		obj.score= score;
		obj.seqname= seqname;
		obj.source= source;
		obj.start= start;
		obj.strand= strand;
		
		return obj;
	}
	
	public boolean equals(Object obj) {
		
		GFFObject anotherGTF;
		try {
			anotherGTF= (GFFObject) obj;
		} catch (ClassCastException e) {
			return false;
		}
		
		if (!anotherGTF.getSeqname().equals(seqname)||
			!anotherGTF.getSource().equals(source)||
			!anotherGTF.getFeature().equals(feature)||
			anotherGTF.getStart()!= start||
			anotherGTF.getEnd()!= end||
			anotherGTF.isStrand()!= isStrand())
			return false;

		// attributes and comments not checked
		
		return true;
	}
	
	public HashMap getAttributes() {
		if (attributes== null)
			attributes= new HashMap<String,String>();
		return attributes;
	}
	
	public String getAttribute(String id) {
		
		if (attributes== null)
			return null;
		
		return (String) attributes.get(id);
	}
	public boolean overlaps(GFFObject anotherObject) {
		if ((!getChromosome().equals(anotherObject.getChromosome()))||
				getStrand()!= anotherObject.getStrand())
				return false;
		if ((getStart()>= anotherObject.getStart()&& getStart()<= anotherObject.getEnd())||
				(anotherObject.getStart()>= getStart()&& anotherObject.getStart()<= getEnd()))
			return true;
		return false;
	}
	
	public String removeAttribute(String id) {
		
		if (attributes== null)
			return null;
		
		return (String) attributes.remove(id);
	}
	public String getExonID() {
		
		if (attributes== null)
			return null;
		
		return getAttribute(EXON_ID_TAG);
	}
	public String getTranscriptID() {
		
		if (attributes== null)
			return null;
		
		return getAttribute(TRANSCRIPT_ID_TAG);
	}
	
	public String getGeneID() {
		
		if (attributes== null)
			return null;
		
		return getAttribute(GENE_ID_TAG);
	}
	
	public String getGeneAlias() {
		
		if (attributes== null)
			return null;
		
		return getAttribute(GENE_ALIAS_TAG);
	}	
	
	public String getIntronID() {
		
		if (attributes== null)
			return null;
		
		return getAttribute(INTRON_ID_TAG);
	}
	
	public String getIntronStatus() {
		
		if (attributes== null)
			return null;
		
		return getAttribute(INTRON_STATUS_TAG);
	}

	
	/**
	 * @return Returns the end.
	 */
	public int getEnd() {
		return end;
	}
	/**
	 * @param end The end to set.
	 */
	public void setEnd(int end) {
		this.end= Math.abs(end);
		if (start>= 0&& start> this.end) {
			int h= start;
			start= this.end;
			this.end= h;
		}
	}
/**
 * @return Returns the feature.
 */
public String getFeature() {
	return feature;
}
/**
 * @param feature The feature to set.
 */
public void setFeature(String feature) {
	
	if (false) {	//!isGff()) {
		//Exception e;
		int i;
		for (i = 0; i < GFFObject.FEATURE_VALID.length; i++) { 
			if (feature.equals(GFFObject.FEATURE_VALID[i]))
				break;
			if (feature.equalsIgnoreCase(GFFObject.FEATURE_VALID[i]))
				System.err.println("check case spelling for "+feature);
		}
		if (i== GFFObject.FEATURE_VALID.length) {
			System.err.println("no valid entry for feature\n\t"+ feature);
			//throw(e);
		}
	}
	this.feature = feature;
}
	/**
	 * @return Returns the score.
	 */
	public float getScore() {
		return score;
	}
	
	public String getScoreString() {
		if (Float.toString(score).equals(Float.toString(Float.NaN)))
			return ".";
		return Float.toString(score);
	}
	/**
	 * @param score The score to set.
	 */
	public void setScore(String scoreStr) throws NumberFormatException {
		
		if (scoreStr.trim().equals(".")) 
			return;
	
		score= Float.parseFloat(scoreStr);
	}

	/**
	 * @param score The score to set.
	 */
	public void setScore(float sc) {
		
		score= sc;
	}
/**
 * @return Returns the seqname.
 */
public String getSeqname() {
	String s= seqname;
	// now always with chr, write own method if you dont want it
	// (better for scaffolds, contigs..
//	if (s.startsWith("chr"))	// not in mart output
//		s= seqname.substring(3);	// "chr..."

	return s;
}
/**
 * @param seqname The seqname to set.
 */
public void setSeqname(String seqname) {
	
//	if (seqname.length()<= 3)
//		seqname= "chr"+ seqname;

//	if (seqname.indexOf('_')>= 0) 	// human known genes
//		seqname= seqname.split("_")[0];
	this.seqname = seqname;
}
	/**
	 * @return Returns the source.
	 */
	public String getSource() {
		return source;
	}
	/**
	 * @param source The source to set.
	 */
	public void setSource(String source) {
		this.source = source;
	}
	/**
	 * @return Returns the start.
	 */
	public int getStart() {
		return start;
	}
	
	public byte getStrand() {
		return strand;
	}
	
	public char getStrandSymbol() {
		if (strand== 1)
			return '+';
		if (strand== -1)
			return '-';
		return '.';
	}
	/**
	 * @param start The start to set.
	 */
	public void setStart(int start) {
		this.start= Math.abs(start);
		if (end>= 0&& this.start> end) {
			int h= this.start;
			this.start= end;
			end= h;
		}
	}
	/**
	 * @return Returns the frame.
	 */
	public int getFrame() {
		return frame;
	}
	
	
	public char getFrameSymbol() {
		if (frame== -1)
			return '.';
		return Integer.toString(frame).charAt(0);
	}
	/**
	 * @param frame The frame to set.
	 */
	public void setFrame(byte frame) {
		if (frame< 0 || frame> 2)
			return;
			//throw new Exception("no valid frame-shift "+frame);
		this.frame = frame;
	}
	
	public void setFrame(String frameStr) throws NumberFormatException {
		
		if (frameStr.trim().equals("."))	// unused
			return;
		else setFrame(Byte.parseByte(frameStr));
	}
/**
 * @return Returns the leadingStrand.
 */
public boolean isStrand() {
	return (strand== 1);
}
/**
 * @param strand The leadingStrand to set.
 */
public void setStrand(String leadingLagging) throws Exception{
	if (leadingLagging.trim().equals("+")|| leadingLagging.trim().equals("1")) {
		strand= 1;
		return;
	} else if (leadingLagging.trim().equals("-")|| leadingLagging.trim().equals("-1")) {
		strand= -1;
		return;
	}
	Exception e= new Exception("no valid mark for orientation! "+leadingLagging);
	throw(e);
}
public void setStrand(byte strandInt) {
	if (strandInt== 1|| strandInt== -1 || strandInt== 0)
		this.strand= strandInt;
	else if (!isGff()) {
		System.err.println("no valid mark for orientation! "+strandInt);
	}
}
	/**
	 * @return Returns the comments.
	 */
	public String getComments() {
		return comments;
	}
	/**
	 * @param comments The comments to set.
	 */
	public void setComments(String comments) {
		this.comments = comments;
	}

	public String getChromosome() {
		String s= getSeqname();
		return s;
		
	}
	public boolean isGff() {
		return gff;
	}
	public void setGff(boolean gff) {
		this.gff = gff;
	}

	public static PositionComparator getDefaultPositionComparator() {
		return defaultPositionComparator;
	}

	public static Hashtable<String, Object> getRedundancyHash() {
		if (redundancyHash == null) {
			redundancyHash = new Hashtable<String, Object>();
		}

		return redundancyHash;
	}
}
