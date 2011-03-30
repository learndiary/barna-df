package fbi.genome.sequencing.rnaseq.simulation;

import java.io.IOException;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Vector;

import fbi.genome.model.Transcript;
import fbi.genome.model.DirectedRegion;
import fbi.genome.model.gff.GFFObject;

public class Fragment implements Comparable {
	
	public static final String GTF_FEATURE= "fragment",
				GTF_ATTRIBUTE_EXONIC_START= "exonic_start",
				GTF_ATTRIBUTE_EXONIC_END= "exonic_end",
				GTF_ATTRIBUTE_EXONIC_LENGTH= "exonic_length",
				GTF_ATTRIBUTE_FRAGMENT_ID= "frag_id";
	
	String ancestorID= null;
	int start= 0, end= 0;
	double weight= 0d, key= 0d;
	public static void main(String[] args) {
		int size= 10000000;
		//memtest
		Vector<Fragment> v= new Vector<Fragment>(size);
		Fragment mother= new Fragment("x",1,1000);
		for (int i = 0; i < size; i++) {
			v.add((Fragment) mother.clone());
		}
		System.out.println("ready.");
		try {
			System.in.read();
		} catch (IOException e) {
			// TODO Auto-generated catch block
		}
	}
	public static LengthComparator defaultLengthComparator= new LengthComparator();
	
	static class LengthComparator implements Comparator<Fragment> {

		//@Override
		public int compare(Fragment o1, Fragment o2) {
			return (o1.length()- o2.length());
		}
		
	}
	
	public Fragment(String newAncestorID, int startPos, int endPos) {
		ancestorID= newAncestorID;
		start= startPos;
		end= endPos;
	}
	public Fragment(String newAncestorID, int startPos, int endPos, int newNb) {
		this(newAncestorID, startPos, endPos);
		nb= newNb;
	}
	
	int nb= 0;
	
	public int length() {
		return end- start+ 1;
	}
	
	public String getAncestorID() {
		return ancestorID;
	}

	public int getStart() {
		return start;
	}

	public int getEnd() {
		return end;
	}

	@Override
	public String toString() {
		
		return "["+start+" "+end+"]"+length();
	}
	
	public String toStringSerialized() {
		
		StringBuffer sb= new StringBuffer(20);
		sb.append(length());
		sb.append("\t");
		sb.append(getAncestorID());
		sb.append("\t");
		sb.append(getStart());
		sb.append("\t");
		sb.append(getEnd());
		
		return sb.toString();
	}
	
	public Fragment(String serialized) {
		int p= serialized.indexOf('\t');
		int p1= serialized.indexOf('\t', p+1);
		ancestorID= serialized.substring(p+1, p1);
		p= p1;
		p1= serialized.indexOf('\t', p+1);
		start= Integer.parseInt(serialized.substring(p+1, p1));
		end= Integer.parseInt(serialized.substring(p1+1, serialized.length()));
		
		// just in case
		nb= 1; key= 0; weight= 0;
	}

	DirectedRegion[] getRegions(Transcript trpt) {
		int genStart= trpt.getGenomicPosition(getStart());
		int genEnd= trpt.getGenomicPosition(getEnd());
		DirectedRegion reg= new DirectedRegion(genStart, genEnd, trpt.getStrand());
		reg.setChromosome(trpt.getChromosome());
		//reg.setSpecies(trpt.getSpecies());
		
		DirectedRegion[] regs= DirectedRegion.intersect(new DirectedRegion[] {reg}, trpt.getExons());
		return regs;
	}
	
	DirectedRegion getGenomicFragRegion(DirectedRegion[] regs) {
		int genStart= 0, genEnd= 0;
		
		int added= 0;
		for (int i = 0; i < regs.length; i++) {
			if (added+ regs[i].getLength()< getStart())
				added+= regs[i].getLength();
			else {
				genStart= regs[i].get5PrimeEdge()+ (getStart()- added);
				break;
			}
		}
		
		added= 0;
		for (int i = 0; i < regs.length; i++) {
			if (added+ regs[i].getLength()< getEnd())
				added+= regs[i].getLength();
			else {
				genEnd= regs[i].get5PrimeEdge()+ (getEnd()- added);
				break;
			}
		}

		DirectedRegion reg= new DirectedRegion(genStart, genEnd, regs[0].getStrand());
		return reg;
	}
	
	DirectedRegion[] getRegions(DirectedRegion[] regs) {
		DirectedRegion reg= getGenomicFragRegion(regs);
		reg.setChromosome(regs[0].getChromosome());
		//reg.setSpecies(trpt.getSpecies());
		
		String id= regs[0].getID();
		regs= DirectedRegion.intersect(new DirectedRegion[] {reg}, regs);
		regs[0].setID(id);
		
		return regs;
	}
	
	public String toStringGTF(Transcript trpt, HashMap<String, String> attributes) {
		
		DirectedRegion[] regs= getRegions(trpt);
//		if (regs.length!= 2)
//			System.currentTimeMillis();
		StringBuffer sb= new StringBuffer(100);
		for (int i = 0; i < regs.length; i++) {
			GFFObject obj= new GFFObject();
			obj.setSeqname(trpt.getChromosome());
			obj.setSource(trpt.getSource());
			obj.setFeature(GTF_FEATURE);
			obj.setStrand(trpt.getStrand());
			obj.setStart(regs[i].getStart());
			obj.setEnd(regs[i].getEnd());
			obj.addAttribute(GFFObject.TRANSCRIPT_ID_TAG, trpt.getTranscriptID());
			obj.addAttribute(GTF_ATTRIBUTE_EXONIC_START, Integer.toString(getStart()));
			obj.addAttribute(GTF_ATTRIBUTE_EXONIC_END, Integer.toString(getEnd()));
			obj.addAttribute(GTF_ATTRIBUTE_EXONIC_LENGTH, Integer.toString(trpt.getExonicLength()));
			obj.addAttribute(GTF_ATTRIBUTE_FRAGMENT_ID, Integer.toHexString(this.hashCode()));
			if (attributes!= null) {
				Iterator<String> keyIter= attributes.keySet().iterator();
				while (keyIter.hasNext()) {
					String key= keyIter.next();
					obj.addAttribute(key, attributes.get(key));
				}
			}
			sb.append(obj.toString());			
			sb.append("\n");
		}
		
		return sb.toString();
	}
	
	/**
	 * @deprecated deactivated due to compilation erros
	 * @param trpt
	 * @param attributes
	 * @return
	 */
	public String toStringBED(Transcript trpt, HashMap<String, String> attributes) {
		/*
		DirectedRegion[] regs= getRegions(trpt);
		int debug= 0;
		if (regs== null) {
			System.currentTimeMillis();
			getRegions(trpt);
		}
		for (int i = 0; i < regs.length; i++) {
			debug+= regs[i].getLength();
		}
		
		BEDobject obj= new BEDobject(regs);
		//String name= Integer.toHexString(this.hashCode());	
		//String name= trpt.getTranscriptID()+"_"+attributes.get(Sequencer.GTF_ATTRIBUTE_READ_NR);
		String name= trpt.getTranscriptID()+"@"+attributes.get(Sequencer.GTF_ATTRIBUTE_RELXPR);

		obj.setName(name);
		//obj.setCol(c);
		
		String bedStr= obj.toString()+ "\n";
		return bedStr;
		*/
		return null;
	}
	/**
	 * @deprecated deactivated due to compilation erros
	 * @param regs
	 * @param tid
	 * @param attributes
	 * @return
	 */
	public String toStringBED(DirectedRegion[] regs, String tid, HashMap<String, String> attributes) {
		/*
		int debug= 0;
			for (int i = 0; i < regs.length; i++) {
			debug+= regs[i].getLength();
		}
		
		BEDobject obj= new BEDobject(regs);
		//String name= Integer.toHexString(this.hashCode());	
		String name= tid+"_"+attributes.get(Sequencer.GTF_ATTRIBUTE_READ_NR);
		obj.setName(name);
		//obj.setCol(c);
		
		String bedStr= obj.toString()+ "\n";
		return bedStr;
		*/
		return null;
	}
	
	public Object clone() {
		Fragment nuFrag= new Fragment(ancestorID, start, end);
		return nuFrag;
	}
	public int getNb() {
		return nb;
	}
	//@Override
	public int compareTo(Object o) {
		//return (int) (key- ((Fragment) o).getKey());	// not for rounding errors
		if (key< ((Fragment) o).getKey())
			return -1;
		else if (key> ((Fragment) o).getKey())
			return 1;		
		return 0;
	}
	public double getWeight() {
		return weight;
	}
	public void setWeight(double weight) {
		this.weight = weight;
	}
	public double getKey() {
		return key;
	}
	public void setKey(double key) {
		this.key = key;
	}
	public void setNb(int nb) {
		this.nb = nb;
	}
	public String toStringGTF(DirectedRegion[] regs, String string,
			HashMap<String, String> h) {
		// TODO Auto-generated method stub
		return null;
	}
}
