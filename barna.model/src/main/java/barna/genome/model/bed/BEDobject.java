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

package barna.genome.model.bed;

//import gphase.gui.ColorMaster;

import barna.commons.ByteArrayCharSequence;
import barna.genome.model.AnnotationObject;
import barna.genome.model.DirectedRegion;
import barna.genome.model.Graph;

import java.awt.*;
import java.util.NoSuchElementException;
import java.util.Vector;

/**
 * 080210 major bugfix in the coordinates of <code>fromGTFObjects()</code>
 * 090811 change coordinate system: internally only positive coordinates
 * @author micha
 *
 */
public class BEDobject extends AnnotationObject {
	
	static int maxRecycle= 100000;
	static Vector<BEDobject> recycleObjV= new Vector<BEDobject>();

	public static BEDobject getRecycleObj(int asize) {
		
		BEDobject o= null;
		
		synchronized (recycleObjV) {
			int z= recycleObjV.size();
			for (int i = 0; i < z; i++) {
				if (recycleObjV.elementAt(z- 1- i).chrom.chars.length>= asize) {
					o= recycleObjV.remove(z- 1- i);
					break;
				}
			}
		}
		
		return o;
	}
	
	public static BEDobject getTryRecycleObj() {	
		
		if (1== 1)
			return null;
		
		BEDobject bed= null;
		try {
			synchronized(recycleObjV) {
				bed= getRecycleObjV().remove(recycleObjV.size()- 1);
			} 
			
		} catch (Exception e) {	// NullPointer, ArrayIndexOutOfBounds
			; // :)
		}
		return bed;
	}
		
	public static BEDobject getRecycleObj() {	
		BEDobject bed= getTryRecycleObj();
		if (bed== null)
			bed= new BEDobject();				

		return bed;
	}
	public static boolean addRecycleObj(BEDobject obj) {
		
		if (1== 1)
			return false;
		
		obj.blockCount= -1;
		obj.blockStarts= null;
		obj.blockSizes= null;
		obj.score= Integer.MIN_VALUE;
		obj.name= null;
		obj.col= null;
		obj.thickStart= -1;
		obj.thickEnd= -1;
		obj.strand= 0;
		obj.start= 0;
		obj.end= 0;
		obj.chrom= null;
		// do NOT set cs to null

		synchronized (recycleObjV) {
			if (recycleObjV.size()>= maxRecycle)
				return false;
			recycleObjV.add(obj);
		}
		
		return true;
	}
	
	public static boolean addRecycleObj_original(BEDobject obj) {
		
		obj.chrom= null;
		obj.start= 0;
		obj.end= 0;
		obj.name= null;
		obj.thickStart= 0;
		obj.thickEnd= 0;
		obj.score= 0;
		obj.col= null;
		obj.blockCount= 0;
		obj.blockSizes= null;
		obj.blockStarts= null;

		synchronized(recycleObjV) {
			if (getRecycleObjV().size()>= maxRecycle)
				return false;
			getRecycleObjV().add(obj);
		}
		return true;
	}
	
	public static Vector<BEDobject> getRecycleObjV() {
		if (recycleObjV == null) {
			recycleObjV = new Vector<BEDobject>(maxRecycle);
		}

		return recycleObjV;
	}
	
	public static String toString(int[] a) {
		StringBuilder sb= new StringBuilder(a.length);
		for (int i = 0; i < a.length; i++) {
			sb.append(Integer.toString(a[i]));
			sb.append(COMMA);
		}
		sb.deleteCharAt(sb.length()- 1);
		return sb.toString();
	}
	
	int start= 0, end= 0, thickStart= -1, thickEnd= -1, blockCount= -1;
	ByteArrayCharSequence cs= null;
	/**
	 * 1 positive, -1 negative, 2 default positive
	 */
	byte strand= 0;	// assume forward, bed line can be incomplete
	int score= Integer.MIN_VALUE;	// could be (-1), defined: 0..1000
	ByteArrayCharSequence name= null, chrom= null, blockSizes= null, blockStarts= null;
	Color col= null;
	public String readSequence() {
		int start= this.start; // (strand> 0)?this.start:Math.abs(end);
		String s= "";
		try {
			int f= getStrand()>= 0? 1: -1;
			if (getBlockSizes()== null)
				s= Graph.readSequence(null, chrom, getStrand()>= 0, f* (start+ 1), f* end);
			else for (int i = 0; getBlockSizes()!= null && i < getBlockCount(); i++) 
				if (getStrand()>= 0)
					s+= Graph.readSequence(null, chrom, getStrand()>= 0, f* (start+ getBlockStart(i)+ 1), f* (start+ getBlockStart(i)+ getBlockSize(i)));
				else
					s= Graph.readSequence(null, chrom, getStrand()>= 0, f* (start+ getBlockStart(i)+ 1), f* (start+ getBlockStart(i)+ getBlockSize(i)))+ s;
		} catch (Exception e) {
			throw new RuntimeException("Problems reading BED object sequence:\n\t"+ toString()+ "\n\t"+ e.getMessage());
		}
		return s;		
	}
	public void readSequence(ByteArrayCharSequence cs) {
		int start= this.start; // (strand> 0)?this.start:Math.abs(end);
		try {
			int f= getStrand()>= 0? 1: -1;
			if (getBlockSizes()== null) {
				int first= f* (start+ 1), 
					last= f* end,
					len= last- first+ 1; 
				cs.ensureLength(cs.start, len);
				Graph.readSequence(null, chrom, getStrand()>= 0, 
						first, last, 
						cs, cs.start, len);
			} else for (int i = 0; getBlockSizes()!= null && i < getBlockCount(); i++) { 
				int first= f* (start+ getBlockStart(i)+ 1), 
					last= f* (start+ getBlockStart(i)+ getBlockSize(i)),
					len= last- first+ 1; 
				if (getStrand()>= 0) {
					cs.ensureLength(cs.end, len);
					Graph.readSequence(null, chrom, getStrand()>= 0, 
							first,	// f* (start+ getBlockStart(i)+ 1) 
							last,	// f* (start+ getBlockStart(i)+ getBlockSize(i)) 
							cs, cs.end, cs.end+ len);
				} else {
					cs.ensureLength(cs.start, len);
					Graph.readSequence(null, chrom, getStrand()>= 0, 
							first, last, 
							cs, cs.start, len);
				}
			}
		} catch (Exception e) {
			throw new RuntimeException("Problems reading BED object sequence:\n\t"+ toString()+ "\n\t"+ e.getMessage());
		}
	}
	
	public BEDobject(ByteArrayCharSequence chromName, byte newStrand) {
		setChrom(chromName);
		setStrand(newStrand);
	}

	public BEDobject(String chromName, byte newStrand) {
		this(new ByteArrayCharSequence(chromName), newStrand);
	}
	public BEDobject() {
	}
	public static BEDobject createBEDobject(ByteArrayCharSequence cs, String chr, int bedStart, int bedEnd) {
		BEDobject o= getTryRecycleObj();
		if (o== null)
			return new BEDobject(cs, chr, bedStart, bedEnd);	// not: cs.cloneCurrentSeq
		
		o.init(cs, o.chrom, bedStart, bedEnd);
		return o;
	}
	
	public BEDobject(ByteArrayCharSequence cs, ByteArrayCharSequence chr, int bedStart, int bedEnd) {
		
		init(cs, chr, bedStart, bedEnd);
	}
	
	public BEDobject(ByteArrayCharSequence cs, String chr, int bedStart, int bedEnd) {
		
		ByteArrayCharSequence chrcs= new ByteArrayCharSequence(chr);
		init(cs, chrcs, bedStart, bedEnd);
	}
	
	private void init(ByteArrayCharSequence cs, ByteArrayCharSequence chromi, int bedStart, int bedEnd) {
		//setChrom(cs.getToken(1, TAB));
		this.cs= ByteArrayCharSequence.cloneSequence(cs, this.cs);

		//chrom= chromi;	// reinit, not sure whether its the cloned sequence
		getChrom();
		start= (bedStart);
		end= (bedEnd);

/*		ByteArrayCharSequence subs= cs.getToken(4, TAB);
		if (subs== null)
			return;
		else
			name= subs;
		subs= cs.getToken(5, TAB);
		if (subs== null) 
			return;
		else {
//			try {
				setScore(BEDobject.encodeInt(subs, 0, subs.length()));
//			} catch (NumberFormatException e) {
//				; // :)
//			}
		}
		
		subs= cs.getToken(6, TAB);
		if (subs== null) 
			return;
		else
			setStrand(BEDobject.parseStrand(subs));
		
						
		subs= cs.getToken(7, TAB);
		if (subs== null) 
			return;
		else {
			//try {
				setThickStart(BEDobject.encodeInt(subs, 0, subs.length()));
//			} catch (NumberFormatException e) {
//				; // :)
//			}
		}
		subs= cs.getToken(8, TAB);
		if (subs== null) 
			return;
		else 
			//try {
				setThickEnd(BEDobject.encodeInt(subs, 0, subs.length()));
//			} catch (NumberFormatException e) {
//				; // :)
//			}
		
		subs= cs.getToken(9, TAB);
		if (subs== null) 
			return;
		else
			setCol(subs);

		subs= cs.getToken(10, TAB);
		if (subs== null)
			return;
		else
			setBlockCount(BEDobject.encodeInt(subs));
		subs= cs.getToken(11, TAB);
		setBlockSizes(subs);
		subs= cs.getToken(12, TAB);
		setBlockStarts(subs);
*/		
	}
	
	public BEDobject(DirectedRegion[] obs) {
		fromGTFObjects(obs);
	}
	
	public BEDobject(String chromName, byte newStrand, int newStart, int newEnd) {
		this(new ByteArrayCharSequence(chromName), newStrand, newStart, newEnd);
	}
	
	public BEDobject(ByteArrayCharSequence chromName, byte newStrand, int newStart, int newEnd) {
		this(chromName, newStrand);
		setStart(newStart);
		setEnd(newEnd);
	}

	private void fromGTFObjects(DirectedRegion[] obs) {
		int minStart= Integer.MAX_VALUE, maxEnd= Integer.MIN_VALUE;
		String chr= null, name= null, id= null;
		byte strand= 0;
		for (int i = 0; i < obs.length; i++) {
			if (Math.abs(obs[i].getStart())< minStart)
				minStart= Math.abs(obs[i].getStart());
			if (Math.abs(obs[i].getEnd())> maxEnd)
				maxEnd= Math.abs(obs[i].getEnd());	
			if (chr== null)
				chr= obs[i].getChromosome();
			else if (!chr.equals(obs[i].getChromosome()))
				System.err.println("Chr not matching "+obs[i].getChromosome()
						+"("+chr+").");
			if (strand== 0)
				strand= (byte) obs[i].getStrand();
			else if (obs[i].getStrand()!= strand)
				System.err.println("Strand not matching "+obs[i].getStrand()
						+"("+strand+").");
			
			String newName= null;
			
			// TODO check whether still needed
//			if (obs[i].getAttribute(LaVista.GTF_ATTRIBUTE_GROUP_ID)== null)
//				newName= obs[i].getID();
//			else
//				newName= (String) obs[i].getAttribute(LaVista.GTF_ATTRIBUTE_GROUP_ID);
//			if (id== null)
//				id= newName;
//			if (obs[i].getAttribute(GTFObject.TRANSCRIPT_ID_TAG)!= null)
//				newName+= "_"+obs[i].getAttribute(GTFObject.TRANSCRIPT_ID_TAG); 
			
			if (name== null) 
				name= newName;
			else if (!newName.equals(name))
				System.err.println("Name/group not matching "+newName
						+"("+name+").");
		}
		--minStart;	// bed is 0-based, browser 1-based
		//++maxEnd;	// AND end position not included in the feature
		
		setChrom(chr);
		setStrand(strand);
		setStart(minStart);
		setEnd(maxEnd); 
		//obj.setScore(score);
		setName(name);
		Color c= null;
//		String id= name;
//		int pos= name.lastIndexOf("_");
//		if (pos>= 0)
//			id= name.substring(0,pos);		
		
		// TODO check whether still needed
//		if (ColorMaster.getDefaultColorMap().get(id)== null)
//			c= ColorMaster.getRandomColor();
//		else
//			c= ColorMaster.getDefaultColorMap().get(id);
		setCol(c);
		
		int[] starts= new int[obs.length], lengths= new int[obs.length];
		if (obs[0].getStrand()< 0) 		
			for (int i = 0; i < obs.length; i++) {
				starts[obs.length-1-i]= Math.abs(obs[i].getStart())- minStart- 1;	// bed is 0-based, browser 1-based
				lengths[obs.length-1-i]= Math.abs(obs[i].getEnd()- obs[i].getStart())+ 1;
			}
		else
			for (int i = 0; i < obs.length; i++) {
				starts[i]= Math.abs(obs[i].getStart())- minStart- 1;	// bed is 0-based, browser 1-based
				lengths[i]= Math.abs(obs[i].getEnd()- obs[i].getStart())+ 1;
			}

		setBlockCount(obs.length);
		setBlockStarts(toString(starts));
		setBlockSizes(toString(lengths));

	}
	
	public int getBlockCount() {
		if (blockCount == -1&& (cs!= null)) {
			ByteArrayCharSequence subs= cs.getToken(9);
			if (subs== null) {
				setBlockCount(0);
			} else
				setBlockCount(BEDobject.encodeInt(subs));
			
		}
		
		return blockCount;
	}

	public void setBlockCount(int blockCount) {
		this.blockCount = blockCount;
	}

	public boolean check() {
		return ((getStart()+ getBlockStart(getBlockCount()- 1)+ getBlockSize(getBlockCount()- 1))
				== getEnd());
	}
	
	public int length() {
		if (getBlockCount()<= 1)
			return (getEnd()- getStart());
		int sum= 0;
		for (int i = 0; i < getBlockCount(); i++) 
			sum+= getBlockSize(i);
		return sum;
	}
	
	private static NumberFormatException defaultNumberFormatException= new NumberFormatException();
	/**
	 * 
	 * @param s
	 * @param from
	 * @param to
	 * @return the decoded number
	 */
	public static final int encodeInt(CharSequence s, int from, int to) {
		int sum= 0, fac= 1;
		for (int i = to- 1; i >= from; --i) {
			int x= ((byte) s.charAt(i))- 48;
			if (x< 0|| x> 9)
				return -1; //throw defaultNumberFormatException;
			sum+= x* fac;
			fac= 10* fac;
		}
		return sum;
	}

	
	public static final int encodeInt(CharSequence s) {
		return encodeInt(s, 0, s.length());
	}
	
	/**
	 * 
	 * @param s
	 * @param sep
	 * @param nr 1-based !!!
	 * @return
	 */
	public static final int getInt(CharSequence s, char sep, int nr) throws NoSuchElementException {
		int last= 0, curr= 0, cnt= 0;		
		for (int i = 0; i < s.length()&& cnt< nr; i++) {
			if (s.charAt(i)== sep) {
				last= curr;
				curr= i;
				++cnt;
			}
		}
		if (cnt< nr- 1)
			throw new NoSuchElementException("Found "+(cnt+1)+" token");
		if (cnt< nr) {
			last= curr;
			curr= s.length();
		}
		if (last!= 0)
			++last;

		int x= encodeInt(s, last, curr);
		return x;
	}
	
	public static final CharSequence getField(CharSequence s, char sep, int nr) throws NoSuchElementException {
		int last= 0, curr= 0, cnt= 0;		
		for (int i = 0; i < s.length()&& cnt< nr; i++) {
			if (s.charAt(i)== sep) {
				last= curr;
				curr= i;
				++cnt;
			}
		}
		if (cnt< nr- 1)
			throw new NoSuchElementException("Found "+(cnt+1)+" token");
		if (curr== 0)
			curr= s.length();
		return s.subSequence(last, curr);
	}
	
	private static final char COMMA= ',';
	public ByteArrayCharSequence getBlockSizes() {
		
		if (blockSizes == null&& cs!= null) {
			ByteArrayCharSequence subs= cs.getToken(10);
			if (subs!= null)
				setBlockSizes(subs);
		}

		return blockSizes;
	}
	public ByteArrayCharSequence getBlockStarts() {
		if (blockStarts == null) {
			ByteArrayCharSequence subs= cs.getToken(11);
			if (subs!= null)
				setBlockStarts(subs);
		}

		return blockStarts;
	}
	public int getBlockSize(int nr) {
		if (getBlockCount()== 0)
			return length();
		return getInt(getBlockSizes(), COMMA, nr+ 1);
	}
	public int getBlockStart(int nr) {
		if (getBlockCount()== 0)
			return 0;
		return getInt(getBlockStarts(), COMMA, nr+ 1);
	}
	
	/**
	 * mem-inefficient!
	 * @param blockSize
	 */
	public void addBlockSize(int blockSize) {
		if (blockSizes== null)
			blockSizes= new ByteArrayCharSequence(Integer.toString(blockSize));
		else 
			blockSizes.append(COMMA+ Integer.toString(blockSize));
	}
	
	public void addBlockStart(int blockStart) {
		if (blockStarts== null)
			blockStarts= new ByteArrayCharSequence(Integer.toString(blockStart));
		else 
			blockStarts.append(COMMA+ Integer.toString(blockStart));
	}
	

	public CharSequence getChrom() {
		if (chrom == null) {
			ByteArrayCharSequence subs= cs.getToken(0);
			if (subs!= null) 
				setChrom(subs);
		}

		return chrom;
	}

	public void setChrom(ByteArrayCharSequence chrom) {
		this.chrom = chrom;
	}

	public void setChrom(String chrom) {
		this.chrom = new ByteArrayCharSequence(chrom);
	}
	
	public boolean isStrandInited() {
		return (getStrand()!= 0);
	}
	
	public Color getCol() {
		if (col == null) {
			if (cs!= null) {
				ByteArrayCharSequence subs= cs.getToken(8);
				if (subs!= null) 
					setCol(subs);
			}
		}

		return col;
	}

	public void setCol(Color col) {
		this.col = col;
	}
	
	public void setCol(int red, int green, int blue) {
		setCol(new Color(red, green, blue));
	}
	
	public void setCol(CharSequence rgbVal) {
		int[] tokens= parseCommaSeparatedInts(rgbVal, null);
		if (tokens.length!= 3) {
			if (!rgbVal.equals("0"))
				System.out.println("WARNING: invalid color "+rgbVal);
			return;
		}
		setCol(tokens[0],tokens[1],tokens[2]);
	}
	
	int[] parseCommaSeparatedInts(CharSequence in, int[] out) {
//		String[] tokens= in.split(",");
//		int[] out= new int[tokens.length];
//		for (int i = 0; i < out.length; i++) 
//			out[i]= Integer.parseInt(tokens[i]);
		
		int cnt= 0;
		for (int i = 0; i < in.length(); i++) 
			if (in.charAt(i)== ',')
				++cnt;
		++cnt;
		if (out== null)
			out= new int[cnt];
		int last= 0;
		cnt= 0;
		for (int i = 0; i < in.length(); i++) { 
			if (in.charAt(i)== ',') {
				out[cnt++]= encodeInt(in, last, i);
				last= i+1;
			}
		}
		if (cnt>= out.length)
			System.currentTimeMillis();
		out[cnt++]= encodeInt(in, last, in.length());
		
		return out;
	}
	
	public void setBlockStarts(ByteArrayCharSequence in) {
		this.blockStarts= in;
	}
	
	public void setBlockStarts(String in) {
		this.blockStarts= new ByteArrayCharSequence(in);
	}
	
	public void setBlockSizes(ByteArrayCharSequence in) {
		this.blockSizes= in;
	}
	
	public void setBlockSizes(String in) {
		this.blockSizes= new ByteArrayCharSequence(in);
	}
	

	public int getEnd() {
		return end;
	}

	public int get5Prime() {
		if (getStrand()>= 0)
			return start;
		return end;
	}
	
	public int get3Prime() {
		if (getStrand()>= 0)
			return end;
		return start;
	}
	
	/**
	 * sets the genomic region, even if fed with neg. coordinates
	 * @param end
	 */
	public void setEnd(int end) {
		
		int abs = Math.abs(end);
		if (start!= 0&& abs< start) {
			this.end= start;
			start= abs;
		} else
			this.end= abs;
		
//		this.end = Math.abs(end);
//		if (start>= 0&& start> this.end) {
//			int h= start;
//			start= this.end;
//			this.end= h;
//		}
	}

	public CharSequence getName() {
		if (name == null) {
			name = cs.getToken(3);
		}
		return name;
	}

	public void setName(ByteArrayCharSequence name) {
		this.name = name;
	}

	public int getScore() {
		if (score == Integer.MIN_VALUE) {
			if (cs!= null) {
				ByteArrayCharSequence subs= cs.getToken(4);
				if (subs!= null) 
					setScore(BEDobject.encodeInt(subs, 0, subs.length()));
			}
		}

		return score;
	}

	public void setScore(int score) {
		this.score = score;
	}

	public int getAbsoluteStart() {
		return start;
	}
	
	public int getAbsoluteEnd() {
		return end;
	}

	/**
	 * sets the genomic region, even if fed with neg. coordinates
	 * @param start
	 */
	public void setStart(int start) {
		int absStart = Math.abs(start);
		if (end!= 0&& absStart> end) {
			this.start= end;
			end= absStart;
		} else
			this.start= absStart;
//		if (end>= 0&& this.start> end) {
//			int h= this.start;
//			this.start= end;
//			end= h;
//		}
			
	}
	
	public int getAbsBlockStart(int nr) {
		return getStart()+ getBlockStart(nr);
	}

	public int getAbsBlockEnd(int nr) {
		return getStart()+ getBlockStart(nr)+ getBlockSize(nr);
	}
	
	public byte getStrand() {
		if (strand == 0) {
			ByteArrayCharSequence subs= cs.getToken(5);
			if (subs!= null)
				setStrand(BEDobject.parseStrand(subs));
		}

		return strand;
	}

	public void setStrand(byte strand) {
		if (strand!= 1&& strand!= -1)
			System.out.println("WARNING: No strand assignment for "+this);

		this.strand = (strand> 0)?(byte) 1:(byte) (-1);
		
//		if (strand< 0) {	// check for swap
//			start= -1* Math.abs(getStart());
//			end= -1* Math.abs(getEnd());
//			if (start> end) {
//				int h= start;
//				start= end;
//				end= h;
//			}
//		}
		
		
		
	}
	
	public void setStrand(String str) {
		setStrand(parseStrand(str));
	}
	private static final char PLUS= '+', MINUS= '-';
	public static byte parseStrand(CharSequence str) {
		if (str.charAt(0)== PLUS)
			return 1;
		else if (str.charAt(0)== MINUS)
			return (-1);
		else
			System.out.println("WARNING: invalid strand tag "+str);
		return 0;
	}

	public int getThickEnd() {
		if (thickEnd == -1) {
			ByteArrayCharSequence subs= cs.getToken(6);
			if (subs!= null) 
				setThickEnd(BEDobject.encodeInt(subs, 0, subs.length()));
		}

		return thickEnd;
	}

	public void setThickEnd(int thickEnd) {
		thickEnd= Math.abs(thickEnd);
		this.thickEnd = thickEnd;
		if (thickStart>= 0&& thickStart> this.thickEnd) {
			int h= thickStart;
			thickStart= this.thickEnd;
			this.thickEnd= h;
		}
	}

	public int getThickStart() {
		if (thickStart == -1) {
			ByteArrayCharSequence subs= cs.getToken(6);
			if (subs!= null) 
				setThickStart(BEDobject.encodeInt(subs, 0, subs.length()));
		}

		return thickStart;
	}

	public void setThickStart(int thickStart) {
		thickStart= Math.abs(thickStart);
		this.thickStart = thickStart;
		if (thickEnd>= 0&& this.thickStart> thickEnd) {
			int h= this.thickStart;
			this.thickStart= thickEnd;
			thickEnd= h;
		}
	}
	
	public void extend(boolean atStart, int diff) {
		
		if (atStart)
			start-= diff;
		else
			end+= diff;
		
		if (blockCount> 0) {
			int[] starts= new int[blockCount], sizes= new int[blockCount];
			for (int i = 0; i < blockCount; i++) {
				starts[i]= getBlockStart(i);
				sizes[i]= getBlockSize(i);
			}
			if (atStart) {
				sizes[0]+= diff;
				for (int i = 1; i < starts.length; i++) 
					starts[i]+= diff;
			} else {
				sizes[sizes.length- 1]+= diff;
			}
			setBlockStarts(toString(starts));
			setBlockSizes(toString(sizes));
		}
	}
	
	public boolean trim(boolean fromStart, int diff) {
		
		if (diff> length())
			return false;
		
		if (getBlockSizes()== null|| getBlockStarts()== null) {
			if (fromStart)
				start+= diff;
			else 
				end-= diff;
			return true;
		}
		int afterStart= 0, afterLen= getBlockCount();
		int bcount= 0;
		StringBuilder sizeb= new StringBuilder(getBlockSizes().length()), startb= new StringBuilder(getBlockStarts().length());
		if (fromStart) {
			start+= diff;
			int left= diff;
			for (int i= 0; i< getBlockCount(); ++i) {
				int x= Math.min(left, getBlockSize(i));
				int y= getBlockSize(i)- x;
				int z= getBlockStart(i); // + x; NO! intrinsic
				if (y== 0) 
					afterStart= i+1;
				else {
					sizeb.append(Integer.toString(y));
					sizeb.append(COMMA);
					startb.append(Integer.toString(z));
					startb.append(COMMA);
					++bcount;
				}
				left-= x;
			}
		} else {
			end-= diff;
			int left= diff;
			int i= 0;
			for (i= getBlockCount()- 1; i>= 0; --i) {	
				int x= Math.min(left, getBlockSize(i));
				left-= x;
				int y= getBlockSize(i)- x;				
				if (y== 0) 
					afterLen= i;
				else {
					sizeb.insert(0, COMMA);
					sizeb.insert(0, Integer.toString(y));
					startb.insert(0, COMMA);
					startb.insert(0, getBlockStart(i));
					++bcount;
				}
			}
		}
		
/*		if (afterStart> 0) {
			int[] news= new int[getBlockSizes().length- afterStart];
			System.arraycopy(getBlockSizes(), afterStart, news, 0, news.length);
			blockSizes= news;
			news= new int[getBlockSizes().length- afterStart];
			System.arraycopy(getBlockStarts(), afterStart, news, 0, news.length);
			blockStarts= news;
			blockCount-= afterStart;
		} 
		if (afterLen!= getBlockSizes().length) {
			int[] news= new int[afterLen];
			System.arraycopy(getBlockSizes(), 0, news, 0, afterLen);
			blockSizes= news;
			news= new int[afterLen];
			System.arraycopy(getBlockStarts(), 0, news, 0, afterLen);
			blockStarts= news;
			blockCount= afterLen;
		}
*/
		if (sizeb.length()== 0) {
			blockStarts= null;
			blockSizes= null;
			blockCount= 0;
		} else {
			sizeb.deleteCharAt(sizeb.length()- 1);
			blockSizes= new ByteArrayCharSequence(sizeb.toString());
			startb.deleteCharAt(startb.length()- 1);
			blockStarts= new ByteArrayCharSequence(startb.toString());
			blockCount= bcount;
		}
		
		return true;
	}
	
	private static final String EMPTY_STRING= "";
	
	private static final char ZERO= '0';
	@Override
	public String toString() {
		
		StringBuilder sb= new StringBuilder();
		sb.append(getChrom());
		sb.append(TAB);
		sb.append(Integer.toString(getStart()));
		sb.append(TAB);
		sb.append(Integer.toString(getEnd()));
		if (getName()!= null) {
			sb.append(TAB);
			sb.append(getName());
		} else {
			sb.append("\tnobody");
		}
		if (getScore()!= Integer.MIN_VALUE) {
			sb.append(TAB);
			sb.append(getScore());
		} else if (getStrand()!= 0|| thickStart>= 0|| getCol()!= null|| blockCount> 0) {
			sb.append(TAB);
			sb.append(ZERO);
		}
		if (getStrand()!= 0) {
			if (getStrand()> 0) {
				sb.append(TAB);
				sb.append("+");
			} else { 
				sb.append(TAB);
				sb.append("-");
			}
		} else if (thickStart>= 0|| getCol()!= null|| blockCount> 0) {
			if (thickStart>= 0|| getCol()!= null|| getBlockCount()> 1)
				sb.append(DOT);
		}
		if (thickStart>= 0&& thickEnd>= 0) {
			sb.append(TAB);
			sb.append(getThickStart());
			sb.append(TAB);
			sb.append(getThickEnd());
		} else if (getCol()!= null|| blockCount> 0) {
			if (thickStart>= 0|| getCol()!= null|| getBlockCount()> 1)
				sb.append(DOT);
		} 
		
		if (getCol()!= null) {
			sb.append(TAB);
			sb.append(getCol().getRed());
			sb.append(COMMA);
			sb.append(getCol().getGreen());
			sb.append(",");
			sb.append(getCol().getBlue());
		} else if (blockCount> 0) {
			if (thickStart>= 0|| getCol()!= null|| getBlockCount()> 1)
				sb.append(DOT);
		} 
		
		if (blockCount> 1) {
			sb.append(TAB);
			sb.append(getBlockCount());
			sb.append(TAB);
			sb.append(getBlockSizes());
			sb.append(TAB);
			sb.append(getBlockStarts());
		}
		
		return sb.toString();
	}
	
	private static final char TAB= '\t', DOT= '.'; 
	
	public int getStart() {
		return start;
	}
	public void setName(String name) {
		if (name== null) {
			this.name= null;
			return;
		}
		this.name = new ByteArrayCharSequence(name);
	}
	
}
