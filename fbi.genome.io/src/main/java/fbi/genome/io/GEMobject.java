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

package fbi.genome.io;

import fbi.commons.ByteArrayCharSequence;
import fbi.commons.Progressable;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.OutputStreamWriter;
import java.util.HashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;


public class GEMobject {
	
	public static final String SEP_MATCHPOS= ",";
	public static final Pattern PATTY_TAB= Pattern.compile("\\t"), PATTY_COMMA= Pattern.compile(SEP_MATCHPOS);
	int[] pos= null;	// = new int[4];
	ByteArrayCharSequence base= null;
	
	/**
	 * SOLEXA:4:1:43:293#0/1   
	 * TTTTTTTTTTTTTTTTTTTTTTTAAATTTTTCCTTT    
	 * aaaaaaaababaabaaa]RBBBBBBBBBBBBBBBBB    
	 * 0:1:6   
	 * chr2-:46494402T24,chr14+:102145117A32T33,chr18-:39752609T32T33,chr2-:73856487T32T33,chr6-:24692626C28G32,chr10-:118490450C25T26,chr11-:30324323T26C36
	 */
	public static final int GEM_NR_SEP= 4;
	
	public static final void GEMtoBED(File inFile, File outFile, Progressable prog) {
		try {
			ThreadedBufferedByteArrayStream buffy= new ThreadedBufferedByteArrayStream(100000, inFile, true);
			BufferedWriter vampire= null;
			if (outFile== null)
				vampire= new BufferedWriter(new OutputStreamWriter(System.out));
			else
				vampire= new BufferedWriter(new FileWriter(outFile));
			
			GEMobject go= new GEMobject();
			int[] err= new int[3];
			byte[] tmpError= null;
			long bytesRead= 0, totBytes= inFile.length();
			int perc= 0;
			if (prog!= null) {
				prog.start("parsing error model from GEM alignment ");
			}
			long readCtr= 0;
			ByteArrayCharSequence cs= new ByteArrayCharSequence(10000);
			char[] maxName= new char[100];
			for (ByteArrayCharSequence line = buffy.readLine(cs); line.end!= 0; line = buffy.readLine(cs)) {
				bytesRead+= line.length()+ 1;	// fsep
				++readCtr;
				if (bytesRead*10d/ totBytes> perc) {
					if (prog!= null)
						prog.progress();
					//System.err.print("*");
					++perc;
				}
				
				if (go.reuse(line)== null)
					continue;
				
				for (int j = go.pos[3]; j < line.length();) {
					
					int l= j;
					for (; j < line.length()&& line.charAt(j)!= ','; j++); 
					String chr= go.getMatchPos(line.subSequence(l, j), go.p);	// find match
					if (chr== null)
						break;
					
					vampire.write(chr);
					vampire.write("\t");
					vampire.write(Integer.toString(go.p[0]));
					vampire.write("\t");
					vampire.write(Integer.toString(go.p[0]+ go.getSequence().length()));
					ByteArrayCharSequence nameCS= go.getName();	// TODO FMRD !!!
					maxName= nameCS.toCharArray(maxName);
					vampire.write(maxName, 0, nameCS.length());
					vampire.write("\t.\t");	// score
					vampire.write(Character.toString(go.p[1]> 0? SYMBOL_POSITIVE: SYMBOL_NEGATIVE));
					vampire.write("\n");
					vampire.flush();
				}
				
			}
			
			//buffy.setStop(true);
			vampire.flush();
			vampire.close();
			
		} catch (Exception e) {
			e.printStackTrace();
		}
		
	}
	
	public static class Match {
		char[] chr;
		int chrLen;
		boolean fwd;
		int pos;
		int[] mmLoc;
		int MMlocLen;
		char[] mm;
		public char[] getChr() {
			return chr;
		}
		public void setChr(ByteArrayCharSequence chr) {
			if (chr== null) {
				this.chr = new char[chr.length()];
				chrLen= chr.length();
			} else if (this.chr.length>= chr.length()) {
				System.arraycopy(chr, 0, this.chr, 0, chr.length());
			}
				
		}
		public int getChrLen() {
			return chrLen;
		}
		public boolean isFwd() {
			return fwd;
		}
		public void setFwd(boolean fwd) {
			this.fwd = fwd;
		}
		public int getPos() {
			return pos;
		}
		public void setPos(int pos) {
			this.pos = pos;
		}
		public int[] getMmLoc() {
			return mmLoc;
		}
		public void setMmLoc(int[] mmLoc) {
			this.mmLoc = mmLoc;
		}
		public int getMMlocLen() {
			return MMlocLen;
		}
		public void setMMlocLen(int mlocLen) {
			MMlocLen = mlocLen;
		}
		public char[] getMm() {
			return mm;
		}
		public void setMm(char[] mm) {
			this.mm = mm;
		}
	}
	
	
	public GEMobject() {
		pos= new int[4];
	}

	public boolean isValid() {
		return (pos!= null);
	}
	
	
	
	public GEMobject reuse(ByteArrayCharSequence s) {
		base= s;
		int lastPos= Integer.MIN_VALUE, c= 0;
		try {
			for (int i = 0; i < s.length(); i++) 
				if (s.charAt(i)== '\t') {	// only tabs now, allow spaces in fasta header tag
					if (i!= lastPos+ 1) {
						pos[c++]= i;
						lastPos= i;
					}
				}
		} catch (ArrayIndexOutOfBoundsException e) {
			return null;
		}
		// too short 
		if(c< pos.length- 1) 
			return null;
		
		if (c< pos.length)	// no qualities
			for (int i = pos.length-1; i>= 2; --i) 
				pos[i]= pos[i-1];
		
		return this;
	}

	/**
	 * 
	 * @param x	1-based
	 * @return
	 */
	private int[] p= new int[2];
	int[] getMatch(int x, int from) {
		p[0]= (from<0)?pos[pos.length-1]:from;
		p[1]= -1;
		for (int i = 0; p[1]< base.length()&& i < x; i++) {
			p[1]= p[0];
			while (base.charAt(++p[1])!= ','&& p[1]< base.length());
		}
		
		return p;
	}
	
	public ByteArrayCharSequence getMatch(int from) {
		
		if (from== -1)
			from= pos[3]+1;
		else if (from== base.length())
			return null;
		
		int p= from+1;
		while (p< base.length()&& base.charAt(p++)!=',');
		
		return base.subSequence(from, p);
	}
	
	Pattern pattMM= Pattern.compile("([A,C,G,T,N,a,c,g,t,n]\\d{1,3}+)");
	public Matcher getMismatch(ByteArrayCharSequence cs, Matcher m, int[] mm) {
		
		if (m== null) {
//			m= pattMM.matcher("chr2-:166306374G19A22T26C27A28C29A32G34C35A36");
//			c= m.groupCount();
			m= pattMM.matcher(cs);
		} 
			
		if (!m.find()) 
			return null;
		mm[0]= base.subSequence(pos[3]+ 1+ m.start(1)+1, pos[3]+ 1+ m.end(1)).parseInt();
		//char c1= base.charAt(pos[0]+ mm[0]);
		mm[2]= base.charAt(pos[0]+ mm[0]);	// subst!
		//char c2= base.charAt(pos[3]+ 1+ m.start(1));
		mm[1]= base.charAt(pos[3]+ 1+ m.start(1));
		
		// DEBUG
		String s= m.group(1);
		char c1= base.charAt(pos[0]+ mm[0]);
		char c2= base.charAt(pos[0]+ mm[0]);
		
		return m;
	}

	/**
	 * 
	 * @param cs
	 * @param pos return by parameter, int[0]= pos, int[1]= -1/1 strand
	 * @return String with chromosome
	 */
	public static final char SYMBOL_POSITIVE= '+', SYMBOL_POSITIVE_GEM= 'F', SYMBOL_NEGATIVE= '-', SYMBOL_NEGATIVE_GEM= 'R';
	public static final Pattern PATTY_PLUS_MINUS= Pattern.compile("(\\+|\\-)");
	private static final HashMap<ByteArrayCharSequence, String> mapChrNames= new HashMap<ByteArrayCharSequence, String>();
	public String getMatchPos(ByteArrayCharSequence cs, int[] pos) {
		
		Matcher m= pattMM.matcher(cs);	// find mismatches
		int end= cs.length();
		if (m.find()) 
			end= m.start(1);
		
		m= PATTY_PLUS_MINUS.matcher(cs.subSequence(0, end));
		int sep= -1;
		if (m.find()) {
			sep= m.start(1);
			pos[1]= cs.charAt(sep)== SYMBOL_POSITIVE_GEM? 1: -1;
		} else
			return null; // no sep found, no hit
		
		ByteArrayCharSequence chrCS= cs.subSequence(1, sep);
		String chr= mapChrNames.get(chrCS);
		if (chr== null) {
			chr= chrCS.toString();
			mapChrNames.put(chrCS, chr);
		}
		
		sep+= 2; // chr1+:123
		pos[0]= cs.subSequence(sep, end).parseInt(); //base.subsequence?
		
		return chr;
	}
	
	public int getReadlength() {
		return pos[1]- pos[0]- 1; // 1st AND last pos excluded
	}

	public boolean hasQualities() {
		if (pos[1]== pos[2])
			return false;
		return true;
	}
	
	public ByteArrayCharSequence getQualities() {
		if (!hasQualities())
			return null;
		return base.subSequence(pos[1]+1, pos[2]);
	}
	
	public ByteArrayCharSequence getSequence() {
		return base.subSequence(pos[0]+1, pos[1]);
	}
	
	public ByteArrayCharSequence getName() {
		return base.subSequence(0, pos[0]);
	}
}
