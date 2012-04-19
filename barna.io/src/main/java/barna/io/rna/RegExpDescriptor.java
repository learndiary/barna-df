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

package barna.io.rna;

import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.regex.PatternSyntaxException;

public class RegExpDescriptor extends Descriptor {

	public static void main(String[] args) {
		String regexp= ".*/([12])/([120])$";
		Pattern paty= Pattern.compile(regexp);
		String[] reads= new String[] {
				"chr1\t10033\t10109\tPAN:2:27:1091:1987/2_strand1",
				"chr1    10033   10109   PRESLEY:3:67:1348:1236/1_strand0",
				"chr1    10033   10109   PRESLEY:3:67:1348:1236/1/0"
		};
		System.err.println("Regexp: "+ regexp);
		for (int i = 0; i < reads.length; i++) {
			Matcher matt= paty.matcher(reads[i]);
			boolean bool= matt.matches();
			System.err.println("Read "+(i+1)+": "+ bool);
			if (!bool)
				continue;
			System.err.println(matt.group(1)+"\t"+matt.group(2));
		}
	}
	
	public static final char CHAR_BR_OPEN= '(', CHAR_BR_CLOSE= ')';
	static final String AMBIGUITY_CHARS= "^$.\\*+{}[]";
	static String attributes= null, regExp= null;
	static Pattern pati= null;
	static Matcher matt= null;
	static String[] anchors= null;
	static char char_mateUnknown= '\n', char_mate1, char_mate2, char_senseUnknown= '\n', char_sense, char_asense;
	static int idPos= -1, pairedPos= -1, strandedPos= -1;
	
	public boolean init(String newAttributes, String newRegExp) {
		attributes= newAttributes;
		regExp= newRegExp;
		if (attributes== null|| regExp== null)
			return false;
		
		try {
			pati= Pattern.compile(regExp);
		} catch (PatternSyntaxException e) {
			return false;
		}
		
		int p= 0;
		//anchors= new String[attributes.length()+ 1];
		for (int i = 0; i < attributes.length(); ++i) {			
			int p1= regExp.indexOf('(', p+ 1);
			if (p1< 0)
				return false;
			//anchors[i]= parseAnchor(regExp.substring(p, p1));
			int p2= regExp.indexOf(')', p1+ 1);
			if (p2< 0)
				return false;
			if (attributes.charAt(i)== CHAR_ID_STRANDED) {
				if (idPos>= 0)
					return false;
				idPos= i;
				continue;
			}
			String s= regExp.substring(p1+ 1, p2);
			if (s.length()< 4|| s.length()> 5|| !(s.startsWith("[")&& s.endsWith("]")))
				return false;
			s= s.substring(1, s.length()- 1);
			for (int j = 0; j < s.length(); j++) {
				if (AMBIGUITY_CHARS.indexOf(s.charAt(j))>= 0)
					return false;
			}
			if (attributes.charAt(i)== CHAR_ID_PAIRED) {
				if (pairedPos>= 0)
					return false;
				int x= 0;
				if (s.length()== 3)
					char_mateUnknown= s.charAt(x++);
				char_mate1= s.charAt(x++);
				char_mate2= s.charAt(x++);
				pairedPos= i;
			} else if (attributes.charAt(i)== CHAR_ID_STRANDED) {
				if (strandedPos>= 0)
					return false;
				int x= 0;
				if (s.length()== 3)
					char_senseUnknown= s.charAt(x++);
				char_sense= s.charAt(x++);
				char_asense= s.charAt(x++);
				strandedPos= i;
			} else 
				return false;
			
			p= p2+ 1;
		}
		if (regExp.indexOf(CHAR_BR_OPEN, p)>= 0|| regExp.indexOf(CHAR_BR_CLOSE)>= 0)
			return false;
		//anchors[anchors.length- 1]= parseAnchor(regExp.substring(p, regExp.length()));
		return true;
	}
	
	/**
	 * trims from end, keeps only non-ambigous chars
	 * @param s
	 * @return
	 * @deprecated 
	 */
	private String parseAnchor(String s) {
		int i= s.length()- 1;
		for (; i>= 0; --i) {
			char c= s.charAt(i);
			int p= AMBIGUITY_CHARS.indexOf(c);
			if (p>= 0)
				break;
		}
		i+= 2; // for '\s';
		if (i>= s.length())
			return null;
		return s.substring(i);
	}
	
	@Override
	public String toString() {		
		return "RegExp "+attributes+" "+regExp;
	}
	
	public int getMode(CharSequence seq, int[] fromTo) {
		
		matt= pati.matcher(seq);
		if (!matt.matches())
			return -1;

		// ID
		fromTo[0]= matt.start(idPos+ 1);
		fromTo[1]= matt.end(idPos+ 1);
		int res= 0;
		if (pairedPos>= 0) {
			int start= matt.start(pairedPos+ 1),
				end= matt.end(pairedPos+ 1);
			if (end!= start+ 1)
				return -1;
			char c= seq.charAt(start); 
			if (c== char_mate1) 
				res|= MODE_MATE1;
			else if (c== char_mate2)
				res|= MODE_MATE2;
			else if (c!= char_mateUnknown)
				return -1;
		}
		if (strandedPos>= 0) {
			int start= matt.start(strandedPos+ 1),
				end= matt.end(strandedPos+ 1);
			if (end!= start+ 1)
				return -1;
			char c= seq.charAt(start); 		
			if (c== char_sense) 
				res|= MODE_MATE1;
			else if (c== char_asense)
				res|= MODE_MATE2;
			else if (c!= char_senseUnknown)
				return -1;
		}
		return res;
	}

	@Override
	public boolean allowsPairs() {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public boolean allowsStranded() {
		// TODO Auto-generated method stub
		return false;
	}
}
