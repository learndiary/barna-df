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

package barna.genome.io.rna;

import barna.genome.model.constants.Constants;

import java.util.HashMap;

public class UniversalReadDescriptor {

	public class Attributes {
		public byte strand= 0, flag= 0;
		public CharSequence id;
		@Override
		public String toString() {
			StringBuilder sb= new StringBuilder(id);
			if (UniversalReadDescriptor.this.isPaired()) {
				sb.append(" mate ");
				sb.append(Byte.toString(flag));
			}
			if (UniversalReadDescriptor.this.isStranded()) {
				sb.append(" sense ");
				sb.append(Byte.toString(strand));
			}
			return sb.toString();
		}
	}
	
	public static String TAG_ID= "ID";
	public static String TAG_PAIR= "MATE"; 
	public static String TAG_STRAND= "STRAND";
	public static String TAG_MATE1SENSE= "MATE1SENSE";
	public static String TAG_MATE2SENSE= "MATE2SENSE";

	public static char SYMBOL_TAG_LEFT= '{';
	public static char SYMBOL_TAG_RIGHT= '}';
	public static char SYMBOL_OPT_LEFT= '(';
	public static char SYMBOL_OPT_RIGHT= ')';
	public static char SYMBOL_SET_LEFT= '[';
	public static char SYMBOL_SET_RIGHT= ']';
	public static char SYMBOL_ESCAPE= '\\';
	public static char SYMBOL_STAR= '*';
	public static char SYMBOL_QUESTION= '?';

	public static String 
		DESCRIPTORID_SIMPLE= "SIMPLE", 
		DESCRIPTORID_PAIRED= "PAIRED", 
		DESCRIPTORID_STRAND_MATE= "STRAND_MATE", 
		DESCRIPTORID_MATE_STRAND_CSHL= "MATE_STRAND_CSHL",
		DESCRIPTORID_MATE1_SENSE= "MATE1_SENSE",
		DESCRIPTORID_MATE2_SENSE= "MATE2_SENSE",
		DESCRIPTORID_SIMULATOR= "SIMULATOR",
		DESCRIPTORID_BARNA= "BARNA";
	
	static HashMap<String, String> mapSimpleDescriptors= new HashMap<String, String>();
	static {
		mapSimpleDescriptors.put(DESCRIPTORID_SIMPLE, 
				SYMBOL_TAG_LEFT+ 
				TAG_ID+
				SYMBOL_TAG_RIGHT);
		mapSimpleDescriptors.put(DESCRIPTORID_PAIRED, 
				SYMBOL_TAG_LEFT+ TAG_ID+ SYMBOL_TAG_RIGHT+
				"/"+
				SYMBOL_TAG_LEFT+ TAG_PAIR+ SYMBOL_TAG_RIGHT+
				SYMBOL_SET_LEFT+ "1,2"+ SYMBOL_SET_RIGHT);
		mapSimpleDescriptors.put(DESCRIPTORID_STRAND_MATE, 
				SYMBOL_TAG_LEFT+ TAG_ID+ SYMBOL_TAG_RIGHT+
				"/"+
				SYMBOL_TAG_LEFT+ TAG_STRAND+ SYMBOL_TAG_RIGHT+
				SYMBOL_SET_LEFT+ "0,1,2"+ SYMBOL_SET_RIGHT+
				"/"+
				SYMBOL_TAG_LEFT+ TAG_PAIR+ SYMBOL_TAG_RIGHT+
				SYMBOL_SET_LEFT+ "1,2"+ SYMBOL_SET_RIGHT);
		mapSimpleDescriptors.put(DESCRIPTORID_MATE_STRAND_CSHL, 
				SYMBOL_TAG_LEFT+ TAG_ID+ SYMBOL_TAG_RIGHT+
				"/"+
				SYMBOL_TAG_LEFT+ TAG_PAIR+ SYMBOL_TAG_RIGHT+
				SYMBOL_SET_LEFT+ "1,2"+ SYMBOL_SET_RIGHT+
				"_strand"+
				SYMBOL_TAG_LEFT+ TAG_STRAND+ SYMBOL_TAG_RIGHT+
				SYMBOL_SET_LEFT+ "0,1,2"+ SYMBOL_SET_RIGHT);
		mapSimpleDescriptors.put(DESCRIPTORID_MATE1_SENSE, 
				SYMBOL_TAG_LEFT+ TAG_ID+ SYMBOL_TAG_RIGHT+
				"/"+
				SYMBOL_TAG_LEFT+ TAG_MATE1SENSE+ SYMBOL_TAG_RIGHT+
				SYMBOL_SET_LEFT+ "1,2"+ SYMBOL_SET_RIGHT);
		mapSimpleDescriptors.put(DESCRIPTORID_MATE2_SENSE, 
				SYMBOL_TAG_LEFT+ TAG_ID+ SYMBOL_TAG_RIGHT+
				"/"+
				SYMBOL_TAG_LEFT+ TAG_MATE2SENSE+ SYMBOL_TAG_RIGHT+
				SYMBOL_SET_LEFT+ "1,2"+ SYMBOL_SET_RIGHT);
		mapSimpleDescriptors.put(DESCRIPTORID_BARNA, 
				SYMBOL_TAG_LEFT+ TAG_ID+ SYMBOL_TAG_RIGHT+
				"/"+
				SYMBOL_TAG_LEFT+ TAG_PAIR+ SYMBOL_TAG_RIGHT+
				SYMBOL_SET_LEFT+ "1,2"+ SYMBOL_SET_RIGHT+
				
				SYMBOL_OPT_LEFT+ TAG_STRAND+ SYMBOL_OPT_RIGHT+
				SYMBOL_SET_LEFT+ "s,a"+ SYMBOL_SET_RIGHT);
		mapSimpleDescriptors.put(DESCRIPTORID_SIMULATOR, 
				SYMBOL_TAG_LEFT+ TAG_ID+ SYMBOL_TAG_RIGHT+
				":"+
				SYMBOL_TAG_LEFT+ TAG_STRAND+ SYMBOL_TAG_RIGHT+
				SYMBOL_SET_LEFT+ "S,A"+ SYMBOL_SET_RIGHT+
				"/"+
				SYMBOL_TAG_LEFT+ TAG_PAIR+ SYMBOL_TAG_RIGHT+
				SYMBOL_SET_LEFT+ "1,2"+ SYMBOL_SET_RIGHT);	
	}
//	static {
//		mapSimpleDescriptors.put(DESCRIPTORID_SIMPLE, "#");
//		mapSimpleDescriptors.put(DESCRIPTORID_PAIRED, "#/@");
//		mapSimpleDescriptors.put(DESCRIPTORID_STRAND_MATE, "#/$/@");
//		mapSimpleDescriptors.put(DESCRIPTORID_MATE_STRAND_CSHL, "#/@_strand$[0,1,2]");
//		mapSimpleDescriptors.put(DESCRIPTORID_MATE1_SENSE, "#/?");
//		mapSimpleDescriptors.put(DESCRIPTORID_MATE2_SENSE, "#/!");
//		mapSimpleDescriptors.put(DESCRIPTORID_BARNA, "#/@~[s,a]");
//		mapSimpleDescriptors.put(DESCRIPTORID_SIMULATOR, "#:*:$[S,A]/@[1,2]");	
//	}
	public static String getDescriptor(CharSequence descriptorID) {
		return mapSimpleDescriptors.get(descriptorID);
	}
	
	public static HashMap<String, String> getMapSimpleDescriptors() {
		return mapSimpleDescriptors;
	}


	String[] separators= null;
	boolean[] mandatory= null;
	int posPair= -1, posStrand= -1, posID= -1;
	
	char symbolMate1= '1', symbolMate2= '2', symbolSense= '1', symbolAsense= '2', symbolNoSense= '0';
	
	public UniversalReadDescriptor() {
	}
	
	private int initSymbolsSense(String descriptor, int i) {
		if (i+5>= descriptor.length()|| descriptor.charAt(i+1)!= '['
			|| descriptor.charAt(i+3)!= ',')
			return i;
		if (descriptor.charAt(i+5)== ']') {
			symbolNoSense= Constants.NL;
			symbolSense= descriptor.charAt(i+ 2);
			symbolAsense= descriptor.charAt(i+ 4);
			return (i+5);
		} else if(i+7< descriptor.length()&& descriptor.charAt(i+5)== ','&& descriptor.charAt(i+7)== ']') {
			symbolNoSense= descriptor.charAt(i+ 2);
			symbolSense= descriptor.charAt(i+ 4);
			symbolAsense= descriptor.charAt(i+ 6);
			return (i+7);
		} else
			return i;
	}

	private int initSymbolsMate(String descriptor, int i) {
		
		if (i+5>= descriptor.length()|| descriptor.charAt(i+1)!= '['
			|| descriptor.charAt(i+3)!= ','|| descriptor.charAt(i+5)!= ']')
			return i;
		symbolMate1= descriptor.charAt(i+ 2);
		symbolMate2= descriptor.charAt(i+ 4);
		return (i+5);
	}
	
	@Override
	public String toString() {

		StringBuilder sb= new StringBuilder();
		for (int i = 0; i < separators.length; i++) {
			if (separators[i]!= null)
				sb.append(separators[i]);
			if (mandatory[i])
				sb.append(SYMBOL_TAG_LEFT);
			else
				sb.append(SYMBOL_OPT_LEFT);
			if (i== posID)
				sb.append(TAG_ID);
			else if (i== posPair&& i== posStrand) {
				if (symbolMate1== symbolSense)
					sb.append(TAG_MATE1SENSE);
				else
					sb.append(TAG_MATE2SENSE);
			} else if (i== posPair) 
				sb.append(TAG_PAIR);
			else if (i== posStrand) 
				sb.append(TAG_STRAND);
			
			if (mandatory[i])
				sb.append(SYMBOL_TAG_RIGHT);
			else
				sb.append(SYMBOL_OPT_RIGHT);

			// append set
			if (i== posPair) {
				sb.append(SYMBOL_SET_LEFT);
				sb.append(symbolMate1);
				sb.append(',');
				sb.append(symbolMate2);
				sb.append(SYMBOL_SET_RIGHT);
			// if pair includes mate info, the following is not executed
			} else if (i== posStrand) {
				sb.append(SYMBOL_SET_LEFT);
				if (symbolNoSense!= Constants.NL) {
					sb.append(symbolNoSense);
					sb.append(',');
				}
				sb.append(symbolSense);
				sb.append(',');
				sb.append(symbolAsense);
				sb.append(SYMBOL_SET_RIGHT);
			}
		}
		
		return sb.toString();
	}
	
	public boolean isStranded() {
		return (posStrand>= 0);
	}
	
	public boolean isPaired() {
		return (posPair>= 0);
	}
	
	public Attributes getAttributes(CharSequence cs, Attributes a) {
		
		if (a== null)
			a= new Attributes();
		
		else {	// re-init
			a.strand= 0;
			a.flag= 0;
			a.id= null;
		}
		
		try {
			
			// parse from end
			int cpos= cs.length()- 1;
			int lastCPos= cpos;
			for (int i = separators.length- 1; i >= 0; --i) {
				if (separators[i]== null) {
					int parsed= 0;
					if (i== posPair) {
						if(!setPair(a, cs, cpos)) {
							if (mandatory[i])
								return null;
						} else
							parsed= 1;
					} else { 
						if (i== posStrand) {
							if (!setStrand(a, cs, cpos)) {
								if (mandatory[i])
									return null;
							} else
								parsed= 1;
						}				
					}
					cpos-= parsed;
					lastCPos= cpos;
					continue;
				}
				
				// scan for separator
				boolean ok= true;
				int left= -1;
				int k = separators[i].length()- 1;
				int wildcards= 0;
				for (; ok&& k>= 0; --k) {
					
					char d= separators[i].charAt(k);
					if (d== SYMBOL_STAR) {
						wildcards= cpos;
						continue;
					}	
					if (d== SYMBOL_QUESTION) {
						++wildcards;
						continue;
					}
					
					for (; cpos>= 0; --cpos) {
						char c= cs.charAt(cpos);
						if (c== d) {
							--cpos;	// set cpos always to next position
							break;
						}
						if (k== separators[i].length()- 1) {
							left= cpos;
							continue; // start of target chars
						}

						assert(k< separators[i].length()- 1);						
						if (wildcards> 0) {
							--wildcards;
							continue;
						}
						ok= false;
						break;
					}
					wildcards= 0;
				}
				
				if (k>= 0|| !ok) {
					if (mandatory[i])
						return null;
					else {
						cpos= lastCPos;
						continue;
					}
				}
				
				if (left!= lastCPos)
					return null;
				if (i== posPair) {
					if(!setPair(a, cs, left)) {
						if (mandatory[i])
							return null;
					}
				} 
				// 20101222, not else for pair-stranded 
				if (i== posStrand) {
					if (!setStrand(a, cs, left)) {
						if (mandatory[i])
							return null;
					}
				}
				//--cpos;
				lastCPos= cpos;
				
			}
			
			a.id= cs.subSequence(0, cpos+ 1);	// 20110808: cpos+ 1 to catch last char
			return a;
			
		} catch (Throwable e) {
			e.printStackTrace();
			return null;
		}
	}

	private boolean setStrand(Attributes a, CharSequence cs, int i) {
		if (cs.charAt(i)== symbolSense)
			a.strand= 1;
		else if (cs.charAt(i)== symbolAsense)
			a.strand= 2;
		else if (symbolNoSense!= Constants.NL&& cs.charAt(i)== symbolNoSense)
			a.strand= 0;
		else
			return false;
		return true;
	}

	private boolean setPair(Attributes a, CharSequence cs, int i) {
		if (cs.charAt(i)== symbolMate1)
			a.flag= 1;
		else if (cs.charAt(i)== symbolMate2)
			a.flag= 2;
		else
			return false;
		return true;
	}

	public Attributes createAttributes() {		
		return new Attributes();
	}

	public void init(String descriptor) {
		
		String uc= descriptor.toUpperCase();
		if (mapSimpleDescriptors.containsKey(uc))
			descriptor= mapSimpleDescriptors.get(uc);
		
		int pos= 0;
		String[] seps= new String[3];
		boolean[] mand= new boolean[3];
		
		String sep= null;
		boolean esc= false;
		for (int i = 0; i < descriptor.length(); ++i) {
			char cc= descriptor.charAt(i);
			if (cc=='\\') 
				esc= !esc;
			else {
				if (cc== SYMBOL_TAG_LEFT|| cc== SYMBOL_OPT_LEFT) {
					boolean opt= cc== SYMBOL_OPT_LEFT; 
					char dd= (opt? SYMBOL_OPT_RIGHT: SYMBOL_TAG_RIGHT);
					int j = i;
					while (++j< descriptor.length()&& descriptor.charAt(j)!= dd);
					if (j>= descriptor.length())
						throw new RuntimeException("No closing symbol found for '"+ cc+ "' at position "+ i+ ".");
					
					String tag= descriptor.substring(i+1, j);
					if (tag.equalsIgnoreCase(TAG_ID)) {
						if (opt)
							throw new RuntimeException("Read ID cannot be optional.");
						posID= pos;
					} else if (tag.equalsIgnoreCase(TAG_STRAND)) {
						posStrand= pos;
						j= initSymbolsSense(descriptor, j);
					} else if (tag.equalsIgnoreCase(TAG_MATE1SENSE)|| tag.equalsIgnoreCase(TAG_MATE2SENSE)) {
						posStrand= pos;
						posPair= pos;
						if (tag.equalsIgnoreCase(TAG_MATE2SENSE)) {
							symbolSense= symbolMate2;
							symbolAsense= symbolMate1;
						}
						// combined coding not allowed to change default symbols
					} else if (tag.equalsIgnoreCase(TAG_PAIR)) {
						if (opt)
							throw new RuntimeException("Mate pair information cannot be optional.");
						posPair= pos;
						j= initSymbolsMate(descriptor, j);
					} else
						throw new RuntimeException("Tag '"+ tag+ "' unknown.");
					seps[pos]= sep;
					mand[pos++]= !opt;
					
					sep= null;
					i= j;
					
				} else  {	// end tag block
				
					// lazily append to separator
					if (!esc)
						if (sep== null)
							sep= descriptor.substring(i, i+1);
						else
							sep+= descriptor.charAt(i);				
				}
			} // end no escape
		} // end iterating descriptor
		
		separators= new String[pos];
		System.arraycopy(seps, 0, separators, 0, pos);
		mandatory= new boolean[pos];
		System.arraycopy(mand, 0, mandatory, 0, pos);
		
		// ID must be present and delimited
		if (posID< 0|| !mandatory[posID])
			throw new RuntimeException("Read descriptor must contain an ID which cannot be optional.");
		else {
			if (posID> 0&& separators[posID- 1]== null)
				throw new RuntimeException("ID (pos "+ posID+ ") lacks left separator.");
			else if (posID+ 1< separators.length&& separators[posID+ 1]== null)
				throw new RuntimeException("ID (pos "+ posID+ ") lacks right separator.");
		}
		
		// everything after last separator and before 
		// the beginning of the next separator must NOT 
		// be optional to avoid ambiguities! 
		for (int i = 0; i < mandatory.length; i++) {
			if (!mandatory[i]) {	
				for (int j = i-1; j >0&& separators[j]== null ; --j) {
					if (!mandatory[j])
						throw new RuntimeException("Ambigous concatenation of optional tags "+ i+ " and "+ j);
				}
				for (int j = i+1; j < mandatory.length&& separators[j]== null; ++j) {
					if (!mandatory[j])
						throw new RuntimeException("Ambigous concatenation of optional tags "+ i+ " and "+ j);
				}
			}
		}
	}
}