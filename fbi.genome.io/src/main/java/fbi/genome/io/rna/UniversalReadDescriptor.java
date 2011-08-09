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

package fbi.genome.io.rna;

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
		mapSimpleDescriptors.put(DESCRIPTORID_SIMPLE, "#");
		mapSimpleDescriptors.put(DESCRIPTORID_PAIRED, "#/@");
		mapSimpleDescriptors.put(DESCRIPTORID_STRAND_MATE, "#/$/@");
		mapSimpleDescriptors.put(DESCRIPTORID_MATE_STRAND_CSHL, "#/@_strand$[0,1,2]");
		mapSimpleDescriptors.put(DESCRIPTORID_MATE1_SENSE, "#/?");
		mapSimpleDescriptors.put(DESCRIPTORID_MATE2_SENSE, "#/!");
		mapSimpleDescriptors.put(DESCRIPTORID_BARNA, "#/@~[s,a]");
		mapSimpleDescriptors.put(DESCRIPTORID_SIMULATOR, "#:*:$[S,A]/@[1,2]");	
	}
	public static String getDescriptor(CharSequence descriptorID) {
		return mapSimpleDescriptors.get(descriptorID);
	}
	
	public static char TAG_ID= '#', TAG_PAIR= '@', TAG_STRAND= '$', TAG_STRAND_OPT= '~', TAG_MATE1_SENSE= '?', TAG_MATE2_SENSE= '!';
	String[] separators= null;
	boolean[] mandatory= null;
	int posPair= -1, posStrand= -1, posID= -1;
	char symbolMate1= '1', symbolMate2= '2', symbolSense= '1', symbolAsense= '2', symbolNoSense= '0';
	
	public UniversalReadDescriptor() {
	}
	
	public boolean init(String descriptor) {
		
		try {
			String uc= descriptor.toUpperCase();
			if (mapSimpleDescriptors.containsKey(uc))
				descriptor= mapSimpleDescriptors.get(uc);
			
			int pos= 0;
			String[] seps= new String[3];
			boolean[] mand= new boolean[3];
			
			String sep= null;
			for (int i = 0; i < descriptor.length(); ++i) {
				char cc= descriptor.charAt(i);
				if (cc== TAG_ID) {
					seps[pos]= sep;
					mand[pos]= true;
					sep= null;
					posID= pos++;
				} else if (cc== TAG_MATE1_SENSE|| cc== TAG_MATE2_SENSE) {
					seps[pos]= sep;
					mand[pos]= true;
					sep= null;
					posStrand= pos;
					posPair= pos++;
					if (cc== TAG_MATE2_SENSE) {
						symbolSense= symbolMate2;
						symbolAsense= symbolMate1;
					}
				} else if (cc== TAG_PAIR) {
					seps[pos]= sep;
					mand[pos]= true;
					sep= null;
					posPair= pos++;
					i= initSymbolsMate(descriptor, i);
				} else if (cc== TAG_STRAND|| cc== TAG_STRAND_OPT) {
					seps[pos]= sep;
					mand[pos]= (cc== TAG_STRAND);
					sep= null;
					posStrand= pos++;
					i= initSymbolsSense(descriptor, i);
				} else {
					if (sep== null)
						sep= descriptor.substring(i, i+1);
					else
						sep+= descriptor.charAt(i);				
				}
			}
			
			separators= new String[pos];
			System.arraycopy(seps, 0, separators, 0, pos);
			mandatory= new boolean[pos];
			System.arraycopy(mand, 0, mandatory, 0, pos);
			
			boolean consistent= true;
			for (int i = 0; i < mandatory.length; i++) {
				if (!mandatory[i]) {	// alles nach letztem sep/anfang und vor naechstem sep/ende muss mandatory sein
					for (int j = i-1; j >0&& separators[j]== null ; --j) {
						if (!mandatory[j]) {
							consistent= false;
							break;
						}
					}
					for (int j = i+1; consistent&& j < mandatory.length&& separators[j]== null; ++j) {
						if (!mandatory[j])
							consistent= false;
					}
				}
			}
			// ID must be present and delimited
			consistent&= posID>= 0&& mandatory[posID];
			if (posID> 0)
				consistent&= separators[posID- 1]!= null;
			else if (posID+ 1< separators.length)
				consistent&= separators[posID+ 1]!= null;
			
			return consistent;
		} catch (Exception e) {
			return false;
		}
	}

	private int initSymbolsSense(String descriptor, int i) {
		if (i+5>= descriptor.length()|| descriptor.charAt(i+1)!= '['
			|| descriptor.charAt(i+3)!= ',')
			return i;
		if (descriptor.charAt(i+5)== ']') {
			symbolNoSense= '\n';
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
			char cc= (posID==i)?TAG_ID:(posPair==i)?TAG_PAIR:TAG_STRAND;
			char c= cc;
			if (!mandatory[i])
				c= Character.toLowerCase(c);
			sb.append(c);
			if (cc== TAG_PAIR) {
				sb.append('[');
				sb.append(symbolMate1);
				sb.append(',');
				sb.append(symbolMate2);
				sb.append(']');
			} else if (cc== TAG_STRAND) {
				sb.append('[');
				sb.append(symbolNoSense);
				sb.append(',');
				sb.append(symbolSense);
				sb.append(',');
				sb.append(symbolAsense);
				sb.append(']');
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
					continue;
				}
				
				// scan for separator
				boolean ok= true;
				int left= -1;
				int k = separators[i].length()- 1;
				for (; ok&& k>= 0; --k) {
					
					char d= separators[i].charAt(k);
					if (d== '*')
						continue;
					
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
						if (separators[i].charAt(k+ 1)== '*')
							continue;
						ok= false;
						break;
					}
				}
				
				if (k>= 0|| !ok) {
					if (mandatory[i])
						return null;
					else {
						cpos= lastCPos;
						continue;
					}
				}
				
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
		else if (symbolNoSense!= '\n'&& cs.charAt(i)== symbolNoSense)
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
}