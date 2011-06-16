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
		public byte strand, flag;
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
	
	static HashMap<String, String> mapSimpleDescriptors= new HashMap<String, String>();
	static {
		mapSimpleDescriptors.put("SIMPLE", "#");
		mapSimpleDescriptors.put("PAIRED", "#/@");
		mapSimpleDescriptors.put("MATE_STRAND", "#/$/@");
		mapSimpleDescriptors.put("MATE_STRAND-CSHL", "#/@_strand$[0,1,2]");
		mapSimpleDescriptors.put("STRAND_MATE", "#/$/@");
		mapSimpleDescriptors.put("MATE1_SENSE", "#/?");
		mapSimpleDescriptors.put("MATE2_SENSE", "#/!");
		mapSimpleDescriptors.put("BARNA", "#/@~[s,a]");
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
	
	public static void main(String[] args) {
		
		UniversalReadDescriptor descriptor= new UniversalReadDescriptor();
		String barnaDesc= "#/@~[s,a]";
		String gingerasDesc= "#/$/@";
		String oldgingerDesc= "#/@_strand$";
		String desc= mapSimpleDescriptors.get("MATE2_SENSE");
		System.out.println(desc);		
		boolean b= descriptor.init(desc);
		System.out.println(desc);
		System.out.println(b);
		if (b)
			System.out.println(descriptor);
		else
			System.exit(-1);
		System.out.println();
		
		String gingerasID1= "BILLIEHOLIDAY:5:100:1000:1190/1/2";
		String gingerasID2= "BILLIEHOLIDAY:5:100:1000:1190/2/2";
		String gingerasID3= "BILLIEHOLIDAY:5:100:1000:1190/2";
		String gingerasID4= "BILLIEHOLIDAY:5:100:1000:1190/2/s";
		
		String barnaID1= "BILLIEHOLIDAY:5:100:1000:1190/1s";
		String barnaID2= "BILLIEHOLIDAY:5:100:1000:1190/1a";
		String barnaID3= "BILLIEHOLIDAY:5:100:1000:1190/2";
		String barnaID4= "BILLIEHOLIDAY:5:100:1000:1190/2a";
		String barnaID5= "BILLIEHOLIDAY:5:100:1000:1190/1";
		
		String oldGingerID1= "BILLIEHOLIDAY:5:100:1000:1190/1_strand2";
		String oldGingerID2= "BILLIEHOLIDAY:5:100:1000:1190/1_strand0";
		
		String newGingerCombined= "MARILYN_0005:7:1:2804:1011#0/1";
		
		String id= newGingerCombined;
		Attributes a= descriptor.getAttributes(id, null);
		System.out.println(id);
		b= (a!= null);
		System.out.println(b);
		if (b) 
			System.out.println(a);
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
		
		try {
		
			// from start
			int cpos= 0, lastPos= 0, lastCPos= 0;
			// TEMPORARILY DEACTIVATED
//			for (int i = 0; i < posID; i++) {
//				if (separators[i]!= null) {
//					int j= cpos, k= 0;
//					for(;j< cs.length()&& k< separators[i].length()
//						&& cs.charAt(j++) == separators[i].charAt(k++););
//					if (j>= cs.length()) 
//						return null;
//					if (k== separators[i].length()) {
//						int dc= cpos- lastCPos, dp= i- lastPos;
//						assert(dc<= dp);
//						boolean opt= (dc< dp)?false:true;
//						for (int h = lastPos, p= lastCPos; h < i; h++) {
//							if (mandatory[h]|| opt) {
//								if (h== posPair){
//									if(!setPair(a, cs, p++))
//										return null;
//								} else if (h== posStrand)
//									setStrand(a, cs, p++);
//							}
//						}
//					} 
//					cpos+= k;
//					if (k== separators[i].length()) { 
//						lastCPos= cpos;
//						lastPos= i;
//					}
//				}
//			}
			int posLeft= lastPos;
				
			// from end
			cpos= cs.length()- 1;
			lastCPos= cpos+ 1;
			lastPos= separators.length;
			for (int i = separators.length- 1; i > posID;) {
				if (separators[i]== null) {
					--i;
					continue;
				}
				int j= cpos, k= separators[i].length()- 1;
				for(;j>= 0&& k>= 0;--j, --k) {
					if (cs.charAt(j)!= separators[i].charAt(k))
						break;
				}
				
				if (j< 0) 
					return null;
				cpos-= separators[i].length()- (k< 0?0:k);
				if (k< 0) {	// sep found
					int p= cpos+ 1+ separators[i].length();
					int dc= lastCPos- p, dp= lastPos- i;
					assert(dc<= dp);
					boolean opt= (dc< dp);
					for (int h = i; h < lastPos; ++h) {
						if (mandatory[h]|| (!opt)) {
							if (h== posPair) {
								if(!setPair(a, cs, p))
									return null;
							} 
							// 20101222, not else for pair-stranded 
							if (h== posStrand) {
								if (!setStrand(a, cs, p))
									return null;
							}
							p++;
						}
					}
					lastCPos= cpos+ 1;
					lastPos= i;
					--i;
				}
			}
			
			a.id= cs.subSequence(posLeft, cpos);
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