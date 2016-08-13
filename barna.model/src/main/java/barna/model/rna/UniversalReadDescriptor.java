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

package barna.model.rna;

import barna.model.constants.Constants;

import java.util.Arrays;
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
    public static String TAG_SENSE= "SENSE";
    public static String TAG_ANTISENSE= "ANTISENSE";
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
        DESCRIPTORID_SENSE= "SENSE",
        DESCRIPTORID_ANTISENSE= "ANTISENSE",
        DESCRIPTORID_PAIRED= "PAIRED",
		DESCRIPTORID_STRAND_MATE= "STRAND_MATE", 
		DESCRIPTORID_MATE_STRAND_CSHL= "MATE_STRAND_CSHL",
		DESCRIPTORID_MATE1_SENSE= "MATE1_SENSE",
		DESCRIPTORID_MATE2_SENSE= "MATE2_SENSE",
		DESCRIPTORID_SIMULATOR= "SIMULATOR",
		DESCRIPTORID_BARNA= "BARNA",
        DESCRIPTORID_CASAVA18 = "CASAVA18";

    public static String DESCRIPTORID_DEFAULT = DESCRIPTORID_SIMPLE;

    static HashMap<String, String> mapSimpleDescriptors= new HashMap<String, String>();
	static {
		mapSimpleDescriptors.put(DESCRIPTORID_SIMPLE, 
				SYMBOL_TAG_LEFT+ 
				TAG_ID+
				SYMBOL_TAG_RIGHT);
        mapSimpleDescriptors.put(DESCRIPTORID_SENSE,
                SYMBOL_TAG_LEFT+
                        TAG_SENSE+
                        SYMBOL_TAG_RIGHT);
        mapSimpleDescriptors.put(DESCRIPTORID_ANTISENSE,
                SYMBOL_TAG_LEFT+
                        TAG_ANTISENSE+
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
        mapSimpleDescriptors.put(DESCRIPTORID_CASAVA18,
                SYMBOL_TAG_LEFT+ TAG_ID+ SYMBOL_TAG_RIGHT+
                        " "+
                        SYMBOL_TAG_LEFT+ TAG_PAIR + SYMBOL_TAG_RIGHT);
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

    public static UniversalReadDescriptor getDefaultDescriptor() {
        UniversalReadDescriptor d = new UniversalReadDescriptor();
        d.init(mapSimpleDescriptors.get(DESCRIPTORID_DEFAULT));
        return d;
    }


    /**
     * Strings that have been provided as separators between the information fields.
     */
	String[] separators= null;

    /**
     * Array of flags that indicates which fields are marked to be required.
     */
	boolean[] mandatory= null;

    /**
     * Position in the field array that specifies the paired-end information.
     */
	int posPair= -1;

    /**
     * Position in the field array that specifies the strand information.
     */
    int posStrand= -1;

    /**
     * Position in the field array that specifies the read identifier.
     */
    int posID= -1;

    /**
     * Flag indicating all mappings of <code>this</code> descriptor are sense
     * with respect to the transcription directionality.
     */
    boolean sense= false;

    /**
     * Flag indicating all mappings of <code>this</code> descriptor are anti-sense
     * with respect to the transcription directionality.
     */
    boolean antisense= false;

    /**
     * Symbol for first mate
     */
    char symbolMate1= '1';
    /**
     * Symbol for second mate
     */
    char symbolMate2= '2';
    /**
     * Symbol for sense directionality
     */
    char symbolSense= '1';
    /**
     * Symbol for anti-sense directionality
     */
    char symbolAsense= '2';
    /**
     * Symbol for cases where the directionality failed to be determined
     */
    char symbolNoSense= '0';

    /**
     * Creates an uninitialized instance.
     */
	private UniversalReadDescriptor() {
	}

    /**
     * Creates an empty descriptor instance for testing.
     * @return empty descriptor instance
     */
    public static UniversalReadDescriptor createTestDescriptor() {
        return new UniversalReadDescriptor();
    }

    /**
     * Creates an instance initialized with the parameters by the provided string,
     * might trigger a <code>RuntimeException</code> during parsing.
     * @param descriptor a string representation of the descriptor attributes
     */
    public UniversalReadDescriptor(String descriptor) {
        init(descriptor);
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

    /**
     * Generates an expression for the rules underlying <code>this</code> descriptor.
     * @return the generated string
     */
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
			if (i== posID) {
                if (sense)
                    sb.append(TAG_SENSE);
                else if (antisense)
                    sb.append(TAG_ANTISENSE);
                else
				    sb.append(TAG_ID);
            } else if (i== posPair&& i== posStrand) {
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

    /**
     * Indicates whether expressions of <code>this</code> descriptor carry strand information.
     * @return <code>true</code> it the descriptor comprises strand information, <code>false</code>
     * otherwise
     */
	public boolean isStranded() {
		return (sense|| antisense|| posStrand>= 0);
	}

    /**
     * Indicates whether expressions of <code>this</code> descriptor carry mate information.
     * @return <code>true</code> it the descriptor comprises mate information, <code>false</code>
     * otherwise
     */
    public boolean isPaired() {
		return (posPair>= 0);
	}
	
	/**
	 * Checks whether <code>this</code> read descriptor is applicable 
	 * to the given character sequence.
	 * @param cs a character sequence
	 * @return <code>true</code> if the character sequence could be 
	 * parsed according to <code>this</code>' rules, <code>false</code>
	 * otherwise
	 */
	public boolean isApplicable(CharSequence cs) {
		try {
			return getAttributes(cs, null)!= null;
		} catch (RuntimeException e) {
			e.printStackTrace();
			return false;
		}
				
	}
	
	/**
	 * Returns the attributes expressed by a given character sequence
	 * with respect to <code>this</code> UniversalReadDescriptor, or
	 * <code>null</code> if the descriptor is not applicable
	 * @param cs the character sequence
	 * @param a an optional <code>Attributes</code> instance for object re-use
	 * @return <code>a</code> initialized with the values of <code>cs</code>,
	 * or <code>null</code> if an parsing error occured
	 */
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
				
				// TODO: security check disabled to ignore 
				// redundant end of Geneva reads
//				if (left!= lastCPos)
//					return null;
				
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

            // set strand for whole readID if necessary
            if (sense)
                a.strand= 1;
            else if (antisense)
                a.strand= 2;

			return a;
			
		} catch (Throwable e) {
			throw new RuntimeException(e);
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
					if (tag.equalsIgnoreCase(TAG_ID)||
                            tag.equalsIgnoreCase(TAG_SENSE)||
                            tag.equalsIgnoreCase(TAG_ANTISENSE)) {
						if (opt)
							throw new RuntimeException("Read ID cannot be optional.");
						posID= pos;
                        if (tag.equalsIgnoreCase(TAG_SENSE))
                            sense= true;
                        else if (tag.equalsIgnoreCase(TAG_ANTISENSE))
                            antisense= true;
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

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        UniversalReadDescriptor that = (UniversalReadDescriptor) o;

        if (antisense != that.antisense) return false;
        if (posID != that.posID) return false;
        if (posPair != that.posPair) return false;
        if (posStrand != that.posStrand) return false;
        if (sense != that.sense) return false;
        if (symbolAsense != that.symbolAsense) return false;
        if (symbolMate1 != that.symbolMate1) return false;
        if (symbolMate2 != that.symbolMate2) return false;
        if (symbolNoSense != that.symbolNoSense) return false;
        if (symbolSense != that.symbolSense) return false;
        if (!Arrays.equals(mandatory, that.mandatory)) return false;
        if (!Arrays.equals(separators, that.separators)) return false;
        if (!this.toString().equals(that.toString())) return false;

        return true;
    }

    @Override
    public int hashCode() {
        int result = separators != null ? Arrays.hashCode(separators) : 0;
        result = 31 * result + (mandatory != null ? Arrays.hashCode(mandatory) : 0);
        result = 31 * result + posPair;
        result = 31 * result + posStrand;
        result = 31 * result + posID;
        result = 31 * result + (sense ? 1 : 0);
        result = 31 * result + (antisense ? 1 : 0);
        result = 31 * result + (int) symbolMate1;
        result = 31 * result + (int) symbolMate2;
        result = 31 * result + (int) symbolSense;
        result = 31 * result + (int) symbolAsense;
        result = 31 * result + (int) symbolNoSense;
        return result;
    }
}