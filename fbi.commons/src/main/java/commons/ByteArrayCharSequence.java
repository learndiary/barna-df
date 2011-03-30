package commons;

import java.util.Arrays;
import java.util.Hashtable;


public class ByteArrayCharSequence implements CharSequence, Comparable<ByteArrayCharSequence> {
	
	public static final byte BYTE_TAB= 9; 
	
	/**
	 * testing
	 * @param args
	 */
	public static void main(String[] args) {
		Hashtable<CharSequence,Object> blash= new Hashtable<CharSequence, Object>();
		ByteArrayCharSequence bas= new ByteArrayCharSequence("Heute");
		blash.put(bas, bas);
		System.out.println(blash.get(bas));
		ByteArrayCharSequence b2= bas.subSequence(1, 3);
		blash.put(b2, b2);
		System.out.println(blash.get(b2));
		ByteArrayCharSequence b3= b2.cloneCurrentSeq();
		blash.put(b3, b3);
		System.out.println(blash.get(b3));

		
/*		String s= "a69";
		for (int i = 0; i < s.length(); i++) 
			System.err.println(s.charAt(i)+ "\t"+ ((byte) s.charAt(i)));
*/		
	}
	
	public static final int defaultSize = 0;
    final static byte[] DigitTens = {
    	48, 48, 48, 48, 48, 48, 48, 48, 48, 48,
    	49, 49, 49, 49, 49, 49, 49, 49, 49, 49,
    	50, 50, 50, 50, 50, 50, 50, 50, 50, 50,
    	51, 51, 51, 51, 51, 51, 51, 51, 51, 51,
    	52, 52, 52, 52, 52, 52, 52, 52, 52, 52,
    	53, 53, 53, 53, 53, 53, 53, 53, 53, 53,
    	54, 54, 54, 54, 54, 54, 54, 54, 54, 54,
    	55, 55, 55, 55, 55, 55, 55, 55, 55, 55,
    	56, 56, 56, 56, 56, 56, 56, 56, 56, 56,
    	57, 57, 57, 57, 57, 57, 57, 57, 57, 57,
    } ; 

    final static byte[] DigitOnes = { 
    	48, 49, 50, 51, 52, 53, 54, 55, 56, 57,
    	48, 49, 50, 51, 52, 53, 54, 55, 56, 57,
    	48, 49, 50, 51, 52, 53, 54, 55, 56, 57,
    	48, 49, 50, 51, 52, 53, 54, 55, 56, 57,
    	48, 49, 50, 51, 52, 53, 54, 55, 56, 57,
    	48, 49, 50, 51, 52, 53, 54, 55, 56, 57,
    	48, 49, 50, 51, 52, 53, 54, 55, 56, 57,
    	48, 49, 50, 51, 52, 53, 54, 55, 56, 57,
    	48, 49, 50, 51, 52, 53, 54, 55, 56, 57,
    	48, 49, 50, 51, 52, 53, 54, 55, 56, 57,
    } ;
    /**
     * All possible chars for representing a number as a String
     */
    final static byte[] digits = {
		48 , 49 , 50 , 51 , 52 , 53 ,
		54 , 55 , 56 , 57 , 97 , 98 ,
		99 , 100 , 101 , 102 , 103 , 104 ,
		105 , 106 , 107 , 108 , 109 , 110 ,
		111 , 112 , 113 , 114 , 115 , 116 ,
		117 , 118 , 119 , 120 , 121 , 122
    };
        
	public byte[] a;

	public int start, end;	// end excluded, heilichs blechle !!

	public ByteArrayCharSequence(byte[] value, String sepString) {
		this(value);
		this.sepString= sepString;
	}
	public ByteArrayCharSequence(byte[] value) {
		a = value;
		this.start = 0;
		this.end = a.length;
		resetFind();
	}
	
	public ByteArrayCharSequence(String s, String sepString) {
		this(s);
		this.sepString= sepString;
	}
	public ByteArrayCharSequence(String s) {
		init(s);
	}
	
	public void init(String s) {
		a= new byte[s.length()];
		for (int i = 0; i < a.length; i++) 
			a[i]= (byte) s.charAt(i);
		start= 0;
		end= s.length();
		resetFind();
	}

	public ByteArrayCharSequence(int capacity, String sepString) {
		this(capacity);
		this.sepString= sepString;
	}
	public ByteArrayCharSequence(int capacity) {
		this(new byte[capacity]);
		this.end= 0;
		resetFind();
	}
		
	public ByteArrayCharSequence(byte[] value, int start, int end, String sepString) {
		this(value, start, end);
		this.sepString= sepString;
	}
	public ByteArrayCharSequence(byte[] value, int start, int end) {
		this(value);
		this.start = start;
		this.end = end;
	}

	//@Override
	public char charAt(int index) {
		return (char) a[start + index];
	}
	
	public void setCharAt(int index, char c) {		
		a[start + index]= (byte) c;
	}
	
	public byte byteAt(int index) {
		return a[start + index];
	}

	//@Override
	public int length() {
		return (end - start);
	}

	/**
	 * The resulting CharSequence operates on the same byte[] !
	 */
	//@Override
	public ByteArrayCharSequence subSequence(int start, int end) {
		// byte[] b= new byte[end-start];
		// System.arraycopy(a, start, b, 0, b.length);
		// return new ByteArrayCharSequence(b);
		return new ByteArrayCharSequence(a, this.start + start, this.start
				+ end);

	}
	
	public int p1, p2, cnt;	// internal find
	private String sepString= "\t";
	/**
	 * 
	 * @param fieldNr 0-based
	 */
	protected void find(int fieldNr) {
		if (fieldNr== cnt)
			return;
		String sepString= this.sepString;
		int sepLen= sepString.length();
		if (fieldNr> cnt) {	// fw
			
			if (p1!= end) {
				for (int i = p2+ 1; cnt < fieldNr&& i< end; ++i) {
					for (int j = 0; j < sepLen; ++j) {
						if (charAt(i)== sepString.charAt(j)) {
							++cnt;
							p1= p2+ 1;
							p2= i;
							// allow only one sep occurrence
							break;
						}
						
					}
				}
				//++p1;
				if (cnt< fieldNr) {
	//				if (cnt== (fieldNr- 1)) {
						p1= (p2== end)? p2: p2+ 1;
						p2= end;
						++cnt;
	//				} else 
	//					p1= p2= end;	// resetFind()
				}
			}
			
		} else {	// rev
			if (p1!= start) {
				//--p1; --p2;
				for (int i = p1- 2; cnt > fieldNr&& i>= start; --i) { 
					for (int j = 0; j < sepLen; ++j) {
						if (charAt(i)== sepString.charAt(j)) {
							--cnt;
							p2= p1- 1;
							p1= i+ 1;
							// allow only one sep occurrence
							break;
						}
					}
				}
				if (cnt> fieldNr) {
	//				if (cnt== (fieldNr+ 1)) {
						p2= (p1== start)? start: p1- 1;
						p1= start;
						--cnt;
	//				} else 
	//					p1= p2= start;	// resetFind();
				} 
			}
		}
		
	}
	
	
	public ByteArrayCharSequence getToken(int fieldNr) {
		if(fieldNr< 0)
			return null;
		
		find(fieldNr);
		if (cnt!= fieldNr)
			return null;
		return subSequence(p1, p2);
	}
	
	public int getTokenInt(int fieldNr) {
		if(fieldNr< 0)
			return Integer.MIN_VALUE;
		
		find(fieldNr);
		if (cnt!= fieldNr) {
			System.err.println(new String(toCharArray()));
			return Integer.MIN_VALUE;
		}
		return parseInt(p1, p2);
	}
	
	/**
	 * @deprecated 1-based
	 * @param nr
	 * @param sep
	 * @return
	 */
	//TODO merge with Tools.CharSequences, duplicated there
	public ByteArrayCharSequence getToken(int nr, char sep) {
		if (nr== 0)
			return null;
		int last= 0, now= 0, cnt= 0;
		for (int i = 0; cnt < nr&& i< length(); i++) 
			if (charAt(i)== sep) {
				++cnt;
				last= now;
				now= i;
			}
		if (last> 0)
			++last;
		if (cnt< nr) {
			if (cnt== nr- 1) {
				if (now+ 1== length())
					return null;
				else
					return subSequence(now+ 1, length());
			}
			else return null;
		}
		if (last== now)
			return null;
		else
			return subSequence(last, now);
	}
	
    final static int [] sizeTable = { 9, 99, 999, 9999, 99999, 999999, 9999999,
        99999999, 999999999, Integer.MAX_VALUE };
    
    protected static int stringSize(int x) {
    	int val= 0;
    	if (x< 0) {
    		val= 1;
    		x= -x;
    	}
		for (int i=0; ; i++)
			if (x <= sizeTable[i])
				return val+ i+ 1;
	}
    
    static int getChars(int i, int endIndex, byte[] buf) {
        int q, r;
        int charPos = endIndex;
        byte sign = 0;

        if (i < 0) { 
            sign = 45;	// '-'
            i = -i;
        }

        // Generate two digits per iteration
        while (i >= 65536) {
            q = i / 100;
        // really: r = i - (q * 100);
            r = i - ((q << 6) + (q << 5) + (q << 2));
            i = q;
            buf [--charPos] = DigitOnes[r];
            buf [--charPos] = DigitTens[r];
        }

        // Fall thru to fast mode for smaller numbers
        // assert(i <= 65536, i);
        for (;;) { 
            q = (i * 52429) >>> (16+3);
            r = i - ((q << 3) + (q << 1));  // r = i-(q*10) ...
            buf [--charPos] = digits [r];
            i = q;
            if (i == 0) break;
        }
        if (sign != 0) {
            buf [--charPos] = sign;
        }
        
        return charPos;
    }

	public boolean replace(int fieldNr, int value) {
		
		if(fieldNr< 0)
			return false;
		
		find(fieldNr);
		if (cnt!= fieldNr)
			return false;
		return replaceCurrField(value);

	}

	public boolean replace(int fieldNr, CharSequence value) {
		
		if(fieldNr< 0)
			return false;
		
		find(fieldNr);
		if (cnt!= fieldNr)
			return false;
		return replaceCurrField(value);

	}
	
	protected boolean replaceCurrField(byte value) {
		if (p1< start|| p1>= end)
			return false;
		
		int diff= 1- (p2- p1);
		if (end+ diff> a.length)
			extend(diff);
		if (diff!=  0) {
			System.arraycopy(a, p2, a, p2+ diff, end- p2);
			end+= diff;
			p2+= diff;
		} 
		assert(p2== p1+ 1);
		a[p1]= value;
		return true;
	}
	
	protected boolean replaceCurrField(int value) {
		if (p1< start|| p1> end)	// for search ends at start/end, p1== p2
			return false;
		
		int digits= value< 0? stringSize(-value)+ 1: stringSize(value);
		int diff= digits- (p2- p1);
		if (end+ diff> a.length)
			extend(diff);
		if (diff!=  0) {
			System.arraycopy(a, p2, a, p2+ diff, end- p2);
			end+= diff;
			p2+= diff;
		} 
		
		int p= getChars(value, p2, a);	// p2 is new end of field
		assert(p== p1);
		return true;
	}
	
	protected boolean replaceCurrField(CharSequence value) {
		if (p1< start|| p1> end)	// for search ends at start/end, p1== p2
			return false;
		
		int digits= value.length();
		int diff= digits- (p2- p1);
		if (end+ diff> a.length)
			extend(diff);
		if (diff!=  0) {
			System.arraycopy(a, p2, a, p2+ diff, end- p2);
			end+= diff;
			p2+= diff;
		} 
		
		int p;
		for (p = p1; p < p1+ digits; ++p) 
			a[p]= (byte) value.charAt(p- p1);
		
		assert(p== p2);
		return true;
	}

	@Override
	public int hashCode() {
		
		int h= 0, len = length();

	    for (int i = 0; i < len; i++) {
	    	h = 31*h + charAt(i);
        }
	    
        return h;
	}
	
	public void append(CharSequence cs) {
		int len= cs.length();
		ensureLength(end, len);
//		int sum= end+ cs.length();
//		if (sum>= a.length)
//			extend(sum- a.length);
		for (int i = 0; i < len; ++i) 
			a[end+ i]= (byte) cs.charAt(i);
		end+= len;
	}
	public void append(byte[] b, int from, int to) {
		int len= to- from;
		ensureLength(end, len);
		System.arraycopy(b, from, a, end, len);
		end+= len;
	}
	public void appendCurrField(byte x) {
		ensureLength(end, 1);
		System.arraycopy(a, p2, a, p2+ 1, end- p2);
		end++;
		a[p2++]= x;
	}
	
	public void appendCurrField(int x) {
		int digits= stringSize(x);
		ensureLength(end, digits);
		System.arraycopy(a, p2, a, p2+ digits, end- p2);
		end+= digits;
		p2+= digits; // end index
		getChars(x, p2, a);
	}
	
	public void appendCurrField(CharSequence value) {
		int digits= value.length();
		ensureLength(end, digits);
		if (digits!=  0) {
			System.arraycopy(a, p2, a, p2+ digits, end- p2);
			end+= digits;
			p2+= digits;
		} 
		
		int p;
		for (p = p1; p < p1+ digits; ++p) 
			a[p]= (byte) value.charAt(p- p1);
		
		assert(p== p2);
	}
	
	public void append(byte x) {
		ensureLength(end, 1);
		a[end++]= x;
	}
	public void append(int x) {
		int len= stringSize(x);
		ensureLength(end, len);
		end+= len;	// end index
		getChars(x, end, a);
	}
	
	public void extend() {
		extend(Math.max(1, (int) (a.length* 0.3)));
	}
	
	public void extend(int x) {
		byte[] b= new byte[a.length+ x];
		System.arraycopy(a, 0, b, 0, a.length);
		a= b;
	}
	
	public ByteArrayCharSequence cloneCurrentSeq() {
		byte[] cloned= getCurrentArray();
		return new ByteArrayCharSequence(cloned, 0, cloned.length, sepString);
	}
	
	public byte[] getCurrentArray() {
		byte[] cloned= new byte[length()];
		System.arraycopy(a, start, cloned, 0, cloned.length);
		return cloned;
	}
	
	@Override
	public String toString() {
		StringBuffer s= new StringBuffer(length());
		for (int i = 0; i < length(); i++) 
			s.append(charAt(i));

		return s.toString();
	}
	
	public char[] toCharArray() {
		char[] cc= new char[length()];
		for (int i = 0; i < cc.length; i++) 
			cc[i]= charAt(i);
		return cc;
	}
	
	public char[] toCharArray(char[] c) {
		if (length()> c.length)
			return toCharArray();
		for (int i = 0; i < length(); i++) 
			c[i]= charAt(i);
		return c;
	}
	
	public boolean startsWith(CharSequence cs) {
		if (length()< cs.length())
			return false;
		for (int i = 0; i < cs.length(); i++) {
			if (cs.charAt(i)!= charAt(i))
				return false;
		}
		return true;
	}
	
	public boolean endsWith(CharSequence cs) {
		if (length()< cs.length())
			return false;
		for (int i = 0; i < cs.length(); i++) {
			if (cs.charAt(cs.length()- i)!= charAt(length()- i))
				return false;
		}
		return true;
	}
	
	public int parseInt() {
		int val= 0;
		for (int i = end-1, pow = 1; i >= start; --i, pow *= 10)
			val += (a[i] - 48) * pow;
		return val;
	}

	public int parseInt(int from, int to) {
		int f= 1;
		if (a[from]== 45) { // '-'
			f= -1;
			++from;
		}
		int val= 0;
		for (int i = to-1, pow = 1; i >= from; --i, pow *= 10)
			val += (a[i] - 48) * pow;
		return (val* f);
	}
	
	/**
	 * Compares two strings lexicographically.
	 * The comparison is based on the Unicode value of each character in
	 * the strings. The character sequence represented by this
	 * <code>String</code> object is compared lexicographically to the
	 * character sequence represented by the argument string. The result is
	 * a negative integer if this <code>String</code> object
	 * lexicographically precedes the argument string. The result is a
	 * positive integer if this <code>String</code> object lexicographically
	 * follows the argument string. The result is zero if the strings
	 * are equal; <code>compareTo</code> returns <code>0</code> exactly when
	 * the {@link #equals(Object)} method would return <code>true</code>.
	 * <p>
	 * This is the definition of lexicographic ordering. If two strings are
	 * different, then either they have different characters at some index
	 * that is a valid index for both strings, or their lengths are different,
	 * or both. If they have different characters at one or more index
	 * positions, let <i>k</i> be the smallest such index; then the string
	 * whose character at position <i>k</i> has the smaller value, as
	 * determined by using the &lt; operator, lexicographically precedes the
	 * other string. In this case, <code>compareTo</code> returns the
	 * difference of the two character values at position <code>k</code> in
	 * the two string -- that is, the value:
	 * <blockquote><pre>
	 * this.charAt(k)-anotherString.charAt(k)
	 * </pre></blockquote>
	 * If there is no index position at which they differ, then the shorter
	 * string lexicographically precedes the longer string. In this case,
	 * <code>compareTo</code> returns the difference of the lengths of the
	 * strings -- that is, the value:
	 * <blockquote><pre>
	 * this.length()-anotherString.length()
	 * </pre></blockquote>
	 *
	 * @param   anotherSeq   the <code>String</code> to be compared.
	 * @return  the value <code>0</code> if the argument string is equal to
	 *          this string; a value less than <code>0</code> if this string
	 *          is lexicographically less than the string argument; and a
	 *          value greater than <code>0</code> if this string is
	 *          lexicographically greater than the string argument.
	 */
	//@Override
	public int compareTo(ByteArrayCharSequence anotherSeq) {
	int len1 = length();
	int len2 = anotherSeq.length();
	int n = Math.min(len1, len2);
	int i = start;
	int j = anotherSeq.start;
	
	for (int m = 0; m < n; m++) {
		if (a[i+m]!= anotherSeq.a[j+m])
			return a[i+m]- anotherSeq.a[j+m];
	}
	return len1-len2;
	
	}
	public int compareTo(CharSequence anotherSeq) {
		int len1 = length();
		int len2 = anotherSeq.length();
		int n = Math.min(len1, len2);
		
		for (int m = 0; m < n; m++) {
			if (charAt(m)!= anotherSeq.charAt(m))
				return charAt(m)- anotherSeq.charAt(m);
		}
		return len1-len2;
		
	}
	
	@Override
	public boolean equals(Object obj) {
		if (obj== null)
			return false;
		CharSequence cs= null;
		try {
			cs= (CharSequence) obj;
		} catch (Exception e) {
			return false;
		}
		if (cs.length()!= length())
			return false;
		for (int i = 0; i < cs.length(); i++) {
			if (cs.charAt(i)!= charAt(i))
				return false;
		}
		return true;
	}

	public static int count(CharSequence cs, char separator) {
		int c= 0, p= 0;
		while (p>= 0) {
			p= indexOf(cs, separator, p, cs.length());
			++c;
		}
		--p;
		return c;
	}
	public static int indexOf(CharSequence cs, char separator, int from, int to) {
		for (int i = from; i < to; ++i) {
			if (cs.charAt(i)== separator)
				return i;
		}
		return -1;
	}
	
	public static ByteArrayCharSequence cloneSequence(ByteArrayCharSequence src, ByteArrayCharSequence target) {
		if (target== null) 
			target= new ByteArrayCharSequence(src.length());
		if (target.a.length< src.length())
			target.a= new byte[src.length()];
		System.arraycopy(src.a, src.start, target.a, 0, src.length());
		target.start= src.start;
		target.end= src.end;
		return target;
	}

	public int countTokens(char tab) {
		int cnt= 0;
		for (int i = start; i < length(); i++) 
			if (charAt(i)== tab)
				++cnt;
		if (length()> 0&& charAt(length()- 1)!= tab)
			++cnt;
		return cnt;
	}

	public void resetFind() {
		p1= start- 1;
		p2= start- 1;
		cnt= -1;
	}
	
	public void reset() {
		start= end= 0;
		resetFind();
	}

	public static void reverseComplement(ByteArrayCharSequence cs) {
		reverseComplement(cs, 0, cs.length());
	}
	public static void reverseComplement(ByteArrayCharSequence cs, int from, int to) {
		complement(cs.a, cs.start+ from, cs.start+ to);
		reverse(cs.a, cs.start+ from, cs.start+ to);
	}

	public static void reverse(ByteArrayCharSequence cs, int from, int to) {
		reverse(cs.a, cs.start+ from, cs.start+ to);
	}
	public static void reverse(byte[] a, int from, int to) {

		// adapted from commons.ArrayUtils
        if (a == null) {
            return;
        }
        int i = from;
        int j = to - 1;
        byte tmp;
        while (j > i) {
            tmp = a[j];
            a[j] = a[i];
            a[i] = tmp;
            --j;
            ++i;
        }

	}

	public static final byte[] CHARS_NORMAL= new byte[]   {45, 65, 67, 71, 75, 77, 78, 82, 84, 85, 88, 89,		// -ACGKMNRTUXY
														   97, 99, 103, 107, 109, 110, 114, 116, 117, 120, 121},// acgkmnrtuxy
							   CHARS_REVERSED= new byte[] {45, 84, 71, 67, 77, 75, 78, 89, 65, 65, 88, 82,		// -TGCMKNYAAXR
														   116, 103, 99, 109, 107, 110, 121, 97, 97, 120, 114};	// tgcmknyaaxr		
	
	public static void complement(ByteArrayCharSequence cs) {
		complement(cs.a, cs.start, cs.end);
	}
	public static void complement(ByteArrayCharSequence cs, int from, int to) {
		complement(cs.a, cs.start+ from, cs.start+ to);
	}
	public static void complement(byte[] a, int from, int to) {
		
		int p;
		for (int i = from; i < to; ++i) {
			p= Arrays.binarySearch(CHARS_NORMAL, a[i]);
			if (p< 0) 
				System.err.println("Complement: unknown symbol "+ ((char) a[i])+ " ("+ a[i]+")");
			a[i]= CHARS_REVERSED[p];
		}
	}

	public void ensureLength(int from, int len) {
		int diff= (from+ len)- a.length;
		if (diff> 0) 
			extend(diff);
	}

	public void toUpperCase() {
		toUpperCase(a, start, end);
	}
	
	public void toUpperCase(int from, int to) {
		toUpperCase(a, start+ from, start+ to);
	}
	
	public void toLowerCase() {
		toLowerCase(a, start, end);
	}
	
	public void toLowerCase(int from, int to) {
		toLowerCase(a, start+ from, start+ to);
	}
	
	public static void toUpperCase(byte[] b, int from, int to) {
		for (int i = from; i < to; i++) {
			if (b[i]>= 97&& b[i]<= 122)	// 'a'..'z'
				b[i]-= 32;
		}
	}

	public static void toLowerCase(byte[] b, int from, int to) {
		for (int i = from; i < to; i++) {
			if (b[i]>= 65&& b[i]<= 90)	// 'A'..'Z'
				b[i]+= 32;
		}
	}
}
