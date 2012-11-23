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

package barna.commons.utils;

import java.text.SimpleDateFormat;
import java.util.*;

/**
 * String constants and utilities
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 * @author Micha Sammeth (gmicha@googlemail.com)
 */
public class StringUtils {
    /**
     * OK string
     */
    public static final String OK = "OK";
    /**
     * Error string
     */
    public static final String ERROR = "ERROR";
    /**
     * Biological characters as spotted in refseq intron annotation of human
     */
    private static final String BIOLOGICAL_CHARS = "ACGTN-acgtnMmKkXxRrYyWwSsBbVvHhDd";
    /**
     * The complements for biological chars
     */
    private static final String COMPLEMENT_CHARS = "TGCAN-tgcanKkMmXxYyRrSsWwVvBbDdHh";


    /**
     * Create a timestamp for the current time
     *
     * @return timestamp timestamp
     */
    public static String getTimestampID() {
        SimpleDateFormat format = new SimpleDateFormat("yyMMddHHmmssSSSS");
        return format.format(new Date());
    }

    /**
     * Reverse the given string
     *
     * @param in the string
     * @return reverse the reverse string
     */
    public static String reverse(String in) {
        StringBuilder buffy = new StringBuilder(in.length());
        for (int i = in.length() - 1; i >= 0; i--) {
            buffy.append(in.charAt(i));
        }
        return buffy.toString();
    }

    /**
     * Returns the reverse complement of a biological sequence. A RuntimeException is thrown
     * of one of the characters could not be translated to its complement.
     *
     * @param sequence the biological sequence
     * @return complement the complement sequence
     * @throws RuntimeException in case a character could not be translated to its complement
     */
    public static String reverseComplement(String sequence) throws RuntimeException {
        return reverse(complement(sequence));
    }

    /**
     * Returns the complement of a biological sequence. A RuntimeException is thrown
     * of one of the characters could not be translated to its complement.
     *
     * @param sequence the biological sequence
     * @return complement the complement sequence
     * @throws RuntimeException in case a character could not be translated to its complement
     */
    public static String complement(String sequence) throws RuntimeException {
        StringBuilder buffy = new StringBuilder(sequence.length());
        for (int i = 0; i < sequence.length(); i++) {
            int p = BIOLOGICAL_CHARS.indexOf(sequence.charAt(i));
            if (p < 0) {
                throw new RuntimeException("Complement: unknown symbol " + sequence.charAt(i) + " in " + sequence);
            }
            buffy.append(COMPLEMENT_CHARS.charAt(p));
        }
        return buffy.toString();
    }

    /**
     * Returns the String representation for a double using the specified number of
     * decimal positions.
     *
     * @param value the value
     * @param dec   decimal
     * @return string string representation
     */
    public static String fprint(double value, int dec) {
        String s = Double.toString(value);
        int p = s.lastIndexOf(".");
        if (p < 0) {
            s += ".";
            for (int i = 0; i < dec; i++) {
                s += "0";
            }
        } else {
            int q = s.indexOf("E");
            String exp = "";
            if (q >= 0) {
                exp = s.substring(q);
            }
            int end = p + dec + 1;
            if (end < s.length()) {
                s = s.substring(0, end);
            } else {
                for (int i = s.length(); i < end; i++) {
                    s += "0";
                }
            }
            s += exp;
        }


        return s;
    }


    /**
     * Append the character c to the given string until given length is reached
     *
     * @param c the character
     * @param s the string
     * @param len length of the final string
     * @param leading append as prefix
     * @return string the final string
     */
    public static String append(char c, String s, int len, boolean leading) {
        if (s.length() >= len) {
            return s;
        }
        StringBuilder sb = new StringBuilder(s);
        for (int i = s.length(); i < len; ++i) {
            if (leading) {
                sb.insert(0, c);
            } else {
                sb.append(c);
            }
        }

        return sb.toString();
    }
    
    /**
     * Creates a <code>String</code> instance with all elements of the 
     * <code>Collection</code> instance provided, separated by the 
     * character specified as &quot;separator&quot;  
     * @param collection set of objects that are to be printed in the
     * string
     * @param separator the separator character
     * @return a string with all elements represented by calling their 
     * <code>toString()</code> method
     */
    public static String toString(Collection collection, char separator) {
    	return toString(Collections.enumeration(collection), separator);
    }	

    /**
     * Creates a <code>String</code> instance with all elements of the 
     * <code>Enumeration</code> instance provided, separated by the 
     * character specified as &quot;separator&quot;  
     * @param enumeration set of objects that are to be printed in the
     * string
     * @param separator the separator character
     * @return a string with all elements represented by calling their 
     * <code>toString()</code> method
     */
    public static String toString(Enumeration enumeration, char separator) {
    	StringBuilder sb= new StringBuilder();
    	while(enumeration.hasMoreElements()) {
    		sb.append(enumeration.nextElement().toString());
    		sb.append(separator);
    	}
    	sb.deleteCharAt(sb.length()- 1);
    	return sb.toString();
    }

    /**
     * Gets the <code>fieldNr</code><i>th</i> token of the provided
     * <code>line</code>. Starts to search from the start of
     * <code>line</code>, inefficient for the last tokens in a line.
     * Trims C/R characters if <code>line</code> is terminated by
     * a newline sequence.
     * @param line the line from which the token is to be extracted
     * @param separator the separator character. Multiple instances are
     * interpreted as separators enclosing empty strings
     * @param fieldNr the number of the token that is to be extracted
     * @return
     */
    public static CharSequence getField(CharSequence line, char separator, int fieldNr) {

        // left separator
        int l= line.length();
        int i= 0, n= 0;
        while(n< fieldNr&& i< l)
            if (line.charAt(i++)== separator)
                ++n;
        if (n< fieldNr)
            return null;
        int j= i;

        // right separator
        ++fieldNr;
        while(n< fieldNr&& j< l)
            if (line.charAt(j++)== separator)
                ++n;

        // trim C/R if EOL spotted
        char c= line.charAt(--j);
        while (c== '\n'|| c== '\r')
            c= line.charAt(--j);
        if (c!= separator)
            ++j;

        return line.subSequence(i, j);
    }

    /**
     * Join a list of elements by calling toString() on the elements.
     * Null elements will throw an exception
     *
     * @param separator the separator
     * @param elements the elements
     * @return joined the joined string
     *
     * @since 1.19
     */
    public static String join(String separator, Collection<?> elements) {
        if(separator == null) throw new NullPointerException("NULL separator not permitted");
        if(elements == null) throw new NullPointerException("NULL element list not permitted");
        StringBuilder b = null;
        for (Object element : elements) {
            if(element == null) throw new NullPointerException("NULL element can not be joined");
            if(b == null){
                b = new StringBuilder();
            }else{
                b.append(separator);
            }
            b.append(element);
        }
        return b== null ? "" : b.toString();
    }
}
