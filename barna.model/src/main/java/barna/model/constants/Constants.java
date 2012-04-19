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

package barna.model.constants;

import java.lang.management.ManagementFactory;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Date;


public class Constants {

	public static final String PROPERTY_KEY_WRAPPER_BASE= "wrapperDir";
	
	public static final byte VERBOSE_DEBUG = 3;
	public static final byte VERBOSE_ERRORS = 2;
	public static final byte VERBOSE_NORMAL = 1;
	public static final byte VERBOSE_SHUTUP = 0;
	public static final String[] VERBOSE_KEYWORDS= new String[] {"SILENT", "VERBOSE", "ERRORS", "DEBUG"};

	public static String globalPfx= null;
	
	public static String getGlobalPfx() {
		if (globalPfx == null) {
			globalPfx = ManagementFactory.getRuntimeMXBean().getName();
			if (globalPfx== null|| globalPfx.length()< 6)
				globalPfx= getTimestampID();
		}

		return globalPfx;
	}

	public static byte verboseLevel= VERBOSE_NORMAL;
	public static final char NL= '\n';
	public static final String TAB= "\t";
	public static final String HASH= "#";
	public static final String EMPTYSTRING= "";
	public static final String STAR= "*";
	public static final String SPACE= " ";
	public static final String DOT= ".";
	public static final String NULL= "0";
	public static final String COLON= ":";
	public static final String SEMICOLON= ";";
	public static final String PAROPEN= "(";
	public static final String PARCLOSE= ")";
	public static final String OK= "OK";
	public static final String ERROR= "ERROR";
	
	public static final String CLI_PAR_LONG= "--", CLI_PAR_SHORT= "-";

	public static String getTime() {		
		return DateFormat.getTimeInstance().format(new Date(System.currentTimeMillis()));
	}
	
	public static String getTimestampID() {
		SimpleDateFormat format= new SimpleDateFormat("yyMMddHHmmssSSSS");
		return format.format(new Date());
	}

	public static final String PROPERTY_TMPDIR= "java.io.tmpdir";
	
	public static final boolean setVerbose(String s) {
		//SILENT, VERBOSE, ERRORS, DEBUG
		s= s.toUpperCase();
		for (int i = 0; i < VERBOSE_KEYWORDS.length; i++) 
			if (s.equals(VERBOSE_KEYWORDS[i])) { 
				verboseLevel= (byte) i;
				return true;
			}
			
		return false;
	}

	public static byte DIR_FORWARD= 0;

	public static byte DIR_BACKWARD= 1;

	public static byte DIR_BOTH= 2;
	
	public static int binarySearch(Object[] a, int from, int to, Object key) {
		
		if (from< 0|| from>= a.length|| to< 0|| to> a.length)
			throw new IllegalArgumentException();
		if (from>= to)
			return -1;
		
		int low = from;
		int high = to- 1;

		while (low <= high) {
		    int mid = (low + high) >> 1;
		    Comparable midVal = (Comparable)a[mid];
		    int cmp = midVal.compareTo(key);

		    if (cmp < 0)
			low = mid + 1;
		    else if (cmp > 0)
			high = mid - 1;
		    else
			return mid; // key found
		}
		return -(low + 1);  // key not found.
	}
}
