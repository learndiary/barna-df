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

import java.util.Arrays;

public class BARNAdescriptor implements ReadDescriptor {

	final public static char SEPARATOR= '/';
	final public static char[] SYMBOLS= new char[] {'a', 's', '1', '2'};
	static int[] BARNA_MODES= Descriptor.MODES.clone();
	// resort, for generality
	{
		char[] symb= SYMBOLS.clone();
		int[] mods= BARNA_MODES.clone();
		Arrays.sort(SYMBOLS);
		for (int i= 0;  i< symb.length; ++i) {
			int p= Arrays.binarySearch(SYMBOLS, symb[i]);
			BARNA_MODES[i]= mods[p];
		}
	}
	
	public static void main(String[] args) {
	}
	/**
	 * 
	 * @param cs
	 * @return -1 for error, 0 for no annotation
	 */
	public int getMode(CharSequence cs, int[] fromTo) {
		assert(fromTo.length>= 2);
		fromTo[0]= 0;
		int i = cs.length()- 1;
		char sep= SEPARATOR;
		for (; i>= 0; --i) 
			if (cs.charAt(i)== sep)
				break;
		if (i< 0) {
			fromTo[1]= cs.length();
			return 0;	// allow for no annotation !
		}
		if (i== cs.length()- 1) {
			fromTo[1]= cs.length()- 1;
			return 0;
		}
		// else
		fromTo[1]= i;
		return parseBarna(cs, i+1);
	}
	
	private int parseBarna(CharSequence cs, int from) {
		
		int res= 0;
		for (int i = from; i < cs.length(); ++i) {
			int p= Arrays.binarySearch(SYMBOLS, cs.charAt(i));
			if (p< 0)
				return 0;	// nothing, or alleged separator
			res|= BARNA_MODES[p];
		}
		return res;
	}
	
	@Override
	public String toString() {
		return "BaRNA ups (.*)/([12])([as])";
	}
	public boolean allowsPend() {
		// TODO Auto-generated method stub
		return false;
	}
	public byte getPairedEndInformation(CharSequence descriptor) {
		// TODO Auto-generated method stub
		return 0;
	}
	public CharSequence getUniqueDescriptor(CharSequence descriptor) {
		// TODO Auto-generated method stub
		return null;
	}
	public boolean isApplicable(CharSequence descriptor) {
		// TODO Auto-generated method stub
		return false;
	}
	public boolean isPairedEnd(CharSequence descriptor) {		
		return true;
	}
	public boolean allowsStranded() {
		return true;
	}
	public byte getStrand(CharSequence descriptor) {
		// TODO Auto-generated method stub
		return 0;
	}
	public boolean isStranded(CharSequence descriptor) {
		// TODO Auto-generated method stub
		return false;
	}
	
}
