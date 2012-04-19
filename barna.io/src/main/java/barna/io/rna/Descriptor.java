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

public abstract class Descriptor {
	final public static byte MATE_UNKNOWN= -1, MATE_1= 0, MATE_2= 1, ORIENTATION_UNKNOWN= -1, ORIENTATION_SENSE= 0, ORIENTATION_ASENSE= 2; 
	final public static int MODE_MATE1= 1, MODE_MATE2= 2, MODE_SENSE= 4, MODE_ASENSE= 8, MAX_ATTRIB= 16;
	final public static int[] MODES= new int[] {MODE_ASENSE, MODE_SENSE, MODE_MATE1, MODE_MATE2};
	final public static char CHAR_ID_ID= 'd', CHAR_ID_PAIRED= 'p', CHAR_ID_STRANDED= 's';
	
	public static byte getMate(int mode) {
		if ((mode& MODE_MATE1)!= 0) {
			assert((mode& MODE_MATE2)== 0);
			return MATE_1;
		}
		if ((mode& MODE_MATE2)!= 0) {
			assert((mode& MODE_MATE1)== 0);
			return MATE_2;
		}
		return MATE_UNKNOWN;
	}

	public static byte getStrand(int mode) {
		if ((mode& MODE_SENSE)!= 0) {
			assert((mode& MODE_ASENSE)== 0);
			return ORIENTATION_SENSE;
		}
		if ((mode& MODE_SENSE)!= 0) {
			assert((mode& MODE_ASENSE)== 0);
			return ORIENTATION_ASENSE;
		}
		return ORIENTATION_UNKNOWN;
	}
	
	public boolean isValid() {
		return true;
	}
	public abstract int getMode(CharSequence cs, int[] fromTo);
	
	public abstract boolean allowsPairs();
	public abstract boolean allowsStranded();
}
