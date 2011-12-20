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
