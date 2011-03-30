package fbi.genome.io.rna;

import java.util.Arrays;

public class CopyOfBARNAdescriptor extends Descriptor {

	final public static char SEPARATOR= '/';
	final public static char[] SYMBOLS= new char[] {'a', 's', '1', '2'};
	static int[] BARNA_MODES= MODES.clone();
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
	@Override
	public boolean allowsPairs() {		
		return true;
	}
	@Override
	public boolean allowsStranded() {
		return true;
	}
	
}
