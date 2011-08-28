package fbi.commons;

import java.util.Comparator;

public class CharsequenceComparator implements Comparator<CharSequence> {

	public static final CharsequenceComparator DEFAULT_CHARSEQUENCE_COMPARATOR= new CharsequenceComparator();
	
	@Override
	public int compare(CharSequence paramT1, CharSequence paramT2) {
		
		int len= Math.min(paramT1.length(), paramT2.length());
		for (int i = 0; i < len; ++i) {
			char c1= paramT1.charAt(i), 
				c2= paramT2.charAt(i);
			if (c1!= c2)
				return (c1-c2);
		}
				
		return paramT1.length()- paramT2.length();
	}
}
