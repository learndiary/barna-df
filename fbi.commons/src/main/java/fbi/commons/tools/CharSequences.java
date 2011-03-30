package fbi.commons.tools;

public class CharSequences {

	/**
	 * maximum one of the fs separates two fields (non-adjacent runs)
	 * @param cs
	 * @param nr
	 * @param sep
	 * @return
	 */
	public static CharSequence getToken(CharSequence cs, int nr, char[] fs) {
		if (nr== 0)
			return null;
		int last= 0, now= 0, cnt= 0;
		for (int i = 0; cnt < nr&& i< cs.length(); i++) {
			for (int j = 0; j < fs.length; j++) {
				if (cs.charAt(i)== fs[j]) {
					++cnt;
					last= now;
					now= i;
					break;
				}
			}
		}
		
		if (last> 0)
			++last;
		if (cnt< nr) {
			if (cnt== nr- 1) {
				if (now+ 1== cs.length())
					return null;
				else
					return cs.subSequence(now+ 1, cs.length());
			}
			else return null;
		}
		if (last== now)
			return null;
		else
			return cs.subSequence(last, now);
	}

	
    public static int compare(CharSequence string, CharSequence anotherString) {
    	
    	int len1 = string.length();
    	int len2 = anotherString.length();
    	int n = Math.min(len1, len2);
    	
	    for (int k = 0; k < n; k++) {
    		char c1 = string.charAt(k);
    		char c2 = anotherString.charAt(k);
    		if (c1 != c2) {
    		    return c1 - c2;
    		}
    		k++;
	    }

    	return len1 - len2;
    }

    public static int parseInt(CharSequence cs) {
		int val= 0;
		for (int i = cs.length()-1, pow = 1; i >= 0; --i, pow *= 10)
			val += (((byte) cs.charAt(i)) - 48) * pow;
		return val;
	}
}
