package fbi.commons;

/**
 * A helper class for tools to work with the 
 * <code>CharSequence</code> interface.
 * 
 * @author Micha Sammeth (gmicha@gmail.com)
 *
 */
public class CharSequences {

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


}
