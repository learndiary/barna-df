package fbi.commons;

import java.io.PrintStream;
import java.text.SimpleDateFormat;
import java.util.Date;

/**
 * String constants and utilities
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class StringUtils {
    public static final String TAB= "\t";
    public static final String HASH= "#";
    public static final String EMPTYSTRING= "";
    public static final String STAR= "*";
    public static final String SPACE= " ";
    public static final String DOT= ".";
    public static final String NULL= "0";
    public static final String COLON= ":";
    public static final String PAROPEN= "(";
    public static final String PARCLOSE= ")";
    public static final String OK= "OK";
    public static final String ERROR= "ERROR";
    /**
     * Biological characters as spotted in refseq intron annotation of human
     */
    private static final String BIOLOGICAL_CHARS = "ACGTN-acgtnMmKkXxRr";
    /**
     * The complements for biological chars
     */
    private static final String COMPLEMENT_CHARS = "TGCAN-tgcanKkMmXxYy";


    /**
     * Create a timestamp for the current time
     *
     * @return timestamp timestamp
     */
    public static String getTimestampID() {
        SimpleDateFormat format= new SimpleDateFormat("yyMMddHHmmssSSSS");
        return format.format(new Date());
    }

    /**
     * Reverse the given string
     *
     * @param in the string
     * @return reverse the reverse string
     */
    public static String reverse(String in) {
        StringBuffer buffy= new StringBuffer(in.length());
        for (int i = in.length()-1; i >=0; i--) {
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
		StringBuffer buffy= new StringBuffer(sequence.length());
		for (int i = 0; i < sequence.length(); i++) {
			int p= BIOLOGICAL_CHARS.indexOf(sequence.charAt(i));
			if (p < 0) {
				throw new RuntimeException("Complement: unknown symbol "+ sequence.charAt(i)+ " in "+ sequence);
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
     * @param dec decimal
     * @return string string representation
     */
    public static String fprint(double value, int dec) {
        String s= Double.toString(value);
        int p= s.lastIndexOf(".");
        if (p< 0) {
            s+= ".";
            for (int i = 0; i < dec; i++)
                s+= "0";
        } else {
            int q= s.indexOf("E");
            String exp= "";
            if (q>= 0)
                exp= s.substring(q);
            int end= p+ dec+ 1;
            if (end< s.length())
                s= s.substring(0, end);
            else
                for (int i = s.length(); i < end; i++)
                    s+= "0";
            s+= exp;
        }


        return s;
    }

    public static int printPercentage(int perc, double val, double base, PrintStream stream) {
        int currPerc= (int) Math.round(val*10d/ base);
        if (currPerc> perc) {
            if (perc== 4|| perc== 9)
                stream.print("+");
            else
                stream.print("*");
            stream.flush();
            ++perc;
        }
        return perc;
    }

    public static String append(char c, String s, int len, boolean leading) {
        if (s.length()>= len)
            return s;
        StringBuilder sb= new StringBuilder(s);
        for (int i = s.length(); i < len; ++i) {
            if (leading)
                sb.insert(0, c);
            else
                sb.append(c);
        }

        return sb.toString();
    }
}
