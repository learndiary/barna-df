package fbi.commons;

import java.io.PrintStream;

public class MyFormatter {
	
	public static String american(long integer) {
		StringBuffer sb= new StringBuffer(Long.toString(integer));
		int offs= 0;
		for (int i = 0; i < sb.length(); ++i) 
			if ((i-offs)%3== 0) {
				sb.insert(sb.length()- offs- i, ',');
				++offs;
			}
		return sb.toString();
	}
	
	public static String add1KSeparator(long x) {
		StringBuilder sb= new StringBuilder(Long.toString(x));
		
		return sb.toString();
	}
	
	public static String fprint(double fp, int dec) {
		String s= java.lang.Double.toString(fp); 
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
