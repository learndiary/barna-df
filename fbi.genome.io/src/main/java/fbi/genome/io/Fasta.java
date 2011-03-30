package fbi.genome.io;

public class Fasta {

	public static byte getReadDir(String s) {
		//chr16:1500431-1600693C;uc002cma.1;9163;1;49
		int p= 0;
		for (int i = 0; i < 3; i++, p= s.indexOf(';', p+1));
		++p;
		byte dir= Byte.parseByte(s.substring(p, s.indexOf(';', p)));
		return dir;
	}
	
	public static String getReadID(String s) {
		//chr16:1500431-1600693C;uc002cma.1;9163;1;49
		int p= 0;
		for (int i = 0; i < 3; i++, p= s.indexOf(';', p+1));
		
		return s.substring(0,p);
	}

	public static final String QFASTA_PLUS = "+";
}
