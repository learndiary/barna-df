package fbi.genome.model.constants;

import java.lang.management.ManagementFactory;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Date;

import commons.Dialogable;
import commons.Progressable;


public class Constants {

	public static final String PROPERTY_KEY_WRAPPER_BASE= "wrapperDir";
	
	public static final byte VERBOSE_DEBUG = 3;
	public static final byte VERBOSE_ERRORS = 2;
	public static final byte VERBOSE_NORMAL = 1;
	public static final byte VERBOSE_SHUTUP = 0;
	public static final String[] VERBOSE_KEYWORDS= new String[] {"SILENT", "VERBOSE", "ERRORS", "DEBUG"};

	public static String globalPfx= null;
	
	public static String getGlobalPfx() {
		if (globalPfx == null) {
			globalPfx = ManagementFactory.getRuntimeMXBean().getName();
			if (globalPfx== null|| globalPfx.length()< 6)
				globalPfx= getTimestampID();
		}

		return globalPfx;
	}

	public static byte verboseLevel= VERBOSE_NORMAL;
	
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
	
	public static final String CLI_PAR_LONG= "--", CLI_PAR_SHORT= "-";
	
	public static Progressable progress= new PrintstreamProgressable(System.err);
	public static Dialogable dialog= new PrintstreamDialogable(System.err, System.in);
	
	public static String getTime() {		
		return DateFormat.getTimeInstance().format(new Date(System.currentTimeMillis()));
	}
	
	public static String getTimestampID() {
		SimpleDateFormat format= new SimpleDateFormat("yyMMddHHmmssSSSS");
		return format.format(new Date());
	}

	public static final String PROPERTY_TMPDIR= "java.io.tmpdir";
	
	public static final boolean setVerbose(String s) {
		//SILENT, VERBOSE, ERRORS, DEBUG
		s= s.toUpperCase();
		for (int i = 0; i < VERBOSE_KEYWORDS.length; i++) 
			if (s.equals(VERBOSE_KEYWORDS[i])) { 
				verboseLevel= (byte) i;
				return true;
			}
			
		return false;
	}

	public static byte DIR_FORWARD= 0;

	public static byte DIR_BACKWARD= 1;

	public static byte DIR_BOTH= 2;
	
	public static int binarySearch(Object[] a, int from, int to, Object key) {
		
		if (from< 0|| from>= a.length|| to< 0|| to> a.length)
			throw new IllegalArgumentException();
		if (from>= to)
			return -1;
		
		int low = from;
		int high = to- 1;

		while (low <= high) {
		    int mid = (low + high) >> 1;
		    Comparable midVal = (Comparable)a[mid];
		    int cmp = midVal.compareTo(key);

		    if (cmp < 0)
			low = mid + 1;
		    else if (cmp > 0)
			high = mid - 1;
		    else
			return mid; // key found
		}
		return -(low + 1);  // key not found.
	}
}
