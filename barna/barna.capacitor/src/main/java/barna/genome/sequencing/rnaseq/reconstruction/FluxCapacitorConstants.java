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

package barna.genome.sequencing.rnaseq.reconstruction;

import java.io.File;
import java.lang.reflect.Method;
import java.util.HashMap;

public class FluxCapacitorConstants {

	public static final String CLI_LONG_PFX= "--";
	public static final String CLI_LONG_BATCH= "batch";
	public static final String CLI_LONG_COMPRESSION= "compress";
	public static final String CLI_LONG_FORCE= "force";
	public static final String CLI_LONG_HELP= "help";
	public static final String CLI_LONG_INSTALL= "install";
	public static final String CLI_LONG_JVM= "jvm";
	public static final String CLI_LONG_LIB= "lib";
	public static final String CLI_LONG_LOCAL= "local";
	public static final String CLI_LONG_FILENAME= "name";
	public static final String CLI_LONG_OUT= "output";
	public static final String CLI_LONG_PAIR= "pair";
	public static final String CLI_LONG_PROFILE= "pro";
	public static final String CLI_LONG_REF= "ref";
	public static final String CLI_LONG_SRA= "sra";
	public static final String CLI_LONG_SSPECIFIC= "sp";
	public static final String CLI_LONG_THREAD= "thread";
	public static final String CLI_LONG_TMP= "tmp";
	public static final String CLI_LONG_TPX= "tpx";
	public static final String CLI_LONG_UNIF= "uniform";
	public static final String CLI_LONG_VERBOSE= "verbose";
	public static final Character  CLI_SHORT_PFX= '-';
	public static final Character CLI_SHORT_BATCH= 'b';
	public static final Character CLI_SHORT_COMPRESSION= 'c';
	public static final Character CLI_SHORT_FORCE= 'f';
	public static final Character CLI_SHORT_LOCAL= 'l';
	public static final Character CLI_SHORT_FILENAME= 'n';
	public static final Character CLI_SHORT_HELP= 'h';
	public static final Character CLI_SHORT_OUT= 'o';
	public static final Character CLI_SHORT_PAIR= 'p';
	public static final Character CLI_SHORT_THREAD= 't';
	public static final Character CLI_SHORT_REF= 'r';
	public static final Character CLI_SHORT_SRA= 's';
	public static final Character CLI_SHORT_UNIF= 'u';
	public static final Character CLI_SHORT_VERBOSE= 'v';
	public static final String GTF_ATTRIBUTE_TOKEN_OBSV= "obs";
	public static final String GTF_ATTRIBUTE_TOKEN_PRED= "pred";
	public static final String GTF_ATTRIBUTE_TOKEN_BALANCED= "bal";
	public static final String GTF_ATTRIBUTE_TOKEN_ALL= "all";
	public static final String GTF_ATTRIBUTE_TOKEN_TID= "split";
	public static final String GTF_ATTRIBUTE_TOKEN_EXC= "uniq";
	public static final String GTF_ATTRIBUTE_TOKEN_READS= "freq";
	public static final String GTF_ATTRIBUTE_TOKEN_RFREQ= "rfreq";
	public static final String GTF_ATTRIBUTE_TOKEN_RPKM= "RPKM"; // rpkm
	public static final String GTF_ATTRIBUTE_TOKEN_COV= "cov";
	public static final String GTF_ATTRIBUTE_TOKEN_FWD= "fwd";
	public static final String GTF_ATTRIBUTE_TOKEN_REV= "rev";
	public static final String GTF_ATTRIBUTE_TOKEN_BID= "bid";
	public static final String GTF_ATTRIBUTE_TOKEN_SEP= "_";
	public static final String GTF_ATTRIBUTE_LENGTH= "slots";
	public static final String GTF_ATTRIBUTE_PROFILE= "profile";
	public static final String GTF_ATTRIBUTE_EXPECT= "expect";
	static final String[] GTF_ATTRIBUTES_BASE= new String[] {GTF_ATTRIBUTE_TOKEN_OBSV, GTF_ATTRIBUTE_TOKEN_PRED, GTF_ATTRIBUTE_TOKEN_BALANCED};
	static final String[] GTF_ATTRIBUTES_RESOLUTION= new String[] {GTF_ATTRIBUTE_TOKEN_ALL, GTF_ATTRIBUTE_TOKEN_TID, GTF_ATTRIBUTE_TOKEN_EXC};
	static final String[] GTF_ATTRIBUTES_MEASUREMENT= new String[] {GTF_ATTRIBUTE_TOKEN_READS, GTF_ATTRIBUTE_TOKEN_RFREQ, GTF_ATTRIBUTE_TOKEN_RPKM};
	public static final String GTF_ATTRIBUTE_PVAL= "falsification";
	public static final String GFF_FEATURE_JUNCTION = "junction";
	public static final String GFF_FEATURE_PAIRED = "paired";
	public static final String GFF_FEATURE_FRAGMENT = "fragment";
	public static final byte STRAND_NONE= 0;
	public static final byte STRAND_ENABLED= 1;
	public static final byte STRAND_SPECIFIC= 2;
	public static byte SHELL_NONE;
	public static byte SHELL_BASH= 1;
	public static byte SHELL_CSH= 2;
	public static byte SHELL_KSH= 3;
	public static final String SUBDIR_NATIVELIBS= "lib"+File.separator+"native";
	public static final String SUBDIR_LPSOLVE= "lpsolve55";
	public static final String SFX_GTF= "gtf";
	public static final String SFX_BED= "bed";
	public static final String SFX_INFLATE= "__inflated";
	public static final String SFX_MAPPED= "_mapped";
	public static final String SFX_NOTMAPPED= "_notmapped";
	public static final String SFX_PROFILES= "_profiles";
	public static final String SFX_LP= "_lp";
	public static final String SFX_INSERTSIZE= "_insertsize";
	static final String PFX_CAPACITOR= "capacitor";
	static final String PFX_MAPPED_READS= "mapped";
	static final String PFX_UNMAPPED_READS= "notmapped";
	static final String obsReadsAllTag= GTF_ATTRIBUTE_TOKEN_OBSV+ GTF_ATTRIBUTE_TOKEN_SEP+ GTF_ATTRIBUTE_TOKEN_ALL+ GTF_ATTRIBUTE_TOKEN_SEP+ GTF_ATTRIBUTE_TOKEN_READS;
	static final String obsReadsSplitTag= GTF_ATTRIBUTE_TOKEN_OBSV+ GTF_ATTRIBUTE_TOKEN_SEP+ GTF_ATTRIBUTE_TOKEN_TID+ GTF_ATTRIBUTE_TOKEN_SEP+ GTF_ATTRIBUTE_TOKEN_READS;
	static final String obsReadsUniqTag= GTF_ATTRIBUTE_TOKEN_OBSV+ GTF_ATTRIBUTE_TOKEN_SEP+ GTF_ATTRIBUTE_TOKEN_EXC+ GTF_ATTRIBUTE_TOKEN_SEP+ GTF_ATTRIBUTE_TOKEN_READS;
	static final String predReadsAllTag= GTF_ATTRIBUTE_TOKEN_PRED+ GTF_ATTRIBUTE_TOKEN_SEP+ GTF_ATTRIBUTE_TOKEN_ALL+ GTF_ATTRIBUTE_TOKEN_SEP+ GTF_ATTRIBUTE_TOKEN_READS;
	static final String predReadsSplitTag= GTF_ATTRIBUTE_TOKEN_PRED+ GTF_ATTRIBUTE_TOKEN_SEP+ GTF_ATTRIBUTE_TOKEN_TID+ GTF_ATTRIBUTE_TOKEN_SEP+ GTF_ATTRIBUTE_TOKEN_READS;
	static final String predReadsUniqTag= GTF_ATTRIBUTE_TOKEN_PRED+ GTF_ATTRIBUTE_TOKEN_SEP+ GTF_ATTRIBUTE_TOKEN_EXC+ GTF_ATTRIBUTE_TOKEN_SEP+ GTF_ATTRIBUTE_TOKEN_READS;
	public static final char CLI_OUT_ALL= 'a';
	public static final char CLI_OUT_BALANCED= 'b';
	public static final char CLI_OUT_COVERAGE= 'c';
	public static final char CLI_OUT_PRED= 'd';
	public static final char CLI_OUT_RFREQ= 'f'; // fraction
	public static final char CLI_OUT_GENE= 'g';
	public static final char CLI_OUT_ISIZE= 'i';
	public static final char CLI_OUT_SJUNCTION= 'j';	// junctions in general.. separate later 
	public static final char CLI_OUT_KEEPSORTED= 'k';
	public static final char CLI_OUT_LP= 'l';
	public static final char CLI_OUT_MAPPED= 'm';
	public static final char CLI_OUT_NOTMAPPED= 'n';
	public static final char CLI_OUT_OBS= 'o';
	public static final char CLI_OUT_PROFILES= 'p';
	public static final char CLI_OUT_UNIQUE= 'q';	// WAS: u
	public static final char CLI_OUT_FREQ= 'r';	// reads 
	public static final char CLI_OUT_SPLIT= 's';
	public static final char CLI_OUT_TRANSCRIPT= 't';
	public static final char CLI_OUT_LOCUS= 'u';	// locus, cluster
	public static final char CLI_OUT_EVENTS= 'v';
	public static final char CLI_OUT_EXON= 'x';	// WAS: e 
	public static final char CLI_OUT_INTRON= 'y';
	public static final String CLI_ABBREV_COST_BOUNDS= "cb";
	public static final String CLI_ABBREV_COST_MODEL= "cm";
	public static final String CLI_ABBREV_COST_SPLIT= "cs";
	public static final String CLI_CMD= "flux";
	protected static HashMap<String, Method> cliLongMap= new HashMap<String, Method>();
	protected static HashMap<Character, Method> cliShortMap= new HashMap<Character, Method>();
	protected static HashMap<String[], String> cliExplMap= new HashMap<String[], String>();
	public static final String VERSION_ID= "1.2";
	public static final String LPSOLVE_LIB_NAME= "lpsolve55";
	public static final String LPSOLVE_JNI_NAME= "lpsolve55j";
	public static final String DEBUG_LP_OUT= "C:\\lp_out";
	static final int BIG= 999999999;
	public static final byte BYTE_0= (byte) 0;
	public static final byte BYTE_1= (byte) 1;
	public static final byte BYTE_MINUS_1= (byte) (-1);
	public static int[] BIN_EXP= new int[] {10, 100};
	public static int[] BIN_LEN= new int[] {1000, 2000};
	static final String normReadsAllTag= GTF_ATTRIBUTE_TOKEN_BALANCED+ GTF_ATTRIBUTE_TOKEN_SEP+ GTF_ATTRIBUTE_TOKEN_ALL+ GTF_ATTRIBUTE_TOKEN_SEP+ GTF_ATTRIBUTE_TOKEN_READS;
	static final String normReadsSplitTag= GTF_ATTRIBUTE_TOKEN_BALANCED+ GTF_ATTRIBUTE_TOKEN_SEP+ GTF_ATTRIBUTE_TOKEN_TID+ GTF_ATTRIBUTE_TOKEN_SEP+ GTF_ATTRIBUTE_TOKEN_READS;
	static final String normReadsUniqTag= GTF_ATTRIBUTE_TOKEN_BALANCED+ GTF_ATTRIBUTE_TOKEN_SEP+ GTF_ATTRIBUTE_TOKEN_EXC+ GTF_ATTRIBUTE_TOKEN_SEP+ GTF_ATTRIBUTE_TOKEN_READS;
	public static final byte MODE_LEARN= 0;
	public static final byte MODE_RECONSTRUCT= 1;
	private static final String NULL_STR= "0";
	static final char UNDERSCORE= '_';
	public static final String PROPERTY_BUILD= "capacitor.build";
	public static final String PROPERTY_JDK= "capacitor.jdk";
	public static final String VALUE_NA= "0";	// NA
	final static String[] L_USER_COMMENTS= new String[] {"[OOOPS] ", "[HEOO] ", "[PLONG] ", "[BAOOO] "};
	static final String lengthAllTag= GTF_ATTRIBUTE_LENGTH+ GTF_ATTRIBUTE_TOKEN_SEP+ GTF_ATTRIBUTE_TOKEN_ALL;
	static final String lengthSplitTag= GTF_ATTRIBUTE_LENGTH+ GTF_ATTRIBUTE_TOKEN_SEP+ GTF_ATTRIBUTE_TOKEN_TID;
	static final String lengthUniqTag= GTF_ATTRIBUTE_LENGTH+ GTF_ATTRIBUTE_TOKEN_SEP+ GTF_ATTRIBUTE_TOKEN_EXC;
	public static byte FORMAT_BED= 0;
	public static byte FORMAT_GTF= 1;
	public static byte FORMAT_SAM= 2;
	static final String FLOAT_STRING_0= "0.0";
	static String FNAME_PROPERTIES= "capacitor.prop";

}
