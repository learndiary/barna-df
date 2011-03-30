package fbi.genome.io.rna;

import fbi.genome.model.bed.BEDobject2;
import fbi.genome.model.Transcript;

import java.util.regex.Pattern;

/**
 * @deprecated
 * @author msammeth
 *
 */
public class FMRD implements ReadDescriptor {

	public static final char ID_PE= 'P';
	public static final char ID_MM= 'M';
	
	// from FluxSimulatorSettings
	public static final byte BYTE_DELIM_FURI= '#';	// flux unique read identifier
	public static final byte BYTE_DELIM_FMOLI= ':';	// flux molecule identifier
	public static final byte BYTE_DELIM_BARNA= '/';

	public static final char[] DELIM_FMRD= new char[] {'|',',',';'};	// flux mapped read descriptor 
	public static final char[] PE_OPT= new char[] {'1','2'};
	public static final Pattern pattPE= Pattern.compile("^[P,p][1,2]");
	public static final Pattern pattMM= Pattern.compile("^M\\d$");	//"^M\\d{}$" 

	// from Sequencer
	public static void addReadName(BEDobject2 obj, Transcript t, Transcript t2, 
			long molNr, byte absDir, int fragStart, int fragEnd, int readStart, int readEnd, boolean pend) {
		
		// FURI
		obj.replace(BEDobject2.FN_NAME, t.getGene().getGeneID());
		obj.appendCurrField(BYTE_DELIM_FMOLI);
		obj.appendCurrField(t.getTranscriptID());
		obj.appendCurrField(BYTE_DELIM_FMOLI);
		obj.appendCurrField((int) (molNr+1));
		obj.appendCurrField(BYTE_DELIM_FMOLI);
		obj.appendCurrField(t.getExonicLength());
		obj.appendCurrField(BYTE_DELIM_FMOLI);
		obj.appendCurrField(fragStart);
		obj.appendCurrField(BYTE_DELIM_FMOLI);
		obj.appendCurrField(fragEnd);
		obj.appendCurrField(BYTE_DELIM_FMOLI);
		obj.appendCurrField(readStart);
		obj.appendCurrField(BYTE_DELIM_FMOLI);
		obj.appendCurrField(readEnd);
		if (pend) {
			obj.appendCurrField(FMRD.DELIM_FMRD[0]);
			obj.appendCurrField(FMRD.ID_PE);
			obj.appendCurrField(absDir== t.getStrand()? PE_OPT[0]: PE_OPT[1]);
		}
		
	}

	public static void appendReadName(BEDobject2 obj, Transcript t, Transcript t2, 
			long molNr, byte absDir, int fragStart, int fragEnd, int readStart, int readEnd, boolean pend) {
		
		// FURI
		obj.append(t.getGene().getGeneID());
		obj.append(BYTE_DELIM_FMOLI);
		obj.append(t.getTranscriptID());
		obj.append(BYTE_DELIM_FMOLI);
		obj.append((int) (molNr+1));
		obj.append(BYTE_DELIM_FMOLI);
		obj.append(t.getExonicLength());
		obj.append(BYTE_DELIM_FMOLI);
		obj.append(fragStart);
		obj.append(BYTE_DELIM_FMOLI);
		obj.append(fragEnd);
		obj.append(BYTE_DELIM_FMOLI);
		obj.append(readStart);
		obj.append(BYTE_DELIM_FMOLI);
		obj.append(readEnd);
		if (pend) {
			obj.append((byte) BYTE_DELIM_BARNA);
			//obj.append((byte) FMRD.ID_PE);
			obj.append((byte) (absDir== t.getStrand()? PE_OPT[0]: PE_OPT[1]));
		}
		
	}
	
	/**
	 * encodes FMRD
	 * @param URI unique read ID
	 * @param pEnd	(-1) not used, 0 fwd, 1 rev
	 * @param cntMM	(-1) not used, ...
	 * @return
	 */
	public static String getFMRD(String URI, int pEnd, int cntMM) {
		return getFMRD(URI, pEnd, cntMM, 0);
	}
	
	public static String getFMRD(String URI, int pEnd, int cntMM, int delim) {
		StringBuilder sb= new StringBuilder(URI);
		if (pEnd>= 0) {
			sb.append(DELIM_FMRD[delim]);
			sb.append(PE_OPT[pEnd]);
		}
		if (cntMM>= 0) {
			sb.append(DELIM_FMRD[delim]);
			sb.append(Integer.toString(cntMM));
		}
		return sb.toString();
	}
	
	/**
	 * @deprecated inefficient
	 * @param FMRD
	 * @return
	 */
	public static String getURIgeneral(String FMRD) {
		int arg= -1, max= 0;
		for (int i = 0; i < DELIM_FMRD.length; i++) {
			String[] st= FMRD.split(Character.toString(DELIM_FMRD[i]));
			if (st.length== 1)
				continue;
			int p= 0;
			for (int j = 0; j < st.length; j++) {
				if (pattPE.matcher(st[j]).matches())
					++p;
			}
			if (st.length> p+ 1&& p> max) {
				arg= i;
				max= p;
			}
		}
		
		if (arg< 0)
			return FMRD;
		String[] st= FMRD.split(Character.toString(DELIM_FMRD[arg]));
		return st[0];
	}
	
	public CharSequence getUniqueDescriptor(CharSequence FMRD) {
		for (int i = 0; i < FMRD.length(); i++) {	// regexp inefficient
			if (FMRD.charAt(i)== '|')	// "\\|" in regexp
				return FMRD.subSequence(0, i);
		}
		return FMRD;
	}
	
	public byte getPairedEndInformation(CharSequence FMRD) {

		int last= 0;
		for (int i = 0; i < FMRD.length(); i++) {	// regexp inefficient
			if (FMRD.charAt(i)== '|') {	// "\\|" in regexp
				if (i- last== 2 && (FMRD.charAt(last)== 'p'|| FMRD.charAt(last)== 'P')) {
					if (FMRD.charAt(last+ 1)== '1')
						return (byte) 1;
					if (FMRD.charAt(last+ 1)== '2')
						return (byte) 2;
				}
				last= i+1;
			}
		}
		if (last!= FMRD.length()) {
			if (FMRD.length()- last== 2 && (FMRD.charAt(last)== 'p'|| FMRD.charAt(last)== 'P')) {
				if (FMRD.charAt(last+ 1)== '1')
					return (byte) 1;
				if (FMRD.charAt(last+ 1)== '2')
					return (byte) 2;
			}
		}
		
		return 0;
	}
	
	public boolean matesByName(String s1, String s2) {
		if (!getUniqueDescriptor(s1).equals(getUniqueDescriptor(s2))) {
			return false;
		}
		byte b1= getPairedEndInformation(s1), b2= getPairedEndInformation(s2);
		if (b1== 0|| b2== 0|| (b1- b2!= 0))
			return true;
		return false;
	}

	public static String getMM(String FMRD) {
		String[] ss= FMRD.split(Character.toString(DELIM_FMRD[0]));
		for (int i = 0; i < ss.length; i++) {
			if (pattMM.matcher(ss[i]).matches())
				return ss[i];
		}
		return null;
	}

	public boolean isApplicable(CharSequence descriptor) {
		// TODO Auto-generated method stub
		return false;
	}

	public boolean isPairedEnd(CharSequence descriptor) {
		if (getPairedEndInformation(descriptor)== 0)
			return false;
		return true;
	}

	public boolean allowsPend() {
		// TODO Auto-generated method stub
		return false;
	}

	public boolean allowsStranded() {
		// TODO Auto-generated method stub
		return false;
	}

	public byte getStrand(CharSequence descriptor) {
		// TODO Auto-generated method stub
		return 0;
	}

	public boolean isStranded(CharSequence descriptor) {
		// TODO Auto-generated method stub
		return false;
	}

}

