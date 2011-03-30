package fbi.genome.sequencing.rnaseq.reconstruction;

import fbi.genome.model.Transcript;
import fbi.genome.model.commons.MyArrays;
import fbi.genome.model.splicegraph.Graph;

import java.io.File;
import java.io.FileOutputStream;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Vector;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

public class Profile {

	static class TProfileByLengthComparator implements Comparator<TProfile> {
		//@Override
		public int compare(TProfile o1, TProfile o2) {
			if (o1.getLength()< o2.getLength())
				return -1;
			if (o1.getLength()> o2.getLength())
				return 1;
			return 0;
		}
	}
	
	static TProfileByLengthComparator defaultTProfileByLengthComparator= null;
	static TProfileByLengthComparator getDefaultTProfileByLengthComaparator() {
		if (defaultTProfileByLengthComparator == null) {
			defaultTProfileByLengthComparator = new TProfileByLengthComparator();
		}

		return defaultTProfileByLengthComparator;
	} 
	
	int binMinReads= 10000, binMinTranscripts= 100;
	float binMaxLengthDistance, binMaxExprDistance; // either abs or factor (when <10)
	Vector<TProfile> profiles;
	TProfile[] profis;
	public final static int LEN_LO= 1000;
	public final static int LEN_UP= 5000;
	public final static int EXP_LO= 10;
	public final static int EXP_UP= 100;
	
	// {1000, 5000}; // 
	public static int[] BIN_LEN= new int[] {500, 1000, 1500, 2000};	// 5 bins, good
	public static int[] BIN_EXP= new int[] {10, 100, 1000};
	
	FluxCapacitor capacitor;
	public Profile(FluxCapacitor cap) {
		profiles= new Vector<TProfile>();
		this.capacitor= cap;
	}
	
	public void finish() {
		System.err.println("[ERROR] fatal call to finish()");
		profis= (TProfile[]) MyArrays.toField(profiles);
		profiles= null;
		System.gc();
	}
	
	public void addProfile(TProfile profile) {
		((Vector<TProfile>) profiles).add(profile);
	}
	
	public TProfile getProfile(int len, boolean strandSpecific, boolean pairedEnd) {
		synchronized (profiles) {
			Vector<TProfile> proV= (Vector<TProfile>) profiles;
			for (int i = 0; i < proV.size(); i++) {
				if (proV.elementAt(i).length()== len)
					return proV.elementAt(i);
			}
			TProfile profile= new TProfile(len, strandSpecific, pairedEnd); 
			profiles.add(profile);
			return profile;
		}
	}
		
	TProfile[] getTProfiles() {
		
		if (profis == null) {
			profis= (TProfile[]) MyArrays.toField(profiles);
			if (profis!= null)
				Arrays.sort(profis, getDefaultTProfileByLengthComaparator());
		}

		return profis;
		
//		if (profiles instanceof Vector) {
//			Vector<TProfile> tpVec= (Vector<TProfile>) profiles;
//			TProfile[] profiles= new TProfile [tpVec.size()];
//			for (int i = 0; i < profiles.length; i++) {
//				profiles[i]= tpVec.elementAt(i);
//			}
//			Arrays.sort(profiles, getDefaultTProfileByLengthComaparator());
//			this.profiles= profiles;
//		}
//		return (TProfile[]) profiles;
	}
	
	public int getNrProfiles() {
		//if (profiles instanceof Vector)
			return ((Vector) profiles).size();
		//return ((TProfile[]) profiles).length;
	}
	
	public long getNrReadsInProfiles() {
		int x= -1;
//		if (profiles instanceof Vector)
			x= ((Vector) profiles).size();
//		else
//			x= ((TProfile[]) profiles).length;

		long cnt= 0;
		for (int i = 0; i < x; i++) {
//			if (profiles instanceof Vector)
				cnt+= ((Vector<TProfile>) profiles).elementAt(i).getReads();
//			else
//				cnt+= ((TProfile[]) profiles)[i].getReads();
		}
		return cnt;
	}
	
	public TProfile getSuperProfile_old(int len, int readLen, int reads, boolean strandSpecific, boolean pairedEnd) {
		
/*		float lenLo= (binMaxLengthDistance<= 10)?len/binMaxLengthDistance:len-binMaxLengthDistance,
				lenHi= (binMaxLengthDistance<= 10)?len*binMaxLengthDistance:len+binMaxLengthDistance,
				redLo= (binMaxExprDistance<= 10)?len/binMaxExprDistance:len-binMaxExprDistance,
				redHi= (binMaxExprDistance<= 10)?len*binMaxExprDistance:len+binMaxExprDistance;
		
		len= Math.max(0,len);
		TProfile dummy= new TProfile(len, strandSpecific, pairedEnd);	// , reads
		if (len<= 0)
			return dummy;
		int p= Arrays.binarySearch(getTProfiles(), dummy, getDefaultTProfileByLengthComaparator());
		if (p< 0)
			p= -(p+1);
		if (p== getTProfiles().length)
			--p;
		
		// try to fulfil both criteria
		if (p< 0)
			System.currentTimeMillis();
		dummy.addProfile(getTProfiles()[p], readLen);
		int counter= recruitProfiles(dummy, readLen, p, lenLo, lenHi, redLo, redHi);
		
		// try only length
		if (counter< binMinTranscripts|| dummy.getReads()< binMinReads) {	// getReadsInBin()
			counter+= recruitProfiles(dummy, readLen, p, lenLo, lenHi, -1, -1);
		}
		// take everything
		if (counter< binMinTranscripts|| dummy.getReads()< binMinReads) {	// getReadsInBin()
			counter+= recruitProfiles(dummy, readLen, p, -1, -1, -1, -1);
		}
		//supa.fill();
		return dummy;
*/
		return null;
	}

	private int recruitProfiles(TSuperProfile supa, int readLen, int[] insertMinMax, int p, float lenLo, float lenHi, float redLo, float redHi) {
		
		int counter= 0;
		for (int i = 1; i < Math.max(getTProfiles().length-p, p); i++) {
			int loIdx= p-i, upIdx= p+i;
			if (loIdx> 0
					&& (lenLo< 0|| getTProfiles()[loIdx].getLength()>= lenLo)
					&& (lenHi< 0|| getTProfiles()[loIdx].getLength()<= lenHi)
					&& (redLo< 0|| getTProfiles()[loIdx].getReads()>= redLo)
					&& (redHi< 0|| getTProfiles()[loIdx].getReads()<= redHi)) {
				supa.addProfile(getTProfiles()[loIdx]);
				++counter;
			}
			if (upIdx< getTProfiles().length
					&& (lenLo< 0|| getTProfiles()[upIdx].getLength()>= lenLo)
					&& (lenHi< 0|| getTProfiles()[upIdx].getLength()<= lenHi)
					&& (redLo< 0|| getTProfiles()[upIdx].getReads()>= redLo)
					&& (redHi< 0|| getTProfiles()[upIdx].getReads()<= redHi)) {
				supa.addProfile(getTProfiles()[upIdx]);
				++counter;
			}
			if (counter>= binMinTranscripts&& supa.getAllReads()>= binMinReads)
				break;
		}
		return counter;
	}

	public int getBinMinReads() {
		return binMinReads;
	}

	public void setBinMinReads(int binMinReads) {
		this.binMinReads = binMinReads;
	}

	public int getBinMinTranscripts() {
		return binMinTranscripts;
	}

	public void setBinMinTranscripts(int binMinTranscripts) {
		this.binMinTranscripts = binMinTranscripts;
	}

	public float getBinMaxLengthDistance() {
		return binMaxLengthDistance;
	}

	public void setBinMaxLengthDistance(float binMaxLengthDistance) {
		this.binMaxLengthDistance = binMaxLengthDistance;
	}

	public float getBinMaxExprDistance() {
		return binMaxExprDistance;
	}

	public void setBinMaxExprDistance(float binMaxExprDistance) {
		this.binMaxExprDistance = binMaxExprDistance;
	}
	
	public void output(File dir) {
		TProfile[] prof= getTProfiles();
		try {
		    ZipOutputStream zos = new ZipOutputStream(new FileOutputStream(dir));
		    HashSet<String> nameSet= new HashSet<String>();
			for (int i = 0; i < prof.length; i++) {
				String fName= Integer.toString(prof[i].getLength())+"_"+Integer.toString((int) prof[i].getReads());
				String name= fName+ ".tpf";
				for (int x= 0;nameSet.contains(name); ++x) {
					if (x> 0)
						System.currentTimeMillis();
					name= fName+"_"+Integer.toString(x)+".tpf";
				}
				nameSet.add(name);
				
			    ZipEntry ze= new ZipEntry(name);
			    zos.putNextEntry(ze);
			    zos.write(prof[i].toByteArray());
			    zos.closeEntry();
			}
			
		    zos.close();
			
		} catch (Exception e) {
			e.printStackTrace();
		}


	}

	/**
	 * @deprecated
	 * @param tlen
	 * @param rpk
	 * @return
	 */
	public UniversalMatrix getMatrix(int tlen, float rpk) {
		int lenBin= 0;
		if (tlen> LEN_LO)
			++lenBin;
		if (tlen> LEN_UP)
			++lenBin;
		
		int expBin= 0;
		if (rpk> EXP_LO)
			++expBin;
		if (rpk> EXP_UP)
			++expBin;
		
		UniversalMatrix m= getMasters()[lenBin]; // [expBin];
		return m;
	}
	
	public UniversalMatrix getMatrix(int tlen) {
		int lenBin= Arrays.binarySearch(BIN_LEN, tlen);
		if (lenBin< 0)
			lenBin= -(lenBin+ 1);
		
		UniversalMatrix m= getMasters()[lenBin];
		return m;
	}
	
	public UniversalMatrix[] masters= null;
	public UniversalMatrix[] getMasters() {
		if (masters == null) {
			masters = new UniversalMatrix[BIN_LEN.length+ 1]; // [3]
			for (int i = 0; i < masters.length; i++) {
				int mlen= i== 0? BIN_LEN[0]/ 2: 
					i>= BIN_LEN.length? BIN_LEN[BIN_LEN.length- 1]: 
						BIN_LEN[i- 1]+ ((BIN_LEN[i]- BIN_LEN[i- 1])/ 2);
				masters[i]= new UniversalMatrix(mlen);
//				for (int j = 0; j < masters[i].length; j++) 
//					masters[i][j]= new UniversalMatrix(mlen);
			}
			
		}

		return masters;
	}
	
	
	
	public TProfile[] getProfis() {
		return profis;
	}

	public int fill() {
		int sum= 0;
		UniversalMatrix[] masters= getMasters();
		for (int i = 0; i < masters.length; i++) {
			sum+= masters[i].fill();
//			for (int j = 0; j < masters[i].length; j++) {
//				sum+= masters[i][j].fill();
//			}
		}
		return sum;
	}
	
	
}
