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

package fbi.genome.sequencing.rnaseq.reconstruction;

import fbi.commons.tools.ArrayUtils;
import fbi.genome.model.Transcript;
import fbi.genome.model.splicegraph.Graph;

import java.io.File;
import java.io.FileOutputStream;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Vector;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

public class TProfileFunction {

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
	
	public static int[] BIN_LEN= new int[] {1000, 3000, 10000};
	public static int[] BIN_EXP= new int[] {10, 100, 1000};
	
	FluxCapacitor capacitor;
	public TProfileFunction(FluxCapacitor cap) {
		profiles= new Vector<TProfile>();
		this.capacitor= cap;
	}
	
	public void finish() {
		System.err.println("[ERROR] fatal call to finish()");
		profis= (TProfile[]) ArrayUtils.toField(profiles);
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
		
	public TProfile[] getTProfiles() {
		
		if (profis == null) {
			profis= (TProfile[]) ArrayUtils.toField(profiles);
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

	public TSuperProfile getMasterProfile(Graph g, Transcript t, int len, 
			int readLen, int[] insertMinMax, int reads, boolean strandSpecific, boolean pairedEnd) {

		int lenBin= 0;
		if (t.getExonicLength()> LEN_LO)
			++lenBin;
		if (t.getExonicLength()> LEN_UP)
			++lenBin;
		
		float rpkm= capacitor.calcRPKM(reads, len);
		int expBin= 0;
		if (rpkm> EXP_LO)
			++expBin;
		if (rpkm> EXP_UP)
			++expBin;
		
		TSuperProfile supi= getMasterProfiles(strandSpecific, pairedEnd, insertMinMax, readLen)[lenBin][expBin];
		//supi.set(g, t);	// RGASP bug
		return supi;
	}
	
	TSuperProfile[][] masters= null;
	public TSuperProfile[][] getMasterProfiles(boolean strandSpec, boolean pairedEnd, int[] instertMinMax, int readLen) {
		if (masters == null) {
			masters = new TSuperProfile[3][3];
			for (int i = 0; i < masters.length; i++) 
				for (int j = 0; j < masters[i].length; j++) 
					masters[i][j]= new TSuperProfile();   
			
			if (getTProfiles()== null|| getTProfiles().length== 0) {
				for (int i = 0; i < masters.length; i++) {
					for (int j = 0; j < masters[i].length; j++) {
						TProfile profi= new TProfile(i==0?LEN_LO:(i==1?LEN_UP:10000), strandSpec, pairedEnd);
						profi.fill(instertMinMax, readLen);
						masters[i][j].addProfile(profi);
					}
				}
			}
			
			for (int i = 0; i < getTProfiles().length; i++) {
				TProfile profi= getTProfiles()[i];
				int len= profi.getLength();
				int reads= profi.getReads();
				int lenBin= 0;
				if (len> LEN_LO)
					++lenBin;
				if (len> LEN_UP)
					++lenBin;
				
				float rpkm= capacitor.calcRPKM(reads, len);
				int expBin= 0;
				if (rpkm> EXP_LO)
					++expBin;
				if (rpkm> EXP_UP)
					++expBin;
				
				masters[lenBin][expBin].addProfile(profi);
				
				getTProfiles()[i]= null;	// remove
			}
			System.gc();
		}

		return masters;
	}
	
	
	
	public TProfile[] getProfis() {
		return profis;
	}
	
	
}
