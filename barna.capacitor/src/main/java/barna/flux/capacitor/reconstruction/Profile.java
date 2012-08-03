/*
 * Copyright (c) 2010, Micha Sammeth
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * The names of its contributors may be not used to endorse or promote
 *       products derived from this software without specific prior written
 *       permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL MICHA SAMMETH BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

package barna.flux.capacitor.reconstruction;

import barna.commons.utils.ArrayUtils;

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
		
	TProfile[] getTProfiles() {
		
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
     * Assigns a bias matrix according to the spliced length
     * and the expression level (i.e., RPK reads per kilobase)
     * detected for that transcript.
     *
	 * @param tlen length of processed transcript
	 * @param rpk reads per kilobase in that transcript
     * @deprecated unused
     * @return unused
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
