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

package barna.flux.capacitor.profile;

import barna.flux.capacitor.matrix.UniversalMatrix;

import java.util.Arrays;

public class Profile {

	int binMinReads= 10000, binMinTranscripts= 100;
	float binMaxLengthDistance, binMaxExprDistance; // either abs or factor (when <10)

    private UniversalMatrix[] masters= null;

	public final static int LEN_LO= 1000;
	public final static int LEN_UP= 5000;
	public final static int EXP_LO= 10;
	public final static int EXP_UP= 100;

	// {1000, 5000}; //
	public static int[] BIN_LEN= new int[] {500, 1000, 1500, 2000};	// 5 bins, good

    public Profile() {
        mappingStats = new MappingStats();
        coverageStats = new CoverageStats();
	}
	
		/**
	 * @deprecated
	 * @param tlen
	 * @param rpk
	 * @return
	 */
	public UniversalMatrix getMatrix(int tlen, double rpk) {
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

        if (lenBin> 0&& lenBin< BIN_LEN.length&& Math.abs(tlen- BIN_LEN[lenBin- 1])< Math.abs(tlen- BIN_LEN[lenBin]))
            --lenBin;

		UniversalMatrix m= getMasters()[lenBin];
		return m;
	}

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

    public void setMasters(UniversalMatrix[] masters) {
        this.masters = masters;
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
	
	private MappingStats mappingStats;

    /**
     * Coverage
     */
    private CoverageStats coverageStats;

    public CoverageStats getCoverageStats() {
        return coverageStats;
    }

    public void setCoverageStats(CoverageStats coverageStats) {
        this.coverageStats = coverageStats;
    }

    public MappingStats getMappingStats() {
        return mappingStats;
    }

    public void setMappingStats(MappingStats mappingStats) {
        this.mappingStats = mappingStats;
    }

    @Override
    public boolean equals(Object other) {
        if (!(other instanceof Profile))
            return false;
        Profile otherProfile = (Profile)other;
        if (this.masters.length != otherProfile.masters.length)
            return false;
        for (int i = 0; i < this.masters.length; i++) {
            if (this.masters[i].getLength() != otherProfile.masters[i].getLength())
                return false;
            if (this.masters[i].sense.length != otherProfile.masters[i].sense.length)
                return false;
            for (int j = 0; j < this.masters[i].sense.length; j++) {
                if (this.masters[i].sense[j] != otherProfile.masters[i].sense[j])
                    return false;
            }
            if (this.masters[i].asense.length != otherProfile.masters[i].asense.length)
                return false;
            for (int j = 0; j < this.masters[i].asense.length; j++) {
                if (this.masters[i].asense[j] != otherProfile.masters[i].asense[j])
                    return false;
            }
            if (this.masters[i].suma != otherProfile.masters[i].suma)
                return false;
            if (this.masters[i].sums != otherProfile.masters[i].sums)
                return false;
        }
        return true;
    }

}
