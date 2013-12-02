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

package barna.model.rna;

import barna.commons.RandomFactory;

import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.Random;

public class SolexaPairedEndDescriptor implements ReadDescriptor {

	public byte getPairedEndInformation(CharSequence descriptor) {
		if (descriptor== null)
			return (byte) 0;
		for (int i = descriptor.length()-2; i >= 0; --i) {	// has to have at least one after
			if (descriptor.charAt(i)== '/') {	// "\\|" in regexp
				if (descriptor.charAt(i+1)== '1')
					return (byte) 1;
				if (descriptor.charAt(i+1)== '2')
					return (byte) 2;
			}
		}
		return (byte) 0;
	}

	public CharSequence getUniqueDescriptor(CharSequence descriptor) {
		if (descriptor== null)
			return null;
		for (int i = descriptor.length()-2; i >= 0; --i)
			if (descriptor.charAt(i)== '/')
				return descriptor.subSequence(0, i);
		return descriptor;
	}

	public boolean isApplicable(CharSequence cs) {
		return true;
	}

	public boolean isPairedEnd(CharSequence cs) {
		if (getPairedEndInformation(cs)== 0)
			return false;
		return true;
	}

	public static final int LANES_PER_FLOWCELL= 8,
		TILES_PER_LANE= 330;
	
	public static final String COLON= ":";
	static int nrReadsAssignedToLanes= 0, nrReadsThisLane= 0, nrLanes= 0;
	private static Random rndSplitter;
	public static final byte NO_MATE= (byte) 0, MATE_1= (byte) 1, MATE_2= (byte) 2;
	/**
	 * 
	 * @param readNr 1-based
	 * @return
	 */
	public static String generateID(int readNr, byte mateInfo) {
/*		>SLXA-B3_604:6:1:533:275
		denotes an unpaired read. 
		It is from lane 6 tile 1 of run 604 on machine SLXA-B3, 
		and the (X,Y) coordinates of the cluster on the tile (in pixel units) 
		are (533,275). 
		Taken together these data are sufficient to uniquely identify a cluster 
		- importantly for multi-run projects, this uniqueness should extend across runs, 
		although this is of course reliant upon machine naming and run numbering being 
		done in a sensible way. 

		>SLXA-B3_604:6:1:533:275/1
		would denote read 1 of a paired read and

		>SLXA-B3_604:6:1:533:275/2
		denotes read 2 of a paired read.
*/
		if (readNr> nrReadsAssignedToLanes) {
			++nrLanes;
			nrReadsThisLane= getRndReadNrLane();
			nrReadsAssignedToLanes+= nrReadsThisLane;
			if (nrLanes> 8) {
				++nrFlowCell;
				runID= null;
				nrLanes%= 8;
			}
		}
		int relReadNr= readNr% (nrReadsAssignedToLanes- nrReadsThisLane);	// nr in current lane
		int tileNr= (int) (relReadNr/ TILES_PER_LANE)+ 1;
		relReadNr-= (relReadNr/ TILES_PER_LANE)* (tileNr- 1);	// nr in current tile
		int x= 0, y= 0;
		if (relReadNr> 0) {
			if (rndSplitter== null)
				rndSplitter= RandomFactory.get();
			x= rndSplitter.nextInt(relReadNr)+ 1;	// [1..relReadNr[
			y= relReadNr- x;
		}
		
		String s= SLXA_PREFIX+ getRunPrefix()+ COLON+ nrLanes+ COLON+ tileNr+ COLON+ x+ COLON+ y;
		if (mateInfo!= NO_MATE)
			s+= mateInfo== MATE_1?"/1":"/2";
		return s;
	}
	
	static Random laneRandom;
	public static final int AVG_NR_READS_PER_LANE= 12345678;
	private static int getRndReadNrLane() {
		if (laneRandom== null)
			laneRandom= RandomFactory.get();
		return AVG_NR_READS_PER_LANE+ (int) (0.1f* AVG_NR_READS_PER_LANE* laneRandom.nextGaussian());
	}
	public static final String SLXA_PREFIX= "SLXA_";
	static String runID= null;
	static int nrFlowCell= 0;
	public static String getRunPrefix() {
		if (runID == null) {
			DateFormat dateFormat = new SimpleDateFormat("MMdd");
            Date date = new Date();
            runID= dateFormat.format(date);
            if (nrFlowCell> 0) {
            	int dd= Integer.parseInt(runID.substring(runID.length()-2, runID.length()));
            	dd+= nrFlowCell;
            	runID= runID.substring(0, runID.length()- 2)+ Integer.toString(dd);
            }
		}

		return runID;
	}

	public boolean allowsPend() {		
		return true;
	}

	public boolean allowsStranded() {
		return false;
	}

	public byte getStrand(CharSequence descriptor) {		
		return 0;
	}

	public boolean isStranded(CharSequence descriptor) {
		return false;
	}
	
}
