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

/**
 * test
 */
package barna.model.splicegraph;

import barna.model.Transcript;
import barna.model.constants.Constants;

/**
 * <code>AbstractEdge</code> implementation for representing 
 * a continuous stretch between two sites.
 * 
 * @author Micha Sammeth (gmicha@gmail.com)
 *
 */
public class SimpleEdge extends AbstractEdge {
	
	/**
	 * Confidence level of the edge's transcript support.
	 */
	protected byte type= Transcript.ID_SRC_UNDEFINED;
	
	/**
	 * Flag indicating optimized edges during AStalavista graph contraction.
	 */
	boolean contracted= false;
	boolean processed= false;
	boolean valid= true;
	boolean exonic= false;
	
	public static String getStringRep(Node v, Node w) {
		return v.getSite().toString()+w.getSite().toString();
	}
	
	String stringRep;
	Partition partition;
	public SimpleEdge(Node newTail, Node newHead) {
		this.tail= newTail;
		this.head= newHead;
		tail.addOutEdge(this);
		head.addInEdge(this);
	}
	
	public String toString() {
		if (stringRep == null) 
			stringRep = getStringRep(getTail(), getHead());

		return stringRep;
	}
	
	/**
	 * Absolute length in Nt, wrapping length() for exons.  
	 */
	public int getNtLength(Transcript tx, int mapLength) {
		return length(); 
	}
	
	/**
	 * can over-/underflow edge position when not possible
	 * @param sense
	 * @param start
	 * @param minMapLen
	 * @param maxMapLen
	 * @return
	 */
	public int getGpos(boolean sense, boolean start, int minMapLen, int maxMapLen) {
		if (sense) {
			if (start)
				return Math.min(getTail().getSite().getPos(), getHead().getSite().getPos()- minMapLen+ 1);
			else
				return getHead().getSite().getPos()- minMapLen+ 1;
		} else {
			if (start)
				return getTail().getSite().getPos()+ minMapLen- 1;
			else
				return Math.max(getHead().getSite().getPos(), getTail().getSite().getPos()+ minMapLen- 1);	
		}
	}
	
	@Override
	public boolean equals(Object obj) {
		SimpleEdge e= (SimpleEdge) obj;
		if (getTail().equals(e.getTail())&& getHead().equals(e.getHead())
				&& SplicingGraph.equalSet(getTranscripts(), e.getTranscripts())
				&& isExonic()== isExonic()		// multiple edges exonic, intronic for eg intron retention
				&& isIntronic()== isIntronic())	
			return true;
		return false;
	}

	public boolean isContracted() {
		return contracted;
	}

	public void setContracted(boolean contracted) {
		this.contracted = contracted;
	}

	public boolean isProcessed() {
		return processed;
	}

	public void setProcessed(boolean processed) {
		this.processed = processed;
	}

	public int length() {
		int len= getHead().getSite().getPos()- getTail().getSite().getPos()+ 1;
		// correct for exon fragments that overlap in exonic positions
		if (getTail().getSite().isRightFlank())
			--len;	
		if (getHead().getSite().isLeftFlank())
			--len;
		if (isExonic())
			assert(len>=0);	// can be 0-edges, see chr1:19,108,320-19,108,335
		return len;
	}
	
	public boolean isIntronic() {
		if (exonic
			|| getTail().getSite().getPos()== Integer.MIN_VALUE
			|| getHead().getSite().getPos()== Integer.MAX_VALUE)
			return false;
		return true;
	}
	
	public boolean isExonic() {		
		return exonic;
	}

	public static PositionComparator getDefaultPositionComparator() {
		return defaultPositionComparator;
	}

	public void setExonic(boolean exonic) {
		this.exonic = exonic;
	}
	
	/**
	 * @deprecated use with directionality
	 * @param t
	 * @param readLen
	 * @return
	 */
	public int[] getFrac(Transcript t, int readLen) {
		int pos= getTail().getSite().getPos();
		int start= t.getExonicPosition(pos);
		if (getTail().getSite().isRightFlank())
			++start;
	//		if (!e.getTail().getSite().isTSS())
	//			start+= readLen- 1;	// narrow for sedge ranges
		pos= getHead().getSite().getPos();
		int end= t.getExonicPosition(pos);
		if (getHead().getSite().isLeftFlank())
			--end;
		
		// TODO: this is questionable, actually we do not accept reads there
		// but due to profile stretching, we can have expectations there
		// so we should allow this region in the last exon
		//if (getHead().getSite().getPos()!= t.get3PrimeEdge())	// !getHead().getSite().isTES()
			end-= readLen- 1;
	
		return new int[] {start,end};
	}
	
	
	public int[] getFrac(Transcript t, int readLen, byte dir) {
		int pos= getTail().getSite().getPos();
		int start= t.getExonicPosition(pos);
		if (getTail().getSite().isRightFlank())
			++start;
	//		if (!e.getTail().getSite().isTSS())
	//			start+= readLen- 1;	// narrow for sedge ranges
		pos= getHead().getSite().getPos();
		int end= t.getExonicPosition(pos);
		if (getHead().getSite().isLeftFlank())
			--end;
		
		if (dir== Constants.DIR_FORWARD)
			end-= readLen- 1;
		else if (dir== Constants.DIR_BACKWARD)
			start+= readLen- 1;
		// in case of both it is ok
		
		return new int[] {start,end};
	}

	/**
	 * Maximum mapping length that can map to the edge. 
	 * @param maxMapLength
	 * @return
	 */
	public int getMapLength(Transcript tx, int maxMapLength) {
		return Math.min(length(), maxMapLength);
	}
	
	/**
	 * Number of different locations with length mapLen.
	 * @param dir
	 * @param mapLenMax
	 * @return
	 */
	public int getEffLength(Transcript tx, byte dir, int mapLenMax) {
		
		return (length()- getMapLength(tx, mapLenMax));
	}
}
