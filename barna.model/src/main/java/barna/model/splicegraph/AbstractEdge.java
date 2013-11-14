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

package barna.model.splicegraph;

import barna.model.Transcript;

import java.util.Comparator;
import java.util.Vector;

/**
 * Abstract class for edge types.
 * 
 * @author Micha Sammeth (gmicha@gmail.com)
 *
 */
public abstract class AbstractEdge {


	/**
	 * A comparator class to compare edges by their <code>tail</code>
	 * and <code>head</code> coordinates.
	 * 
	 * @author Micha Sammeth (gmicha@gmail.com)
	 *
	 */
	public static class PositionComparator implements Comparator<AbstractEdge> {
		
		/**
		 * Returns <code>-1</code> respectively <code>1</code> if the upstream
		 * site of the first edge is left respectively right of the upstream 
		 * site of the second edge, or, if both edges are delimited by the 
		 * same upstream site, if the downstream site of the first edge is 
		 * left respectively right of the second edge. If both edges are delimited
		 * by the same sites, <code>0</code> is returned.
		 * @return <code>(-1)</code>, <code>1</code> or <code>0</code> depending
		 * on whether the first edge is left, right or at the same coordinates as
		 * the second edge. 
		 */
		public int compare(AbstractEdge arg0, AbstractEdge arg1) {
			int s11= arg0.getTail().getSite().getPos(),
				s12= arg0.getHead().getSite().getPos(),
				s21= arg1.getTail().getSite().getPos(),
				s22= arg1.getHead().getSite().getPos();
			
			if (s11< s21)
				return -1;
			if (s21< s11)
				return 1;
			if (s12< s22)
				return -1;
			if (s22< s12)
				return 1;
			return 0;
		}
	}
	
	/**
	 * Default comparator for edges by positions of the delimiting sites.
	 */
	static PositionComparator defaultPositionComparator= new PositionComparator();
	

	/**
	 * Transcripts supporting the edge.
	 */
	protected long[] transcripts;
	
	/**
	 * Downstream site delimiting the edge.
	 */
	protected Node head;
	
	/**
	 * Upstream site delimiting the edge.
	 */
	protected Node tail;

	/**
	 * Superedges that include <code>this</code> edge.
	 */
	protected Vector<SuperEdge> superEdges= null;


	/**
	 * Sets the transcript support of <code>this</code> edge.
	 * 
	 * @param transcripts new transcript support of the edge
	 */
	public void setTranscripts(long[] transcripts) {
		this.transcripts = transcripts;
	}

	/**
	 * Returns the transcript support of <code>this</code> edge.
	 * @return transcript support of <code>this</code> edge
	 */
	public long[] getTranscripts() {
		return transcripts;
	}

	/**
	 * Returns the downstream site delimiting the edge.
	 * @return the downstream site delimiting the edge.
	 */
	public Node getHead() {
		return head;
	}

	/**
	 * Sets the downstream site delimiting the edge.
	 * @param head the downstream site delimiting the edge
	 */
	public void setHead(Node head) {
		this.head = head;
	}

	/**
	 * Returns the upstream site delimiting the edge.
	 * @return the upstream site delimiting the edge.
	 */
	public Node getTail() {
		return tail;
	}

	/**
	 * Sets the upstream site delimiting the edge.
	 * @param tail the upstream site delimiting the edge
	 */
	public void setTail(Node tail) {
		this.tail = tail;
	}

	/**
	 * Retrieves the last included position on the left
	 * respectively the right flank of the edge.
	 * @param tail <code>true</code> if left flank, 
	 * <code>false</code> otherwise
	 * @return genomic position of the last included 
	 * flank position
	 */
	public int getDelimitingPos(boolean tail) {

		if (tail) {
			int pos= getTail().getSite().getPos();
			if (getTail().getSite().isRightFlank()&&
                    ((isExonic()&& getHead().getSite().isRightFlank())          // AD
                        || (isIntronic()&& getHead().getSite().isLeftFlank()))) // intron
				++pos;
			return pos;
		} else {
			int pos= getHead().getSite().getPos();
			if (getHead().getSite().isLeftFlank()&&
                    ((isExonic()&& getTail().getSite().isLeftFlank())               // AA
                        || (isIntronic()&& getTail().getSite().isRightFlank()))) {  // intron
				--pos;
            }
			return pos;
		}
	}

	/**
	 * Returns leftmost and rightmost transcript coordinate covered 
	 * by <code>this</code> edge.
	 * @param t the transcript to which the coordinates refer
	 * @param readLen the read length
	 * @return a tuple {start,end} between which the edge extends
	 * on the transcript
	 */
	public abstract int[] getFrac(Transcript t, int readLen);

	/**
	 * Returns leftmost and rightmost transcript coordinate covered 
	 * by <code>this</code> edge.
	 * @param t the transcript to which the coordinates refer
	 * @param readLen the read length
	 * @param dir the directionality
	 * @return a tuple {start,end} between which the edge extends
	 * on the transcript
	 */
	public abstract int[] getFrac(Transcript t, int readLen, byte dir);
	
	public abstract int getNtLength(Transcript tx, int mapLength);
	
	/**
	 * Returns whether an intronic stretch is represented.
	 * @return <code>true</code> if all parts of this edge are intronic,
	 * <code>false</code> otherwise.
	 */
	public abstract boolean isIntronic();
	
	/**
	 * Returns whether an exonic stretch is represented.
	 * @return <code>true</code> if all parts of this edge are exonic,
	 * <code>false</code> otherwise.
	 */
	public abstract boolean isExonic();

    /**
     * Returns whether an all-intronic stretch is represented.
     * @return <code>true</code> if all parts of this edge are all-intronic,
     * <code>false</code> otherwise.
     */
    public abstract boolean isAllIntronic();

    /**
	 * Returns super-edges that include <code>this</code> edge.
	 * @return super-edges that include <code>this</code> edge
	 */
	public Vector<SuperEdge> getSuperEdges() {
		return superEdges;
	}

	/**
	 * Add another super-edge to the set of super-edges that 
	 * include <code>this</code> edge.
	 * @param superEdge the super-edge to be added
	 */
	public void addSuperEdge(SuperEdge superEdge) {
		if (superEdges== null)
			superEdges= new Vector<SuperEdge>(1,1);
		for (int i = 0; i < superEdges.size(); i++) {
			if (superEdges.elementAt(i)== superEdge)
				return;
		}
		superEdges.add(superEdge);
	}

	/**
	 * Return the length described by the current edge instance.
	 * @return length of the stretch described by <code>this</code> edge
	 */
	public abstract int length();
	
	@Override
	/**
	 * Calculate hash codes by string conversion.
	 */
	public int hashCode() {		
		return toString().hashCode();
	}

	/**
	 * Compute effective length of edge
	 * @param t a base transcript
	 * @param dir directionality
	 * @param mapLenMax maximum length of mappings
	 * @return the effective length of the edge
	 */
	public abstract int getEffLength(Transcript t, byte dir, int mapLenMax);
	
	/**
	 * Returns the genomic position.
	 * @param sense directionality
	 * @param start	exonic position
	 * @param minMapLen minimum length of mappings
	 * @param maxMapLen maximum length of mappings
	 * @return the corresponding genomic position
	 */
	public abstract int getGpos(boolean sense, boolean start, int minMapLen, int maxMapLen);
}
