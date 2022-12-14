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

import barna.model.SpliceSite;
import barna.model.Transcript;
import barna.model.constants.Constants;

import java.util.Comparator;

public class SuperEdge extends AbstractEdge {
	public static class EdgeByTailComparator implements Comparator<SimpleEdge> {
		public int compare(SimpleEdge arg0, SimpleEdge arg1) {
			if (arg0.getTail().getSite().getPos()< arg1.getTail().getSite().getPos())
				return -1;
			if (arg0.getTail().getSite().getPos()> arg1.getTail().getSite().getPos())
				return 1;
			return 0;
		}
	}
	static EdgeByTailComparator defaultEdgeByTailComparator= new EdgeByTailComparator();
	public static class EdgeTuple {
		public SimpleEdge[] edges;
		public EdgeTuple(SimpleEdge[] edges) {
			this.edges= edges;
		}
		@Override
		public boolean equals(Object obj) {
			EdgeTuple other= (EdgeTuple) obj;
			//Arrays.equals(a, a2);	// take oID 4 speed
			if (edges.length!= other.edges.length)
				return false;
			for (int i = 0; i < edges.length; i++) 
				if (edges[i]!= other.edges[i])
					return false;			
			return true;
		}
		
		@Override
		public int hashCode() {
			int val= 0;
			for (int i = 0; i < edges.length; i++) {
				val+= edges[i].hashCode();
			}
			return val; //(edges[0].toString()+ edges[1].toString()).hashCode();
		}
	}

	/**
	 * Edge set comprised by <code>this</code> instance of
	 * Superedge.
	 */
	AbstractEdge[] edges= null;
	
	/**
	 * Flag marking non-contiguous set of edges.
	 */
	boolean pairedEnd= false;

    /**
     * The type <code>this</code> edge is representing.
     */
    byte type= TYPE_NA;

    /**
     * Uninitialized type
     */
    public final static byte TYPE_NA= -1;

    /**
     * Simple exon or segment of an exon
     */
    public final static byte TYPE_SIMPLE= 0;

    /**
     * Superedge describign a continuous
     * sequence of exonic segments
     */
    public final static byte TYPE_COMPOSED= 1;

    /**
     * Superedge describing a tuple of edges
     * which may be identical
     */
    public final static byte TYPE_PAIRED= 2;

    /**
     * Set of any type of (super-)edges.
     */
    public final static byte TYPE_SET= 10;


    /**
     * Constructor to create anything but an edgeset.
	 * @param edges vector of edges
	 * @param supp vector of transcripts supporting the edges
	 * @param pend flag to indicate paired-end reads
	 */
	public SuperEdge(AbstractEdge[] edges, long[] supp, boolean pend) {

		setEdges(edges);// this(edges);	// also sets transcript set
		setTranscripts(supp);

        if(edges.length== 1)
            this.type= TYPE_SIMPLE;
        else {
            if (pend)
                this.type= TYPE_PAIRED;
            else
                this.type= TYPE_COMPOSED;
        }
	}

    /**
     * Constructor for an edgeset.
     * @param edges edges that the set comprises
     */
    public SuperEdge(AbstractEdge[] edges) {
        setEdges(edges);
        this.type= TYPE_SET;
    }

        public AbstractEdge[] getEdges() {
		return edges;
	}
	
	public static int getFirstEJ(AbstractEdge[] edges) {
		int pos= edges[0].getHead().getSite().getPos();
		if (edges[0].getHead().getSite().isLeftFlank())
			--pos;
		return pos;
	}
	
	public int getFirstEJ() {
		return getFirstEJ(edges);
	}
		
	
	public static int getLastEJ(AbstractEdge[] edges) {
        int pos= edges[edges.length-1].getHead().getSite().getPos();
        if (edges[edges.length-1].getTail().getSite().isRightFlank())
            ++pos;
        return pos;
	}
	
	public int getLastEJ() {
		return getLastEJ(edges);
	}

    public int countEJ() {
        return countEJ(edges);
    }

    private int countEJ(AbstractEdge[] edges) {
        int count = 0;
        for (int i = 0; i<edges.length-1;i++) {
            if (edges[i].getHead().getSite().isDonor() && edges[i+1].getTail().getSite().isAcceptor()) {
                    count++;
            }
        }
        return count;
    }

    public boolean isSet() {
        return type== TYPE_SET;
    }

    public boolean isSpliceJunction() {
        if (isSet())
            return false;

        if (this.countEJ() > 0)
            return (this.type== TYPE_PAIRED? false: true);
        return false;
    }

    public boolean isIntronic() {
        if (isSet())
            return false;
        for(AbstractEdge e : edges)
            if(e.isIntronic())
                return true;
		return false;
	}
	
	public boolean isExonic() {

        if (isSet())
            return false;

		for (int i = 0; i < edges.length; i++) 
			if(!edges[i].isExonic())
				return false;
		return true;
	}

    public boolean isAllIntronic() {
        if (isSet())
            return false;

        boolean b= true;
        for(AbstractEdge e : edges)
            if(e.isAllIntronic())
                b&= true;
            else
                b&=false;
        return b;
    }
	
	public void setEdges(AbstractEdge[] edges) {
//		long[] trpts= edges[0].getTranscripts();
		for (int i = 0; i < edges.length; i++) 
			edges[i].addSuperEdge(this);
//			if (i> 0) {
//				trpts= Graph.intersect(trpts, edges[i].getTranscripts());
				
				// !!!super_edges
				// EJ only exist with epos(before)= epos(after)- 1
//			}
			// yes, intronic edges can interrupt an superedge, .. and also pair-ends
	//			if (i> 0&& edges[i-1].getHead()!= edges[i].getTail())
	//				System.err.println("ERROR: super-edge with non-continuous edge-set.");			
//		}
		this.edges = edges;
//		if (edges[0].toString().equals("22409639-22409798^22410529-22410563[22410563[22410608^")
//				&& edges[1].toString().equals("22410529-22410563[22410563[22410608^"))
//			System.currentTimeMillis();
		tail= edges[0].getTail();
		head= edges[edges.length-1].getHead();
//		setTranscripts(trpts);
	}

	public static EdgeByTailComparator getDefaultEdgeByTailComparator() {
		return defaultEdgeByTailComparator;
	}

	/**
	 * length that spans from last nt before 1st splice junction
	 * to nt after last splice junction
	 */
	public int length() {
//		int len= edges[edges.length-1].getTail().getSite().getPos()- 
//			edges[0].getHead().getSite().getPos()+ 1;
		int len= 2;	// border nts
		for (int i = 1; i < edges.length-1; i++) 
			len+= edges[i].length();		
		
		return len;
	}
	
	@Override
	public int getNtLength(Transcript tx, int mapLength) {
		if (isPend()) {
			System.err.println("what about pends");
			return -1;
			
		} else {
			int first= tx.getExonicPosition(edges[0].getHead().getSite().getPos()),
				last= tx.getExonicPosition(edges[edges.length- 1].getTail().getSite().getPos());
			return ((first+ mapLength)- (last- mapLength));
		}
	}
	
	public int getEffLength(byte dir, int mapLenMax) {

		//System.err.println("check effLen");
		assert(dir== Constants.DIR_FORWARD^ dir== Constants.DIR_BACKWARD);
		int effLen= -1;
		if (isPend()) {	
			
			effLen= dir== Constants.DIR_FORWARD? 
					edges[0].getEffLength(dir, mapLenMax):
					edges[edges.length- 1].getEffLength(dir, mapLenMax);
		
		} else {	// ej / sj
			int interNt= 1;
			for (int i = 1; i < edges.length- 1; i++) 
				interNt+= edges[i].length();
            // delimit possible read positions (slots)
			effLen= Math.min(mapLenMax- interNt- 1, // max slots with read overlapping at least 1nt of first/last edge
                    Math.min(edges[0].length(), edges[edges.length- 1].length()));  // or the size of first/last edge
		}
		
		return effLen;
	}
	
	@Override
	/**
	 * compare additionally the set of edges
	 */
	public boolean equals(Object obj) {
		//if ((!(obj instanceof SuperEdge))|| (!super.equals(obj)))
		if (!(obj instanceof SuperEdge))
			return false;
		SuperEdge se= ((SuperEdge) obj);
		if (isPend()!= se.isPend())
			return false;
		AbstractEdge[] e= se.getEdges();
		if (e.length!= edges.length)
			return false;
		for (int i = 0; i < e.length; i++) {
			if (e[i]!= edges[i])
				return false;
		}
		return true;
	}

	public String toString() {
		StringBuffer sb= new StringBuffer();	// getTail().getSite().toString()
		for (int i = 0; i < edges.length; i++) {
//			if (edges[i].getTail()!= n)
				sb.append(edges[i].getTail().getSite().toString());
			sb.append(edges[i].getHead().getSite().toString());
		}
		
		// dont change at running time, hashCode uses it
		if(isPend())
			sb.append("PE");
		
		return sb.toString();
	}

	public boolean isPend() {
		return type== TYPE_PAIRED;
	}

	public static int[] getPEfrac(Transcript t, int readLen, AbstractEdge[] edges, byte dir) {
		assert(edges.length== 2);
		int[] left= edges[0].getFrac(t, readLen,dir);
		int[] right= edges[1].getFrac(t, readLen,dir);
		
		return new int[] {left[0], left[1], right[0], right[1]};
	}
	
	public static int[] getFrac(Transcript t, int readLen, AbstractEdge[] edges, byte dir) {
		
		int[] res= null;
		int start= t.getExonicPosition(getFirstEJ(edges));
		int end= t.getExonicPosition(getLastEJ(edges));
		
		// TODO changed from (start-rDelta)
		// int eDelta= end- start+ 1;	
		// int rDelta= readLen- eDelta;
		// start= Math.max(0, start-rDelta);

		int genFirst= edges[0].getTail().getSite().getPos();
		int firstPos= t.getExonicPosition(genFirst);
		if (edges[0].getTail().getSite().isRightFlank())
			++firstPos;
		
		SpliceSite endSite= edges[edges.length-1].getHead().getSite();
		int genLast= endSite.getPos();
		int lastPos= t.getExonicPosition(genLast);
		if (endSite.isLeftFlank())
			--lastPos;
		
		int first= firstPos, last= lastPos;
		if (dir== Constants.DIR_FORWARD) {
			last= Math.min(start, lastPos-(readLen-1)); 
		} else if (dir== Constants.DIR_BACKWARD) {
			first= Math.max(firstPos, end-(readLen-1));
		}
		// in case of both its ok

		res= new int[] {first,last};	// TODO start,end

		return res;

	}
	
	public int getGpos(boolean sense, boolean start, int minMapLen, int maxMapLen) {

		assert(edges.length>= 2);
		
		// paired-end
		if (this.type== TYPE_PAIRED) {
			if (sense)
				return edges[0].getGpos(sense, start, minMapLen, maxMapLen);
			else
				return edges[edges.length- 1].getGpos(sense, start, minMapLen, maxMapLen);
		}

		// junctions
		int p= -1;
		int interSum= 0;
		for (int i = 1; i < edges.length- 1; i++) 
			interSum+= edges[i].length();
		// int maxOverhang= maxMapLen- interSum- 2, minOverhang= minMapLen- interSum- 2; 
		if (sense) {
			if (start)
				p= Math.max(edges[0].getTail().getSite().getPos(), 
						edges[0].getHead().getSite().getPos()- (maxMapLen- (interSum+ 2)));	// 1 overlap + 1 length->position
			else // end
				p= Math.min(edges[0].getHead().getSite().getPos(),
						edges[0].getHead().getSite().getPos()- (minMapLen- (interSum+ edges[edges.length- 1].length()+ 2))); 
		} else { // asense
			if (start)
				p= Math.max(edges[edges.length- 1].getTail().getSite().getPos(), 
						edges[edges.length- 1].getTail().getSite().getPos()+ minMapLen- (interSum+ edges[0].length()+ 2));
			else
				p= Math.min(edges[edges.length- 1].getHead().getSite().getPos(),
						edges[edges.length- 1].getTail().getSite().getPos()+ maxMapLen- (interSum+ 2));
						
		}
		return p;
	}
	
	
	public Node getTail() {
		return edges[0].getTail();
	}

	
	public Node getHead() {
		return edges[edges.length- 1].getHead();
	}
	
	
	@Override
	/**
	 * @deprecated avoid []s
	 * should not be called anymore, as only atomic edges are
	 * interrogated!
	 * returns int[0]=int[1]+1 in case of nonexisting fragments
	 * 
	 */
	public int[] getFrac(Transcript t, int readLen) {
		if (isPend()) 
			return getPEfrac(t, readLen, edges);
		else 
			return getFrac(t, readLen, edges);
	}
	
	@Override
	public int[] getFrac(Transcript t, int readLen, byte dir) {
		if (isPend()) 
			return getPEfrac(t, readLen, edges, dir);
		else 
			return getFrac(t, readLen, edges, dir);
	}

	/**
	 * @deprecated
	 * @param t a transcript
	 * @param readLen read length
	 * @param edges vector of edges
	 * @return vector of four positions (left, left, right, right)
	 */
	public static int[] getPEfrac(Transcript t, int readLen, AbstractEdge[] edges) {
		assert(edges.length== 2);
		int[] left= edges[0].getFrac(t, readLen);
		int[] right= edges[1].getFrac(t, readLen);
		
		return new int[] {left[0], left[1], right[0], right[1]};
	}

	/**
	 * @deprecated
	 * @param t a transcript
	 * @param readLen read length
	 * @param edges vector of edges
	 * @return vector of two positions (first, last)
	 */
	public static int[] getFrac(Transcript t, int readLen, AbstractEdge[] edges) {
			int[] res= null;
			int start= t.getExonicPosition(getFirstEJ(edges));
			int end= t.getExonicPosition(getLastEJ(edges));
			// TODO changed from (start-rDelta)
			// int eDelta= end- start+ 1;	
			// int rDelta= readLen- eDelta;
			// start= Math.max(0, start-rDelta);
	
			int genFirst= edges[0].getTail().getSite().getPos();
			int firstPos= t.getExonicPosition(genFirst);
			if (edges[0].getTail().getSite().isRightFlank())
				++firstPos;
			int first= Math.max(firstPos, end-(readLen-1));
			
			SpliceSite endSite= edges[edges.length-1].getHead().getSite();
			int genLast= endSite.getPos();
			int lastPos= t.getExonicPosition(genLast);
			if (endSite.isLeftFlank())
				--lastPos;
			int last= Math.min(start, lastPos-(readLen-1)); 
	
			//end+= rDelta;
	//		if (first> last)
	//			System.currentTimeMillis();
			
			// TODO the method returns the theoretical coordinates
			// check downstream whether there is an existing length
			// assert(first<= last);
			res= new int[] {first,last};	// TODO start,end
	
			return res;
	
		}
}
