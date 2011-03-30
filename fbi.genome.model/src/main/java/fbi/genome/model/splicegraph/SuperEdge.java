package fbi.genome.model.splicegraph;

import fbi.genome.model.SpliceSite;
import fbi.genome.model.Transcript;
import fbi.genome.model.constants.Constants;
//import genome.lp.Function;

import java.util.Comparator;

public class SuperEdge extends Edge {
	public static class EdgeByTailComparator implements Comparator<Edge> {
		public int compare(Edge arg0, Edge arg1) {
			if (arg0.getTail().getSite().getPos()< arg1.getTail().getSite().getPos())
				return -1;
			if (arg0.getTail().getSite().getPos()> arg1.getTail().getSite().getPos())
				return 1;
			return 0;
		}
	}
	static EdgeByTailComparator defaultEdgeByTailComparator= new EdgeByTailComparator();
	static class EdgeTuple {
		public Edge[] edges;
		public EdgeTuple(Edge[] edges) {
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
	
	static EdgeTuple first= null;
	Edge[] edges= null;
	Edge[] phantomEdges= null;
	boolean pend= false;
		
	
	/**
	 * @param edges
	 * @param g
	 * @param pend
	 */
	public SuperEdge(Edge[] edges, long[] supp, boolean pend) {

		setEdges(edges);// this(edges);	// also sets transcript set
		setTranscripts(supp);
		setPend(pend);
//		Transcript[] b4= g.decodeTset(transcripts);
		
//		if (!pend) {
////			
////		} else {
//			tail.addOutEdge(this);
//			head.addInEdge(this);
//			
//			Node[] n= g.getNodesInGenomicOrder();
//			for (int i = 0; i < edges.length-1; i++) {
//				int p= Arrays.binarySearch(n, edges[i].getHead(), Node.getDefaultPositionTypeComparator());
//				assert(p>= 0);
//				
//				Iterator<Edge> iter= n[p].getOutEdges().iterator();
//				while (iter.hasNext()) {
//					Edge e= iter.next();
//					if (e.isExonic())
//						continue;
//					if (e.getHead()== edges[i+1].getTail()) {
//						transcripts= Graph.intersect(transcripts, e.getTranscripts());
//						break;
//					}
//				}
//			}
//		}		
//		Transcript[] after= g.decodeTset(transcripts);
		System.currentTimeMillis();
	}

	public Edge[] getEdges() {
		return edges;
	}
	
	public static int getFirstEJ(Edge[] edges) {
		int pos= edges[0].getHead().getSite().getPos();
		if (edges[0].getHead().getSite().isLeftFlank())
			--pos;
		return pos;
	}
	
	public int getFirstEJ() {
		return getFirstEJ(edges);
	}
		
	
	public static int getLastEJ(Edge[] edges) {
		int pos= edges[edges.length-1].getTail().getSite().getPos();
		if (edges[edges.length-1].getTail().getSite().isRightFlank())
			++pos;
		return pos;
	}
	
	public int getLastEJ() {
		return getLastEJ(edges);
	}
	
	public boolean isIntronic() {
		return false;
	}
	
	public boolean isExonic() {
		for (int i = 0; i < edges.length; i++) 
			if(!edges[i].isExonic())
				return false;
		return true;
	}
	
	public void setEdges(Edge[] edges) {
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
	
	public Edge[] getPhantomEdges() {
		if (phantomEdges == null) {
			phantomEdges = new Edge[edges.length];
			for (int i = 0; i < phantomEdges.length; i++) {
				phantomEdges[i]= new Edge();
				phantomEdges[i].tail= edges[i].tail;
				phantomEdges[i].head= edges[i].head;
			}
		}

		return phantomEdges;
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
	public int getMapLength(Transcript tx, int maxMapLength) {
		if (isPend()) {
			System.err.println("what about pends");
			return -1;
			
		} else {
			int start= tx.getExonicPosition(edges[0].getTail().getSite().getPos()),
				first= tx.getExonicPosition(edges[0].getHead().getSite().getPos()),
				last= tx.getExonicPosition(edges[edges.length- 1].getTail().getSite().getPos()),
				end= tx.getExonicPosition(edges[edges.length- 1].getHead().getSite().getPos());
			return Math.min(Math.min(last- start, end- first)+ 1, maxMapLength);
		}
	}
	public int getEffLength(byte dir, int mapLenMax) {

		System.err.println("check effLen");
		assert(dir== Constants.DIR_FORWARD^ dir== Constants.DIR_BACKWARD);
		int effLen= -1;
		if (isPend()) {	
			
			effLen= dir== Constants.DIR_FORWARD? 
					edges[0].getEffLength(null, dir, mapLenMax):
					edges[edges.length- 1].getEffLength(null, dir, mapLenMax);
		
		} else {	// ej / sj
			int interNt= 1;
			for (int i = 1; i < edges.length- 1; i++) 
				interNt+= edges[i].length();
			effLen= Math.min(mapLenMax, 
					dir== Constants.DIR_FORWARD? edges[0].length(): edges[edges.length- 1].length())
					- interNt;
		}
		
		return effLen;
	}
	
	public double getNtCoverage(Transcript tx, byte dir, int mapLenMax) {
		
		if (isPend()) {
			System.err.println("check pend");
			double cov= getReadNr();
			if (dir== Constants.DIR_BOTH) {
				cov*= 2;
				cov/= getEffLength(Constants.DIR_FORWARD, mapLenMax)+ 
						getEffLength(Constants.DIR_BACKWARD, mapLenMax);
			} else if (dir== Constants.DIR_FORWARD)
				cov/= getEffLength(Constants.DIR_FORWARD, mapLenMax);
			else if (dir== Constants.DIR_BACKWARD)
				cov/= getEffLength(Constants.DIR_FORWARD, mapLenMax);
			else
				assert(false);
			return cov;
		} else
			return super.getNtCoverage(tx, dir, mapLenMax);
		
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
		Edge[] e= se.getEdges();
		if (e.length!= edges.length)
			return false;
		for (int i = 0; i < e.length; i++) {
			if (e[i]!= edges[i])
				return false;
		}
		return true;
	}

	/**
	 * finds (exonic) edge containing <code>genomicPos</code>
	 * @param genomicPos
	 * @return
	 */
	private Edge findEdge(int genomicPos) {
		if (genomicPos< edges[0].getTail().getSite().getPos()
				|| genomicPos> edges[edges.length- 1].getHead().getSite().getPos())
			return null;
		int i;
		for (i = 0; i < edges.length; ++i) 
			if (edges[i].getTail().getSite().getPos()> genomicPos)
				break;
		assert(i>= 0&& i<= edges.length);
		if (genomicPos> edges[i- 1].getHead().getSite().getPos())
			return null;
		return edges[i- 1];		
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
		return pend;
	}

	public void setPend(boolean pend) {
		this.pend = pend;
	}
	

	public static int[] getPEfrac(Transcript t, int readLen, Edge[] edges, byte dir) {
		assert(edges.length== 2);
		int[] left= edges[0].getFrac(t, readLen,dir);
		int[] right= edges[1].getFrac(t, readLen,dir);
		
		return new int[] {left[0], left[1], right[0], right[1]};
	}
	
	public static int[] getFrac(Transcript t, int readLen, Edge[] edges, byte dir) {
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
		if (pend) {
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
		int maxOverhang= maxMapLen- interSum- 2, minOverhang= minMapLen- interSum- 2; 
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
	
	
	@Override
	public float[] getCoverage(boolean sense, int minMapLen, int maxMapLen) {
//		if (isPend()) {
//			if (sense)
//				return edges[0].getCoverage(sense, minMapLen, maxMapLen);
//			else
//				return edges[edges.length- 1].getCoverage(sense, minMapLen, maxMapLen);
//		} else
		return super.getCoverage(sense, minMapLen, maxMapLen);
	}
	
	public float getCoverage(int p, boolean sense, int minMapLen, int maxMapLen) {
		
		assert(isExonic());
		
		int start= getGpos(sense, true, minMapLen, maxMapLen);
		int end= getGpos(sense, false, minMapLen, maxMapLen);
		if(p< start|| p> end)
			return -1;
		
		float[] a= getCoverage(sense, minMapLen, maxMapLen);
		if (a== null)
			return -1;
		
		int sStart= getGpos(sense, true, minMapLen, maxMapLen), sEnd= getGpos(sense, false, minMapLen, maxMapLen);
		if (p< sStart|| p> sEnd)
			return -1;
		if ((sense&& coverage== null)|| ((!sense)&& coverageRev== null))
			return 0;
		return sense?coverage[p- sStart]:coverageRev[p- sStart]; 
		
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
	 * @param t
	 * @param readLen
	 * @param edges
	 * @return
	 */
	public static int[] getPEfrac(Transcript t, int readLen, Edge[] edges) {
		assert(edges.length== 2);
		int[] left= edges[0].getFrac(t, readLen);
		int[] right= edges[1].getFrac(t, readLen);
		
		return new int[] {left[0], left[1], right[0], right[1]};
	}

	/**
	 * @deprecated
	 * @param t
	 * @param readLen
	 * @param edges
	 * @return
	 */
	public static int[] getFrac(Transcript t, int readLen, Edge[] edges) {
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
