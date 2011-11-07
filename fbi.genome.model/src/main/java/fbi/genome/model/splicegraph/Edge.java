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

/**
 * test
 */
package fbi.genome.model.splicegraph;

import java.io.PrintStream;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Vector;

import fbi.genome.model.Transcript;
import fbi.genome.model.constants.Constants;

public class Edge {
	
	public static class PositionComparator implements Comparator<Edge> {
		public int compare(Edge arg0, Edge arg1) {
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
	
	static PositionComparator defaultPositionComparator= new PositionComparator();
	
	byte type= Transcript.ID_SRC_UNDEFINED;
	boolean contracted= false, processed= false, valid= true, exonic= false;	
	
	public static String getStringRep(Node v, Node w) {
		return v.getSite().toString()+w.getSite().toString();
	}
	
	Node tail, head;
	String stringRep;
	long[] transcripts;
	Partition partition;
	int readNr= 0, revReadNr= 0;
	Vector<SuperEdge> superEdges= null;
	float[] coverage= null, coverageRev= null;	// array over all the edge (NOT: -readlength)
	
	public void incrReadNr() {
		++readNr;
	}
	
	public boolean containsPos(boolean sense, int minMapLen, int maxMapLen, int genomicPos) {
		float[] a= getCoverage(sense, minMapLen, maxMapLen);
		if (a== null)
			return false;
		int start= getGpos(sense, true, minMapLen, maxMapLen); // getGpos(sense, true, minMapLen);
		int rPos= genomicPos- start;
		if (rPos< 0|| rPos>= a.length)
			return false;	// reads that overhang start/end
		return true;
	}
	
	public void incrReadNr(boolean sense, int minMapLen, int maxMapLen, int genomicPos, float w) {
		float[] a= getCoverage(sense, minMapLen, maxMapLen);
		int start= getGpos(sense, true, minMapLen, maxMapLen); // getGpos(sense, true, minMapLen);
		int rPos= genomicPos- start;
		if (rPos< 0|| rPos>= a.length)
			return;	// reads that overhang start/end
		++a[rPos];
		incrReadNr();
	}
	
	public float[] getCoverage(boolean sense, int minMapLen, int maxMapLen) {
		
		float[] a= sense? coverage: coverageRev;
		if (a == null) {
			int x= getGpos(sense, true, minMapLen, maxMapLen), // getTail().getSite().getPos(),
				y= getGpos(sense, false, minMapLen, maxMapLen);// getHead().getSite().getPos(); 	
			if (x< getTail().getSite().getPos()|| y> getHead().getSite().getPos()|| x> y) {
				getGpos(sense, true, minMapLen, maxMapLen);
				getGpos(sense, false, minMapLen, maxMapLen);
				return null;
			}
			int len= y- x+ 1;
			a = new float[len];
			Arrays.fill(a, 0f);
			if (sense)
				coverage= a;
			else
				coverageRev= a;
		}

		return a;
	}
	
	
	public void decrReadNr() {
		--readNr;
	}
	public void decrRevReadNr() {
		--revReadNr;
	}
		
	public Edge(Node newTail, Node newHead) {
		this.tail= newTail;
		this.head= newHead;
		tail.addOutEdge(this);
		head.addInEdge(this);
	}
	
	protected Edge() {
		
	}

	public Node getHead() {
		return head;
	}

	public Node getTail() {
		return tail;
	}
	
	@Override
	public int hashCode() {		
		return toString().hashCode();
	}
	
	public String toString() {
		if (stringRep == null) 
			stringRep = getStringRep(getTail(), getHead());

		return stringRep;
	}
	
	/**
	 * 
	 * @param p genomic position (sign stranded)
	 * @param sense sense or antisense
	 * @param minMapLen minimum mapping length
	 * @param maxMapLen maximum mapping length
	 * @return
	 */
	public float getCoverage(int p, boolean sense, int minMapLen, int maxMapLen) {
		
		assert(isExonic());
		
		int start= tail.getSite().getPos();
		int end= head.getSite().getPos();
		if(p< start|| p> end)
			return -1;
		
		int sStart= getGpos(sense, true, minMapLen, maxMapLen),
			sEnd= getGpos(sense, false, minMapLen, maxMapLen);
		if (p< sStart|| p> sEnd)
			return -1;
		float[] a= sense? coverage: coverageRev;
		if (a== null)
			return 0;
		return a[p- sStart]; 
		
	}
	
	public float[] getCoverageFlat(boolean sense, int minMapLen, int maxMapLen) {
		int start= getTail().getSite().getPos(), end= getHead().getSite().getPos(); 
		float[] res= new float[end - start+ 1];
		Arrays.fill(res, 0);
		for (int i = 0; i < res.length; i++) {
			float c= getCoverage(i+ start, sense, minMapLen, maxMapLen);
			if (c>= 0)
				res[i]+= c;
		}
		
		Vector<SuperEdge> v= getSuperEdges();
		for (int i = 0; i < res.length; i++) {			
			for (int j = 0; v!= null&& j < v.size(); j++) {
				SuperEdge se= v.elementAt(j);				
				if ((sense&& se.edges[0]== this)|| (se.edges[se.edges.length- 1]== this&& !sense)) {
					float c= se.getCoverage(i+ start, sense, minMapLen, maxMapLen);	// sense
					if (c>= 0) {	// (-1) for not allowed pos
						res[i]+= c;
					}  
				} 
			}
			
		}
		
		return res;
	}
	
	public void printCoverage(PrintStream p, SpliceGraph g, int minMapLen, int maxMapLen) {
		
		// assume exonic
		int start= getTail().getSite().getPos(), end= getHead().getSite().getPos();
		int sSense= Math.max(getGpos(true, true, minMapLen, maxMapLen),start),
			eSense= Math.min(getGpos(true, false, minMapLen, maxMapLen),end),
			sASense= Math.max(getGpos(false, true, minMapLen, maxMapLen), start),
			eASense= Math.min(getGpos(false, false, minMapLen, maxMapLen), end);
		Transcript[] tt= g.decodeTset(transcripts);
		StringBuilder txString= new StringBuilder(tt.length* 10);
		for (int i = 0; i < tt.length; i++) 
			txString.append(tt[i].getTranscriptID()+",");
		txString.deleteCharAt(txString.length()- 1);
		
		float[] a= getCoverageFlat(true, minMapLen, maxMapLen), 
			b= getCoverageFlat(false, minMapLen, maxMapLen);
		Vector<SuperEdge> v= getSuperEdges();
		HashMap<SuperEdge, String> mapSEID= new HashMap<SuperEdge, String>();
		
		if (a.length!= b.length)
			System.currentTimeMillis();
		// 
/*		for (int i = sSense; i <= eASense; ++i) {
			
			p.print(Integer.toString(i)+ "\t");	// genomic position

			if (i<= eSense) {
				p.print(a[i- sSense]+ ",");
			} else
				p.print("NA,");
			if (i>= sASense)
				p.print(b[i- sASense]+",");	// flat all
			else
				p.print("NA,");
			p.print(txString);
			
			if (1==1) {
				p.println();
				continue;
			}
*/
		for (int i= 0; i< b.length; ++i) {
			p.print(Float.toString(a[i])+",");		// cov fwd
			p.print(Float.toString(b[i])+",");	// cov rev
			p.print(txString);								// tid
			if (1==1) {
				p.println();
				continue;
			}
			
			/* splice junctions */
			for (int j = 0; j < v.size(); ++j) {
				SuperEdge se= v.elementAt(j);
				if (se.edges[0]== this) {
					float c= se.getCoverage(i, true, minMapLen, maxMapLen);	// sense
					if (c>= 0) {	// (-1) for not allowed pos
						String id= mapSEID.get(se);
						if (id== null) {
							Transcript[] ttt= g.decodeTset(se.transcripts);
							StringBuilder sb= new StringBuilder(tt.length* 10);
							for (int k = 0; k < ttt.length; ++k) 
								sb.append(ttt[k].getTranscriptID()+",");
							sb.deleteCharAt(sb.length()- 1);
							id= sb.toString();
							mapSEID.put(se, id);
						}
						p.print("\t"+ c+ ","+ "NA"+ id);
					}  
				} else // mutex by definition of a sedge (!)
				if (se.edges[se.edges.length- 1]== this) {
					float c= se.getCoverage(i, false, minMapLen, maxMapLen);	// asense
					if (c>= 0) {	// (-1) for not allowed pos
						String id= mapSEID.get(se);
						if (id== null) {
							Transcript[] ttt= g.decodeTset(se.transcripts);
							StringBuilder sb= new StringBuilder(tt.length* 10);
							for (int k = 0; k < ttt.length; ++k) 
								sb.append(ttt[k].getTranscriptID()+",");
							sb.deleteCharAt(sb.length()- 1);
							id= sb.toString();
							mapSEID.put(se, id);
						}
						p.print("\tNA,"+ c+ ","+ id);
					}
				}
				
				// pend?
				Vector<SuperEdge> w= se.getSuperEdges();
				for (int m = 0; w!= null&& m < w.size(); ++m) {
					se= w.elementAt(m);
					if (se.edges[0]== this) {
						float c= se.getCoverage(i, true, minMapLen, maxMapLen);	// sense
						if (c>= 0) {	// (-1) for not allowed pos
							String id= mapSEID.get(se);
							if (id== null) {
								Transcript[] ttt= g.decodeTset(se.transcripts);
								StringBuilder sb= new StringBuilder(tt.length* 10);
								for (int k = 0; k < ttt.length; ++k) 
									sb.append(ttt[k].getTranscriptID()+",");
								sb.deleteCharAt(sb.length()- 1);
								id= sb.toString();
								mapSEID.put(se, id);
							}
							p.print("\t"+ c+ ","+ "NA"+ id);
						}  
					} else // mutex by definition of a sedge (!)
					if (se.edges[se.edges.length- 1]== this) {
						float c= se.getCoverage(i, false, minMapLen, maxMapLen);	// asense
						if (c>= 0) {	// (-1) for not allowed pos
							String id= mapSEID.get(se);
							if (id== null) {
								Transcript[] ttt= g.decodeTset(se.transcripts);
								StringBuilder sb= new StringBuilder(tt.length* 10);
								for (int k = 0; k < ttt.length; ++k) 
									sb.append(ttt[k].getTranscriptID()+",");
								sb.deleteCharAt(sb.length()- 1);
								id= sb.toString();
								mapSEID.put(se, id);
							}
							p.print("\tNA,"+ c+ ","+ id);
						}
					}
				}
			}
		}
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
		Edge e= (Edge) obj;
		if (getTail().equals(e.getTail())&& getHead().equals(e.getHead())
				&& SpliceGraph.equalSet(getTranscripts(), e.getTranscripts())
				&& isExonic()== isExonic()		// multiple edges exonic, intronic for eg intron retention
				&& isIntronic()== isIntronic())	
			return true;
		return false;
	}

	public long[] getTranscripts() {
		return transcripts;
	}

	public void setTranscripts(long[] transcripts) {
		this.transcripts = transcripts;
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

	public int getReadNr() {
		return readNr;
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

	public Vector<SuperEdge> getSuperEdges() {
		return superEdges;
	}

	public void addSuperEdge(SuperEdge superEdge) {
		if (superEdges== null)
			superEdges= new Vector<SuperEdge>(1,1);
		for (int i = 0; i < superEdges.size(); i++) {
			if (superEdges.elementAt(i)== superEdge)
				return;
		}
		superEdges.add(superEdge);
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
	public int getFrac(boolean tail) {
		if (tail) {
			int pos= getTail().getSite().getPos();
			if (getTail().getSite().isRightFlank())
				++pos;
			return pos;
		} else {
			int pos= getHead().getSite().getPos();
			if (getHead().getSite().isLeftFlank())
				--pos;
			return pos;
		}
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

	public void incrRevReadNr() {
		++revReadNr;
	}

	public int getRevReadNr() {
		return revReadNr;
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
	 * Absolute length in Nt, wrapping length() for exons.  
	 */
	public int getNtLength(Transcript tx, int mapLength) {
		return length();
	}
	
	
	/**
	 * Coverage of nucleotide positions by mapping nucleotides.
	 * @param dir
	 * @param mapLenMax
	 * @return
	 */
	public double getNtCoverage(Transcript tx, byte dir, int mapLenMax) {
		
		if (dir== Constants.DIR_BOTH) {
			double sense= getNtCoverage(tx, Constants.DIR_FORWARD, mapLenMax);
			double asense= getNtCoverage(tx, Constants.DIR_BACKWARD, mapLenMax);
			return (sense+ asense);
		}
		
		assert(dir== Constants.DIR_FORWARD^ dir== Constants.DIR_BACKWARD);
		double cov= dir== Constants.DIR_FORWARD? getReadNr(): getRevReadNr();
		int mapLen= getMapLength(tx, mapLenMax);
		cov*= mapLen;
		cov/= getNtLength(tx, mapLen);
		return cov;
	}

	public double getBpCoverage(Transcript tx, byte dir, int mapLenMax) {
		if (dir== Constants.DIR_BOTH) {
			double sense= getBpCoverage(tx, Constants.DIR_FORWARD, mapLenMax);
			double asense= getBpCoverage(tx, Constants.DIR_BACKWARD, mapLenMax);
			return (sense+ asense);
		}
		
		assert(dir== Constants.DIR_FORWARD^ dir== Constants.DIR_BACKWARD);
		double cov= dir== Constants.DIR_FORWARD? getReadNr(): getRevReadNr();
		cov/= getEffLength(tx, dir, mapLenMax);
		return cov;
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
