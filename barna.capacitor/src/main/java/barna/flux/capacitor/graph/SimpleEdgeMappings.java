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

package barna.flux.capacitor.graph;

import barna.model.Transcript;
import barna.model.constants.Constants;
import barna.model.splicegraph.Node;
import barna.model.splicegraph.SimpleEdge;
import barna.model.splicegraph.SplicingGraph;
import barna.model.splicegraph.SuperEdge;

import java.io.PrintStream;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Vector;

public class SimpleEdgeMappings extends SimpleEdge implements MappingsInterface {

	Mappings mappings= null;
	
	@Override
	public Mappings getMappings() {
		return mappings;
	}

	
	public SimpleEdgeMappings(Node newTail, Node newHead) {
		super(newTail, newHead);
		mappings= new Mappings();
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

	public double getBpCoverage(Transcript tx, byte dir, int mapLenMax) {
		if (dir== Constants.DIR_BOTH) {
			double sense= getBpCoverage(tx, Constants.DIR_FORWARD, mapLenMax);
			double asense= getBpCoverage(tx, Constants.DIR_BACKWARD, mapLenMax);
			return (sense+ asense);
		}
		
		assert(dir== Constants.DIR_FORWARD^ dir== Constants.DIR_BACKWARD);
		double cov= dir== Constants.DIR_FORWARD? mappings.getReadNr(): mappings.getRevReadNr();
		cov/= getEffLength(tx, dir, mapLenMax);
		return cov;
	}

	public float[] getCoverage(boolean sense, int minMapLen, int maxMapLen) {
		
		float[] a= sense? mappings.coverage: mappings.coverageRev;
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
				mappings.coverage= a;
			else
				mappings.coverageRev= a;
		}
	
		return a;
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
		
		int start= getTail().getSite().getPos();
		int end= getHead().getSite().getPos();
		if(p< start|| p> end)
			return -1;
		
		int sStart= getGpos(sense, true, minMapLen, maxMapLen),
			sEnd= getGpos(sense, false, minMapLen, maxMapLen);
		if (p< sStart|| p> sEnd)
			return -1;
		float[] a= sense? mappings.coverage: mappings.coverageRev;
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
				if ((sense&& se.getEdges()[0]== this)|| (se.getEdges()[se.getEdges().length- 1]== this&& !sense)) {
					float c= ((SuperEdgeMappings) se).getCoverage(i+ start, sense, minMapLen, maxMapLen);	// sense
					if (c>= 0) {	// (-1) for not allowed pos
						res[i]+= c;
					}  
				} 
			}
			
		}
		
		return res;
	}

	@Override
	public int getMapLength(Transcript tx, int maxMapLength) {
		// TODO Auto-generated method stub
		return 0;
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
		double cov= dir== Constants.DIR_FORWARD? mappings.getReadNr(): mappings.getRevReadNr();
		int mapLen= getMapLength(tx, mapLenMax);
		cov*= mapLen;
		cov/= getNtLength(tx, mapLen);
		return cov;
	}

	public void incrReadNr(boolean sense, int minMapLen, int maxMapLen, int genomicPos, float w) {
		float[] a= getCoverage(sense, minMapLen, maxMapLen);
		int start= getGpos(sense, true, minMapLen, maxMapLen); // getGpos(sense, true, minMapLen);
		int rPos= genomicPos- start;
		if (rPos< 0|| rPos>= a.length)
			return;	// reads that overhang start/end
		++a[rPos];
		mappings.incrReadNr();
	}

	public void printCoverage(PrintStream p, SplicingGraph g, int minMapLen, int maxMapLen) {
			
			// assume exonic
			int start= getTail().getSite().getPos(), end= getHead().getSite().getPos();
			int sSense= Math.max(getGpos(true, true, minMapLen, maxMapLen),start),
				eSense= Math.min(getGpos(true, false, minMapLen, maxMapLen),end),
				sASense= Math.max(getGpos(false, true, minMapLen, maxMapLen), start),
				eASense= Math.min(getGpos(false, false, minMapLen, maxMapLen), end);
			Transcript[] tt= g.decodeTset(getTranscripts());
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
					if (se.getEdges()[0]== this) {
						float c= ((SuperEdgeMappings) se).getCoverage(i, true, minMapLen, maxMapLen);	// sense
						if (c>= 0) {	// (-1) for not allowed pos
							String id= mapSEID.get(se);
							if (id== null) {
								Transcript[] ttt= g.decodeTset(se.getTranscripts());
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
					if (se.getEdges()[se.getEdges().length- 1]== this) {
						float c= ((SuperEdgeMappings) se).getCoverage(i, false, minMapLen, maxMapLen);	// asense
						if (c>= 0) {	// (-1) for not allowed pos
							String id= mapSEID.get(se);
							if (id== null) {
								Transcript[] ttt= g.decodeTset(se.getTranscripts());
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
						if (se.getEdges()[0]== this) {
							float c= ((SuperEdgeMappings) se).getCoverage(i, true, minMapLen, maxMapLen);	// sense
							if (c>= 0) {	// (-1) for not allowed pos
								String id= mapSEID.get(se);
								if (id== null) {
									Transcript[] ttt= g.decodeTset(se.getTranscripts());
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
						if (se.getEdges()[se.getEdges().length- 1]== this) {
							float c= ((SuperEdgeMappings) se).getCoverage(i, false, minMapLen, maxMapLen);	// asense
							if (c>= 0) {	// (-1) for not allowed pos
								String id= mapSEID.get(se);
								if (id== null) {
									Transcript[] ttt= g.decodeTset(se.getTranscripts());
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

}
