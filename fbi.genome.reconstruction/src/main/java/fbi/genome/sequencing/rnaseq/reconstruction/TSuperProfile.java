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

import java.util.Vector;

import fbi.genome.model.Transcript;
import fbi.genome.model.constants.Constants;
import fbi.genome.model.splicegraph.AbstractEdge;
import fbi.genome.model.splicegraph.SplicingGraph;
import fbi.genome.model.splicegraph.SuperEdge;
import fbi.genome.sequencing.rnaseq.graph.AnnotationMapper;

public class TSuperProfile {

	Vector<TProfile> v= new Vector<TProfile>(2,2);
	public TSuperProfile() {
		
	}

public void addProfile(TProfile pro) {
		v.add(pro);
	}
	
	public Vector<TProfile> getProfiles() {
		return v;
	}
	
	public static int[] rel2absolute(double[] relPos, int len) {
		int[] coords= new int[relPos.length];
		for (int j = 0; j < coords.length; j++) {
			// 090820 floor, floor <-> 
			// exclusive regions for back-normalization needed
			// otherwise fracs> transcriptcount for gene (and reads also)
//			coords[j]= (j%2== 0)? (int) Math.floor(relPos[j]* len)
//							   : (int) Math.floor(relPos[j]* len);
			coords[j]= (int) (relPos[j]* len); // too slow: Math.floor()
			if (coords[j]> len)
				coords[j]= len;
		}
		return coords;
	}
	
	int[] countedReads= null;
	public int getReads(AnnotationMapper g, Transcript t, int x, int readLen, int[] insertMinMax) {
		
		if (countedReads == null) {
			countedReads = new int[v.size()];
			
			Vector<Vector<AbstractEdge>> vv= new Vector<Vector<AbstractEdge>>(1);
			vv.add(new Vector<AbstractEdge>());
			g.getRPK(t, false, SplicingGraph.ETYPE_AL, vv);
			long[] sig= g.encodeTset(t);
			
			for (int i = 0; i < countedReads.length; i++) 
				countedReads[i]= 0;
			
			int elen= t.getExonicLength();
			if (elen>= readLen) {
				for (int j = 0; j < vv.elementAt(0).size(); j++) {
					AbstractEdge e= vv.elementAt(0).elementAt(j);
					long[] inter= SplicingGraph.intersect(sig, e.getTranscripts());
					if (SplicingGraph.isNull(inter))
						continue;
					
					if (insertMinMax!= null) {
						for (int k = 0; e.getSuperEdges()!= null&& k < e.getSuperEdges().size(); k++) {
							SuperEdge se= e.getSuperEdges().elementAt(k);
							inter= SplicingGraph.intersect(se.getTranscripts(), sig);
							if (!se.isPend()|| se.getEdges()[0]!= e|| SplicingGraph.isNull(inter))
								continue;
							int[] bounds= se.getFrac(t, readLen);
							double[] relBounds= GraphLPsolver.bounds2rel(bounds, elen- readLen);
							for (int i = 0; i < countedReads.length; i++) {
								bounds= rel2absolute(relBounds, v.elementAt(i).length());
								countedReads[i]+= v.elementAt(i).getArea(bounds, readLen, insertMinMax, FluxCapacitorConstants.BYTE_0);
							}
						}
					} else {
						int[] bounds= e.getFrac(t, readLen);
						double[] relBounds= GraphLPsolver.bounds2rel(bounds, elen- readLen);
						for (int i = 0; i < countedReads.length; i++) {
							bounds= rel2absolute(relBounds, v.elementAt(i).length());
							countedReads[i]+= v.elementAt(i).getArea(bounds, readLen, insertMinMax, FluxCapacitorConstants.BYTE_0);
						}
					}
				}
			}
			
			// TODO kill test
//			if (t.getTranscriptID().contains("Gnf1-RA"))
//				debug= true;
			double test= 0d;
			if (elen>= readLen) {
				for (int j = 0; j < vv.elementAt(0).size(); j++) {
					AbstractEdge e= vv.elementAt(0).elementAt(j);
					long[] inter= SplicingGraph.intersect(sig, e.getTranscripts());
					if (SplicingGraph.isNull(inter))
						continue;
					
					if (insertMinMax!= null) {
						for (int k = 0; e.getSuperEdges()!= null&& k < e.getSuperEdges().size(); k++) {
							SuperEdge se= e.getSuperEdges().elementAt(k);
							inter= SplicingGraph.intersect(se.getTranscripts(), sig);
							if (!se.isPend()|| se.getEdges()[0]!= e|| SplicingGraph.isNull(inter))
								continue;
							int[] bounds= se.getFrac(t, readLen);
							double[] relBounds= GraphLPsolver.bounds2rel(bounds, elen- readLen);
							double frac= getAreaFrac(g, t, relBounds, readLen, insertMinMax, FluxCapacitorConstants.BYTE_0);
							test+= frac;
		//					if (debug) {
		//						System.out.print(frac+"\t"+se);
		//						for (int i = 0; i < se.getEdges().length; i++) 
		//							System.out.print("\t"+ se.getEdges()[i]);
		//						System.out.println();
		//					}
						}
					} else {
						int[] bounds= e.getFrac(t, readLen);
						double[] relBounds= GraphLPsolver.bounds2rel(bounds, elen- readLen);
						double frac= getAreaFrac(g, t, relBounds, readLen, insertMinMax, FluxCapacitorConstants.BYTE_0);
						test+= frac;
					}
				}
			}
			
			if (insertMinMax!= null)
				try {assert(elen< (2* readLen)+ insertMinMax[0]|| (test!= 0&& Math.abs(1d- test)< 0.1));}catch(AssertionError e){
					if (Constants.verboseLevel> Constants.VERBOSE_NORMAL)
						System.err.println("[ASSERTION] in "+ getClass().getName()+".getReads():\n\ttest= "+test+", readLen= "+readLen+", insertMin "+insertMinMax[0]);
				}
			else
				try {assert(elen<= readLen|| (test!= 0&& Math.abs(1d- test)< 0.1));}catch(AssertionError e){
					if (Constants.verboseLevel> Constants.VERBOSE_NORMAL)
						System.err.println("[ASSERTION] in "+ getClass().getName()+".getReads():\n\ttest= "+test+", readLen= "+readLen+", elen "+elen);
				}
		}

		return countedReads[x];
	}
	
	public int getAllReads() {
		int sum= 0;
		for (int i = 0; i < v.size(); i++) 
			sum+= v.elementAt(i).getReads();
		return sum;
	}
	
	public int getReads(AnnotationMapper g, Transcript t, int readLen, int[] insertMinMax) {
		int sum= 0;
		for (int i = 0; i < v.size(); i++) 
			sum+= getReads(g, t, i, readLen, insertMinMax);
		return sum;
	}
	public double getArea(AnnotationMapper g, Transcript t, double[] relPos, int readLen, int[] insertMinMax, byte dir) {
		return getAreaFrac(g, t, relPos, readLen, insertMinMax, dir, false);
	}
	public double getAreaFrac(AnnotationMapper g, Transcript t, double[] relPos, int readLen, int[] insertMinMax, byte dir) {
		return getAreaFrac(g, t, relPos, readLen, insertMinMax, dir, true);
	}
	double getAreaFrac(AnnotationMapper g, Transcript t, double[] relPos, int readLen, int[] insertMinMax, byte dir, boolean divide) {
		
		int[] coords= new int[relPos.length];
		double sum= 0, total= 0;
		for (int i = 0; i < v.size(); i++) {
			coords= rel2absolute(relPos, v.elementAt(i).length());
			sum+= v.elementAt(i).getArea(coords, readLen, insertMinMax, dir);
			total+= getReads(g, t, i, readLen, insertMinMax);
//			if (insertMinMax!= null&& relPos.length== 2) // hack for pend<>single discrepancy
//				total+= getReads(i, readLen, insertMinMax);
		}
	
		
		if (sum== 0|| (!divide))
			return sum;
		// else
		return (sum/ total);
	}
	
	public int project(int[][] b) {
		if (b== null|| b.length== 0)
			return -1;
		int sum= 0;
		for (int i = 0; i < v.size(); i++) {
			sum+= v.elementAt(i).m.project(b);
		}
		return sum;
	}
	
}
