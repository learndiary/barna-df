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
 * invested a lot in new graph contraction (now 2 edges between two vertexes can exist)
 * mergePartitions not correct -- cannot merge partitions has to save n-combinations realized
 * time benchmark ok -> 15 sec RefSeq, <5min mRNA
 */

package fbi.genome.model.splicegraph;

import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Set;
import java.util.Vector;

import fbi.genome.model.SpliceSite;

public class Node {
	
	public static class PositionTypeComparator extends SpliceSite.PositionTypeComparator {
		public int compare(Object o1, Object o2) {
			return super.compare(((Node) o1).getSite(),(((Node) o2).getSite()));			
		}
	}
	static PositionTypeComparator defaultPositionTypeComparator= new PositionTypeComparator();
	
	SpliceSite site;
	long[] transcripts= null;
	Vector<Edge> outEdges= new Vector<Edge>(2);	// 2,1f
	Vector<Edge> inEdges= new Vector<Edge>(2);	// 2,1f
	boolean processed= false;	// coloring for contracting graph
	HashMap<Node, Vector<Path>> fromNodeMap= new HashMap<Node, Vector<Path>>();
	HashMap<Edge, Vector<Edge>> outPartitionMap= null;
	Vector<Vector<long[]>> forbiddenCombinations= null;
	int outPartitionSize= 0;
	Vector<Partition> partitions= new Vector();
	
	private HashMap<Vector<Edge>, Vector<Path>> mapVecEdgeVecPath;
	private Vector<Vector<Path>> vecVecPath;
	private HashMap<Edge, Edge> mapEdgeEdge;
	private Iterator<Edge> iterEdge;
	private HashMap<Vector<Edge>, Vector<Edge>> mapVecEdgeVecEdge;
	private Iterator<Vector<Edge>> iterVecEdge;
	//private mapVecEdgeVecPath;
	
	public Node(SpliceSite newSite, long[] newTranscripts) {
		this.site= newSite;
		this.transcripts= newTranscripts;
	}
	
	public void splitPathes(int k, Vector<Path> pathes) {
			// collect buckets
//		HashMap<Edge, Vector<Path>> bucketHash= new HashMap<Edge, Vector<Path>>(pathes.size(),1f);
//		for (int i = 0; i < pathes.size(); i++) {
//			Vector<Path> v= bucketHash.get(pathes.elementAt(i).getSourceEdge());			
//			if (v== null)
//				v= bucketHash.get(pathes.elementAt(i).getSinkEdge());	// necessary
//			// if both exist, bucket for outEdge and inEdge, then they are the same bucket (hopefully)
//			if (v== null)
//				v= new Vector<Path>(pathes.size()/ 2);
//			v.add(pathes.elementAt(i));
//			bucketHash.put(pathes.elementAt(i).getSourceEdge(), v);
//			bucketHash.put(pathes.elementAt(i).getSinkEdge(), v);
//		}
//		
//		if ()
	}
	
	
	public HashMap<Edge, Vector<Edge>> getOutPartitionMap() {
		if (outPartitionMap == null) {
			outPartitionMap = new HashMap<Edge, Vector<Edge>>(outEdges.size());
			Iterator<Edge> iter= getOutEdges().iterator();
			while (iter.hasNext()) {
				Edge e= iter.next();
				Vector<Edge> v= new Vector<Edge>(1);
				v.add(e);
				outPartitionMap.put(e,v);
			}
			outPartitionSize= getOutEdges().size();
		}

		return outPartitionMap;
	}
	
	public int getOutPartitionSize() {
		if (outPartitionSize== 0)
			getOutPartitionMap();
		return outPartitionSize;
	}
	
	@Override
	public int hashCode() {
		return getSite().getPos();
	}
	
	@Override
	public boolean equals(Object obj) {		
		return getSite().equals(((Node) obj).getSite());
	}
	
	public void addOutEdge(Edge e) {
		outEdges.add(e);
	}
	
	public void addInEdge(Edge e) {
		inEdges.add(e);
	}

	public Vector<Edge> getInEdges() {
		return inEdges;
	}

	public Vector<Edge> getOutEdges() {
		return outEdges;
	}
	
	public Vector<long[]> getOutPartitions() {
		return getPartitions(outEdges);
	}
	
	public Vector<long[]> getInPartitions() {
		return getPartitions(inEdges);
	}
	
	Vector<long[]> getPartitions(Collection<Edge> e) {
		Vector<long[]> v= new Vector<long[]>(e.size());
		Iterator<Edge> edgeIter= e.iterator();
		while (edgeIter.hasNext())
			v.add(edgeIter.next().getTranscripts());
		return v;
	}
	
	public void removeInEdge(Edge e) {
		for (int i = 0; i < inEdges.size(); i++) {
			if (inEdges.elementAt(i)== e) {
				inEdges.remove(i);	// do not use equals()
				return;
			}
		}
	}
	
	public void removeOutEdge(Edge e) {
		for (int i = 0; i < outEdges.size(); i++) {
			if (outEdges.elementAt(i)== e) {
				outEdges.remove(i);	// do not use equals()
				return;
			}
		}
	}
	
	public long[] getTranscripts() {
		return transcripts;
	}

	public SpliceSite getSite() {
		return site;
	}
	
	public String toString() {
		return getSite().toString();
	}
	
	public boolean isProcessed() {
		return processed;
	}

	public void setProcessed(boolean processed) {
		this.processed = processed;
	}

	public Set<Node> getFromNodes() {
		return fromNodeMap.keySet();
	}
	
	public HashMap<Node, Vector<Path>> getFromNodeMap() {
		return fromNodeMap;
	}

	public boolean checkTuple(Path[] tuple) {
		if (forbiddenCombinations== null|| forbiddenCombinations.size()== 0)
			return true;
		for (int i = 0; i < forbiddenCombinations.size(); i++) {
			int j= 0; 
			for (; j < tuple.length; j++) {
				int k= 0;
				for (; k < forbiddenCombinations.elementAt(j).size(); k++) {
					if (Graph.isNull(Graph.intersect(tuple[j].getTranscripts(), forbiddenCombinations.elementAt(i).elementAt(k))))
						break;
				}
				if (k< forbiddenCombinations.elementAt(j).size())
					break;	// at least one is not in the forbidden partition
			}
			if (j== tuple.length)
				return false;
		}
		return true;	// combination valid
	}
	
	void addForbiddenCombination(Vector<long[]> partition) {
		if (forbiddenCombinations== null)
			forbiddenCombinations= new Vector(1,1);
		forbiddenCombinations.add(partition);
	}
	
	public void mergePartitions_test(Vector<Edge> splits) {
		
//		Iterator<Edge> iterEdge= splits.iterator();
//		while (iterEdge.hasNext()) {
//			long[] inT= iterEdge.next().getTranscripts();
//			Iterator<Edge> iterEdge2= outPartitionMap.keySet().iterator();
//			while (iterEdge2.hasNext()) {
//				Edge e= iterEdge2.next();
//				long[] outT= e.getTranscripts();
//				if (Graph.isNull(Graph.intersect(inT, outT)))
//					continue;
//				Vector<long[]> partPathes= outPartitionMap.get(e);
//				for (int i = 0; i < partPathes.size(); i++) {
//					long[] inter= Graph.intersect(partPathes.elementAt(i), inT);
//					if (!Graph.isNull(inter)) {
//						long[] without= Graph.without(partPathes.elementAt(i), inT);
//						if (!Graph.isNull(without)) {
//							// update partPathes and links
//							
//						}
//					}
//						
//				}
//			}
//		}
//		
//		
//		mapEdgeEdge= new HashMap<Edge, Edge>(splitPathes.size(),1f);
//		for (int i = 0; i < splitPathes.size(); i++) {
//				mapEdgeEdge.put(splitPathes.elementAt(i).elementAt(0).getSourceEdge(),
//						splitPathes.elementAt(i).elementAt(0).getSourceEdge());
//		}
//		
//		iterEdge= mapEdgeEdge.keySet().iterator();
//		Vector<Edge> v= new Vector<Edge>(mapEdgeEdge.values());
//		while (iterEdge.hasNext()) {
//			Edge e= iterEdge.next();
//			outPartitionMap.remove(e);
//			outPartitionMap.put(e, v);
//		}
//		
//		mapVecEdgeVecEdge= new HashMap<Vector<Edge>, Vector<Edge>>(outPartitionMap.size(),1f);
//		iterVecEdge= outPartitionMap.values().iterator();
//		while (iterVecEdge.hasNext())
//			mapVecEdgeVecEdge.put(iterVecEdge.next(), null);
//		outPartitionSize= mapVecEdgeVecEdge.size();
	}

	public void setFromNodeMap(HashMap<Node, Vector<Path>> fromNodeMap) {
		this.fromNodeMap = fromNodeMap;
	}

	public static PositionTypeComparator getDefaultPositionTypeComparator() {
		return defaultPositionTypeComparator;
	}
}
