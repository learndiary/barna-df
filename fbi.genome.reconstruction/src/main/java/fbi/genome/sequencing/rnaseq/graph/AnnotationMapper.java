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

package fbi.genome.sequencing.rnaseq.graph;

import java.util.Arrays;
import java.util.Iterator;
import java.util.Vector;

import fbi.commons.Log;
import fbi.genome.io.BufferedIterator;
import fbi.genome.io.rna.UniversalReadDescriptor;
import fbi.genome.io.rna.UniversalReadDescriptor.Attributes;
import fbi.genome.model.Gene;
import fbi.genome.model.SpliceSite;
import fbi.genome.model.Transcript;
import fbi.genome.model.bed.BEDobject;
import fbi.genome.model.bed.BEDobject2;
import fbi.genome.model.splicegraph.Edge;
import fbi.genome.model.splicegraph.Node;
import fbi.genome.model.splicegraph.SplicingGraph;
import fbi.genome.model.splicegraph.SuperEdge;

public class AnnotationMapper extends SplicingGraph {

	public long nrMappingsLocus= 0;
	public long nrMappingsMapped= 0;
	public long nrMappingsNotMapped= 0;
	public long nrMappingsNotMappedAsPair= 0;
	public long nrMappingsWrongStrand= 0;
	public long nrMappingsWrongPairOrientation= 0;
	public long nrMappingsLocusMultiMaps= 0;
	
	long mappedReads= 0;

	public AnnotationMapper(Gene gene) {
		super(gene);
		constructGraph();
		getNodesInGenomicOrder();	// important ??!
		transformToFragmentGraph();

	}
	
	public Edge getEdge(BEDobject2 obj) {
			
		Vector<Edge> v= edgeVector; //new Vector<Edge>();
		v.removeAllElements();
		
		Node[] nodes= getNodesInGenomicOrder();
		boolean parallel= true;
		if (obj.getStrand()!= gene.getStrand())
			parallel= false;
		int bMax= obj.getBlockCount();
		if (bMax== 0)
			++bMax;
		SpliceSite dummySite= new SpliceSite(0, SpliceSite.TYPE_ACCEPTOR, null);
		Node dummyNode= new Node(dummySite, null);
		int[] starts= new int[bMax], ends= new int[bMax];
		if (bMax== 1) {
			starts[0]= obj.getStart();
			ends[0]= obj.getEnd();
		} else { 
			obj.resetBlocks();
			for (int i = 0; i < ends.length; i++) {
				starts[i]= obj.getNextBlockStart();
				ends[i]= starts[i]+ obj.getNextBlockSize();
				++starts[i];
			}
		}
		for (int i = parallel?0:(bMax- 1); parallel?i < bMax:i>=0; i+=parallel?1:(-1)) {
			
			// convert it to directionality of trpt, 
			// now with reads on antisense strand
			int prime5= (gene.getStrand()>= 0?
					//obj.getAbsBlockStart(i):obj.getAbsBlockEnd(i))* gene.getStrand(),
					starts[i]: ends[i])* gene.getStrand(),
				prime3= (gene.getStrand()>= 0?
					//obj.getAbsBlockEnd(i):obj.getAbsBlockStart(i))* gene.getStrand();
					ends[i]:starts[i])* gene.getStrand();
			
			dummyNode.getSite().setPos(prime5);
			int pStart= Arrays.binarySearch(nodes, dummyNode, Node.getDefaultPositionTypeComparator());
			if (pStart< 0) 
				pStart= (-(pStart+1))-1; // +1, -1 for the node before ins point

				// collect all normal edges the regions align to
				// regs are exonic regions from 1 read, they are each contained in exactly one edge
				// NO !!! a read can span more than one exonic stretch !
			Node n= nodes[pStart];
			
			while (n.getSite().getPos()< prime3||
					(n.getSite().isLeftFlank()&& n.getSite().getPos()<= prime3)) {
				Vector<Edge> outEdges= n.getOutEdges();
				int sizeB4= v.size();
				Edge iEdge= null;
				for (int j = 0; j < outEdges.size(); j++) {
					if (outEdges.elementAt(j).isExonic()) {							
						v.add(outEdges.elementAt(j));	// there can only be one exonic outedge
						n= outEdges.elementAt(j).getHead();
						break;	// inner loop
					} else
						iEdge= outEdges.elementAt(j);
				}
				if (v.size()== sizeB4) {
					// lemma: there can only be one
//					if (iEdge!= null) {
//						if (iEdge.getTail().getSite().getPos()< rr.get3PrimeEdge())	// < corrects for exon flank pos
//							v.add(iEdge);
//					}
					break;
				}
			}
			
		}
		
		// with introns cannot longer be
		if (v.size()== 0)
			return null;	// can be now, due to antisense mapping !
		
		if (v.size()== 1) {
			int a= Math.abs(v.elementAt(0).getTail().getSite().getPos()), b= Math.abs(v.elementAt(0).getHead().getSite().getPos());
			int edgeLeft= Math.min(a,b), edgeRight= Math.max(a,b);
			a= obj.getStart(); // Math.abs(obj.getAbsBlockStart(0)); 
			b= obj.getEnd(); // Math.abs(obj.getAbsBlockEnd(bMax- 1)- 1);
			int regLeft= Math.min(a,b), regRight= Math.max(a,b);
			if (edgeLeft> regLeft|| edgeRight< regRight)
				return null; // exceeding gene
			return v.elementAt(0);	// finished, read spans only one edge
		}
		
		assert(v.size()!= 0);
		Edge[] ve= new Edge[v.size()];
		for (int i = 0; i < ve.length; i++) 
			ve[i]= v.elementAt(i);
		Arrays.sort(ve, SuperEdge.getDefaultEdgeByTailComparator());	// sort for comparison
		
		// else, check if there is already a superedge with this edge-set
		Vector<SuperEdge> seV= ve[0].getSuperEdges();		
		SuperEdge se= null;
		for (int j = 0; seV!= null&& j < seV.size(); j++) {
			if (seV.elementAt(j).isPend())
				continue;
			Edge[] e= seV.elementAt(j).getEdges();	// these are sorted
			if (e.length!= ve.length)
				continue;
			int k;
			for (k = 0; k < e.length; k++) {
				if (e[k]!= ve[k])	// neg strand not necesarily sorted..
					break;
			}
			if (k== e.length) {
				se= seV.elementAt(j);
				break;	// se found
			}
		}
		if (se== null) {
			// TODO can be now due to antisense overlap
			return null;
		}
		int a= Math.abs(se.getTail().getSite().getPos()), b= Math.abs(se.getHead().getSite().getPos());
		int edgeLeft= Math.min(a,b), edgeRight= Math.max(a,b);
		a= obj.getStart(); // Math.abs(obj.getAbsBlockStart(0)- 1); 
		b= obj.getEnd(); // Math.abs(obj.getAbsBlockEnd(bMax- 1)- 1);
		int regLeft= Math.min(a,b), regRight= Math.max(a,b);
		if (edgeLeft> regLeft|| edgeRight< regRight)
			return null;
//		if (se== null) {
//			Edge[] e= new Edge[v.size()];
//			for (int j = 0; j < e.length; j++) 
//				e[j]= v.elementAt(j);
//			se= new SuperEdge(e);
//		}
		return se;

	}

	Attributes getAttributes(BEDobject2 o, UniversalReadDescriptor d, Attributes attributes) {
		
		CharSequence tag= o.getName();
		attributes= d.getAttributes(tag, attributes);
		if (attributes== null) {
			Log.warn("Error in read ID: could not parse read identifier "+ tag);
			return null;
		}
		if (d.isPaired()&& attributes.flag<= 0) {
			Log.warn("Error in read ID: could not find mate in "+ tag);
			return null;
		}
		if (d.isStranded()&& attributes.strand< 0) {
			Log.warn("Error in read ID: could not find strand in "+ tag);
			return null;
		}
		return attributes;
	}
	
	public void map(BufferedIterator beds, UniversalReadDescriptor descriptor2) {

			if (beds== null) 
				return;
		
			// init
			BEDobject2 dobject, dobject2;
			CharSequence lastName= null;
			UniversalReadDescriptor.Attributes 
				attributes= descriptor2.createAttributes(), 
				attributes2= descriptor2.createAttributes();
			boolean paired= descriptor2.isPaired();
			boolean stranded= descriptor2.isStranded();
			nrMappingsLocus= 0;
			nrMappingsLocusMultiMaps= 0;
			nrMappingsMapped= 0;
			nrMappingsNotMapped= 0;
			nrMappingsNotMappedAsPair= 0;
			nrMappingsWrongPairOrientation= 0;
			nrMappingsWrongStrand= 0;
		
			// map read pairs
			while (beds.hasNext()) {
				
				dobject= new BEDobject2(beds.next());
				++nrMappingsLocus;
				CharSequence name= dobject.getName();
				if (name.equals(lastName))
					++nrMappingsLocusMultiMaps;
				lastName= name; 
				
				attributes= getAttributes(dobject, descriptor2, attributes);
				if (paired&& attributes.flag== 2)	// don't iterate twice, for counters
					continue;
				Edge target= getEdge2(dobject);
				if (target== null) {
					++nrMappingsNotMapped;
					continue;	// couldn't map
				}

				byte refStrand= trpts[0].getStrand();	// TODO get from edge
				if (stranded) {
					boolean sense= dobject.getStrand()== refStrand;
					byte dir= attributes.strand;
					if ((dir== 2&& sense)|| (dir== 1&& !sense)) {
						++nrMappingsWrongStrand;
						continue;
					}
				}

				if (paired) {
					
					// scan for mates
					beds.mark();
					while (beds.hasNext()) {
						dobject2= new BEDobject2(beds.next());
						attributes2= getAttributes(dobject2, descriptor2, attributes2);
						if (!attributes.id.equals(attributes2.id))
							break;						
						if (attributes2== null|| attributes2.flag== 1)
							continue;

						Edge target2= getEdge2(dobject2);
						if (target2== null) {
							++nrMappingsNotMapped;
							continue;
						}

						// check again strand in case one strand-info had been lost
						if (stranded) {
							boolean sense= dobject2.getStrand()== refStrand;
							byte dir= attributes2.strand;
							if ((dir== 2&& sense)|| (dir== 1&& !sense)) {
								++nrMappingsWrongStrand;
								continue;
							}
						}
						
						// check directionality (sequencing-by-synthesis)
						// 20101222: check also that the leftmost (in genomic direction) 
						// is sense (in genomic direction)
						if (dobject.getStrand()== dobject2.getStrand()
								|| (dobject.getStart()< dobject2.getStart()&& dobject.getStrand()!= Transcript.STRAND_POS)
								|| (dobject2.getStart()< dobject.getStart()&& dobject2.getStrand()!= Transcript.STRAND_POS)) {
							nrMappingsWrongPairOrientation+= 2;
							continue;
						}
						
						// find common super-edge
						Vector<Edge> w= new Vector<Edge>();
						if (target.getFrac(true)< target2.getFrac(true)) {
							w.add(target);
							w.add(target2);
						} else {
							w.add(target2);
							w.add(target);
						}
						SuperEdge se= getSuperEdge(w, true, null);
						if (se== null) {
							nrMappingsNotMappedAsPair+= 2;
							continue;	
						}
						se.incrReadNr();
						nrMappingsMapped+= 2;	
					}
					beds.reset();

				} else {	// single reads, strand already checked
					boolean sense= trpts[0].getStrand()== dobject.getStrand();	// TODO get from edge
					if (sense)
						target.incrReadNr();
					else
						target.incrRevReadNr();
					++nrMappingsMapped;
				}
			}
			setMappedReads(nrMappingsMapped);	// TODO 


		}

	public long getNrMappingsMapped() {
		return nrMappingsMapped;
	}

	public long getNrMappingsNotMappedAsPair() {
		return nrMappingsNotMappedAsPair;
	}

	public long getNrMappingsWrongPairOrientation() {
		return nrMappingsWrongPairOrientation;
	}

	public Edge getEdge_old(BEDobject obj) {
				
				Vector<Edge> v= edgeVector; //new Vector<Edge>();
				v.removeAllElements();
				
				Node[] nodes= getNodesInGenomicOrder();
				boolean parallel= true;
				if (obj.getStrand()!= gene.getStrand())
					parallel= false;
				int bMax= obj.getBlockCount();
				if (bMax== 0)
					++bMax;
				SpliceSite dummySite= new SpliceSite(0, SpliceSite.TYPE_ACCEPTOR, null);
				Node dummyNode= new Node(dummySite, null);
				for (int i = parallel?0:(bMax- 1); parallel?i < bMax:i>=0; i+=parallel?1:(-1)) {
					
					// convert it to directionality of trpt, 
					// now with reads on antisense strand
					int prime5= (gene.getStrand()>= 0?obj.getAbsBlockStart(i):obj.getAbsBlockEnd(i))* gene.getStrand(),
						prime3= (gene.getStrand()>= 0?obj.getAbsBlockEnd(i):obj.getAbsBlockStart(i))* gene.getStrand();
					
					dummyNode.getSite().setPos(prime5);
					int pStart= Arrays.binarySearch(nodes, dummyNode, Node.getDefaultPositionTypeComparator());
					if (pStart< 0) 
						pStart= (-(pStart+1))-1; // +1, -1 for the node before ins point
		
						// collect all normal edges the regions align to
						// regs are exonic regions from 1 read, they are each contained in exactly one edge
						// NO !!! a read can span more than one exonic stretch !
					Node n= nodes[pStart];
					
					while (n.getSite().getPos()< prime3||
							(n.getSite().isLeftFlank()&& n.getSite().getPos()<= prime3)) {
						Vector<Edge> outEdges= n.getOutEdges();
						int sizeB4= v.size();
						Edge iEdge= null;
						for (int j = 0; j < outEdges.size(); j++) {
							if (outEdges.elementAt(j).isExonic()) {							
								v.add(outEdges.elementAt(j));	// there can only be one exonic outedge
								n= outEdges.elementAt(j).getHead();
								break;	// inner loop
							} else
								iEdge= outEdges.elementAt(j);
						}
						if (v.size()== sizeB4) {
							// lemma: there can only be one
	//						if (iEdge!= null) {
	//							if (iEdge.getTail().getSite().getPos()< rr.get3PrimeEdge())	// < corrects for exon flank pos
	//								v.add(iEdge);
	//						}
							break;
						}
					}
					
				}
				
				// with introns cannot longer be
				if (v.size()== 0)
					return null;	// can be now, due to antisense mapping !
				
				if (v.size()== 1) {
					int a= Math.abs(v.elementAt(0).getTail().getSite().getPos()), b= Math.abs(v.elementAt(0).getHead().getSite().getPos());
					int edgeLeft= Math.min(a,b), edgeRight= Math.max(a,b);
					a= Math.abs(obj.getAbsBlockStart(0)); 
					b= Math.abs(obj.getAbsBlockEnd(bMax- 1)- 1);
					int regLeft= Math.min(a,b), regRight= Math.max(a,b);
					if (edgeLeft> regLeft|| edgeRight< regRight)
						return null; // exceeding gene
					return v.elementAt(0);	// finished, read spans only one edge
				}
				
				assert(v.size()!= 0);
				Edge[] ve= new Edge[v.size()];
				for (int i = 0; i < ve.length; i++) 
					ve[i]= v.elementAt(i);
				Arrays.sort(ve, SuperEdge.getDefaultEdgeByTailComparator());	// sort for comparison
				
				// else, check if there is already a superedge with this edge-set
				Vector<SuperEdge> seV= ve[0].getSuperEdges();		
				SuperEdge se= null;
				for (int j = 0; seV!= null&& j < seV.size(); j++) {
					Edge[] e= seV.elementAt(j).getEdges();	// these are sorted
					if (e.length!= ve.length)
						continue;
					int k;
					for (k = 0; k < e.length; k++) {
						if (e[k]!= ve[k])	// neg strand not necesarily sorted..
							break;
					}
					if (k== e.length) {
						se= seV.elementAt(j);
						break;	// se found
					}
				}
				if (se== null) {
					// TODO can be now due to antisense overlap
					return null;
				}
				int a= Math.abs(se.getTail().getSite().getPos()), b= Math.abs(se.getHead().getSite().getPos());
				int edgeLeft= Math.min(a,b), edgeRight= Math.max(a,b);
				a= Math.abs(obj.getAbsBlockStart(0)- 1); 
				b= Math.abs(obj.getAbsBlockEnd(bMax- 1)- 1);
				int regLeft= Math.min(a,b), regRight= Math.max(a,b);
				if (edgeLeft> regLeft|| edgeRight< regRight)
					return null;
		//		if (se== null) {
		//			Edge[] e= new Edge[v.size()];
		//			for (int j = 0; j < e.length; j++) 
		//				e[j]= v.elementAt(j);
		//			se= new SuperEdge(e);
		//		}
				return se;
		
			}

	/**
		 * maps read to atomic edge or to SJ
		 * @param obj
		 * @return
		 */
		public Edge getEdge2(BEDobject2 obj) {
			
			int bcount= obj.getBlockCount();
			int bstart= obj.getStart()+ 1, bend= obj.getEnd();	// to normal space
	
	
			if (bcount< 2) {
				
				Edge e= getEdge2(bstart, bend);
				return e;
				
			} else {	// split-maps
				assert(bcount== 2);	// paolo only maps 1 split
				
	//			if (obj.getName().equals("HWUSI-EAS626_1:5:92:163:105/2"))
	//				System.currentTimeMillis();
				
				int size= obj.getNextBlockSize(), 			
					size2= obj.getNextBlockSize();	// ask for all
				
				Edge e= getEdge2(bstart, bstart+ size- 1);
				if (e== null)
					return null;
				else {
					int epos= trpts[0].getStrand()>= 0? e.getHead().getSite().getPos():
						Math.abs(e.getTail().getSite().getPos());
					if (epos!= bstart+ size- 1)
						return null;
				}
				
				Edge f= getEdge2(bend- size2+ 1, bend);
				if (f== null)
					return null;
				else {
					int epos= trpts[0].getStrand()>= 0? f.getTail().getSite().getPos():
						Math.abs(f.getHead().getSite().getPos());
					if (epos!= bend- size2+ 1)
						return null;
				}
				
				// either of them can be a superEdge, but no pe
				if (e== f)
					return null;	// intron maps to exonic stretch
				// chr19	1612377	1615303	HWUSI-EAS627_1:6:84:296:216/1	1	-	0	0	0,0,0	2	55,20	0,2906
				// new transcript evidence
				if (isNull(intersect(e.getTranscripts(), f.getTranscripts())))
					return null;
				// get intermediate intron
				Node g= null, h= null;
				if (e.getHead().getSite().getPos()< f.getHead().getSite().getPos()) { 
					g= e.getHead();
					h= f.getTail();
				} else {
					g= f.getHead();
					h= e.getTail();
				}
				long[] part= null;
				for (int i = 0; i < g.getOutEdges().size(); i++) {
					if (g.getOutEdges().elementAt(i).isExonic())
						continue;
					if (g.getOutEdges().elementAt(i).getHead()== h) {
						part= g.getOutEdges().elementAt(i).getTranscripts();
						break;
					}
				}
				if (part== null)
					return null;
				
				Vector<Edge> v= new Vector<Edge>();
				if (e instanceof SuperEdge) {
					SuperEdge se= (SuperEdge) e;
					assert(!se.isPend());
					Edge[] ee= se.getEdges();
					for (int i = 0; i < ee.length; i++) 
						v.add(ee[i]);
				} else
					v.add(e);
				if (f instanceof SuperEdge) {
					SuperEdge se= (SuperEdge) f;
					assert(!se.isPend());
					Edge[] ee= se.getEdges();
					for (int i = 0; i < ee.length; i++) 
						v.add(ee[i]);
				} else
					v.add(f);
				
				return getSuperEdge(v, false, part);
			}
		}

	/**
		 * maps a continous stretch of positions between
		 * gstart and gend to an atomic edge, or, a series
		 * of edges corresponding to an EEJ
		 *  
		 * @param bstart start on positive strand
		 * @param bend end on positive strand
		 * @return
		 */
		public Edge getEdge2(int bstart, int bend) {
			
			Vector<Edge> v= new Vector<Edge>();
			int[] su= getSpliceUniverse();
			Node[] nodes= getNodesInGenomicOrder();
			
			int gstart= bstart, gend= bend;
			byte strand= trpts[0].getStrand();
			if (strand< 0) {	// neg strand
				gstart= -bend;
				gend= -bstart;
			}
			
	//		if (gstart== 26755100&& gend== 26755174)
	//			System.currentTimeMillis();
			
			// must be 2 or more adjacent exon-exon junctions
			int p= Arrays.binarySearch(su, gstart);	// anchor<= 5'end of read
			if (p< 0) {
				p= -(p+ 1);	// falls before
				if (nodes[p].getSite().isLeftFlank()
					&& nodes[p- 1].getSite().isRightFlank()) {	// p!= 0, for src
					
					boolean found= checkExonicEdge(nodes[p- 1], nodes[p]);
					if (!found)
						return null;	// in intron
				}
				//else, in both cases
				--p;	// p> 0, src
			} else {	// hits exact
				if (nodes[p].getSite().isRightFlank())	// always
					--p; // p+= strand>= 0? -1: 1;
			}
			
			int q= Arrays.binarySearch(su, gend);	// anchor>= 3'end of read
			if (q< 0) {
				q= -(q+ 1);	// falls before
				if (nodes[q].getSite().isLeftFlank()&&
						nodes[q- 1].getSite().isRightFlank()) {
	
					boolean found= checkExonicEdge(nodes[q- 1], nodes[q]);
					if (!found)
						return null;	// in intron
				}// else nothing, falls before q marks end
			} else {	// hits exact			
				if (nodes[q].getSite().isLeftFlank())
					++q;
			}
				
			
			// get chain of edges, if exists
			Node head= nodes[p];
			while (head!= nodes[q]) {
				if (head.getOutEdges().size()<= 0)
					System.currentTimeMillis();
				assert(head.getOutEdges().size()> 0);
				Iterator<Edge> iter= head.getOutEdges().iterator();
				boolean found= false;
				while (iter.hasNext()) {
					Edge e= iter.next();
					if (!e.isExonic())
						continue;
					v.add(e);
					head= e.getHead();
					found= true;
					break;
				}
				if (!found)
					return null;	// intron in between
			}
			if (v.size()== 0)
				System.currentTimeMillis();
			if (v.size()== 0)
				System.currentTimeMillis();
			assert(v.size()> 0);
			
			// trivial case
			if (v.size()== 1)
				return v.elementAt(0);
			
			// else it is a EEJ
			return getSuperEdge(v, false, null);
	
		}

	public long getMappedReads() {
		return mappedReads;
	}

	public void setMappedReads(long mappedReads) {
		this.mappedReads = mappedReads;
	}

}
