package fbi.genome.sequencing.rnaseq.reconstruction;

import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Vector;

import org.jfree.util.Log;

import fbi.genome.io.Bufferediterator;
import fbi.genome.io.BufferedIteratorMemory;
import fbi.genome.io.rna.UniversalReadDescriptor;
import fbi.genome.io.rna.UniversalReadDescriptor.Attributes;
import fbi.genome.model.DirectedRegion;
import fbi.genome.model.Gene;
import fbi.genome.model.SpliceSite;
import fbi.genome.model.Transcript;
import fbi.genome.model.bed.BEDobject2;
import fbi.genome.model.constants.Constants;
import fbi.genome.model.splicegraph.Edge;
import fbi.genome.model.splicegraph.Graph;
import fbi.genome.model.splicegraph.Node;
import fbi.genome.model.splicegraph.SuperEdge;

public class AnnotationMapper extends Graph {

	long nrMappingsLocus= 0;
	long nrMappingsMapped= 0;
	long nrMappingsNotMapped= 0;
	long nrMappingsNotMappedAsPair= 0;
	long nrMappingsWrongStrand= 0;
	long nrMappingsWrongPairOrientation= 0;
	long nrMappingsLocusMultiMaps= 0;
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
	
	public void map(Bufferediterator beds, UniversalReadDescriptor descriptor2) {

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
				
				dobject= beds.next();
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
						dobject2= beds.next();
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

}
