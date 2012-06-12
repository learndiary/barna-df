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

import barna.commons.log.Log;
import barna.flux.capacitor.reconstruction.FluxCapacitorSettings;
import barna.io.BufferedIterator;
import barna.io.rna.UniversalReadDescriptor;
import barna.io.rna.UniversalReadDescriptor.Attributes;
import barna.model.*;
import barna.model.bed.BEDMapping;
import barna.model.splicegraph.*;
import barna.model.splicegraph.Node;
import org.w3c.dom.*;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

/**
 * Class to handle the <i>annotation</i> mapping of <i>genomic</i> mappings
 * in a locus (gene).
 *
 * @author Michael Sammeth
 */
public class AnnotationMapper extends SplicingGraph {

    /**
     * The number of mappings in the current locus.
     */
	public long nrMappingsLocus= 0;

    /**
     * The number of mappings in the current locus that map to the annotation.
     */
	public long nrMappingsMapped= 0;

    /**
     * The number of mappings in the current locus that do not map to the annotation
     */
	public long nrMappingsNotMapped= 0;

    /**
     * The number of mappings that mapped to the annotation but no mate was found/mapped.
     */
	public long nrMappingsNotMappedAsPair= 0;

    /**
     * The number of mappings that did not match the expected strand.
     */
	public long nrMappingsWrongStrand= 0;

    /**
     * The number of paired and annotation-mapped reads that did not exhibit the expected orientation.
     */
	public long nrMappingsWrongPairOrientation= 0;

    /**
     * The number of mappings that are mapped multiple times in the locus.
     */
	public long nrMappingsLocusMultiMaps= 0;
	
	public AnnotationMapper(Gene gene) {
		super(gene);
		constructGraph();
		getNodesInGenomicOrder();	// important ??!
		transformToFragmentGraph();
	}
	
	public AbstractEdge getEdge(BEDMapping obj) {
			
		Vector<SimpleEdge> v= edgeVector; //new Vector<Edge>();
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
				Vector<SimpleEdge> outEdges= n.getOutEdges();
				int sizeB4= v.size();
				SimpleEdge iEdge= null;
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
		SimpleEdge[] ve= new SimpleEdge[v.size()];
		for (int i = 0; i < ve.length; i++) 
			ve[i]= v.elementAt(i);
		Arrays.sort(ve, SuperEdge.getDefaultEdgeByTailComparator());	// sort for comparison
		
		// else, check if there is already a superedge with this edge-set
		Vector<SuperEdge> seV= ve[0].getSuperEdges();		
		SuperEdge se= null;
		for (int j = 0; seV!= null&& j < seV.size(); j++) {
			if (seV.elementAt(j).isPend())
				continue;
			AbstractEdge[] e= seV.elementAt(j).getEdges();	// these are sorted
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

	Attributes getAttributes(BEDMapping o, UniversalReadDescriptor d, Attributes attributes) {
		
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
			Log.warn("Error in read ID: could not find strand in " + tag);
			return null;
		}
		return attributes;
	}

    /*@Override
    protected int createKeyElements(Document doc, Element graph) {
        int i = super.createKeyElements(doc, graph);
        graph.appendChild(createKeyElement(doc, ++i, "edge", "readNr", "int"));
        graph.appendChild(createKeyElement(doc, ++i, "edge", "revReadNr", "int"));
        return i;
    }

    @Override
    protected Element createEdgeElement(SimpleEdge se, Document doc, int id, String target, String source) {
        Element edge = super.createEdgeElement(se,doc,id,target,source);
        Element data = doc.createElement("data");
        data.setAttribute("key","d4");
        data.appendChild(doc.createTextNode(Integer.toString(((SimpleEdgeMappings) se).getMappings().getReadNr())));
        edge.appendChild(data);
        data = doc.createElement("data");
        data.setAttribute("key","d5");
        data.appendChild(doc.createTextNode(Integer.toString(((SimpleEdgeMappings) se).getMappings().getRevReadNr())));
        edge.appendChild(data);
        return edge;
    }

    @Override
    protected Boolean findEdgeElement(org.w3c.dom.Node edge) {
        String head = edge.getChildNodes().item(0).getChildNodes().item(0).getNodeValue();
        String tail = edge.getChildNodes().item(1).getChildNodes().item(0).getNodeValue();
        String transcripts = edge.getChildNodes().item(2).getChildNodes().item(0).getNodeValue();
        int readNr = Integer.parseInt(edge.getChildNodes().item(3).getChildNodes().item(0).getNodeValue());
        int revReadNr = Integer.parseInt(edge.getChildNodes().item(4).getChildNodes().item(0).getNodeValue());
        for (SimpleEdge se : edgeHash.values())  {
            if (se.getHead().getSite().toString().equals(head)  && se.getTail().getSite().toString().equals(tail) && Arrays.toString(decodeTset(se.getTranscripts())).equals(transcripts) && ((SimpleEdgeMappings)se).getMappings().getReadNr() == readNr && ((SimpleEdgeMappings)se).getMappings().getRevReadNr() == revReadNr)
                return true;
        }
        return false;
    }   */

    /**
	 * Maps genome-mapped reads into the graph.
	 * 
	 * @param lineIterator iterator of input lines
	 * @param settings
	 */
	public void map(BufferedIterator lineIterator, FluxCapacitorSettings settings) {

			if (lineIterator== null) 
				return;
		
			UniversalReadDescriptor descriptor= settings.get(FluxCapacitorSettings.READ_DESCRIPTOR);
			File insertFile= settings.get(FluxCapacitorSettings.INSERT_FILE);
			BufferedWriter buffy= null;
			if (insertFile!= null) 
				try {
					buffy= new BufferedWriter(new FileWriter(insertFile, true));
				} catch (Exception e) {
					throw new RuntimeException(e);
				}
			
			// init
			BEDMapping dobject, dobject2;
			CharSequence lastName= null;
			UniversalReadDescriptor.Attributes 
				attributes= descriptor.createAttributes(), 
				attributes2= descriptor.createAttributes();
			boolean paired= descriptor.isPaired();
			boolean stranded= descriptor.isStranded();
			nrMappingsLocus= 0;
			nrMappingsLocusMultiMaps= 0;
			nrMappingsMapped= 0;
			nrMappingsNotMapped= 0;
			nrMappingsNotMappedAsPair= 0;
			nrMappingsWrongPairOrientation= 0;
			nrMappingsWrongStrand= 0;
		
			// map read pairs
			while (lineIterator.hasNext()) {

				dobject= new BEDMapping(lineIterator.next());
				++nrMappingsLocus;
				CharSequence name= dobject.getName();
				if (name.equals(lastName))
					++nrMappingsLocusMultiMaps;
				lastName= name; 
				
				attributes= getAttributes(dobject, descriptor, attributes);
				if (paired&& attributes.flag== 2)	// don't iterate twice, for counters
					continue;
				AbstractEdge target= getEdge2(dobject);
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
					lineIterator.mark();
					while (lineIterator.hasNext()) {
						dobject2= new BEDMapping(lineIterator.next());
						attributes2= getAttributes(dobject2, descriptor, attributes2);
						if (!attributes.id.equals(attributes2.id))
							break;						
						if (attributes2== null|| attributes2.flag== 1)
							continue;

						AbstractEdge target2= getEdge2(dobject2);
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
						Vector<AbstractEdge> w= new Vector<AbstractEdge>();
						if (target.getDelimitingPos(true)< target2.getDelimitingPos(true)) {
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
						((SuperEdgeMappings) se).getMappings().incrReadNr();
						nrMappingsMapped+= 2;
						if (buffy!= null) 
							writeInsert(buffy, se, dobject, dobject2, attributes2.id);
					}
					lineIterator.reset();

				} else {	// single reads, strand already checked
                    boolean sense= trpts[0].getStrand()== dobject.getStrand();	// TODO get from edge
					if (target.isAllIntronic()) {
                        if (sense)
                            ((SimpleEdgeIntronMappings)target).incrReadNr(dobject.getStart(), dobject.getEnd());
                    }   else {
                        if (sense)
                            ((MappingsInterface) target).getMappings().incrReadNr();
                        else
                            ((MappingsInterface) target).getMappings().incrRevReadNr();
                    }
					++nrMappingsMapped;
				}
			} // end: while(iter.hasNext())
			
			// close insert writer
			if (buffy!= null)
				try {
					buffy.close();
				} catch (IOException e) {
					throw new RuntimeException(e);
				}
			
		}

	/**
	 * Writes the insert described by 2 paired mates to the provided writer.
	 * @param buffy writer for output
	 * @param se a super-edge
	 * @param dobject 1st bed object
	 * @param dobject2 2nd bed object
	 */
	protected void writeInsert(BufferedWriter buffy, SuperEdge se,
			BEDMapping dobject, BEDMapping dobject2, CharSequence id) {
		
		Transcript[] tt= decodeTset(se.getTranscripts());
		int[] isizes= new int[tt.length];
		int startMin, endMax;
		if (dobject.getStart()< dobject2.getStart()) {
			startMin= dobject.getStart();
			endMax= dobject2.getEnd();
		} else {
			startMin= dobject2.getStart();
			endMax= dobject.getEnd();
		}
		for (int i = 0; i < tt.length; i++) {
			int epos1= tt[i].getExonicPosition(startMin),
					epos2= tt[i].getExonicPosition(endMax);
			if (epos1== Integer.MIN_VALUE|| epos2== Integer.MIN_VALUE)
				isizes[i]= -1;
			else
				isizes[i]= epos2- epos1+ 1;
			
		}
		
		// make non-redundant
		Arrays.sort(isizes);
		int v= -1, c= 0;
		for (int i = 0; i < isizes.length; i++) {
			if (isizes[i]>= 0&& isizes[i]!= v)
				++c;
			v= isizes[i];
		}
		
		// output
		String pfx= dobject.getChromosome()+ "\t"+ startMin+ "\t"+ endMax+ "\t"+ id+ "\t",
				sfx= "\t"+ (tt[0].isForward()? "+": "-")+ "\t0\t0\t";
		v= -1;
		for (int i = 0, cc= 0; i < isizes.length; i++) {
			if (isizes[i]< 0)
				continue;
			if (isizes[i]!= v) {
				++cc;
				try {
					buffy.write(pfx+ Integer.toString(isizes[i])+ sfx+ cc+ ","+ Integer.toString(c)+ ","+ tt.length+ "\n");
				} catch (IOException e) {
					throw new RuntimeException(e);
				}
			}
			v= isizes[i];
		}
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

	/**
	 * maps read to atomic edge or to SJ
	 * @param obj
	 * @return
	 */
	public AbstractEdge getEdge2(BEDMapping obj) {
		
		int bcount= obj.getBlockCount();
		int bstart= obj.getStart()+ 1, bend= obj.getEnd();	// to normal space


		if (bcount< 2) {
			
			AbstractEdge e= getEdge2(bstart, bend);				
			return e;
			
		} else {	// split-maps
			
			// paolo only maps 1 split
			//assert(bcount== 2);	
			// yes, Paolo, but the simulator can map more..
			
//			if (obj.getName().equals("HWUSI-EAS626_1:5:92:163:105/2"))
//				System.currentTimeMillis();
			
			int size= obj.getNextBlockSize(), 			
				size2= obj.getNextBlockSize();	// ask for all
			
			AbstractEdge e= getEdge2(bstart, bstart+ size- 1);
			if (e== null)
				return null;
			else {
				int epos= trpts[0].getStrand()>= 0? e.getHead().getSite().getPos():
					Math.abs(e.getTail().getSite().getPos());
				if (epos!= bstart+ size- 1)
					return null;
			}
			
			AbstractEdge f= getEdge2(bend- size2+ 1, bend);
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
			
			Vector<AbstractEdge> v= new Vector<AbstractEdge>();
			if (e instanceof SuperEdge) {
				SuperEdge se= (SuperEdge) e;
				assert(!se.isPend());
				AbstractEdge[] ee= se.getEdges();
				for (int i = 0; i < ee.length; i++) 
					v.add(ee[i]);
			} else
				v.add(e);
			if (f instanceof SuperEdge) {
				SuperEdge se= (SuperEdge) f;
				assert(!se.isPend());
				AbstractEdge[] ee= se.getEdges();
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
	public AbstractEdge getEdge2(int bstart, int bend) {
		
		Vector<AbstractEdge> v= new Vector<AbstractEdge>();
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
			/*if (nodes[p].getSite().isLeftFlank()
				&& nodes[p- 1].getSite().isRightFlank()) {	// p!= 0, for src
				
				boolean found= checkExonicEdge(nodes[p- 1], nodes[p]);
				if (!found)
					return null;	// in intron
			} */
			//else, in both cases
			--p;	// p> 0, src
		} else {	// hits exact
			if (nodes[p].getSite().isRightFlank())	// always
				--p; // p+= strand>= 0? -1: 1;
		}
		
		int q= Arrays.binarySearch(su, gend);	// anchor>= 3'end of read
		if (q< 0) {
			q= -(q+ 1);	// falls before
			/*if (nodes[q].getSite().isLeftFlank()&&
					nodes[q- 1].getSite().isRightFlank()) {

				boolean found= checkExonicEdge(nodes[q- 1], nodes[q]);
				if (!found)
					return null;	// in intron
			}// else nothing, falls before q marks end*/
		} else {	// hits exact			
			if (nodes[q].getSite().isLeftFlank())
				++q;
		}
			
		
		// get chain of edges, if exists
		Node head= nodes[p];
		while (head!=nodes[q]) {

			assert(head.getOutEdges().size()> 0);
			Iterator<SimpleEdge> iter= head.getOutEdges().iterator();
			boolean found= false;
			while (iter.hasNext()) {

				SimpleEdge e= iter.next();
                // ... && ((!e.isExonic())|| (!e.isAllIntronic()))
				if ((nodes[p].getSite().isLeftFlank()&&!e.isExonic())||(nodes[p].getSite().isRightFlank()&&!e.isAllIntronic()))
					continue;
				v.add(e);
				head= e.getHead();
				found= true;
                // only one exonic outedge that leads to next node
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

	@Override
	protected SimpleEdge createSimpleEdge(Node v, Node w, long[] newTset) {
        SimpleEdgeMappings e;
        if (v.getSite().isDonor() && w.getSite().isAcceptor())
            e = new SimpleEdgeIntronMappings(v,w);
        else
            e= new SimpleEdgeMappings(v, w);
		e.setTranscripts(newTset);
		return e;
	}
	
	@Override
	protected SuperEdge createSuperEdge(AbstractEdge[] newEdges,
			long[] newTset, boolean isPairedEnd) {
		SuperEdgeMappings se= new SuperEdgeMappings(newEdges, newTset, isPairedEnd);
		return se;
	}
	
	public void getRPK(ASEvent event, boolean pend, byte etype, Vector<Vector<AbstractEdge>> v) {
	
		SpliceSite[][] c= event.getSpliceChains();
		int[] maxLen= new int[2];
		for (int i = 0; i < c.length; i++) {
			int lastPos= event.getSrc().getPos();
			for (int j = 0; j < c[i].length; lastPos= c[i][j++].getPos()) {
				if (c[i][j].isRightFlank())
					maxLen[i]+= c[i][j].getPos()- lastPos+ 1;
			}
			if (event.getSnk().isRightFlank())
				maxLen[i]+= event.getSnk().getPos()- lastPos+ 1; 
		}
		Transcript[][] t= event.getTranscripts();
		long[][] sig= new long[event.getDimension()][];
		for (int i = 0; i < sig.length; i++) 
			sig[i]= encodeTset(t[i]);
		Node srcNode= getNode(event.getSrc());
		for (int i = 0; i < event.getDimension(); i++) 
			getRPK(srcNode, sig[i], null, maxLen[i], pend, etype, v.elementAt(i));
		
		// add junctions
		Node snkNode= getNode(event.getSnk());
		SimpleEdge[] e= new SimpleEdge[2];
		for (int i = 0; i < srcNode.getInEdges().size(); i++) {
			if (srcNode.getInEdges().elementAt(i).isExonic()) {
				e[0]= srcNode.getInEdges().elementAt(i);
				break;
			}
		}
		for (int i = 0; i < snkNode.getOutEdges().size(); i++) {
			if (snkNode.getOutEdges().elementAt(i).isExonic()) {
				e[1]= snkNode.getOutEdges().elementAt(i);
				break;
			}
		}
		
		for (int x = 0; x < e.length; x++) {
			if (e[x]== null)
				continue;
			for (int i = 0; e[x].getSuperEdges()!= null&& i < e[x].getSuperEdges().size(); i++) {
				SuperEdge se= e[x].getSuperEdges().elementAt(i);
				if ((!pend)&& se.isPend())
					continue;
				add(se, sig, v);
				if (pend)
					for (int j = 0; se.getSuperEdges()!= null&& j < se.getSuperEdges().size(); j++) 
						add(se.getSuperEdges().elementAt(j), sig, v);
			}
		}
	}

	public void getRPK(Exon exon, Transcript tx, boolean pend, byte etype, Vector<Vector<AbstractEdge>> v) {
		
		long[] sig= encodeTset(new Transcript[] {tx});
		int maxLen= exon.getLength();
		Node n= getNode(exon.getAcceptor());
		
		getRPK(n, sig, null, maxLen, pend, etype, v.elementAt(0));
		
	}

	public void getRPK(Gene g, boolean pend, byte edgeType, Vector<Vector<AbstractEdge>> v) {
			
			Iterator<SimpleEdge> i= getEdgeHash().values().iterator();
			while (i.hasNext()) {
				SimpleEdge e= i.next();
				if (!(e.isExonic()))	// || e.isIntronic()
					continue;
				//if (!pend) {
					if (checkEtype(edgeType, e))
						v.elementAt(0).add(e);
				//}
				for (int j = 0; e.getSuperEdges()!= null
								&& j < e.getSuperEdges().size(); j++) {
					SuperEdge se= (SuperEdge) e.getSuperEdges().elementAt(j);
					if (se.getEdges()[0]!= e)
						continue;
	//				if (se.isPend()) {
	//					if (pend&& checkEtype(edgeType, se))
	//						v.elementAt(0).add(se);
	//				} else {
	//					if (pend) 
	//						for (int k = 0; se.getSuperEdges()!= null
	//									&& k < se.getSuperEdges().size(); k++) {
	//							if (checkEtype(edgeType, se))
	//								v.elementAt(0).add(se.getSuperEdges().elementAt(k));
	//						}
	//					else { // EE- or SJ
	//						if (checkEtype(edgeType, se))
	//							v.elementAt(0).add(se);
	//					}
	//				}
					if (se.isPend()) 
						continue;
					if (checkEtype(edgeType, se))
						v.elementAt(0).add(se);
				}
			}
		}

	/*
		 * 							int sf= decodeCount(e.getSuperEdges().elementAt(j).getTranscripts());
								if (forward>= 0) {
									int v= e.getSuperEdges().elementAt(j).getReadNr();
									val[0]+= v;	
									val[1]+= v/ (float) sf;
								}
			assert(len== val[2]);	// must be composed of complete blocks
			assert(val[1]% 1== 0);
			if (!exon.getDonor().isTES())
				val[2]-= readLen- 1;
	
		 */
		public void getRPK(Node n, long[] sig, long[] noSig, int maxLen, boolean pend, byte etype, Vector<AbstractEdge> v) {
			
			int len= 0;
			while(len< maxLen) {
				SimpleEdge e= null;
				int i = 0;
				for (; len< maxLen&& i < n.getOutEdges().size(); i++) {
					if ((!n.getOutEdges().elementAt(i).isExonic())
							|| isNull(intersect(sig, n.getOutEdges().elementAt(i).getTranscripts())))
						continue;
					e= n.getOutEdges().elementAt(i);
					len+= e.length();
					// don't add paired-ends, always single
	//				if ((!pend)&&
					if (
							(noSig== null|| isNull(intersect(e.getTranscripts(), noSig)))
							&& (checkEtype(etype, e)))
						v.add(e);
	//				if (len> rpk[1])
	//					System.currentTimeMillis();
					for (int j = 0; e.getSuperEdges()!= null&& j < e.getSuperEdges().size(); j++) {
						SuperEdge se= e.getSuperEdges().elementAt(j);					
						if (se.getEdges()[0]!= e|| 
								(noSig!= null&& !isNull(intersect(e.getTranscripts(), noSig))))
							continue;
						if (se.length()+ len> maxLen&& !se.isPend())
							continue;
						// don't add super-edges
						//if (((se.isPend()&& pend)|| ((!se.isPend())&& (!pend)))
						if ((!se.isPend())
								&& checkEtype(etype, se))
							v.add(se);
						if (se.isPend())	// && !pend
							continue;
						
						// dont add super-edges
	//					if (pend&& (!se.isPend())) {
	//						for (int k = 0; se.getSuperEdges()!= null&& k < se.getSuperEdges().size(); k++) {
	//							SuperEdge sse= se.getSuperEdges().elementAt(k); 
	//							if (sse.getEdges()[0]!= se)
	//								continue;
	//							assert(se.getSuperEdges().elementAt(k).isPend());
	//							if ((noSig== null
	//									|| isNull(intersect(e.getTranscripts(), noSig)))
	//									&& checkEtype(etype, sse))
	//								v.add(sse);
	//						}
	//					}
					}
					n= e.getHead();
					i= -1;
				}
				if (len< maxLen&& i== n.getOutEdges().size()) {
					Node o= null;
					for (int j = 0; j < n.getOutEdges().size(); j++) {
						SimpleEdge f= n.getOutEdges().elementAt(j);
						if (f.isExonic()|| (isNull(intersect(sig, f.getTranscripts()))))
							continue;
						if (o== null|| f.getHead().getSite().getPos()< o.getSite().getPos())
							o= f.getHead();
					}
					assert(o!= null);
					n= o;
				}
					
			}
			assert(len== maxLen);
		}

	public void getRPK(Transcript tx, boolean pend, byte etype, Vector<Vector<AbstractEdge>> v) {
		
		long[] sig= encodeTset(new Transcript[] {tx});
		int maxLen= tx.getExonicLength();
		Node n= getNode(tx.getExons()[0].getAcceptor());
		
		getRPK(n, sig, null, maxLen, pend, etype, v.elementAt(0));
	}

    /**
     * Given a transcript the method returns the number of reads belonging to the SuperEdges which span
     * across one splice junction.
     * @param t the transcript to follow
     * @return a <code>Map</code> with the <code>SuperEdge</code> string as key and the number of reads as value.
     */
    public Map<String,Integer> getSJReads(Transcript t, boolean paired) {
        Map<SuperEdge,Integer> nodesReads = new HashMap<SuperEdge,Integer>();
        long[] tsupp = encodeTset(t);
        Node n = null;
        for (int i = 1; i < getNodesInGenomicOrder().length-1; i++) {
            n = getNodesInGenomicOrder()[i];
            if (n.getSite().isLeftFlank() && !isNull(intersect(n.getTranscripts(),tsupp))) {
                Vector<SimpleEdge> ev = n.getOutEdges();
                for (SimpleEdge e : ev) {
                    if (!isNull(intersect(e.getTranscripts(),tsupp)) && e.getSuperEdges() != null) {
                        for (SuperEdge se : e.getSuperEdges()) {
                            if (!isNull(intersect(se.getTranscripts(),tsupp)) && se.isSpliceJunction()) {
                                /*if (nodesReads.get(se)!=null)
                                    continue;*/
                                if (paired) {
                                    if (se.getSuperEdges()!=null)
                                        nodesReads.put(se,countReads(se.getSuperEdges(),tsupp));
                                } else {
                                    nodesReads.put(se, ((SuperEdgeMappings) se).getMappings().getReadNr()+((SuperEdgeMappings) se).getMappings().getRevReadNr());
                                }
                            }
                        }
                    }
                }
            }
        }
        TreeMap<String,Integer> sjReads = new TreeMap<String, Integer>();
        for (SuperEdge edge : nodesReads.keySet()) {
            int donor = edge.getEdges()[0].getHead().getSite().isDonor()?edge.getEdges()[0].getHead().getSite().getPos() : -1;
            int acceptor = edge.getEdges()[1].getTail().getSite().isAcceptor()?edge.getEdges()[1].getTail().getSite().getPos() : -1;
            if (donor ==-1 || acceptor ==-1)
                continue;
            if (sjReads.containsKey(donor+"-"+acceptor))
                sjReads.put(donor+"-"+acceptor,sjReads.get(donor+"-"+acceptor)+ nodesReads.get(edge));
            else
                sjReads.put(donor+"-"+acceptor, nodesReads.get(edge));
        }
        return sjReads;
    }

    public Map<String,Integer> getSJReads(boolean paired) {
        Map<SuperEdge,Integer> nodesReads = new HashMap<SuperEdge,Integer>();
        Node n = null;
        for (int i = 1; i < getNodesInGenomicOrder().length-1; i++) {
            n = getNodesInGenomicOrder()[i];
            if (n.getSite().isLeftFlank()) {
                Vector<SimpleEdge> ev = n.getOutEdges();
                for (SimpleEdge e : ev) {
                    if (e.getSuperEdges() != null) {
                        for (SuperEdge se : e.getSuperEdges()) {
                            if (se.isSpliceJunction()) {
                                if (nodesReads.get(se)!=null)
                                    continue;
                                if (paired) {
                                    if (se.getSuperEdges()!=null)
                                        nodesReads.put(se,countReads(se.getSuperEdges()));
                                } else {
                                    nodesReads.put(se, ((SuperEdgeMappings) se).getMappings().getReadNr()+((SuperEdgeMappings) se).getMappings().getRevReadNr());
                                }
                            }
                        }
                    }
                }
            }
        }
        TreeMap<String,Integer> sjReads = new TreeMap<String, Integer>();
        for (SuperEdge edge : nodesReads.keySet()) {
            int donor = edge.getEdges()[0].getHead().getSite().isDonor()?edge.getEdges()[0].getHead().getSite().getPos() : -1;
            int acceptor = edge.getEdges()[1].getTail().getSite().isAcceptor()?edge.getEdges()[1].getTail().getSite().getPos() : -1;
            if (donor ==-1 || acceptor ==-1)
                continue;
            if (sjReads.containsKey(donor+"-"+acceptor))
                sjReads.put(donor+"-"+acceptor,sjReads.get(donor+"-"+acceptor)+ nodesReads.get(edge));
            else
                sjReads.put(donor+"-"+acceptor, nodesReads.get(edge));
        }
        return sjReads;
    }

    protected int countReads(Vector<SuperEdge> seVector) {
        int reads = 0;
        for (SuperEdge se : seVector) {
            /*if (se.getSuperEdges()!=null)
                reads+=countReads(se.getSuperEdges());  */
            int n = 1;
            if (se.getEdges()[0].getClass().isAssignableFrom(SuperEdgeMappings.class) && se.getEdges()[0].equals(se.getEdges()[1]))
                n = 2;
            reads+=n*((SuperEdgeMappings)se).getMappings().getReadNr();
        }
        return reads;
    }

    protected int countReads(Vector<SuperEdge> seVector, long[] tsupp) {
        int reads = 0;
        for (SuperEdge se : seVector) {
            if (!isNull(intersect(se.getTranscripts(),tsupp))) {
                /*if (se.getSuperEdges()!=null)
                    reads+=countReads(se.getSuperEdges(),tsupp); */
                reads+=((SuperEdgeMappings)se).getMappings().getReadNr();
            }
        }
        return reads;
    }

    public void writeSJReads(ZipOutputStream out, boolean paired) throws IOException {
        try {
            out.putNextEntry(new ZipEntry(gene.getGeneID()+"/"));
            for (Transcript t : gene.getTranscripts()) {
                String s = "";
                Map<String,Integer> reads = getSJReads(t,paired);
                out.putNextEntry(new ZipEntry(gene.getGeneID()+"/" + t.getTranscriptID()));
                for (String k : reads.keySet()) {
                    s = new String(k + "\t" + reads.get(k) + "\n");
                    out.write(s.getBytes());
                }
                out.closeEntry();
            }
            out.closeEntry();
        } finally {
            out.flush();
        }
    }


    public Map<String,Integer[]> getAllIntronicReads() {
        Map<String,Integer[]> nodesReads = new HashMap<String,Integer[]>();
        Node n = null;
        for (int i = 1; i < getNodesInGenomicOrder().length-1; i++) {
            n = getNodesInGenomicOrder()[i];
            if (n.getSite().isRightFlank()) {
                Vector<SimpleEdge> ev = n.getOutEdges();
                for (SimpleEdge e : ev) {
                    if (e.isAllIntronic()) {
                        SimpleEdgeIntronMappings e1 = (SimpleEdgeIntronMappings)e;
                        nodesReads.put(e.toString(), new Integer[]{e1.getMappings().getReadNr()+e1.getMappings().getRevReadNr(),e1.getBinCoverage()});
                    }
                }
            }
        }
        return nodesReads;
    }

    public Map<String,Integer[]> getAllIntronicReads(Transcript t) {
        Map<String,Integer[]> nodesReads = new HashMap<String,Integer[]>();
        long[] tsupp = encodeTset(t);
        Node n = null;
        for (int i = 1; i < getNodesInGenomicOrder().length-1; i++) {
            n = getNodesInGenomicOrder()[i];
            if (n.getSite().isRightFlank() && !isNull(intersect(n.getTranscripts(),tsupp))) {
                Vector<SimpleEdge> ev = n.getOutEdges();
                for (SimpleEdge e : ev) {
                    if (e.isAllIntronic() && !isNull(intersect(e.getTranscripts(),tsupp))) {
                        SimpleEdgeIntronMappings e1 = (SimpleEdgeIntronMappings)e;
                        nodesReads.put(e.toString(), new Integer[]{e1.getMappings().getReadNr()+e1.getMappings().getRevReadNr(),e1.getBinCoverage()});
                    }
                }
            }
        }
        return nodesReads;
    }

    public void writeAllIntronicReads(ZipOutputStream out) throws IOException {
        try {
            out.putNextEntry(new ZipEntry(gene.getGeneID()+"/"));
            for (Transcript t : gene.getTranscripts()) {
                String s = "";
                Map<String,Integer[]> reads = getAllIntronicReads();
                out.putNextEntry(new ZipEntry(gene.getGeneID()+"/" + t.getTranscriptID()));
                for (String k : reads.keySet()) {
                    s = new String(k + "\t" + reads.get(k) + "\n");
                    out.write(s.getBytes());
                }
                out.closeEntry();
            }
            out.closeEntry();
        } finally {
            out.flush();
        }
    }

}
