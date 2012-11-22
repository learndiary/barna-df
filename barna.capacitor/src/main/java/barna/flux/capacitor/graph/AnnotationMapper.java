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
import barna.flux.capacitor.graph.ComplexCounter.CounterType;
import barna.flux.capacitor.reconstruction.FluxCapacitorSettings;
import barna.io.MSIterator;
import barna.io.rna.UniversalReadDescriptor;
import barna.io.rna.UniversalReadDescriptor.Attributes;
import barna.model.*;
import barna.model.bed.BEDMapping;
import barna.model.splicegraph.*;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

/**
 * Class to handle the <i>annotation</i> mapping of <i>genomic</i> mappings
 * in a locus (gene).
 *
 * @author Michael Sammeth
 */
public class AnnotationMapper extends SplicingGraph {

    // Size of the window to search for splicing variants next to existing splice sites
    static final int VARIANT_WINDOW= 60;

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

    /**
     * For counting various stuff
     */
    ComplexCounter cc= null;

    /**
     * Default type(s) for counter
     */
    static final EnumSet<CounterType> DEFAULT_COUNTER_TYPES = EnumSet.of(CounterType.SIMPLE);

    public AnnotationMapper(Gene gene) {
        this(gene, DEFAULT_COUNTER_TYPES);
    }

    public AnnotationMapper(Gene gene, EnumSet<CounterType> counterTypes) {
		super(gene);
		constructGraph();
        getNodesInGenomicOrder();    //TODO important ??!
		transformToFragmentGraph();
        if (!counterTypes.isEmpty())
            cc = new ComplexCounter(counterTypes);
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
            starts[0] = obj.getStart();
            ends[0] = obj.getEnd();
        } else {
            obj.resetBlocks();
            for (int i = 0; i < ends.length; i++) {
                starts[i] = obj.getNextBlockStart();
                ends[i] = starts[i] + obj.getNextBlockSize();
                ++starts[i];
            }
        }
        for (int i = parallel ? 0 : (bMax - 1); parallel ? i < bMax : i >= 0; i += parallel ? 1 : (-1)) {

            // convert it to directionality of trpt,
            // now with reads on antisense strand
            int prime5 = (gene.getStrand() >= 0 ?
                    //obj.getAbsBlockStart(i):obj.getAbsBlockEnd(i))* gene.getStrand(),
                    starts[i] : ends[i]) * gene.getStrand(),
                    prime3 = (gene.getStrand() >= 0 ?
                            //obj.getAbsBlockEnd(i):obj.getAbsBlockStart(i))* gene.getStrand();
                            ends[i] : starts[i]) * gene.getStrand();

            dummyNode.getSite().setPos(prime5);
            int pStart = Arrays.binarySearch(nodes, dummyNode, Node.getDefaultPositionTypeComparator());
            if (pStart < 0)
                pStart = (-(pStart + 1)) - 1; // +1, -1 for the node before ins point

            // collect all normal edges the regions align to
            // regs are exonic regions from 1 read, they are each contained in exactly one edge
            // NO !!! a read can span more than one exonic stretch !
            Node n = nodes[pStart];

            while (n.getSite().getPos() < prime3 ||
                    (n.getSite().isLeftFlank() && n.getSite().getPos() <= prime3)) {
                Vector<SimpleEdge> outEdges = n.getOutEdges();
                int sizeB4 = v.size();
                SimpleEdge iEdge = null;
                for (int j = 0; j < outEdges.size(); j++) {
                    if (outEdges.elementAt(j).isExonic()) {
                        v.add(outEdges.elementAt(j));    // there can only be one exonic outedge
                        n = outEdges.elementAt(j).getHead();
                        break;    // inner loop
                    } else
                        iEdge = outEdges.elementAt(j);
                }
                if (v.size() == sizeB4) {
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
        if (v.size() == 0)
            return null;    // can be now, due to antisense mapping !

        if (v.size() == 1) {
            int a = Math.abs(v.elementAt(0).getTail().getSite().getPos()), b = Math.abs(v.elementAt(0).getHead().getSite().getPos());
            int edgeLeft = Math.min(a, b), edgeRight = Math.max(a, b);
            a = obj.getStart(); // Math.abs(obj.getAbsBlockStart(0));
            b = obj.getEnd(); // Math.abs(obj.getAbsBlockEnd(bMax- 1)- 1);
            int regLeft = Math.min(a, b), regRight = Math.max(a, b);
            if (edgeLeft > regLeft || edgeRight < regRight)
                return null; // exceeding gene
            return v.elementAt(0);    // finished, read spans only one edge
        }

        assert (v.size() != 0);
        SimpleEdge[] ve = new SimpleEdge[v.size()];
        for (int i = 0; i < ve.length; i++)
            ve[i] = v.elementAt(i);
        Arrays.sort(ve, SuperEdge.getDefaultEdgeByTailComparator());    // sort for comparison

        // else, check if there is already a superedge with this edge-set
        Vector<SuperEdge> seV = ve[0].getSuperEdges();
        SuperEdge se = null;
        for (int j = 0; seV != null && j < seV.size(); j++) {
            if (seV.elementAt(j).isPend())
                continue;
            AbstractEdge[] e = seV.elementAt(j).getEdges();    // these are sorted
            if (e.length != ve.length)
                continue;
            int k;
            for (k = 0; k < e.length; k++) {
                if (e[k] != ve[k])    // neg strand not necesarily sorted..
                    break;
            }
            if (k == e.length) {
                se = seV.elementAt(j);
                break;    // se found
            }
        }
        if (se == null) {
            // TODO can be now due to antisense overlap
            return null;
        }
        int a = Math.abs(se.getTail().getSite().getPos()), b = Math.abs(se.getHead().getSite().getPos());
        int edgeLeft = Math.min(a, b), edgeRight = Math.max(a, b);
        a = obj.getStart(); // Math.abs(obj.getAbsBlockStart(0)- 1);
        b = obj.getEnd(); // Math.abs(obj.getAbsBlockEnd(bMax- 1)- 1);
        int regLeft = Math.min(a, b), regRight = Math.max(a, b);
        if (edgeLeft > regLeft || edgeRight < regRight)
            return null;
//		if (se== null) {
//			Edge[] e= new Edge[v.size()];
//			for (int j = 0; j < e.length; j++) 
//				e[j]= v.elementAt(j);
//			se= new SuperEdge(e);
//		}
        return se;

    }

	Attributes getAttributes(Mapping mapping, UniversalReadDescriptor desc, Attributes attributes) {

		CharSequence tag= mapping.getName();
		attributes= desc.getAttributes(tag, attributes);
        if (attributes == null) {
            Log.warn("Error in read ID: could not parse read identifier " + tag);
            return null;
        }
		if (desc.isPaired()&& attributes.flag<= 0) {
            Log.warn("Error in read ID: could not find mate in " + tag);
            return null;
        }
		if (desc.isStranded()&& attributes.strand< 0) {
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
     * @param mappings iterator of input lines
     * @param settings
     */
	public void map(MSIterator<Mapping> mappings, FluxCapacitorSettings settings) {

        if (mappings == null)
            return;

        UniversalReadDescriptor descriptor = settings.get(FluxCapacitorSettings.READ_DESCRIPTOR);
        File insertFile = settings.get(FluxCapacitorSettings.INSERT_FILE);
        BufferedWriter buffy = null;
        if (insertFile != null)
            try {
                buffy = new BufferedWriter(new FileWriter(insertFile, true));
            } catch (Exception e) {
                throw new RuntimeException(e);
            }

        // init
			Mapping mapping, otherMapping;
        CharSequence lastName = null;
        UniversalReadDescriptor.Attributes
                attributes = descriptor.createAttributes(),
                attributes2 = descriptor.createAttributes();
        boolean paired = descriptor.isPaired();
        boolean stranded = descriptor.isStranded();
        nrMappingsLocus = 0;
        nrMappingsLocusMultiMaps = 0;
        nrMappingsMapped = 0;
        nrMappingsNotMapped = 0;
        nrMappingsNotMappedAsPair = 0;
        nrMappingsWrongPairOrientation = 0;
        nrMappingsWrongStrand = 0;

        // map read pairs
        while (mappings.hasNext()) {

				mapping= mappings.next();
            ++nrMappingsLocus;
				CharSequence name= mapping.getName();
            if (name.equals(lastName))
                ++nrMappingsLocusMultiMaps;
            lastName = name;

				attributes= getAttributes(mapping, descriptor, attributes);
            if (paired && attributes.flag == 2)    // don't iterate twice, for counters
                continue;
				AbstractEdge target= getEdge2(mapping);
            if (target == null) {
                ++nrMappingsNotMapped;
                continue;    // couldn't map
            }

            byte refStrand = trpts[0].getStrand();    // TODO get from edge
            if (stranded) {
					boolean sense= mapping.getStrand()== refStrand;
                byte dir = attributes.strand;
                if ((dir == 2 && sense) || (dir == 1 && !sense)) {
                    ++nrMappingsWrongStrand;
                    continue;
                }
            }

            if (paired) {

                // scan for mates
//                mappings.mark();
                Iterator<Mapping> mates = mappings.getMates(mapping, descriptor);

                while (mates.hasNext()) {
						otherMapping= mates.next();
//						attributes2= getAttributes(otherMapping, descriptor, attributes2);
//                    if (!attributes.id.equals(attributes2.id))
//                        break;
//                    if (attributes2 == null || attributes2.flag == 1)
//                        continue;

						AbstractEdge target2= getEdge2(otherMapping);
                    if (target2 == null) {
                        ++nrMappingsNotMapped;
                        continue;
                    }

                    // check again strand in case one strand-info had been lost
                    if (stranded) {
							boolean sense= otherMapping.getStrand()== refStrand;
                        byte dir = attributes2.strand;
                        if ((dir == 2 && sense) || (dir == 1 && !sense)) {
                            ++nrMappingsWrongStrand;
                            continue;
                        }
                    }

                    // check directionality (sequencing-by-synthesis)
                    // 20101222: check also that the leftmost (in genomic direction)
                    // is sense (in genomic direction)
						if (mapping.getStrand()== otherMapping.getStrand()
								|| (mapping.getStart()< otherMapping.getStart()&& mapping.getStrand()!= Transcript.STRAND_POS)
								|| (otherMapping.getStart()< mapping.getStart()&& otherMapping.getStrand()!= Transcript.STRAND_POS)) {
                        nrMappingsWrongPairOrientation += 2;
                        continue;
                    }

                    // find common super-edge
                    Vector<AbstractEdge> w = new Vector<AbstractEdge>();
                    if (target.getDelimitingPos(true) < target2.getDelimitingPos(true)) {
                        w.add(target);
                        w.add(target2);
                    } else {
                        w.add(target2);
                        w.add(target);
                    }
                    SuperEdge se = getSuperEdge(w, true, null);
                    if (se == null) {
                        nrMappingsNotMappedAsPair += 2;
                        continue;
                    }

                    if (target.getClass().isAssignableFrom(SimpleEdgeIntronMappings.class) && target.isAllIntronic()) {
                        ((SimpleEdgeIntronMappings) target).incrReadNr(mapping.getStart(), mapping.getEnd(), false);
                    }
                    if (target2.getClass().isAssignableFrom(SimpleEdgeIntronMappings.class) && target2.isAllIntronic() && !target2.equals(target)) {
                        ((SimpleEdgeIntronMappings) target2).incrReadNr(otherMapping.getStart(), otherMapping.getEnd(), false);
                    }
                    ((SuperEdgeMappings) se).getMappings().incrReadNr();
                    if (se.isExonic()) {
                        nrMappingsMapped += 2;
                    }
                    if (buffy != null)
							writeInsert(buffy, se, mapping, otherMapping, attributes2.id);
                }
//                mappings.reset();

            } else {    // single reads, strand already checked
					boolean sense= trpts[0].getStrand()== mapping.getStrand();	// TODO get from edge
                    if (target.isAllIntronic()) {
                        if (sense)
                                ((SimpleEdgeIntronMappings) target).incrReadNr(mapping.getStart(), mapping.getEnd(), true);
                        else
                                ((SimpleEdgeIntronMappings) target).incrRevReadNr(mapping.getStart(), mapping.getEnd(), true);
                    } else {
                        if (sense)
                            ((MappingsInterface) target).getMappings().incrReadNr();
                        else
                            ((MappingsInterface) target).getMappings().incrRevReadNr();
                        ++nrMappingsMapped;
                }
            }
        } // end: while(iter.hasNext())

        // close insert writer
        if (buffy != null)
            try {
                buffy.close();
            } catch (IOException e) {
                throw new RuntimeException(e);
            }

    }

    /**
     * Writes the insert described by 2 paired mates to the provided writer.
     *
     * @param buffy    writer for output
     * @param se       a super-edge
	 * @param mapping 1st bed object
	 * @param otherMapping 2nd bed object
     */
    protected void writeInsert(BufferedWriter buffy, SuperEdge se,
			Mapping mapping, Mapping otherMapping, CharSequence id) {

        Transcript[] tt = decodeTset(se.getTranscripts());
        int[] isizes = new int[tt.length];
        int startMin, endMax;
		if (mapping.getStart()< otherMapping.getStart()) {
			startMin= mapping.getStart();
			endMax= otherMapping.getEnd();
        } else {
			startMin= otherMapping.getStart();
			endMax= mapping.getEnd();
        }
        for (int i = 0; i < tt.length; i++) {
            int epos1 = tt[i].getExonicPosition(startMin),
                    epos2 = tt[i].getExonicPosition(endMax);
            if (epos1 == Integer.MIN_VALUE || epos2 == Integer.MIN_VALUE)
                isizes[i] = -1;
            else
                isizes[i] = epos2 - epos1 + 1;

        }

        // make non-redundant
        Arrays.sort(isizes);
        int v = -1, c = 0;
        for (int i = 0; i < isizes.length; i++) {
            if (isizes[i] >= 0 && isizes[i] != v)
                ++c;
            v = isizes[i];
        }

        // output
		String pfx= mapping.getChromosome()+ "\t"+ startMin+ "\t"+ endMax+ "\t"+ id+ "\t",
                sfx = "\t" + (tt[0].isForward() ? "+" : "-") + "\t0\t0\t";
        v = -1;
        for (int i = 0, cc = 0; i < isizes.length; i++) {
            if (isizes[i] < 0)
                continue;
            if (isizes[i] != v) {
                ++cc;
                try {
                    buffy.write(pfx + Integer.toString(isizes[i]) + sfx + cc + "," + Integer.toString(c) + "," + tt.length + "\n");
                } catch (IOException e) {
                    throw new RuntimeException(e);
                }
            }
            v = isizes[i];
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
     *
     * @param obj
     * @return
     */
	public AbstractEdge getEdge2(Mapping obj) {

        int bcount = obj.getBlockCount();
        int     bstart = obj.getStart() + 1,
                bend = obj.getEnd();    // to normal space


        if (bcount < 2) {

            AbstractEdge e = getEdge2(bstart, bend);
            return e;

        } else {    // split-maps

            // paolo only maps 1 split
            // assert(bcount== 2);
            // yes, Paolo, but the simulator can map more..
            // now also Paolo can map more than 1 split

            int[] bstarts= new int[bcount], bsizes= new int[bcount];
            for (int i = 0; i < bcount; ++i) {
                bstarts[i]= bstart+ obj.getNextBlockStart();
                bsizes[i]= obj.getNextBlockSize();
            }

            // count 5'- and 3'-variant split-mappings
            if (bcount > 1 && cc != null) { // TODO condition for counting
                int[] su= getSpliceUniverse();
                for (int i = 0; i < bcount; ++i) {

                    if (i> 0) { // check left flank
                        boolean valid= false;
                        int p= Arrays.binarySearch(su, bstarts[i]);
                        if (p< 0)
                            p= -(p+ 1);
                        if (p> 0) { // check distance to next upstream site
                            if (Math.abs(bstarts[i]- su[p- 1])<= (VARIANT_WINDOW/ 2)) {
                                valid= checkSpliceSiteType(su[p- 1], SpliceSite.TYPE_ACCEPTOR);
                            }
                        }
                        if ((!valid)&& p< su.length) { // check distance to next downstream site
                            if (Math.abs(bstarts[i]- su[p])<= (VARIANT_WINDOW/ 2)) {
                                valid= checkSpliceSiteType(su[p], SpliceSite.TYPE_ACCEPTOR);
                            }
                        }
                        if (valid) {
                            String id= bstart+ bstarts[i- 1]+ bsizes[i- 1]- 1+ "-"+ bstart+ bstarts[i];
                            cc.increment(id, CounterType.SIMPLE);
                        }
                    }
                    if (i< bcount) { // check right flank
                        boolean valid= false;
                        int p= Arrays.binarySearch(su, bstarts[i]+ bsizes[i]- 1);
                        if (p< 0)
                            p= -(p+ 1);
                        else
                            continue;
                        if (p> 0) { // check distance to next upstream site
                            if (Math.abs(bstarts[i]- su[p- 1])<= (VARIANT_WINDOW/ 2)) {
                                valid= checkSpliceSiteType(su[p- 1], SpliceSite.TYPE_DONOR);
                                if (valid)
                                    break;
                            }
                        }
                        if ((!valid)&& p< su.length) { // check distance to next downstream site
                            if (Math.abs(bstarts[i]- su[p])<= (VARIANT_WINDOW/ 2))
                                valid= checkSpliceSiteType(su[p], SpliceSite.TYPE_DONOR);
                        }
                        if (valid) {
                            String id= bstart+ bstarts[i- 1]+ bsizes[i- 1]- 1+ "^"+ bstart+ bstarts[i];
                            cc.increment(id, CounterType.SIMPLE);
                        }
                    }

                } /* end for all blocks */


            }

                Vector<AbstractEdge> v = new Vector<AbstractEdge>(bcount, 1);
            long[] inter= null, part= null;
            for (int i = 0; i < bcount; ++i) {
                // next coordinates
                int pos= bstarts[i];
                int size = bsizes[i];

                // try to retrieve exonic segment
                AbstractEdge e = getEdge2(pos, pos + size - 1);

                // abort if no segment found
                if (e==null)
                    return null;

                // update common tx support
                inter = (inter== null? e.getTranscripts(): intersect(inter, e.getTranscripts()));

                // abort
                if (isNull(inter)       // no common tx support
                   || v.contains(e))    // novel intron: splice-junction within a single exonic segment
                    return null;

                // check left flank for same coordinate as alignment
                if ((i> 0)&& (Math.abs(trpts[0].getStrand() >= 0 ? e.getTail().getSite().getPos() : e.getHead().getSite().getPos())
                        != pos))
                    return null;
                // check right flank for same coordinate as alignment
                if ((i< bcount- 1)&& (Math.abs(trpts[0].getStrand() >= 0 ? e.getHead().getSite().getPos() : e.getTail().getSite().getPos())
                        != pos+ size- 1))
                    return null;

                // get tx support of intermediate introns, not the same of exon segment intersection!
                if (v.size()> 0) {
                    AbstractEdge f= v.elementAt(v.size()- 1);
                    Node g = null, h = null;
                    if (f.getHead().getSite().getPos() < e.getHead().getSite().getPos()) {
                        g = f.getHead();
                        h = e.getTail();
                    } else {
                        g = e.getHead();
                        h = f.getTail();
                    }
                    for (int j = 0; j < g.getOutEdges().size(); ++j) {
                        AbstractEdge d= g.getOutEdges().elementAt(j);
                        if (d.isExonic()|| d.isAllIntronic())
                            continue;
                        if (d.getHead() == h) {
                            part = (part== null? d.getTranscripts(): intersect(part, d.getTranscripts()));
                            break;
                        }
                    }
                    if (part == null)
                        return null;
                }

                // else, successfully add next segment(s)
                if (e instanceof SuperEdge) {
                    SuperEdge se = (SuperEdge) e;
                    assert (!se.isPend());
                    AbstractEdge[] ee = se.getEdges();
                    for (int j = 0; j < ee.length; ++j)
                        v.add(ee[j]);   // decompose multi-segment hits
                } else
                    v.add(e);

            }

            return getSuperEdge(v, false, part);
        }
    }

    private boolean checkSpliceSiteType(int position, byte type) {

        Vector<SpliceSite> v= gene.getSpliceSites(position);
        for (SpliceSite spliceSite : v) {
            if (spliceSite.isAcceptor()) {
                return true;
            }

        }
        return false;
    }


    /**
     * Maps a continous stretch of positions between
     * <code>bstart</code> and <code>end</code to an atomic edge,
     * or, a series of edges corresponding to an EEJ.
     *
     * @param bstart start on positive strand
     * @param bend   end on positive strand
     * @return the edge mapped
     */
    public AbstractEdge getEdge2(int bstart, int bend) {

        Vector<AbstractEdge> v = new Vector<AbstractEdge>();
        int[] su = getSpliceUniverse(); // copies site positions, TODO comparator for binarySearch() on nodes[]
        Node[] nodes = getNodesInGenomicOrder();

        // get genomic start/end position
        int gstart = bstart, gend = bend;
        byte strand = trpts[0].getStrand();
        if (strand < 0) {    // neg strand
            gstart = -bend;
            gend = -bstart;
        }

        // get the nr of the splice sites (a,b) delimiting the genomic area of the mapping:
        // a is the closest splice site upstream of gstart
        // b is the closest splice site downstream of gend

        // find closest site upstream of genomic start of mapping
        int p = Arrays.binarySearch(su, gstart);
        if (p < 0) {
            p = -(p + 1);   // index after the start of mapping..
            --p;            // ..correct down
        } else {    // mapping starts exactly at the position of a site
            if (nodes[p].getSite().isRightFlank())    // always
                --p; // correct down for the first site upstream of the mapping
        }

        // find closest site downstream of genomic end of mapping
        int q = Arrays.binarySearch(su, gend);
        if (q < 0) {
            q = -(q + 1);    // index after end of mapping, do not correct down
        } else {    // mapping ends exactly at the position of a site
            if (nodes[q].getSite().isLeftFlank())
                ++q; // correct up for the first site downstream of the mapping
        }

        // genomic start after genomic end (should not occur)
        // first site upstream corresponds to/is downstream of first site downstream of mapping
        if (p >= q)
            return null;

        // adjacent sites, exactly one segment between a and b
        // => check whether this segment is an exon or (all-)intronic
        if (p == q - 1 && ((!checkExonicEdge(nodes[p], nodes[q])) && !checkAllIntronicEdge(nodes[p], nodes[q])))
            return null;

        // non-adjacent sites, connected by more than one segment: either consecutive exonic segments, or split-map
        // require that the first segment (p,p+1) and the last segment (q-1,q) are exonic
        // OBS: mappings to all-intronic regions map to ONE single edge
        else {
            if (p < q - 1 && !(checkExonicEdge(nodes[p], nodes[p + 1]) && checkExonicEdge(nodes[q - 1], nodes[q])))
                return null;
        }

        // get the chain of segments between a and b
        Node head = nodes[p];
        while (head != nodes[q]) {

            assert (head.getOutEdges().size() > 0);
            Iterator<SimpleEdge> iter = head.getOutEdges().iterator();

            // iterate out-edges, find the one to continue the segment-chain
            boolean found = false;
            while (iter.hasNext()) {

                SimpleEdge e = iter.next();
                // the backbone of genomic segments is exclusively
                // composed by exonic and all-intronic edges.
                // OBS: intersection of transcript support CANNOT be empty
                if (!e.isExonic() && !e.isAllIntronic()) {
                    continue;
                }
                v.add(e);
                head = e.getHead(); // continue building chain with this edge
                found = true;

                break;
            }

            // continue only if there is at least 1 exonic or all-intronic
            if (!found)
                return null;
        }

        assert (v.size() > 0);  // at least one edge connecting p with q

        // trivial case, when p= q-1
        if (v.size() == 1)
            return v.elementAt(0);

        // otherwise create a corresponding super-edge, method returns null if
        // there is no common transcript support
        return getSuperEdge(v, false, null);

    }

    @Override
    protected SimpleEdge createSimpleEdge(Node v, Node w, long[] newTset, byte type) {
        SimpleEdgeMappings e;
        if (type == SimpleEdge.ALL_INTRONIC) {
            e = new SimpleEdgeIntronMappings(v, w);
        } else
            e = new SimpleEdgeMappings(v, w);
        e.setTranscripts(newTset);
        return e;
    }

    @Override
    protected SuperEdge createSuperEdge(AbstractEdge[] newEdges,
                                        long[] newTset, boolean isPairedEnd) {
        SuperEdgeMappings se = new SuperEdgeMappings(newEdges, newTset, isPairedEnd);
        return se;
    }

    public void getRPK(ASEvent event, boolean pend, byte etype, Vector<Vector<AbstractEdge>> v) {

        SpliceSite[][] c = event.getSpliceChains();
        int[] maxLen = new int[2];
        for (int i = 0; i < c.length; i++) {
            int lastPos = event.getSrc().getPos();
            for (int j = 0; j < c[i].length; lastPos = c[i][j++].getPos()) {
                if (c[i][j].isRightFlank())
                    maxLen[i] += c[i][j].getPos() - lastPos + 1;
            }
            if (event.getSnk().isRightFlank())
                maxLen[i] += event.getSnk().getPos() - lastPos + 1;
        }
        Transcript[][] t = event.getTranscripts();
        long[][] sig = new long[event.getDimension()][];
        for (int i = 0; i < sig.length; i++)
            sig[i] = encodeTset(t[i]);
        Node srcNode = getNode(event.getSrc());
        for (int i = 0; i < event.getDimension(); i++)
            getRPK(srcNode, sig[i], null, maxLen[i], pend, etype, v.elementAt(i));

        // add junctions
        Node snkNode = getNode(event.getSnk());
        SimpleEdge[] e = new SimpleEdge[2];
        for (int i = 0; i < srcNode.getInEdges().size(); i++) {
            if (srcNode.getInEdges().elementAt(i).isExonic()) {
                e[0] = srcNode.getInEdges().elementAt(i);
                break;
            }
        }
        for (int i = 0; i < snkNode.getOutEdges().size(); i++) {
            if (snkNode.getOutEdges().elementAt(i).isExonic()) {
                e[1] = snkNode.getOutEdges().elementAt(i);
                break;
            }
        }

        for (int x = 0; x < e.length; x++) {
            if (e[x] == null)
                continue;
            for (int i = 0; e[x].getSuperEdges() != null && i < e[x].getSuperEdges().size(); i++) {
                SuperEdge se = e[x].getSuperEdges().elementAt(i);
                if ((!pend) && se.isPend())
                    continue;
                add(se, sig, v);
                if (pend)
                    for (int j = 0; se.getSuperEdges() != null && j < se.getSuperEdges().size(); j++)
                        add(se.getSuperEdges().elementAt(j), sig, v);
            }
        }
    }

    public void getRPK(Exon exon, Transcript tx, boolean pend, byte etype, Vector<Vector<AbstractEdge>> v) {

        long[] sig = encodeTset(new Transcript[]{tx});
        int maxLen = exon.getLength();
        Node n = getNode(exon.getAcceptor());

        getRPK(n, sig, null, maxLen, pend, etype, v.elementAt(0));

    }

    public void getRPK(Gene g, boolean pend, byte edgeType, Vector<Vector<AbstractEdge>> v) {

        Iterator<SimpleEdge> i = getEdgeHash().values().iterator();
        while (i.hasNext()) {
            SimpleEdge e = i.next();
            if (!(e.isExonic()))    // || e.isIntronic()
                continue;
            //if (!pend) {
            if (checkEtype(edgeType, e))
                v.elementAt(0).add(e);
            //}
            for (int j = 0; e.getSuperEdges() != null
                    && j < e.getSuperEdges().size(); j++) {
                SuperEdge se = (SuperEdge) e.getSuperEdges().elementAt(j);
                if (se.getEdges()[0] != e)
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

        int len = 0;
        while (len < maxLen) {
            SimpleEdge e = null;
            int i = 0;
            for (; len < maxLen && i < n.getOutEdges().size(); i++) {
                if ((!n.getOutEdges().elementAt(i).isExonic())
                        || isNull(intersect(sig, n.getOutEdges().elementAt(i).getTranscripts())))
                    continue;
                e = n.getOutEdges().elementAt(i);
                len += e.length();
                // don't add paired-ends, always single
                //				if ((!pend)&&
                if (
                        (noSig == null || isNull(intersect(e.getTranscripts(), noSig)))
                                && (checkEtype(etype, e)))
                    v.add(e);
                //				if (len> rpk[1])
                //					System.currentTimeMillis();
                for (int j = 0; e.getSuperEdges() != null && j < e.getSuperEdges().size(); j++) {
                    SuperEdge se = e.getSuperEdges().elementAt(j);
                    if (se.getEdges()[0] != e ||
                            (noSig != null && !isNull(intersect(e.getTranscripts(), noSig))))
                        continue;
                    if (se.length() + len > maxLen && !se.isPend())
                        continue;
                    // don't add super-edges
                    //if (((se.isPend()&& pend)|| ((!se.isPend())&& (!pend)))
                    if ((!se.isPend())
                            && checkEtype(etype, se))
                        v.add(se);
                    if (se.isPend())    // && !pend
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
                n = e.getHead();
                i = -1;
            }
            if (len < maxLen && i == n.getOutEdges().size()) {
                Node o = null;
                for (int j = 0; j < n.getOutEdges().size(); j++) {
                    SimpleEdge f = n.getOutEdges().elementAt(j);
                    if (f.isExonic() || (isNull(intersect(sig, f.getTranscripts()))))
                        continue;
                    if (o == null || f.getHead().getSite().getPos() < o.getSite().getPos())
                        o = f.getHead();
                }
                assert (o != null);
                n = o;
            }

        }
        assert (len == maxLen);
    }

    public void getRPK(Transcript tx, boolean pend, byte etype, Vector<Vector<AbstractEdge>> v) {

        long[] sig = encodeTset(new Transcript[]{tx});
        int maxLen = tx.getExonicLength();
        Node n = getNode(tx.getExons()[0].getAcceptor());

        getRPK(n, sig, null, maxLen, pend, etype, v.elementAt(0));
    }

    /**
     * Count reads which map to splice junction.
     *
     * @param paired whether the annotation mapping consider pairing information or not
     * @return a <code>Map</code> with the junction id string as key and the number of reads as value.
     */
    public Map<String, Integer> getSJReads(boolean paired) {
        Map<SuperEdge, Integer> nodesReads = new HashMap<SuperEdge, Integer>();
        Node n = null;
        for (int i = 1; i < getNodesInGenomicOrder().length - 1; i++) {
            n = getNodesInGenomicOrder()[i];
            if (n.getSite().isLeftFlank()) {
                Vector<SimpleEdge> ev = n.getOutEdges();
                for (SimpleEdge e : ev) {
                    if (e.getSuperEdges() != null) {
                        for (SuperEdge se : e.getSuperEdges()) {
                            if (se.isExonic() && se.isSpliceJunction()) {
                                if (nodesReads.get(se) != null)
                                    continue;
                                if (paired) {
                                    if (se.getSuperEdges() != null)
                                        nodesReads.put(se, countSJReads(se.getSuperEdges()));
                                } else {
                                    nodesReads.put(se, ((SuperEdgeMappings) se).getMappings().getReadNr() + ((SuperEdgeMappings) se).getMappings().getRevReadNr());
                                }
                            }
                        }
                    }
                }
            }
        }
        TreeMap<String, Integer> sjReads = new TreeMap<String, Integer>();
        for (SuperEdge edge : nodesReads.keySet()) {
            for (int i = 0; i < edge.getEdges().length - 1; i++) {
                int donor = -1, acceptor = -1;
                SimpleEdge e = (SimpleEdge) edge.getEdges()[i];
                SimpleEdge e1 = (SimpleEdge) edge.getEdges()[i + 1];
                if (e.equals(e1))
                    continue;
                if (e.getHead().getSite().isDonor() && e1.getTail().getSite().isAcceptor()) {
                    donor = e.getHead().getSite().getPos();
                    acceptor = e1.getTail().getSite().getPos();
                }
                if (donor == -1 || acceptor == -1)
                    continue;
                String sjID = gene.getGeneID(decodeTset(edge.getTranscripts())) + "^" + donor + "^" + acceptor;
                if (sjReads.containsKey(sjID))
                    sjReads.put(sjID, sjReads.get(sjID) + nodesReads.get(edge));
                else
                    sjReads.put(sjID, nodesReads.get(edge));
            }
        }
        return sjReads;
    }

    protected int countSJReads(Vector<SuperEdge> seVector) {
        int reads = 0;
        for (SuperEdge se : seVector) {
            int n = 1;
            if (se.getEdges()[0].getClass().isAssignableFrom(SuperEdgeMappings.class) && se.getEdges()[0].equals(se.getEdges()[se.getEdges().length - 1]))
                n = 2;
            reads += n * ((SuperEdgeMappings) se).getMappings().getReadNr();
        }
        return reads;
    }


    /**
     * Count reads which map to all-intronic regions.
     *
     * @param paired whether the annotation mapping consider pairing information or not
     * @return a <code>Map</code> with the intron id string as key and an array with the count and the
     *         distribution of the reads over the intron as the value.
     */
    public Map<String, Float[]> getAllIntronicReads(boolean paired) {
        Map<String, Float[]> nodesReads = new TreeMap<String, Float[]>();
        Node n = null;
        for (int i = 1; i < getNodesInGenomicOrder().length - 1; i++) {
            n = getNodesInGenomicOrder()[i];
            if (n.getSite().isRightFlank()) {
                Vector<SimpleEdge> ev = n.getOutEdges();
                for (SimpleEdge e : ev) {
                    if (e.isAllIntronic()) {
                        SimpleEdgeIntronMappings e1 = (SimpleEdgeIntronMappings) e;
                        String intronID = gene.getGeneID(decodeTset(e1.getTranscripts())) + "^" + e1.getTail().getSite().getPos() + "^" + e1.getHead().getSite().getPos();
                        if (paired) {
                            if (e1.getSuperEdges() != null)
                                nodesReads.put(intronID, new Float[]{(float) countAllIntronicReads(e1.getSuperEdges()), e1.getBinCoverage()});
                        } else
                            nodesReads.put(intronID, new Float[]{(float) e1.getMappings().getReadNr() + e1.getMappings().getRevReadNr(), e1.getBinCoverage()});
                    }
                }
            }
        }
        return nodesReads;
    }

    protected int countAllIntronicReads(Vector<SuperEdge> seVector) {
        int reads = 0;
        for (SuperEdge se : seVector) {
            if (se.isPend()) {
                int n = 1;
                if (se.getEdges()[0].getClass().isAssignableFrom(SimpleEdgeIntronMappings.class) && se.getEdges()[0].equals(se.getEdges()[se.getEdges().length - 1]))
                    n = 2;
                reads += n * ((SuperEdgeMappings) se).getMappings().getReadNr();
            }
        }
        return reads;
    }
}
