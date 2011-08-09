package fbi.genome.sequencing.rnaseq.reconstruction;

import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Vector;

import fbi.genome.io.rna.UniversalReadDescriptor;
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

	HashSet<CharSequence> mapReadOrPairIDs;
	UniversalReadDescriptor.Attributes attributes= null;
	UniversalReadDescriptor descriptor2;
	HashMap<CharSequence, Vector<BEDobject2>[]> mapEndsOfPairs;
	long nrReadsMapped= 0;
	long nrReadsLoci= 0;
	long nrMappingsForced= 0;
	long nrMappingsReadsOrPairs= 0;
	long nrMappingsWrongStrand= 0;
	long nrPairsWrongOrientation= 0;
	long nrPairsNoTxEvidence= 0;
	long nrLocusMultimaps= 0;
	public AnnotationMapper(Gene gene) {
		super(gene);
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

	/**
	 * add a SINGLE read
	 * @param regs
	 */
	int mapRead(BEDobject2 dobject, boolean force) {
		// find the edge(s) where the regions align
		// if these do not form a continous chain, create a new edge
		
//			GFFObject[] gtfs= GFFObject.fromBed(dobject);	// TODO kill intermediary GTFs !!!
//			DirectedRegion[] regs= new DirectedRegion[gtfs.length];
//			for (int i = 0; i < regs.length; i++) 
//				regs[i]= new DirectedRegion(gtfs[i]);
		
		boolean paired= descriptor2.isPaired();
		boolean stranded= descriptor2.isStranded();
		
		if (force&& mapReadOrPairIDs.contains(dobject.getName())) {
			return 0;
		}
		
		CharSequence tag= dobject.getName();
		attributes= descriptor2.getAttributes(tag, attributes);
		if (attributes== null) {
			if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) 
				System.err.println("Invalid read identifier "+ tag);
			return 0;	
		}
		byte flag= 0; //getFlag(dobject); 
		if (paired) {
			flag= attributes.flag;
			if (flag<= 0) {
				if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) 
					System.err.println("Error in readID:\n"+ dobject.getName());
				return 0;
			}
		}
		CharSequence ID= tag;	//getID(dobject);
		if (paired|| stranded) 
			ID= attributes.id;

//				if (ID.equals("HWUSI-EAS626_1:5:82:1446:1847"))
//					System.currentTimeMillis();
		
		Edge target= getEdge2(dobject);
		
		if (target== null)
			return 0;
		if (force) {
			boolean sense= trpts[0].getStrand()== dobject.getStrand();
			if (sense)
				target.incrReadNr();
			else
				target.incrRevReadNr();
			mapReadOrPairIDs.add(dobject.getName());
			return 1;
		}
		
		
		byte refStrand= trpts[0].getStrand();
		if (stranded) {
			boolean sense= dobject.getStrand()== refStrand;
			byte dir= attributes.strand;
			if ((dir== 2&& sense)|| (dir== 1&& !sense)) {
				++nrMappingsWrongStrand;
				return 0;
			}
		}
		byte antiflag= (byte) ((flag==1)?2:1);
		int mapCtr= 0;
		
		
		// add first/single read
		if (paired) { /* PAIRED END */

			//int mappedIDsBefore= mapReadOrPairIDs.size();
			Vector<BEDobject2>[] vv= mapEndsOfPairs.get(ID);
			Vector<BEDobject2> v= null;
			if (vv!= null)
				v= vv[antiflag- 1];
			for (int i = 0; v!= null
				&&i < v.size(); i++) {
				
				BEDobject2 dobject2= v.elementAt(i);
				if (dobject.getStrand()== dobject2.getStrand()
						// 20101222: check also that the leftmost (in genomic direction) is sense (in genomic direction) 
						|| (dobject.getStart()< dobject2.getStart()&& dobject.getStrand()!= Transcript.STRAND_POS)
						|| (dobject2.getStart()< dobject.getStart()&& dobject2.getStrand()!= Transcript.STRAND_POS)) {
					++nrPairsWrongOrientation;
					continue;
				}
				
				Edge target2= getEdge2(dobject2);
				if (target2== null)
					continue;

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
					++nrPairsNoTxEvidence;
					continue;	// no tx evidence
				}
				se.incrReadNr();
				++mapCtr;
				
				
//						if (gene.getGeneID().equals("chr12:58213712-58240747C")) 
//							try {
//								testWriter.write(dobject.toString()+ "\n");
//								testWriter.write(dobject2.toString()+ "\n");
//							} catch (Exception e) {
//								e.printStackTrace();
//							}
					

				mapReadOrPairIDs.add(dobject.getName());
				mapReadOrPairIDs.add(dobject2.getName()); // !!! must have same id as bed object

//					if (outputMapped) {
//						writeMappedRead(dobject);
//						writeMappedRead(dobject2);
//					}
			}
			
			//Vector<DirectedRegion[]>[] vv= null;
			if (vv== null) {
				vv= new Vector[] {new Vector<DirectedRegion>(5,5),
						new Vector<DirectedRegion>(5,5)};
				mapEndsOfPairs.put(ID, vv);
			} 
			vv[flag- 1].add(dobject);
			
			return mapCtr; 	// (mapReadOrPairIDs.size()> mappedIDsBefore);
			
			
		} else { /* SINGLE READS */
			
			//incrementProfile(g, target, dobject, sense);
			mapCtr= 1;
			if (!mapReadOrPairIDs.add(dobject.getName()))
				++nrLocusMultimaps;
//				if (outputMapped)
//					writeMappedRead(dobject);
			return mapCtr;
		}
		
	}

	void map(BEDobject2[] beds) {

			// init
			nrMappingsReadsOrPairs= 0;
			mapReadOrPairIDs= new HashSet<CharSequence>(beds== null?0:beds.length, 1f);
			if (descriptor2.isPaired())
				mapEndsOfPairs = new HashMap<CharSequence, Vector<BEDobject2>[]>(beds== null?0:beds.length/ 2, 1f);
			attributes= descriptor2.createAttributes();
		
			if (beds== null|| beds.length== 0)
				return;
			
			boolean output= false;
			long t0 = System.currentTimeMillis();
			
			// map read pairs
			for (int j = 0; beds!= null&& j< beds.length; ++j) {

				int xyxx= mapRead(beds[j], false);
				//nrMappingsReadsOrPairs+= xxx;
			}
			nrMappingsReadsOrPairs+= mapReadOrPairIDs.size()/ 2;

			// map single reads, mapping on single edges and inc count
			if (false) {
				for (int j = 0; beds!= null&& j< beds.length; ++j) {
					int xxx= mapRead(beds[j], true);
					nrMappingsForced+= xxx;
				}
			}
			
			// increase pair number
/*				int gstart= gene.get5PrimeEdge(), gend= gene.get3PrimeEdge();
				if (gene.isReverse()) {
					int h= Math.abs(gstart);
					gstart= Math.abs(gend);
					gend= h;
				}
				for (int j = 0; beds!= null&& j< beds.length; ++j) {
					if (beds[j].getStart()>= gstart&& beds[j].getEnd()<= gend)
						continue;
					ByteArrayCharSequence id= beds[j].getName();
					if (mapReadOrPairIDs.contains(id))
						continue;
					byte flag= getFlag(beds[j]);  
					char antiflag= ((flag==1)?'2':'1');
					id.setCharAt(id.length()- 1, antiflag);
					if (mapReadOrPairIDs.contains(id)) 
						++nrMappingsReadsOrPairs;
				}
*/				
			
			int cntNomapped= 0;
			//HashSet<CharSequence> tmpMap= new HashSet<CharSequence>();
			for (int i = 0; i < beds.length; i++) {
				if (mapReadOrPairIDs.contains(beds[i].getName())) 
					continue;
//				if (outputNotmapped)
//					writeNotmappedRead(beds[i]);
//				BEDobject.addRecycleObj(beds[i]);

				//tmpMap.add(beds[i].getName());
				++cntNomapped;
				
				beds[i]= null;
			}

			
/*			Iterator<CharSequence> iter= tmpMap.iterator();
			while(iter.hasNext()) {
				CharSequence cseq= iter.next();
				if (mapReadOrPairIDs.contains(cseq)) {
					System.err.println(cseq.hashCode());
					System.currentTimeMillis();
				}
			}*/
			
			if (output) {
				System.err.println(", mapping "+((System.currentTimeMillis()- t0)/ 1000)+" sec.");
				System.err.flush();
			}
			
			long notMapped= 0;
			if (beds!= null)  {	// && mapOnly
				if (descriptor2.isPaired()) {
					//assert(mappedReads== myGraph.mappedIDSet.size()); // no, multimaps
					//mappedReads*= 2;
				}
				notMapped= beds.length- nrMappingsReadsOrPairs;
				setMappedReads(nrMappingsReadsOrPairs); 
				nrReadsMapped+= nrMappingsReadsOrPairs;
				nrReadsLoci+= beds.length;
			}
	//		if (mapOnly)
	//			return;
			if (notMapped> 0) { 
				
				if (Constants.verboseLevel>= Constants.VERBOSE_DEBUG)
					System.err.println("[WARNING] locus "+gene.getReferenceTranscript().getTranscriptID()
						+" couldnt map "+notMapped+" of "+beds.length+" mappings.");
			}

		}

}
