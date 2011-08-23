package fbi.genome.sequencing.rnaseq.reconstruction;

import fbi.commons.ByteArrayCharSequence;
import fbi.commons.Execute;
import fbi.commons.Log;
import fbi.commons.StringUtils;
import fbi.commons.file.FileHelper;
import fbi.commons.system.SystemInspector;
import fbi.commons.thread.SyncIOHandler2;
import fbi.commons.thread.ThreadedQWriter;
import fbi.genome.io.bed.BEDiteratorArray;
import fbi.genome.io.bed.BEDwrapper;
import fbi.genome.io.bed.BufferedBEDiterator;
import fbi.genome.io.gff.GFFReader;
import fbi.genome.io.rna.ReadDescriptor;
import fbi.genome.io.rna.UniversalReadDescriptor;
import fbi.genome.model.*;
import fbi.genome.model.bed.BEDobject;
import fbi.genome.model.bed.BEDobject2;
import fbi.genome.model.commons.MyFile;
import fbi.genome.model.constants.Constants;
import fbi.genome.model.gff.GFFObject;
import fbi.genome.model.splicegraph.Edge;
import fbi.genome.model.splicegraph.Graph;
import fbi.genome.model.splicegraph.Node;
import fbi.genome.model.splicegraph.SuperEdge;
import lpsolve.LpSolve;
import lpsolve.VersionInfo;

import java.io.*;
import java.lang.reflect.Method;
import java.text.SimpleDateFormat;
import java.util.*;
import java.util.concurrent.atomic.AtomicLong;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;
import java.util.zip.ZipOutputStream;


public class FluxCapacitor implements ReadStatCalculator {

	Object lock= new Object();

	ReadDescriptor descriptor;
	UniversalReadDescriptor descriptor2;
	
	class GTFreaderThread extends Thread {
		
		public GTFreaderThread() {
			super("GTF_reader");
		}
		
		@Override
		public void run() {
			try {
				getGTFreader().read();
			} catch (Exception e) {				
				e.printStackTrace();
			}		
		}
	}
	
	class LocusSolver extends Thread {
		Gene gene= null;
		ASEvent[] events= null;
		BEDobject2[] beds= null;
		boolean decompose= false;
		Thread threadBefore= null;
		int nrMappingsReadsOrPairs;
		HashSet<CharSequence> mapReadOrPairIDs;
		HashMap<CharSequence, Vector<BEDobject2>[]> mapEndsOfPairs;
		long[] sigall= null;

		private float invariantTestObsSplitFreq= 0, invariantTestPredSplitFreq= 0;
		
		public LocusSolver(Gene newGene, BEDobject2[] newBeds, boolean decompose) {
			//super(newGene.getGeneID());
			
			this.gene= newGene; 
			this.beds= newBeds;
			this.decompose= decompose;
			
			nrMappingsReadsOrPairs= 0;
			mapReadOrPairIDs= new HashSet<CharSequence>(newBeds== null?0:newBeds.length, 1f);
			if (pairedEnd)
				mapEndsOfPairs = new HashMap<CharSequence, Vector<BEDobject2>[]>(newBeds== null?0:newBeds.length/ 2, 1f);
		}
		
		public void run() {

			try {
				//int mapped= mapReadOrPairIDs.size();
				if (decompose) {
					Graph myGraph= getGraph(this.gene);
					sigall= myGraph.encodeTset(myGraph.trpts);
					map(myGraph, this.gene, this.beds);
					
					GraphLPsolver mySolver= null;
					if (mapReadOrPairIDs.size()> 0&& this.gene.getTranscriptCount()> 1) {	// OPTIMIZE			
						mySolver= getSolver(myGraph, nrMappingsReadsOrPairs); // not: getMappedReadcount()
						mySolver.run();
					}
	//				if (this.gene.getTranscriptCount()== 1)
	//					System.currentTimeMillis();
					outputGFF(myGraph, events, mySolver);
				} else {
					// map all reads
					if (this.gene.getTranscriptCount()== 1) {
						++nrSingleTranscriptLearn;
						learn(this.gene.getTranscripts()[0], beds);	
					}
				}
			} catch (Throwable e) {
				System.err.println("\n[ERROR] in locus "+ gene.getGeneID());
				e.printStackTrace();
				System.err.print("\n\tcontinuing ");
				System.err.flush();
			}
			
			// cleanup
//			for (int i = 0; (beds!= null)&& i < beds.length; i++) 
//				if (beds[i]!= null)
//					BEDobject.addRecycleObj(beds[i]);
			
			beds= null;
			gene= null;
			// makes it terribly slow
			//System.gc();
			
//			synchronized(FluxCapacitor.this.threadPool) {
				FluxCapacitor.this.threadPool.remove(this);
//			}
		}

		Graph getGraph(Gene gene) {
				boolean output= false;
				
				// construct graph
			long t0= System.currentTimeMillis();
			
			Graph myGraph= new Graph(gene);
			myGraph.createDefaultCoordComparator(readLenMin);
			myGraph.constructGraph();
			
			//if (outputLocus) {
				myGraph.setRetrieveDSEvents(true);
				myGraph.setRetrieveVSEvents(true);
				events= myGraph.getEvents(eventDim);
			//}
			
			myGraph.getNodesInGenomicOrder();	// important ??!
			myGraph.transformToFragmentGraph();
			if (output) {
				System.err.print(", transform "+((System.currentTimeMillis()- t0)/ 1000)+" sec, ");
				System.err.flush();
			}
			
			int nrSJ= myGraph.addEJ(readLenMin);
			if (output) {
				System.err.print(", EJ "+((System.currentTimeMillis()- t0)/ 1000)+" sec, ");
				System.err.flush();
			}
		
			if (pairedEnd) {
				int nrPE= addPE(myGraph, insertMinMax, readLenMin);
				if (output) {
					System.err.print(", PE "+((System.currentTimeMillis()- t0)/ 1000)+" sec, ");
					System.err.flush();
				}
			}
		
			return myGraph;
		}

		Graph map(Graph myGraph, Gene g, BEDobject2[] beds) {
			
			
			if (beds== null|| beds.length== 0)
				return null;
			
			boolean output= false;
			long t0 = System.currentTimeMillis();
			
//			if (g.getTranscripts()[0].getTranscriptID().equals("ENST00000407072"))
//				System.currentTimeMillis();
			nrMappingsReadsOrPairs= 0;
			for (int j = 0; beds!= null&& j< beds.length; ++j) {

				int xxx= mapRead(myGraph, beds[j]);
				nrMappingsReadsOrPairs+= xxx;
			}
			
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
			
			int notMapped= 0;
			if (beds!= null)  {	// && mapOnly
				if (pairedEnd) {
					//assert(mappedReads== myGraph.mappedIDSet.size()); // no, multimaps
					//mappedReads*= 2;
				}
				notMapped= beds.length- nrMappingsReadsOrPairs;
				myGraph.setMappedReads(nrMappingsReadsOrPairs); 
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

			return myGraph;
		}
				
		private int nrLocusMultimaps= 0;
		/**
		 * add a SINGLE read
		 * @param regs
		 */
		int mapRead(Graph g, BEDobject2 dobject) {
			// find the edge(s) where the regions align
			// if these do not form a continous chain, create a new edge
			
//			GFFObject[] gtfs= GFFObject.fromBed(dobject);	// TODO kill intermediary GTFs !!!
//			DirectedRegion[] regs= new DirectedRegion[gtfs.length];
//			for (int i = 0; i < regs.length; i++) 
//				regs[i]= new DirectedRegion(gtfs[i]);
			
			byte flag= getFlag(dobject);  	
			CharSequence ID= getID(dobject); 	

			Edge target= g.getEdge(dobject);
			
			if ((target== null)|| ((!pairedEnd)&&  (!(target instanceof SuperEdge))
					&& target.length()< readLenMin))
				return 0;

			byte refStrand= g.trpts[0].getStrand();
			boolean sense= dobject.getStrand()== refStrand;
			byte antiflag= (byte) ((flag==1)?2:1);
			int mapCtr= 0;
			
			
			// add first/single read
			if (pairedEnd) { /* PAIRED END */

				//int mappedIDsBefore= mapReadOrPairIDs.size();
				Vector<BEDobject2>[] vv= mapEndsOfPairs.get(ID);
				Vector<BEDobject2> v= null;
				if (vv!= null)
					v= vv[antiflag- 1];
				for (int i = 0; v!= null
					&&i < v.size(); i++) {
					
					BEDobject2 dobject2= v.elementAt(i);
					Edge target2= g.getEdge(dobject2);
					if (target2== null)
						continue;

					// check whether they map within isize constraints
					int j = 0;
					for (; target.getSuperEdges()!= null&& j < target.getSuperEdges().size(); j++) {
						if (!target.getSuperEdges().elementAt(j).isPend())
							continue;
						Edge[] ee= target.getSuperEdges().elementAt(j).getEdges();
						int k = 0;
						//for (; k < ee.length&& ee[k]!= target2; k++);	// TODO binarySearch
						// for the case target== target2, better check that there is no other edge
						for (; k < ee.length; k++) {
							if (ee[k]!= target&& ee[k]!= target2)
								break;
						}
							
						if (k== ee.length)
							break;	// common superedge found@ j
					}
					SuperEdge se= null;
					if (target.getSuperEdges()!= null&& j< target.getSuperEdges().size()) 
						se= target.getSuperEdges().elementAt(j);
					else {
						continue;	// not possible paired-end
					}

					se.incrReadNr();
					++mapCtr;
					
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

				if (sense|| (strand!= FluxCapacitorConstants.STRAND_ENABLED)) {
					target.incrReadNr();
					mapCtr= 1;
				} else if (strand!= FluxCapacitorConstants.STRAND_SPECIFIC) {
					target.incrRevReadNr();
					mapCtr= 1;
				} else {
					++nrMappingsWrongStrand;
					mapCtr= 0;
				}
				
				
				
				if (!mapReadOrPairIDs.add(dobject.getName()))
					++nrLocusMultimaps;
//				if (outputMapped)
//					writeMappedRead(dobject);
				return mapCtr;
			}
			
		}
		
		private CharSequence getID(BEDobject2 dobject) {
			CharSequence id= null; 	
			if (pairedEnd) {
				id= descriptor.getUniqueDescriptor(dobject.getName());
			} else
				id= dobject.getName();

			return id;
		}
		
		private byte getFlag(BEDobject2 dobject) {
			byte flag= 0;  	
			if (pairedEnd) 
				flag= (byte) (descriptor.getPairedEndInformation(dobject.getName()));
			return flag;
		}
		
		public int getMappedReadcount() {
			if (pairedEnd)
				return mapReadOrPairIDs.size()/ 2;
			return mapReadOrPairIDs.size();
		}
		

		private boolean writeMappedRead(BEDobject o) {
			if (getFileMappedReads()== null)
				return false;
			try {
				getWriterMappedReads().write(o.toString()+ "\n");
				return true;
			} catch (IOException e) {			
				e.printStackTrace();
				return false;
			}
		}
		
		private boolean writeNotmappedRead(BEDobject o) {
			if (getFileNotMappedReads()== null)
				return false;
			try {
				getWriterNotmappedReads().write(o.toString()+ "\n");
				return true;
			} catch (IOException e) {
				e.printStackTrace();
				return false;
			}
		}
		
		
		private void outputGFF(Graph g, ASEvent[] events, GraphLPsolver solver) {
			++nrLoci;
			if (solver!= null) 
				++nrLociExp;
			double perM= nrReadsAll/ 1000000d;
			// deprecated
			boolean unsolvedSystem= false;	
			double valOF= solver== null?0: solver.getValObjFunc();
			if (valOF> FluxCapacitorConstants.BIG) { 
				++nrUnsolved;
				unsolvedSystem= true;
			}
			String pv= getAttributeOF(valOF, solver, getMappedReadcount());

			// prebuild rpkm hash
			HashMap<String, Float> rpkmMap= null;
			if (outputBalanced) {
				rpkmMap= new HashMap<String, Float>(g.trpts.length, 1f);
				for (int i = 0; i < g.trpts.length; i++) {
					float val= 0f, rpkm= 0f;
					if (solver== null) {
						// no reads
						val= nrMappingsReadsOrPairs;
					} else {
//						if (solver.getTrptExprHash().get(g.trpts[i].getTranscriptID())!= 0)
//							System.currentTimeMillis();
						val= (float) (solver.getNFactor()* solver.getTrptExprHash().get(g.trpts[i].getTranscriptID()));
						if (val< 1- costBounds[0]) // 1- 0.95
							val= 0;
						try {
							assert(g.trpts.length> 1|| val== nrMappingsReadsOrPairs);
						} catch (AssertionError e) {
							System.err.println(val+ " x "+ nrMappingsReadsOrPairs);
							solver.getNFactor();
						}
					}
					if (pairedEnd)
						val*= 2;	// count both ends for paired-end reads
					if (val> 0&& !(outputObs|| outputPred))
						++nrTxExp;
					rpkm= calcRPKM(val, g.trpts[i].getExonicLength());
					rpkmMap.put(g.trpts[i].getTranscriptID(), rpkm);
				}
			}
			
			
			// reproduce original
			boolean foundTranscripts= false, foundExons= false;
			if (getGTFreader().isKeepOriginalLines()&& origLines!= null) {
				foundExons= true;
				for (int i = 0; i < origLines.size(); i++) {
					String s= origLines.elementAt(i);
					String feat= GFFObject.getField(3, s);
					String tid= GFFObject.getTranscriptID(s);
					int tx= 0;
					if ((feat.equals(feat.equals(Transcript.GFF_FEATURE_TRANSCRIPT))
							|| feat.equals(Exon.GFF_FEATURE_EXON))
							&& (outputObs|| outputPred))
						for (tx = 0; tx < g.trpts.length; tx++) 
							if (g.trpts[tx].getTranscriptID().equals(tid))
								break;
					if (tx>= g.trpts.length) {
						System.err.println("\nTranscript "+ tid+ " not found in: ");
						for (int j = 0; j < g.trpts.length; j++) 
							System.err.println("\t"+ g.trpts[j].getTranscriptID());
						System.err.println();
					}
					
					if (feat.equals(Transcript.GFF_FEATURE_TRANSCRIPT)&& outputTranscript) {
						foundTranscripts= true;
						StringBuilder sb= new StringBuilder(s);
						int x= sb.length();	// trim
						while(Character.isWhitespace(sb.charAt(--x)))
							sb.delete(x, x+1);
						if (sb.charAt(x)!= ';')
							sb.append("; ");
						else
							sb.append(Constants.SPACE);
						
						if ((outputObs|| outputPred)&& tx< g.trpts.length)
							getGTF(sb, g.trpts[tx], solver, g, perM, pv, true);
						else if (outputBalanced) {
							sb.append(FluxCapacitorConstants.GTF_ATTRIBUTE_TOKEN_RPKM);
							sb.append(Constants.SPACE);
							if (rpkmMap.containsKey(tid))
								sb.append(String.format("%1$f", rpkmMap.get(tid).floatValue()));	// rgasp parser does not like scientific notation
							else
								sb.append(Constants.NULL);
							sb.append(";\n");
						}
						
						try {
							getWriter().write(sb.toString());
						} catch (IOException e) {
							e.printStackTrace();
						} 
					} else if (feat.equals(Exon.GFF_FEATURE_EXON)&& outputExon) {
						
						StringBuilder sb= new StringBuilder(s); 
						int x= sb.length();
						while(Character.isWhitespace(sb.charAt(--x)))
							sb.delete(x, x+1);
						if (sb.charAt(x)!= ';')
							sb.append("; ");
						else
							sb.append(' ');
						
												
						if ((outputObs|| outputPred)&& tx< g.trpts.length) {
							int start= Integer.parseInt(GFFObject.getField(4, s));
							int end= Integer.parseInt(GFFObject.getField(5, s));
							int j = 0;
							for (; j < g.trpts[x].getExons().length; j++) {
								int begin= Math.abs(g.trpts[x].getExons()[j].getStart()),
									ende= Math.abs(g.trpts[x].getExons()[j].getEnd());
								if (begin> ende) {
									int h= begin;
									begin= ende;
									ende= h;
								}
								if (begin== start&& ende== end)
									break;
							}
							getGTF(sb, g.trpts[x].getExons()[j], g.trpts[i], g, solver, unsolvedSystem, perM, pv, true);
						} else if (outputBalanced) {
							sb.append(FluxCapacitorConstants.GTF_ATTRIBUTE_TOKEN_RPKM);
							sb.append(Constants.SPACE);
							if (rpkmMap.containsKey(tid))
								sb.append(String.format("%1$f", rpkmMap.get(tid).floatValue()));	// rgasp parser does not like scientific notation
							else
								sb.append(Constants.NULL);
							sb.append(";\n");
						}
						
						
						try {
							getWriter().write(sb.toString());
						} catch (IOException e) {
							e.printStackTrace();
						}
					} else {
						try {
							getWriter().write(s+ System.getProperty("line.separator"));
						} catch (IOException e) {
							e.printStackTrace();
						}
					}
				}
			}
			
				
			
			StringBuilder sb= new StringBuilder();
			// LOCUS TODO genes 
			if (outputGene) {
				if (outputObs|| outputPred) {
					getGTF(sb, g.trpts[0].getGene(), g, solver, perM, pv);	
					try {assert(testInvariant(invariantTestObsSplitFreq, 
							pairedEnd?nrMappingsReadsOrPairs*2:nrMappingsReadsOrPairs, 0.05));}	// min: 0.01
					catch (AssertionError e) {
						if (Constants.verboseLevel> Constants.VERBOSE_NORMAL)
							System.err.println("[ASSERTION] "+getClass().getName()+".outputGFF():\n\tinvariantTestObsSplitFreq= "
									+ invariantTestObsSplitFreq+ ", nrMappingsReadsOrPairs= "+ (pairedEnd?nrMappingsReadsOrPairs*2:nrMappingsReadsOrPairs)
									+ "\n\tlocus: "+ g.trpts[0].getTranscriptID());
					};
					try {assert(testInvariant(invariantTestPredSplitFreq, 
							pairedEnd?nrMappingsReadsOrPairs*2:nrMappingsReadsOrPairs, 0.1));}
					catch (AssertionError e) {
						if (Constants.verboseLevel> Constants.VERBOSE_NORMAL)
							System.err.println("[ASSERTION] "+getClass().getName()+".outputGFF():\n\tinvariantTestPredSplitFreq= "
									+ invariantTestPredSplitFreq+ ", nrMappingsReadsOrPairs= "+ (pairedEnd?nrMappingsReadsOrPairs*2:nrMappingsReadsOrPairs)
									+ "\n\tlocus: "+ g.trpts[0].getTranscriptID());
					};
					
				} else if (outputBalanced) {
				}
			}

			
			// TRANSCRIPTS
			if (outputTranscript|| outputExon|| outputSJunction) {
				float invariantObsAllTx= 0, invariantPredAllTx= 0,
					invariantObsAllEx= 0, invariantPredAllEx= 0;
				for (int i = 0; i < g.trpts.length; i++) {
					++nrTx;
//					float invariantObsTx= invariantTestObsSplitFreq,
//					invariantPredTx= invariantTestPredSplitFreq;
					float invariantObsTx= 0, invariantPredTx= 0;
					if (outputTranscript&& !foundTranscripts) {
						if (outputObs|| outputPred) {
							getGTF(sb, g.trpts[i], solver, g, perM, pv, false);	// writer.write
							invariantObsAllTx+= invariantTestObsSplitFreq; //invariantObsTx;
							invariantPredAllTx+= invariantTestPredSplitFreq; // invariantPredTx;
							invariantObsTx= invariantTestObsSplitFreq;
							invariantPredTx= invariantTestPredSplitFreq;
							if (invariantPredTx> 0)
								++nrTxExp;

						} else if (outputBalanced) {
							GFFObject obj= GFFObject.createGFFObject(g.trpts[i]);
							sb.append(obj.toString());
							int x= sb.length();
							while(Character.isWhitespace(sb.charAt(--x)))
								sb.delete(x, x+1);
							if (sb.charAt(x)!= ';')
								sb.append("; ");
							else
								sb.append(Constants.SPACE);
							
							sb.append(FluxCapacitorConstants.GTF_ATTRIBUTE_TOKEN_RPKM);
							sb.append(Constants.SPACE);							
							//sb.append(rpkmMap.get(g.trpts[i].getTranscriptID()));
							sb.append(String.format("%1$f", rpkmMap.get(g.trpts[i].getTranscriptID()).floatValue()));	// rgasp parser does not like scientific notation
							sb.append(";\n");
						}
					}
					// EXONS
					float invariantObsEx= 0, invariantPredEx= 0;
					if (outputExon&& !foundExons) {
						Exon[] exons=  g.trpts[i].getExons();
						for (int j = 0; j < exons.length; j++) {
							getGTF(sb, exons[j], g.trpts[i], g, solver, unsolvedSystem, perM, pv, false);
							invariantObsEx+= invariantTestObsSplitFreq;
							invariantPredEx+= invariantTestPredSplitFreq;
						}
					}
					
					// SJ
					if (outputSJunction) {
						Vector<Vector<Edge>> eeV= new Vector<Vector<Edge>>(5,5);
						eeV.add(new Vector<Edge>());
						g.getRPK(g.trpts[i], pairedEnd, Graph.ETYPE_SJ, eeV);
						long[][] sig= new long[][]{g.encodeTset(g.trpts[i])};
						for (int j = 0; j < eeV.elementAt(0).size(); j++) { 
							getGTF(sb, eeV.elementAt(0).elementAt(j), sig, g, solver, perM);
							invariantObsEx+= invariantTestObsSplitFreq;
							invariantPredEx+= invariantTestPredSplitFreq;
						}
					}
					invariantObsAllEx+= invariantObsEx;
					invariantPredAllEx+= invariantPredEx;
					
					if (outputExon&& outputSJunction&& outputTranscript) {
						try {assert(testInvariant(invariantObsEx, invariantObsTx, 0.05));}	// min: 0.02
						catch (AssertionError e) {
							if (Constants.verboseLevel> Constants.VERBOSE_NORMAL)
								System.err.println("[ASSERTION] "+getClass().getName()+".outputGFF():\n\tinvariantObsEx= "
										+ invariantObsEx+ ", invariantObsTx= "+ invariantObsTx
										+ "\n\tlocus: "+ g.trpts[0].getTranscriptID());
						};
						try {assert(testInvariant(invariantPredEx, invariantPredTx, 0.1));}
						catch (AssertionError e) {
							if (Constants.verboseLevel> Constants.VERBOSE_NORMAL)
								System.err.println("[ASSERTION] "+getClass().getName()+".outputGFF():\n\tinvariantPredEx= "
										+ invariantPredEx+ ", invariantPredTx= "+ invariantPredTx
										+ "\n\tlocus: "+ g.trpts[0].getTranscriptID());
						};
					}
				}
				if (outputTranscript) {
					try {assert(testInvariant(invariantObsAllTx, 
							pairedEnd?nrMappingsReadsOrPairs*2:nrMappingsReadsOrPairs, 0.05));}	// min: 0.01
					catch (AssertionError e) {
						if (Constants.verboseLevel> Constants.VERBOSE_NORMAL)
							System.err.println("[ASSERTION] "+getClass().getName()+".outputGFF():\n\tinvariantObsAllTx= "
									+ invariantObsAllTx+ ", nrMappingsReadsOrPairs= "+ (pairedEnd?nrMappingsReadsOrPairs*2:nrMappingsReadsOrPairs)
									+ "\n\tlocus: "+ g.trpts[0].getTranscriptID());
					};
					try {assert(testInvariant(invariantPredAllTx, 
							pairedEnd?nrMappingsReadsOrPairs*2:nrMappingsReadsOrPairs, 0.1));}
					catch (AssertionError e) {
						if (Constants.verboseLevel> Constants.VERBOSE_NORMAL)
							System.err.println("[ASSERTION] "+getClass().getName()+".outputGFF():\n\tinvariantPredAllTx= "
									+ invariantPredAllTx+ ", nrMappingsReadsOrPairs= "+ (pairedEnd?nrMappingsReadsOrPairs*2:nrMappingsReadsOrPairs)
									+ "\n\tlocus: "+ g.trpts[0].getTranscriptID());
					};
				}
				if (outputExon&& outputSJunction) {
					try {assert(testInvariant(invariantObsAllEx, 
							pairedEnd?nrMappingsReadsOrPairs*2:nrMappingsReadsOrPairs, 0.05));}	// min: 0.02
					catch (AssertionError e) {
						if (Constants.verboseLevel> Constants.VERBOSE_NORMAL)
							System.err.println("[ASSERTION] "+getClass().getName()+".outputGFF():\n\tinvariantObsAllEx= "
									+ invariantObsAllEx+ ", nrMappingsReadsOrPairs= "+ (pairedEnd?nrMappingsReadsOrPairs*2:nrMappingsReadsOrPairs)
									+ "\n\tlocus: "+ g.trpts[0].getTranscriptID());
					};
					try {assert(testInvariant(invariantPredAllEx, 
							pairedEnd?nrMappingsReadsOrPairs*2:nrMappingsReadsOrPairs, 0.1));}
					catch (AssertionError e) {
						if (Constants.verboseLevel> Constants.VERBOSE_NORMAL)
							System.err.println("[ASSERTION] "+getClass().getName()+".outputGFF():\n\tinvariantPredAllEx= "
									+ invariantPredAllEx+ ", nrMappingsReadsOrPairs= "+ (pairedEnd?nrMappingsReadsOrPairs*2:nrMappingsReadsOrPairs)
									+ "\n\tlocus: "+ g.trpts[0].getTranscriptID());
					};
				}
			}
			
			// EVENTS
			if (outputEvent) {
				HashMap<Object,Double> tExpMap= null;
				if (solver!= null) {
					tExpMap= solver.getTrptExprHash();
					Object[] keys= tExpMap.keySet().toArray();	
					for (int i = 0; i < keys.length; i++) {
						if (!(keys[i] instanceof String))
							continue;
						if (tExpMap.get(keys[i])<0)
							tExpMap.put((String) keys[i], 0d);	// TODO ugly
					}
				}
				for (int i = 0; events!= null&& i < events.length; i++) {
					if (outputObs|| outputPred)
						getGTF(sb, events[i], g, solver, unsolvedSystem, perM, pv, tExpMap);
					else
						++nrEvents;
					if (outputBalanced) {
						sb.append(events[i].toStringGTF());
						sb.append(Constants.SPACE);
						sb.append("\"");
						boolean allPos= true;
						for (int j = 0; j < events[i].getTranscripts().length; j++) {
							float sum= 0;
							for (int k = 0; k < events[i].getTranscripts()[j].length; k++) 
								sum+= rpkmMap.get(events[i].getTranscripts()[j][k].getTranscriptID());
							sb.append(sum);
							sb.append(",");
							allPos&= (sum> 0);
						}
						if (allPos&& !(outputObs|| outputPred))
							++nrEventsExp;
						sb.replace(sb.length()- 1, sb.length(), "\";\n");
					}

				}
			}

			// FRAGMENTS and XJUNCTIONS
			if (false&& solver!= null) {
				ArrayList<Edge> cc= new ArrayList<Edge>();
				if (solver!= null) {
					Iterator<Object> iter= solver.getConstraintHash().keySet().iterator();
					while (iter.hasNext()) {
						Object o= iter.next();
						if (o instanceof Edge)
							cc.add((Edge) o);
					}
				}
				Collections.sort(cc, Edge.getDefaultPositionComparator());

				Iterator<Edge> iter= cc.iterator();
				while (iter.hasNext()) {
					Edge e= iter.next();
					// no INTRONS
					if ((!(e instanceof SuperEdge))&& (!e.isExonic()))
						continue;
					getGTF(sb, e, new long[][]{e.getTranscripts()}, g, solver, perM);
				}
			}
				
			try {
				write(sb);
			} catch (Exception e) {
				e.printStackTrace();
			}
			
		}

		private boolean testInvariant(double invariant, double reference, double stringency) {
			double delta= Math.abs(reference== 0?invariant: (invariant- reference)/ reference);
			if (delta> stringency) {
				if (invariant<= 0)
					return true; 	// catch 0-predictions
				return false;
			}
			
			return true;
		}

		public void setThreadBefore(Thread threadBefore) {
			this.threadBefore = threadBefore;
		}

		private void learn_old(Transcript tx, BEDobject2[] beds) {
				
				if (beds== null|| beds.length== 0)
					return;
				
				int elen= tx.getExonicLength();	// TODO this is the effective length
				if (elen< readLenMin)
					return;	// TODO discards reads
				
				int plen= elen- (readLenMin- 1);
				TProfile profile= null;
				if (!uniform) {
					// not any more, now expr x length unique, and IDs
					// profile= func.getProfile(plen, (strand== STRAND_ENABLED), pairedEnd);
					profile= new TProfile(tx.getGene().getGeneID(), plen, (strand== FluxCapacitorConstants.STRAND_ENABLED), pairedEnd);
					func.profiles.add(profile);
				}
				HashMap<CharSequence, BEDobject2[][]> mapPends= new HashMap<CharSequence, BEDobject2[][]>();
				nrReadsSingleLoci+= beds.length;
				// possible read pairs
				// slow
/*				for (int i = 0; i < beds.length; i++) {
					for (int j = i+1; j < beds.length; j++) {
						if (((beds[i].getStrand()== 0&& beds[j].getStrand()== 0)
								|| (beds[i].getStrand()== -beds[j].getStrand()))
								&& FMRD.matesByName(beds[i].getName(), beds[j].getName()))
							++nrReadsSingleLociPotentialPairs;
					}
				}
*/				
				for (int i = 0; i < beds.length; i++) {
					
					if (beds[i]== null)
						continue;
					
					int len= beds[i].getLength();
					if (readLenMin< 0|| len< readLenMin)
						readLenMin= len;
					if (len> readLenMax)	// readLenMax< 0|| 
						readLenMax= len;
					int[] ePos= getExonicPos(tx, beds[i], elen);
					
					if (ePos== null) {
						//ePos= getExonicPos(tx, beds[i]);
						++nrReadsSingleLociNoAnnotation;
						continue; 	// doesnt align to the transcript
					}

					if (strand== FluxCapacitorConstants.STRAND_SPECIFIC&& beds[i].getStrand()!= tx.getStrand()) {
						++nrMappingsWrongStrand;
						continue;
					}

					assert(ePos[1]> ePos[0]);
					if(ePos[0]< 0|| ePos[1]>= elen) {
						//ePos= getExonicPos(tx, beds[i]);
						System.err.println("[ahuuuaaa] learn() does not report an exonic position for read\n\t"+ beds[i]);
						continue;
					}
					
					++nrReadsSingleLociMapped;
					if (pairedEnd) {

						byte flag=  (byte) (descriptor.getPairedEndInformation(beds[i].getName())- 1);	// (Fasta.getReadDir(beds[i].getName())- 1))
						CharSequence id= descriptor.getUniqueDescriptor(beds[i].getName());	//Fasta.getReadID(beds[i].getName())
						//CharSequence id= beds[i].getName();
						BEDobject2[][] oo= mapPends.get(id);
						if (oo== null) {
							oo= new BEDobject2[2][];
							mapPends.put(id, oo);
						} 
						if (oo[flag]== null) 
							oo[flag]= new BEDobject2[] {beds[i]};
						else {
							BEDobject2[] op= new BEDobject2[oo[flag].length+ 1];
							System.arraycopy(oo[flag], 0, op, 0, oo[flag].length);
							op[op.length- 1]= beds[i];
							oo[flag]= op;
						}
						
						for (int j = 0; j < oo.length; j++) {
							if (j==flag|| oo[j]== null)
								continue;
							for (int k = 0; k < oo[j].length; k++) {
								int[] ePos2= getExonicPos(tx, oo[j][k], elen);
								if (ePos[0]> ePos2[0]) {
									int[] h= ePos;
									ePos= ePos2;
									ePos2= h;
								}
								if (!uniform)
									profile.addReadPair(ePos[0], ePos2[0], readLenMin);
								addInsertSize(ePos2[1]- ePos[0]+ 1);
								++nrReadsSingleLociPairsMapped;	
//								if (nrReadsSingleLociPairsMapped> nrReadsSingleLoci)
//									System.currentTimeMillis();
		//						int delta= ePos2[0]- ePos[0];
							}
						}
							
					} else {
						byte dir= Constants.DIR_FORWARD;
						if (strand== FluxCapacitorConstants.STRAND_ENABLED&& beds[i].getStrand()!= tx.getStrand())
							dir= Constants.DIR_BACKWARD;
						
						if (!uniform)
							profile.addRead(ePos[0],readLenMin,dir);
					}
				}
				
				if ((!uniform)&& profile.getReads()== 0) {
					//profile.eliminate0regs();
					func.profiles.remove(profile);
				}
				//assert(func.getNrReadsInProfiles()== nrReadsLoci);
			}

		private GraphLPsolver getSolver(Graph g, int mappedReads) {
		
			GraphLPsolver solver= new GraphLPsolver(g, readLenMin, 
					pairedEnd?insertMinMax:null, mappedReads, 
					(strand== FluxCapacitorConstants.STRAND_ENABLED), 
					pairedEnd);
			if (outputLP)
				solver.setFileLPdir(getFileLP());
			solver.costModel= costModel;	// COSTS_LINEAR
			solver.setCostSplit(costSplit);
			solver.setFunc(func);
			solver.setReadLen(readLenMin);
			solver.setFlow(true);
			solver.costBounds= costBounds;
		
			return solver;
		}

		private String getGTF(StringBuilder sb, ASEvent event, Graph g, GraphLPsolver solver, boolean unsolvedSystem, 
					double perM, String pv, HashMap<Object,Double> tExpMap) {
				
		//		for (int i = 0; i < eeV.size(); i++) 
		//			eeV.elementAt(i).removeAllElements();
				Vector<Vector<Edge>> eeV= new Vector<Vector<Edge>>(5,5);
				while (eeV.size()< event.getDimension()) 
					eeV.add(new Vector<Edge>());
				
				g.getRPK(event, pairedEnd, Graph.ETYPE_AL, eeV);
				sb.append(event.toStringGTF());
		
				long[][] sig= new long[event.getDimension()][];
				for (int i = 0; i < sig.length; i++) { 
					sig[i]= g.createAllArray();		
					for (int j = 0; j < sig.length; j++) 
						sig[i]= j== i? sig[i]: 
							Graph.unite(sig[i], g.encodeTset(event.getTranscripts()[j]));
				}
				
				double splitReads= getGTFappend(sb, g, solver, eeV, perM, sig);
				++nrEvents;
				if (splitReads> 0)
					++nrEventsExp;
		
				for (int i = 0; i < eeV.size(); i++) 
					eeV.elementAt(i).removeAllElements();
				
				return sb.toString();
			}

		private String getGTF(StringBuilder sb, Edge e, long[][] sig, Graph g, GraphLPsolver solver, 
					double perM) {
				
				Vector<Vector<Edge>> eeV= new Vector<Vector<Edge>>(5,5);
				eeV.add(new Vector<Edge>(1));
				eeV.elementAt(0).add(e);
				
				sb.append(g.trpts[0].getChromosome());
				sb.append("\t");
				Transcript[] tt= g.decodeTset(e.getTranscripts());
		//		if (tt[0].getTranscriptID().equals("ENST00000407980"))
		//			System.currentTimeMillis();
				sb.append(Transcript.getSource(tt));		
				sb.append("\t");
				if (e instanceof SuperEdge) {
					if (((SuperEdge) e).isPend())
						sb.append(FluxCapacitorConstants.GFF_FEATURE_PAIRED);
					else 
						sb.append(FluxCapacitorConstants.GFF_FEATURE_JUNCTION);
				} else
					sb.append(FluxCapacitorConstants.GFF_FEATURE_FRAGMENT);
				sb.append("\t");
				
				int[] frac= e.getFrac(tt[0], readLenMin);
				int start= Math.abs(tt[0].getGenomicPosition(frac[0]));
				int end= Math.abs(tt[0].getGenomicPosition(frac[1]));
				if (start>end) {
					int h= start;
					start= end;
					end= h;
				}
				
				sb.append(Integer.toString(start));
				sb.append("\t");
				sb.append(Integer.toString(end));
				sb.append("\t.\t");
				sb.append(GFFObject.getStrandSymbol(tt[0].getStrand()));
				sb.append("\t.\t");
				
				sb.append(GFFObject.TRANSCRIPT_ID_TAG+" \"");
				for (int j = 0; j < tt.length; j++) {
					sb.append(tt[j].getTranscriptID());
					if (j< tt.length-1)
						sb.append(",");
				}
				sb.append("\";");
				
				sb.append(" ");
				sb.append(GFFObject.GENE_ID_TAG);
				sb.append(" \"");
				sb.append(tt[0].getGene().getGeneID());
				sb.append("\";");
		
				sb.append(" ");
				sb.append("edges");
				sb.append(" \"");
				sb.append(e.toString());
				sb.append("\";");			
				
//				sb.append(" ");
//				sb.append(GTF_ATTRIBUTE_LENGTH);
//				sb.append(" ");	// \"
//				//sb.append(frac[1]- frac[0]+ 1);
//				sb.append(getLength(g, eeV.elementAt(0), e.getTranscripts(), false));
//				sb.append(";");	// \"			
				
				// TODO why not?
				getGTFappend(sb, g, solver, eeV, perM, sig);
	
				// here was the alternative:
/*				sb.append(" ");
				sb.append(GTF_ATTRIBUTE_TOKEN_OBSV);
				sb.append(" \"");
				sb.append(e.getReadNr());
				sb.append("\";");			
						
				sb.append(" ");
				sb.append(GTF_ATTRIBUTE_TOKEN_PRED);
				sb.append(" \"");
				for (int j = 0; j < tt.length; j++) {
					sb.append(solver.getNFactor()* solver.getTrptExprHash().get(tt[j].getTranscriptID())* 
							solver.getSuperProfileMap().get(tt[j].getTranscriptID()).getAreaFrac(
							GraphLPsolver.bounds2rel(frac, tt[j].getExonicLength()- readLenMin), readLenMin, 
							strandSpecific?TProfile.DIR_BOTH:TProfile.DIR_FORWARD));
					if (j< tt.length-1)
						sb.append(",");
				}
				sb.append("\";");			
		
				sb.append(" ");
				sb.append(GTF_ATTRIBUTE_TOKEN_PRED+"_total");
				sb.append(" \"");
				sb.append(solver.getReads(eeV.elementAt(0), BYTE_0, e.getTranscripts()));
				sb.append("\";");			
		
				
				// expectations
				sb.append(" ");
				sb.append(GTF_ATTRIBUTE_EXPECT);
				sb.append(" \"");
				for (int j = 0; j < tt.length; j++) {
					sb.append(solver.getSuperProfileMap().get(tt[j].getTranscriptID()).getAreaFrac(
							GraphLPsolver.bounds2rel(frac, tt[j].getExonicLength()- readLenMin), readLenMin, 
							strandSpecific?TProfile.DIR_BOTH:TProfile.DIR_FORWARD));
					if (j< tt.length-1)
						sb.append(",");
				}
				sb.append("\";");
				
				// profiles
				sb.append(" ");
				sb.append(GTF_ATTRIBUTE_PROFILE);
				sb.append(" \"");
				for (int j = 0; j < tt.length; j++) {
					Vector<TProfile> v= solver.getSuperProfileMap().get(tt[j].getTranscriptID()).getProfiles();
					for (int i = 0; i < v.size(); i++) {
						sb.append(v.elementAt(i).length());
						if (i< v.size()-1)
							sb.append(":");
					}
					if (j< tt.length-1)
						sb.append(",");
				}
				sb.append("\";");
				
				sb.append("\n");
*/				
				return sb.toString();
			}

		private String getGTF(StringBuilder sb, Exon exon, Transcript t, Graph g, GraphLPsolver solver, boolean unsolvedSystem, 
					double perM, String pv, boolean attributesOnly) {

				if (!attributesOnly) {
					GFFObject obj= GFFObject.createGTFObjects(exon, t)[0];
					sb.append(obj.toString());
				}
			
		//		if (eeV.size()< 1)
		//			eeV.add(new Vector<Edge>());
		//		else
		//			eeV.elementAt(0).removeAllElements();
				Vector<Vector<Edge>> eeV= new Vector<Vector<Edge>>(5,5);
				eeV.add(new Vector<Edge>());
				
				//if (g.readCount> 0) // get lengths
				g.getRPK(exon, t, pairedEnd, Graph.ETYPE_AL, eeV);
		
				//containerIntA1A1[0][0]= g.readCount> 0? getLength(eeV.elementAt(0), null, true, false): (exon.getLength()- readLen);
				long[][] sig= new long[][]{g.encodeTset(t)};
				getGTFappend(sb, g, solver, eeV, perM, sig);
				eeV.elementAt(0).removeAllElements(); 
				
				return sb.toString();
			}

		private String getGTF(StringBuilder sb, Gene gene, Graph g, GraphLPsolver solver, double perM, String pv) {
			
			//clearEdgeContainer(1);
			Vector<Vector<Edge>> eeV= new Vector<Vector<Edge>>(5,5);
			eeV.add(new Vector<Edge>());
			
			GFFObject obj= GFFObject.createGFFObject(gene);
			sb.append(obj.toString());
			//if (g.readCount> 0) // for getting lengths 
			g.getRPK(gene, pairedEnd, Graph.ETYPE_AL, eeV);
			
			
			//lenExon[0][0]= g.readCount> 0? getLength(eeV.elementAt(0), null): (t.getExonicLength()- readLen);
			//containerIntA1A1[0][0]= getLength(eeV.elementAt(0), null, true, false);
			//debug= true;
			long[][] sig= new long[][]{g.createAllArray()};
			getGTFappend(sb, g, solver, eeV, perM, sig);
			//debug= false;
			return sb.toString();
		}	
		
		private String getGTF(StringBuilder sb, Transcript t, GraphLPsolver solver, Graph g, double perM, String pv, boolean attributesOnly) {
				
				GFFObject obj= GFFObject.createGFFObject(t);
				sb.append(obj.toString());

				Vector<Vector<Edge>> eeV= new Vector<Vector<Edge>>(5,5);
				eeV.add(new Vector<Edge>());
				//if (g.readCount> 0) // get lengths 
				g.getRPK(t, pairedEnd, Graph.ETYPE_AL, eeV);
		
				//lenExon[0][0]= g.readCount> 0? getLength(eeV.elementAt(0), null): (t.getExonicLength()- readLen);
				//containerIntA1A1[0][0]= getLength(eeV.elementAt(0), null, true, false);
				long[][] others= new long[1][];
				others[0]= Graph.without(g.createAllArray(), g.encodeTset(new Transcript[] {t}));
				
				long[][] sig= new long[][]{g.encodeTset(t)};	// containerLongA1A[0]
				getGTFappend(sb, g, solver, eeV, perM, sig);
				
				eeV.elementAt(0).removeAllElements();
				
				return sb.toString();
			}

		private double getGTFappend(StringBuilder sb, Graph g, GraphLPsolver solver, Vector<Vector<Edge>> eeV, double perM, long[][] tid) {
				
				invariantTestObsSplitFreq= 0; 
				invariantTestPredSplitFreq= 0;
				sb.append(" ");
				Vector<double[]> containerVecLen= new Vector<double[]>();
				for (int x = 0; x < 3; x++) {	// virtual length sum, split, uniq
					boolean output= true;
					if ((x== 0&& !outputAll)|| (x== 1&& !outputSplit)|| (x==2&& !outputUnique))
						output= false;
					if (output) {
						sb.append(FluxCapacitorConstants.GTF_ATTRIBUTE_LENGTH);
						sb.append(FluxCapacitorConstants.GTF_ATTRIBUTE_TOKEN_SEP);
						sb.append(x== 0? FluxCapacitorConstants.GTF_ATTRIBUTE_TOKEN_ALL: (x==1? FluxCapacitorConstants.GTF_ATTRIBUTE_TOKEN_TID: FluxCapacitorConstants.GTF_ATTRIBUTE_TOKEN_EXC));
						sb.append(" ");	// \"
					}
					boolean excl= x==2? true: false;
					for (int i = 0; i < tid.length; i++) {
						long[] sig= x== 0? sigall: tid[i];
						if (i>= containerVecLen.size())
							containerVecLen.add(new double[3]);
						containerVecLen.elementAt(i)[x]= getLength(g, eeV.elementAt(i), sig, excl);
						if (output) {
							sb.append(Double.toString(containerVecLen.elementAt(i)[x]));
							sb.append(",");
						}
					}
					if (output) {
						sb.deleteCharAt(sb.length()- 1);
						sb.append("; ");	// \"
					}
				}
		
				double ms= miss?factor():1d;
				for (int i = 0; i < 3; i++) { // obs, pred, norm
					boolean output1= true;
					if ((i== 0&& !outputObs)|| (i== 1&& !outputPred)|| (i== 2&& !outputBalanced))
						output1= false;
					ms= (i== 0)?Math.round(ms):ms;
					
					ReadStatCalculator calc= (i== 0? FluxCapacitor.this: solver);
					for (int j = 0; j < 3; j++) {	// sum, split, uniq
						boolean output2= true;
						if ((!output1)|| (j== 0&& !outputAll)|| (j== 1&& !outputSplit)|| (j==2&& !outputUnique))
							output2= false;
						boolean excl= j==2? true: false;
						for (int x = 0; x < 2; x++) {	// reads, coverage
							boolean output= true;
							if ((!output2)|| (x== 0&& !outputFreq)
									|| (x== 1&& !outputRfreq)
									|| (x== 2&& !outputRcov))
								output= false;
							if (output) {
								sb.append(i==0? FluxCapacitorConstants.GTF_ATTRIBUTE_TOKEN_OBSV: (i==1? FluxCapacitorConstants.GTF_ATTRIBUTE_TOKEN_PRED:FluxCapacitorConstants.GTF_ATTRIBUTE_TOKEN_BALANCED));
								sb.append(FluxCapacitorConstants.GTF_ATTRIBUTE_TOKEN_SEP);
								sb.append(j== 0?FluxCapacitorConstants.GTF_ATTRIBUTE_TOKEN_ALL:(j==1? FluxCapacitorConstants.GTF_ATTRIBUTE_TOKEN_TID: FluxCapacitorConstants.GTF_ATTRIBUTE_TOKEN_EXC));
								sb.append(FluxCapacitorConstants.GTF_ATTRIBUTE_TOKEN_SEP);				
								sb.append(x== 0?FluxCapacitorConstants.GTF_ATTRIBUTE_TOKEN_READS: FluxCapacitorConstants.GTF_ATTRIBUTE_TOKEN_COV);	// pairedEnd?GTF_ATTRIBUTE_TOKEN_COV:GTF_ATTRIBUTE_TOKEN_RPKM
								sb.append(" ");	// no \"
							}
							for (int m = 0; m < tid.length; m++) {
								long[] sig= j== 0? sigall: tid[m];
		
								if (calc== null) {
									if (output)
										sb.append(FluxCapacitorConstants.VALUE_NA);
								} else {
									
									boolean normalized= i== 2? true: false;
									double val= ms* (j== 0? calc.getReads(eeV.elementAt(m), FluxCapacitorConstants.BYTE_0, sig, normalized):  // ALL								
										calc.getReadsAvg(eeV.elementAt(m), FluxCapacitorConstants.BYTE_0, g, sig, excl, normalized));	// TID, excl 							
									
//									if (val!= 0&& i>= 1&& j== 2&& containerVecLen.elementAt(m)[j]== 0)
//										System.currentTimeMillis();
//									
//									val= ms* (j== 0? calc.getReads(eeV.elementAt(m), BYTE_0, sig, normalized):  // ALL								
//										calc.getReadsAvg(eeV.elementAt(m), BYTE_0, g, sig, excl, normalized));	// TID, excl
//									System.currentTimeMillis();
									
									// TODO only negatives?
									try{assert(true|| val== 0|| val>= 0.0000000001);}catch(AssertionError err){
										if (Constants.verboseLevel> Constants.VERBOSE_NORMAL)
											System.err.println("Encountered value < e-10: "+ val);
										val= 0;	
									}
										

									if (x== 0) { //READS
										if (j== 1) {
											if (i== 0)
												invariantTestObsSplitFreq+= val;
											else if (i== 1)
												invariantTestPredSplitFreq+= val;
										}
										// why they should be int
//										if (i== 0)
//											sb.append(Integer.toString((int) Math.round(val)));
//										else
										if (output)
											sb.append(Float.toString((float) val));
									} else {	// RFREQ, COVERAGE
										double length= containerVecLen.elementAt(m)[j];
		//								if (length== 0&& val!= 0) {
		//									System.currentTimeMillis();
		//									getLength(g, eeV.elementAt(0), sig, excl);
		//									val= j== 0? calc.getReads(eeV.elementAt(m), BYTE_0, sig):  // ALL								
		//										calc.getReadsAvg(eeV.elementAt(m), BYTE_0, g, sig, excl);	// TID, excl
		//								}
										if (length== 0) {
											try {assert(val== 0);} catch (AssertionError e){
												if (Constants.verboseLevel> Constants.VERBOSE_NORMAL)
													System.err.println("Value found at 0 length: "+ val);
												val= 0;	
											};
										} else {
											val/= length;
										}
										if (output)
											sb.append(length== 0? FluxCapacitorConstants.FLOAT_STRING_0: Float.toString((float) val));
									} 
								}
								if (output)
									sb.append(",");
							}
							if (output) {
								sb.deleteCharAt(sb.length()- 1);
								sb.append("; ");	// no \"
							}
						}
					}
				}
				
				sb.append("\n"); // here?
				
				return invariantTestPredSplitFreq;
			}

		private void outputGFF_save(Graph g, ASEvent[] events, GraphLPsolver solver) {
					++nrLoci;
					if (solver!= null) 
						++nrLociExp;
					double perM= nrReadsAll/ 1000000d;
					// deprecated
						boolean unsolvedSystem= false;	
						double valOF= solver== null?0: solver.getValObjFunc();
						if (valOF> FluxCapacitorConstants.BIG) { 
							++nrUnsolved;
							unsolvedSystem= true;
						}
						String pv= getAttributeOF(valOF, solver, getMappedReadcount());
		
					StringBuilder sb= new StringBuilder();
					// LOCUS TODO genes 
					if (outputGene) {
						getGTF(sb, g.trpts[0].getGene(), g, solver, perM, pv);	
						try {assert(testInvariant(invariantTestObsSplitFreq, 
								pairedEnd?nrMappingsReadsOrPairs*2:nrMappingsReadsOrPairs, 0.05));}	// min: 0.01
						catch (AssertionError e) {
							if (Constants.verboseLevel> Constants.VERBOSE_NORMAL)
								System.err.println("[ASSERTION] "+getClass().getName()+".outputGFF():\n\tinvariantTestObsSplitFreq= "
										+ invariantTestObsSplitFreq+ ", nrMappingsReadsOrPairs= "+ (pairedEnd?nrMappingsReadsOrPairs*2:nrMappingsReadsOrPairs)
										+ "\n\tlocus: "+ g.trpts[0].getTranscriptID());
						};
						try {assert(testInvariant(invariantTestPredSplitFreq, 
								pairedEnd?nrMappingsReadsOrPairs*2:nrMappingsReadsOrPairs, 0.1));}
						catch (AssertionError e) {
							if (Constants.verboseLevel> Constants.VERBOSE_NORMAL)
								System.err.println("[ASSERTION] "+getClass().getName()+".outputGFF():\n\tinvariantTestPredSplitFreq= "
										+ invariantTestPredSplitFreq+ ", nrMappingsReadsOrPairs= "+ (pairedEnd?nrMappingsReadsOrPairs*2:nrMappingsReadsOrPairs)
										+ "\n\tlocus: "+ g.trpts[0].getTranscriptID());
						};
					}
		
					
					// TRANSCRIPTS
					if (outputTranscript|| outputExon|| outputSJunction) {
						float invariantObsAllTx= 0, invariantPredAllTx= 0,
							invariantObsAllEx= 0, invariantPredAllEx= 0;
						for (int i = 0; i < g.trpts.length; i++) {
							++nrTx;
		//					float invariantObsTx= invariantTestObsSplitFreq,
		//					invariantPredTx= invariantTestPredSplitFreq;
							float invariantObsTx= 0, invariantPredTx= 0;
							if (outputTranscript) {
								getGTF(sb, g.trpts[i], solver, g, perM, pv, false);	// writer.write
								invariantObsAllTx+= invariantTestObsSplitFreq; //invariantObsTx;
								invariantPredAllTx+= invariantTestPredSplitFreq; // invariantPredTx;
								invariantObsTx= invariantTestObsSplitFreq;
								invariantPredTx= invariantTestPredSplitFreq;
								if (invariantPredTx> 0)
									++nrTxExp;
							}
							// EXONS
							float invariantObsEx= 0, invariantPredEx= 0;
							if (outputExon) {
								Exon[] exons=  g.trpts[i].getExons();
								for (int j = 0; j < exons.length; j++) {
									getGTF(sb, exons[j], g.trpts[i], g, solver, unsolvedSystem, perM, pv, false);
									invariantObsEx+= invariantTestObsSplitFreq;
									invariantPredEx+= invariantTestPredSplitFreq;
								}
							}
							
							// SJ
							if (outputSJunction) {
								Vector<Vector<Edge>> eeV= new Vector<Vector<Edge>>(5,5);
								eeV.add(new Vector<Edge>());
								g.getRPK(g.trpts[i], pairedEnd, Graph.ETYPE_SJ, eeV);
								long[][] sig= new long[][]{g.encodeTset(g.trpts[i])};
								for (int j = 0; j < eeV.elementAt(0).size(); j++) { 
									getGTF(sb, eeV.elementAt(0).elementAt(j), sig, g, solver, perM);
									invariantObsEx+= invariantTestObsSplitFreq;
									invariantPredEx+= invariantTestPredSplitFreq;
								}
							}
							invariantObsAllEx+= invariantObsEx;
							invariantPredAllEx+= invariantPredEx;
							
							if (outputExon&& outputSJunction&& outputTranscript) {
								try {assert(testInvariant(invariantObsEx, invariantObsTx, 0.05));}	// min: 0.02
								catch (AssertionError e) {
									if (Constants.verboseLevel> Constants.VERBOSE_NORMAL)
										System.err.println("[ASSERTION] "+getClass().getName()+".outputGFF():\n\tinvariantObsEx= "
												+ invariantObsEx+ ", invariantObsTx= "+ invariantObsTx
												+ "\n\tlocus: "+ g.trpts[0].getTranscriptID());
								};
								try {assert(testInvariant(invariantPredEx, invariantPredTx, 0.1));}
								catch (AssertionError e) {
									if (Constants.verboseLevel> Constants.VERBOSE_NORMAL)
										System.err.println("[ASSERTION] "+getClass().getName()+".outputGFF():\n\tinvariantPredEx= "
												+ invariantPredEx+ ", invariantPredTx= "+ invariantPredTx
												+ "\n\tlocus: "+ g.trpts[0].getTranscriptID());
								};
							}
						}
						if (outputTranscript) {
							try {assert(testInvariant(invariantObsAllTx, 
									pairedEnd?nrMappingsReadsOrPairs*2:nrMappingsReadsOrPairs, 0.05));}	// min: 0.01
							catch (AssertionError e) {
								if (Constants.verboseLevel> Constants.VERBOSE_NORMAL)
									System.err.println("[ASSERTION] "+getClass().getName()+".outputGFF():\n\tinvariantObsAllTx= "
											+ invariantObsAllTx+ ", nrMappingsReadsOrPairs= "+ (pairedEnd?nrMappingsReadsOrPairs*2:nrMappingsReadsOrPairs)
											+ "\n\tlocus: "+ g.trpts[0].getTranscriptID());
							};
							try {assert(testInvariant(invariantPredAllTx, 
									pairedEnd?nrMappingsReadsOrPairs*2:nrMappingsReadsOrPairs, 0.1));}
							catch (AssertionError e) {
								if (Constants.verboseLevel> Constants.VERBOSE_NORMAL)
									System.err.println("[ASSERTION] "+getClass().getName()+".outputGFF():\n\tinvariantPredAllTx= "
											+ invariantPredAllTx+ ", nrMappingsReadsOrPairs= "+ (pairedEnd?nrMappingsReadsOrPairs*2:nrMappingsReadsOrPairs)
											+ "\n\tlocus: "+ g.trpts[0].getTranscriptID());
							};
						}
						if (outputExon&& outputSJunction) {
							try {assert(testInvariant(invariantObsAllEx, 
									pairedEnd?nrMappingsReadsOrPairs*2:nrMappingsReadsOrPairs, 0.05));}	// min: 0.02
							catch (AssertionError e) {
								if (Constants.verboseLevel> Constants.VERBOSE_NORMAL)
									System.err.println("[ASSERTION] "+getClass().getName()+".outputGFF():\n\tinvariantObsAllEx= "
											+ invariantObsAllEx+ ", nrMappingsReadsOrPairs= "+ (pairedEnd?nrMappingsReadsOrPairs*2:nrMappingsReadsOrPairs)
											+ "\n\tlocus: "+ g.trpts[0].getTranscriptID());
							};
							try {assert(testInvariant(invariantPredAllEx, 
									pairedEnd?nrMappingsReadsOrPairs*2:nrMappingsReadsOrPairs, 0.1));}
							catch (AssertionError e) {
								if (Constants.verboseLevel> Constants.VERBOSE_NORMAL)
									System.err.println("[ASSERTION] "+getClass().getName()+".outputGFF():\n\tinvariantPredAllEx= "
											+ invariantPredAllEx+ ", nrMappingsReadsOrPairs= "+ (pairedEnd?nrMappingsReadsOrPairs*2:nrMappingsReadsOrPairs)
											+ "\n\tlocus: "+ g.trpts[0].getTranscriptID());
							};
						}
					}
					
					// EVENTS
					if (outputEvent) {
						HashMap<Object,Double> tExpMap= null;
						if (solver!= null) {
							tExpMap= solver.getTrptExprHash();
							Object[] keys= tExpMap.keySet().toArray();	
							for (int i = 0; i < keys.length; i++) {
								if (!(keys[i] instanceof String))
									continue;
								if (tExpMap.get(keys[i])<0)
									tExpMap.put((String) keys[i], 0d);	// TODO ugly
							}
						}
						for (int i = 0; events!= null&& i < events.length; i++) {
							getGTF(sb, events[i], g, solver, unsolvedSystem, perM, pv, tExpMap);
						}
					}
		
					// FRAGMENTS and XJUNCTIONS
					if (false&& solver!= null) {
						ArrayList<Edge> cc= new ArrayList<Edge>();
						if (solver!= null) {
							Iterator<Object> iter= solver.getConstraintHash().keySet().iterator();
							while (iter.hasNext()) {
								Object o= iter.next();
								if (o instanceof Edge)
									cc.add((Edge) o);
							}
						}
						Collections.sort(cc, Edge.getDefaultPositionComparator());
		
						Iterator<Edge> iter= cc.iterator();
						while (iter.hasNext()) {
							Edge e= iter.next();
							// no INTRONS
							if ((!(e instanceof SuperEdge))&& (!e.isExonic()))
								continue;
							getGTF(sb, e, new long[][]{e.getTranscripts()}, g, solver, perM);
						}
					}
						
					try {
						write(sb);
					} catch (Exception e) {
						e.printStackTrace();
					}
					
				}

		private void learn(Transcript tx, BEDobject2[] beds) {
						
						if (beds== null|| beds.length== 0)
							return;
						
						int elen= tx.getExonicLength();	// TODO this is the effective length
						if (elen< readLenMin)
							return;	// TODO discards reads
						
						int binLen= Arrays.binarySearch(FluxCapacitorConstants.BIN_LEN, elen);
						if (binLen< 0)
							binLen= -(binLen+ 1);
						
						//int xLen= binLen< BIN_LEN.length? profileStub[binLen].length: 10000;
						float f= 20f/ elen;
						
						int plen= elen- (readLenMin- 1);
						TProfile profile= null;
						if (!uniform) {
							// not any more, now expr x length unique, and IDs
							// profile= func.getProfile(plen, (strand== STRAND_ENABLED), pairedEnd);
							profile= new TProfile(tx.getGene().getGeneID(), plen, (strand== FluxCapacitorConstants.STRAND_ENABLED), pairedEnd);
							func.profiles.add(profile);
						}
						HashMap<CharSequence, BEDobject2[][]> mapPends= new HashMap<CharSequence, BEDobject2[][]>();
						nrReadsSingleLoci+= beds.length;
						// possible read pairs
						// slow
		/*				for (int i = 0; i < beds.length; i++) {
							for (int j = i+1; j < beds.length; j++) {
								if (((beds[i].getStrand()== 0&& beds[j].getStrand()== 0)
										|| (beds[i].getStrand()== -beds[j].getStrand()))
										&& FMRD.matesByName(beds[i].getName(), beds[j].getName()))
									++nrReadsSingleLociPotentialPairs;
							}
						}
		*/				
						int cntSusp= 0;
						for (int i = 0; i < beds.length; i++) {
							
							if (beds[i]== null)
								continue;
							
							int len= beds[i].getLength();
							if (readLenMin< 0|| len< readLenMin)
								readLenMin= len;
							if (len> readLenMax)	// readLenMax< 0|| 
								readLenMax= len;
							int[] ePos= getExonicPos(tx, beds[i], elen);
							
							if (ePos== null) {
								//ePos= getExonicPos(tx, beds[i]);
								++nrReadsSingleLociNoAnnotation;
								continue; 	// doesnt align to the transcript
							}
		
							if (strand== FluxCapacitorConstants.STRAND_SPECIFIC&& beds[i].getStrand()!= tx.getStrand()) {
								++nrMappingsWrongStrand;
								continue;
							}
		
							assert(ePos[1]> ePos[0]);
							if(ePos[0]< 0|| ePos[1]>= elen) {
								//ePos= getExonicPos(tx, beds[i]);
								System.err.println("[ahuuuaaa] learn() does not report an exonic position for read\n\t"+ beds[i]);
								continue;
							}
							
							++nrReadsSingleLociMapped;
							if (pairedEnd) {
		
								byte flag=  (byte) (descriptor.getPairedEndInformation(beds[i].getName())- 1);	// (Fasta.getReadDir(beds[i].getName())- 1))
								CharSequence id= descriptor.getUniqueDescriptor(beds[i].getName());	//Fasta.getReadID(beds[i].getName())
								//CharSequence id= beds[i].getName();
								BEDobject2[][] oo= mapPends.get(id);
								if (oo== null) {
									oo= new BEDobject2[2][];
									mapPends.put(id, oo);
								} 
								if (oo[flag]== null) 
									oo[flag]= new BEDobject2[] {beds[i]};
								else {
									BEDobject2[] op= new BEDobject2[oo[flag].length+ 1];
									System.arraycopy(oo[flag], 0, op, 0, oo[flag].length);
									op[op.length- 1]= beds[i];
									oo[flag]= op;
								}
								
								for (int j = 0; j < oo.length; j++) {
									if (j==flag|| oo[j]== null)
										continue;
									for (int k = 0; k < oo[j].length; k++) {
										int[] ePos2= getExonicPos(tx, oo[j][k], elen);
										if (ePos[0]> ePos2[0]) {
											int[] h= ePos;
											ePos= ePos2;
											ePos2= h;
										}
										if (!uniform)
											profile.addReadPair(ePos[0], ePos2[0], readLenMin);
										addInsertSize(ePos2[1]- ePos[0]+ 1);
										++nrReadsSingleLociPairsMapped;	
		//								if (nrReadsSingleLociPairsMapped> nrReadsSingleLoci)
		//									System.currentTimeMillis();
				//						int delta= ePos2[0]- ePos[0];
									}
								}
									
							} else {
								
								
								byte dir= Constants.DIR_FORWARD;
								if (strand== FluxCapacitorConstants.STRAND_ENABLED) {
									//assert(ePos[0]> 0);	// can be, out of range
									if (!tx.getTranscriptID().startsWith("AT1G04050")) {
										if (ePos[0]> 0&& ePos[1]<= elen) {
											if (beds[i].getStrand()!= tx.getStrand()) {
												dir= Constants.DIR_BACKWARD;
												if (ePos[1]<= elen) {
													int p= (int) ((ePos[1]- 1)* f);
													if (binLen== 2&& p>=30&& p<= 40) {
														if (tx.getTranscriptID().startsWith("AT1G04050"))
															System.currentTimeMillis();
														++cntSusp;
													}
													++profileStubRev[binLen][p];
												}
											} else {
												if (ePos[0]> 0) {
													int p= (int) ((ePos[0]- 1)* f);
													if (binLen== 2&& p>=30&& p<= 40) {
														if (tx.getTranscriptID().startsWith("AT1G04050"))
															System.currentTimeMillis();
														++cntSusp;
													}
													++profileStub[binLen][p];
												}
											}
										}
									}
								} else if (!uniform)
									profile.addRead(ePos[0],readLenMin,dir);
							}
						}
						
						if (cntSusp> 100)
							System.currentTimeMillis();
						
						if ((!uniform)&& profile.getReads()== 0) {
							//profile.eliminate0regs();
							func.profiles.remove(profile);
						}
						//assert(func.getNrReadsInProfiles()== nrReadsLoci);
					}
	}
	
	public static boolean debug= false;
	public static boolean outputPbClusters= false;
	public boolean pairedEnd= false, stranded= false, force= false;
	int tolerance= 1000;	// +/- gene-near region
	public static boolean 
		cheatDoNotExit= false,
		cheatLearn= false, 
		cheatDisableFCheck= true,
		cheatDisableCleanup= true,
		cheatCopyFile= true,
		doUseLocusNormalization= false;
	
	static void exit(int code) {
		String pfx= "[ASTALAVISTA] ";
		if (code< 0)
			pfx= "[CIAO] ";
		if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) 
			System.err.println(pfx+ "I'm exiting.");
		System.exit(code);
	}
	
	/**
	 * @deprecated
	 * @return
	 */
	boolean checkPreliminaries() {
		if (outputProfiles&& uniform) {
			outputProfiles= false;
//			if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
//				System.err.println("[HEY] there are no profiles if you specify uniformal read distribution!\n"
//						+ "\t(parameter ["+CLI_SHORT_PFX+ CLI_SHORT_UNIF+"|"+ CLI_LONG_PFX+ CLI_LONG_UNIF+"])");
//			return false;
		}
		if (!(outputExon|| outputSJunction|| outputTranscript|| outputGene|| outputEvent)) {
			if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
				System.err.println("[WARNING] no features for output selected (["
						+ FluxCapacitorConstants.CLI_SHORT_PFX+ FluxCapacitorConstants.CLI_SHORT_OUT+ "|"+ FluxCapacitorConstants.CLI_LONG_PFX+ FluxCapacitorConstants.CLI_LONG_OUT+ "] [ejtgv]");
		}
		if (outputExon|| outputSJunction|| outputTranscript|| outputGene|| outputEvent) {
			if (!(outputObs|| outputPred|| outputBalanced)) {
				if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
					System.err.println("[UPS] no base specified for feature abundances");
				return false;
			}
			if (!(outputAll|| outputSplit|| outputUnique)) {
				if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
					System.err.println("[UHLALA] no scope specified for feature abundances");
				return false;
			}
			if (!(outputFreq|| outputRfreq|| outputRcov)) {
				if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
					System.err.println("[UHLALA] no measurement specified for feature abundances");
				return false;
			}
		}
		
		if (fileOutDir== null&& (outputMapped|| outputNotmapped|| outputISize|| outputProfiles|| outputLP)) {
			if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
				System.err.println("[UPS] no filename specified, cannot output additional information");
			return false;
		}
			
		
		if (outputISize&& !pairedEnd) {
			outputISize= false;
//			if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
//				System.err.println("[AIAI] insert sizes are only available in paired end mode, parameter [" 
//						+ Character.toString(CLI_SHORT_PFX)+ CLI_SHORT_PAIR+ "|"+ CLI_LONG_PFX+ CLI_LONG_PAIR+ "]");
//			return false;
		}
		
		if (strand== FluxCapacitorConstants.STRAND_SPECIFIC&& pairedEnd) {
			if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
				System.err.println("[NONO] strand specific reads can (currently) not be mate pairs, no?!");
			return false;
		}
		
		return true;
	}
	
	private static class StringArrayByFirstComparator implements Comparator<String[]> {
		public int compare(String[] o1, String[] o2) {
			return o1[0].compareTo(o2[0]);
		}
	}
	
	class LocusSolver2 extends Thread {
			Gene gene= null;
			ASEvent[] events= null;
			BufferedBEDiterator beds= null;
			boolean decompose= false;
			Thread threadBefore= null;
			int nrMappingsReadsOrPairs;
			HashSet<CharSequence> mapReadOrPairIDs;
			HashMap<CharSequence, Vector<BEDobject2>[]> mapEndsOfPairs;
			long[] sigall= null;
			UniversalReadDescriptor.Attributes attributes= null;
	
			private float invariantTestObsSplitFreq= 0, invariantTestPredSplitFreq= 0;
			
			public LocusSolver2(Gene newGene, BufferedBEDiterator newBeds, boolean decompose) {
				//super(newGene.getGeneID());
				
				this.gene= newGene; 
				this.beds= newBeds;
				this.decompose= decompose;
				
				nrMappingsReadsOrPairs= 0;
				mapReadOrPairIDs= new HashSet<CharSequence>();
				if (pairedEnd)
					mapEndsOfPairs = new HashMap<CharSequence, Vector<BEDobject2>[]>();
				attributes= descriptor2.createAttributes();
			}
			
			
			public void run() {
					
				try {
					//int mapped= mapReadOrPairIDs.size();
					if (decompose) {
						
						// BUG 110301: do not use mapTrivial, problems with split-maps
						
//						if (this.gene.getTranscriptCount()== 1) {
//							mapTrivial(gene.getTranscripts()[0], beds);
//							outputGFF(null, null, null);
//						} else {
							//Graph myGraph= getGraph(this.gene);
							AnnotationMapper mapper= new AnnotationMapper(this.gene);
							//map(myGraph, this.gene, this.beds); 
							mapper.map(this.beds, descriptor2);
							nrReadsLoci+= mapper.nrReadsLoci;
							nrReadsMapped+= mapper.getNrMappingsReadsOrPairs();
							nrMappingsReadsOrPairs+= mapper.getMapReadOrPairIDs().size()/ 2;
							nrPairsNoTxEvidence+= mapper.getNrPairsNoTxEvidence();
							nrPairsWrongOrientation+= mapper.getNrPairsWrongOrientation();
							
							GraphLPsolver2 mySolver= null;
							// != mapReadOrPairIDs.size()> 0, does also count singles
//							if (nrMappingsReadsOrPairs> 0&& this.gene.getTranscriptCount()> 1) {	// OPTIMIZE			
//								mySolver= getSolver(myGraph, nrMappingsReadsOrPairs* 2); // not: getMappedReadcount()
//								mySolver.run();
//							}
							if (mapper.nrMappingsReadsOrPairs> 0&& this.gene.getTranscriptCount()> 1) {	// OPTIMIZE			
								mySolver= getSolver(mapper, (int) (mapper.nrMappingsReadsOrPairs* 2)); // not: getMappedReadcount()
								mySolver.run();
							}
			//				if (this.gene.getTranscriptCount()== 1)
			//					System.currentTimeMillis();
							//outputGFF(myGraph, events, mySolver);
							outputGFF(mapper, events, mySolver);
							
//						}
					} else {
						// map all reads
						if (this.gene.getTranscriptCount()== 1) {
							++nrSingleTranscriptLearn;
							learn(this.gene.getTranscripts()[0], beds);	
						}
					}
				} catch (Throwable e) {
					System.err.println("\n[ERROR] in locus "+ gene.getGeneID());
					e.printStackTrace();
					System.err.print("\n\tcontinuing ");
					System.err.flush();
				}
				
				// cleanup
	//			for (int i = 0; (beds!= null)&& i < beds.length; i++) 
	//				if (beds[i]!= null)
	//					BEDobject.addRecycleObj(beds[i]);
				
				beds= null;
				gene= null;
				// makes it terribly slow
				//System.gc();
				
	//			synchronized(FluxCapacitor.this.threadPool) {
					FluxCapacitor.this.threadPool.remove(this);
	//			}
			}
	
			
			/**
			 * @deprecated problems with split-maps !
			 * @param tx
			 * @param beds
			 */
			private void mapTrivial(Transcript tx, BEDobject2[] beds) {
				
				nrMappingsReadsOrPairs= 0;
				if (beds== null|| beds.length== 0)
					return;
				
				// read pairing
				Arrays.sort(beds, BEDobject2.DEFAULT_ID_COMPARATOR);
				for (int i = 0; i < beds.length; i++) {
					CharSequence tag= beds[i].getName();
					if (descriptor2.getAttributes(tag, attributes)== null) {
						if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
							System.err.println("[IGNORED] Descriptor does not match readID "+ tag);
						continue;
					}
					
					// filter stranded
					if (stranded) {
						if ((tx.getStrand()== beds[i].getStrand()&& attributes.strand== 2)
								|| (tx.getStrand()!= beds[i].getStrand()&& attributes.strand== 1)) {
							beds[i]= null;
							continue;
						}
					}
					
					// filter pairs
					if (pairedEnd) {
						byte flag=  attributes.flag; //(byte) (descriptor.getPairedEndInformation(beds[i].getName())- 1);	// (Fasta.getReadDir(beds[i].getName())- 1))
						if (flag< 1) {
							if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
								System.err.println("Error in readID:\n"+ tag);
							continue;
						}
						CharSequence id= attributes.id; // descriptor.getUniqueDescriptor(tag);	//Fasta.getReadID(beds[i].getName())
						int sep= i, end= i+1;
						for (; end < beds.length; ++end) {
							if (!beds[end].getName().startsWith(id))
								break;
							
							if (sep== i) {
								tag= beds[end].getName();
								if (descriptor2.getAttributes(tag, attributes)== null) {
									if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
										System.err.println("[IGNORED] Descriptor does not match readID "+ tag);
									continue;
								}
								byte flag2= attributes.flag; // (byte) (descriptor.getPairedEndInformation(beds[end].getName())- 1);
								if (flag2< 1) {
									System.err.println("Error in readID:\n"+ beds[end].getName());
									continue;
								}
								if (flag2!= flag)
									sep= end;
							}
								
						}
						
						// no pair
						if (sep== i) {
							for (int j = i; j < end; j++) 
								beds[j]= null;
							i= end- 1;
							continue;
						}
						
						boolean first= false, second= false;
						for (int j = i; j < sep; j++) {
							if (contains(tx, beds[j])) {
								first= true;
								break;
							}
						}
						for (int j = sep; j < end; j++) {
							if (contains(tx, beds[j])) {
								second= true;
								break;
							}
						}
						// none of the pair
						if (!(first&& second)) {	// BUG 110301: changed from !(first|| second)
							for (int j = i; j < end; j++) 
								beds[j]= null;
						} else {
							++nrMappingsReadsOrPairs;
						}
						i= end- 1;
					
					} // end pend
					
				}
				

				if (!doUseLocusNormalization)
					return;
				
				// build individual matrix
				int tlen= tx.getExonicLength();
				byte tstrand= tx.getStrand();
				UniversalMatrix m= new UniversalMatrix(Math.max(tlen/ 100, 10));
				for (int i = 0; i < beds.length; i++) {
					if (beds[i]== null)
						continue;
					int p= getBpoint(tx, beds[i]);
					if (p< 0|| p>= tlen) 
						continue;		
					boolean sense= beds[i].getStrand()== tstrand;
					m.add(p, beds[i].getLength(), tlen, sense?Constants.DIR_FORWARD:Constants.DIR_BACKWARD); 
				}
				
				
				// normalize biases out
				//UniversalMatrix m= profile.getMatrix(tlen);	// no good idea
				if (nrMappingsReadsOrPairs> 100) {
					//better individual matrix, also for saturation biases
					//UniversalMatrix m= profile.getMatrix(tlen);
					
					/*if (tx.getTranscriptID().contains("ENST00000262241")||
							tx.getTranscriptID().contains("ENST00000252818"))
						System.currentTimeMillis();
					*/
					//System.err.println("\n"+ tx.getTranscriptID());
					int w= m.sense.length/ 5;
					int sumsOrig= m.sums, sumaOrig= m.suma;
					m.sums= Kernel.smoothen(Kernel.KERNEL_EPANECHNIKOV, 
							w, m.sense);
					m.suma= Kernel.smoothen(Kernel.KERNEL_EPANECHNIKOV, 
							w, m.asense);
					if (sumsOrig!= m.sums|| sumaOrig!= m.suma) {
						System.currentTimeMillis(); // can happen, rounding errors, eg 1->0
					}
					if (m.sums!= 0&& m.suma!= 0) {
						double f= m.getNfactor(0.2d);
						nrMappingsReadsOrPairs*= f;
					}
				}


			}


			Graph getGraph(Gene gene) {
					boolean output= false;
					
					// construct graph
				long t0= System.currentTimeMillis();
				
				Graph myGraph= new Graph(gene);
				//myGraph.createDefaultCoordComparator(readLenMin);
				myGraph.constructGraph();
				
				//if (outputLocus) {
//					myGraph.setRetrieveDSEvents(true);
//					myGraph.setRetrieveVSEvents(true);
//					if (myGraph.trpts.length> 1)
//						events= myGraph.getEvents(eventDim);
				//}
				
				myGraph.getNodesInGenomicOrder();	// important ??!
				myGraph.transformToFragmentGraph();
				if (output) {
					System.err.print(", transform "+((System.currentTimeMillis()- t0)/ 1000)+" sec, ");
					System.err.flush();
				}
				
/*				int nrSJ= myGraph.addEJ(readLenMin);
				if (output) {
					System.err.print(", EJ "+((System.currentTimeMillis()- t0)/ 1000)+" sec, ");
					System.err.flush();
				}

				insertMinMax= new int[] {0,1000};
				if (pairedEnd) {
					int nrPE= addPE(myGraph, insertMinMax, readLenMin);
					if (output) {
						System.err.print(", PE "+((System.currentTimeMillis()- t0)/ 1000)+" sec, ");
						System.err.flush();
					}
				}
*/				
			
				return myGraph;
			}
	
			Graph map(Graph myGraph, Gene g, BEDobject2[] beds) {
				
				
				if (beds== null|| beds.length== 0)
					return null;
				
				boolean output= false;
				long t0 = System.currentTimeMillis();
				
				// map read pairs
				for (int j = 0; beds!= null&& j< beds.length; ++j) {
	
					int xyxx= mapRead2(myGraph, beds[j], false);
					//nrMappingsReadsOrPairs+= xxx;
				}
				nrMappingsReadsOrPairs+= mapReadOrPairIDs.size()/ 2;

				// map single reads, mapping on single edges and inc count
				if (false) {
					for (int j = 0; beds!= null&& j< beds.length; ++j) {
						int xxx= mapRead2(myGraph, beds[j], true);
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
				
				int notMapped= 0;
				if (beds!= null)  {	// && mapOnly
					if (pairedEnd) {
						//assert(mappedReads== myGraph.mappedIDSet.size()); // no, multimaps
						//mappedReads*= 2;
					}
					notMapped= beds.length- nrMappingsReadsOrPairs;
					myGraph.setMappedReads(nrMappingsReadsOrPairs); 
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
	
				return myGraph;
			}
					
			private int nrLocusMultimaps= 0;
			/**
			 * add a SINGLE read
			 * @param regs
			 */
			int mapRead2(Graph g, BEDobject2 dobject, boolean force) {
				// find the edge(s) where the regions align
				// if these do not form a continous chain, create a new edge
				
	//			GFFObject[] gtfs= GFFObject.fromBed(dobject);	// TODO kill intermediary GTFs !!!
	//			DirectedRegion[] regs= new DirectedRegion[gtfs.length];
	//			for (int i = 0; i < regs.length; i++) 
	//				regs[i]= new DirectedRegion(gtfs[i]);
				
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
				if (pairedEnd) {
					flag= attributes.flag;
					if (flag<= 0) {
						if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) 
							System.err.println("Error in readID:\n"+ dobject.getName());
						return 0;
					}
				}
				CharSequence ID= tag;	//getID(dobject);
				if (pairedEnd|| stranded) 
					ID= attributes.id;
	
//				if (ID.equals("HWUSI-EAS626_1:5:82:1446:1847"))
//					System.currentTimeMillis();
				
				Edge target= g.getEdge2(dobject);
				
				if (target== null)
					return 0;
				if (force) {
					boolean sense= g.trpts[0].getStrand()== dobject.getStrand();
					if (sense)
						target.incrReadNr();
					else
						target.incrRevReadNr();
					mapReadOrPairIDs.add(dobject.getName());
					return 1;
				}
				
				
				byte refStrand= g.trpts[0].getStrand();
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
				if (pairedEnd) { /* PAIRED END */
	
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
						
						Edge target2= g.getEdge2(dobject2);
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
						SuperEdge se= g.getSuperEdge(w, true, null);
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
			
			private CharSequence getID(BEDobject2 dobject) {
				CharSequence id= null; 	
				if (pairedEnd) {
					id= descriptor.getUniqueDescriptor(dobject.getName());
				} else
					id= dobject.getName();
	
				return id;
			}
			
			private byte getFlag(BEDobject2 dobject) {
				byte flag= 0;  	
				if (pairedEnd) 
					flag= (byte) (descriptor.getPairedEndInformation(dobject.getName()));
						
				return flag;
			}
			
			public int getMappedReadcount() {
				if (pairedEnd)
					return mapReadOrPairIDs.size()/ 2;
				return mapReadOrPairIDs.size();
			}
			
	
			private boolean writeMappedRead(BEDobject o) {
				if (getFileMappedReads()== null)
					return false;
				try {
					getWriterMappedReads().write(o.toString()+ "\n");
					return true;
				} catch (IOException e) {			
					e.printStackTrace();
					return false;
				}
			}
			
			private boolean writeNotmappedRead(BEDobject o) {
				if (getFileNotMappedReads()== null)
					return false;
				try {
					getWriterNotmappedReads().write(o.toString()+ "\n");
					return true;
				} catch (IOException e) {
					e.printStackTrace();
					return false;
				}
			}
			
			
			private void outputGFF(Graph g, ASEvent[] events, GraphLPsolver2 solver) {
				++nrLoci;
				if (solver!= null) 
					++nrLociExp;
				double perM= nrReadsAll/ 1000000d;
				// deprecated
				boolean unsolvedSystem= false;	
				double valOF= solver== null?0: solver.getValObjFunc();
				if (valOF> FluxCapacitorConstants.BIG) { 
					++nrUnsolved;
					unsolvedSystem= true;
				}
				//String pv= getAttributeOF(valOF, solver, getMappedReadcount());
				Transcript[] tt= gene.getTranscripts();
	
				// prebuild rpkm hash
				HashMap<String, Float> rpkmMap= null;
				if (outputBalanced) {
					rpkmMap= new HashMap<String, Float>(tt.length, 1f);
					for (int i = 0; i < tt.length; i++) {
						Transcript tx= tt[i];
						float val= 0f, rpkm= 0f;
						if (solver== null) {
							// no reads
							val= nrMappingsReadsOrPairs;
						} else {
	//						if (solver.getTrptExprHash().get(g.trpts[i].getTranscriptID())!= 0)
	//							System.currentTimeMillis();
							val= (float) (solver.getTrptExprHash().get(tx.getTranscriptID()).doubleValue());
							if (val< 1- costBounds[0]) // 1- 0.95
								val= 0;
							try {
								assert(tt.length> 1|| val== nrMappingsReadsOrPairs);
							} catch (AssertionError e) {
								System.err.println(val+ " x "+ nrMappingsReadsOrPairs);
								solver.getNFactor();
							}
						}
//						if (pairedEnd)
//							val*= 2;	// count both ends for paired-end reads
						if (val> 0&& !(outputObs|| outputPred))
							++nrTxExp;
						rpkm= calcRPKM(val, tx.getExonicLength());
						rpkmMap.put(tx.getTranscriptID(), rpkm);
						// TODO chk
						if (Float.isNaN(rpkmMap.get(tt[i].getTranscriptID()).floatValue()))
							System.currentTimeMillis();
					}
				}
				
				
				// reproduce original
				boolean foundTranscripts= false, foundExons= false;
				if (getGTFreader().isKeepOriginalLines()&& origLines!= null) {
					foundExons= true;
					for (int i = 0; i < origLines.size(); i++) {
						String s= origLines.elementAt(i);
						String feat= GFFObject.getField(3, s);
						String tid= GFFObject.getTranscriptID(s);
						int tx= 0;
						if ((feat.equals(feat.equals(Transcript.GFF_FEATURE_TRANSCRIPT))
								|| feat.equals(Exon.GFF_FEATURE_EXON))
								&& (outputObs|| outputPred))
							for (tx = 0; tx < tt.length; tx++) 
								if (tt[tx].getTranscriptID().equals(tid))
									break;
						if (tx>= tt.length) {
							System.err.println("\nTranscript "+ tid+ " not found in: ");
							for (int j = 0; j < tt.length; j++) 
								System.err.println("\t"+ tt[j].getTranscriptID());
							System.err.println();
						}
						
						if (feat.equals(Transcript.GFF_FEATURE_TRANSCRIPT)&& outputTranscript) {
							foundTranscripts= true;
							StringBuilder sb= new StringBuilder(s);
							int x= sb.length();	// trim
							while(Character.isWhitespace(sb.charAt(--x)))
								sb.delete(x, x+1);
							if (sb.charAt(x)!= ';')
								sb.append("; ");
							else
								sb.append(Constants.SPACE);
							
							if ((outputObs|| outputPred)&& tx< tt.length)
								; //getGTF(sb, tt[tx], solver, g, perM, pv, true);
							else if (outputBalanced) {
								sb.append(FluxCapacitorConstants.GTF_ATTRIBUTE_TOKEN_RPKM);
								sb.append(Constants.SPACE);
								if (rpkmMap.containsKey(tid))
									sb.append(String.format("%1$f", rpkmMap.get(tid).floatValue()));	// rgasp parser does not like scientific notation
								else
									sb.append(Constants.NULL);
								sb.append(";\n");
							}
							
							try {
								getWriter().write(sb.toString());
							} catch (IOException e) {
								e.printStackTrace();
							} 
						} else if (feat.equals(Exon.GFF_FEATURE_EXON)&& outputExon) {
							
							StringBuilder sb= new StringBuilder(s); 
							int x= sb.length();
							while(Character.isWhitespace(sb.charAt(--x)))
								sb.delete(x, x+1);
							if (sb.charAt(x)!= ';')
								sb.append("; ");
							else
								sb.append(' ');
							
													
							if ((outputObs|| outputPred)&& tx< tt.length) {
								int start= Integer.parseInt(GFFObject.getField(4, s));
								int end= Integer.parseInt(GFFObject.getField(5, s));
								int j = 0;
								for (; j < tt[x].getExons().length; j++) {
									int begin= Math.abs(tt[x].getExons()[j].getStart()),
										ende= Math.abs(tt[x].getExons()[j].getEnd());
									if (begin> ende) {
										int h= begin;
										begin= ende;
										ende= h;
									}
									if (begin== start&& ende== end)
										break;
								}
								//getGTF(sb, tt[x].getExons()[j], tt[i], g, solver, unsolvedSystem, perM, pv, true);
							} else if (outputBalanced) {
								sb.append(FluxCapacitorConstants.GTF_ATTRIBUTE_TOKEN_RPKM);
								sb.append(Constants.SPACE);
								if (rpkmMap.containsKey(tid))
									sb.append(String.format("%1$f", rpkmMap.get(tid).floatValue()));	// rgasp parser does not like scientific notation
								else
									sb.append(Constants.NULL);
								sb.append(";\n");
							}
							
							
							try {
								getWriter().write(sb.toString());
							} catch (IOException e) {
								e.printStackTrace();
							}
						} else if (outputUnknown) {
							try {
								getWriter().write(s+ System.getProperty("line.separator"));
							} catch (IOException e) {
								e.printStackTrace();
							}
						}
					}
				}
				
					
				
				StringBuilder sb= new StringBuilder();
				// LOCUS TODO genes 
				if (outputGene) {
					if (outputObs|| outputPred) {
						//getGTF(sb, g.trpts[0].getGene(), g, solver, perM, pv);	
						try {assert(testInvariant(invariantTestObsSplitFreq, 
								pairedEnd?nrMappingsReadsOrPairs*2:nrMappingsReadsOrPairs, 0.05));}	// min: 0.01
						catch (AssertionError e) {
							if (Constants.verboseLevel> Constants.VERBOSE_NORMAL)
								System.err.println("[ASSERTION] "+getClass().getName()+".outputGFF():\n\tinvariantTestObsSplitFreq= "
										+ invariantTestObsSplitFreq+ ", nrMappingsReadsOrPairs= "+ (pairedEnd?nrMappingsReadsOrPairs*2:nrMappingsReadsOrPairs)
										+ "\n\tlocus: "+ g.trpts[0].getTranscriptID());
						};
						try {assert(testInvariant(invariantTestPredSplitFreq, 
								pairedEnd?nrMappingsReadsOrPairs*2:nrMappingsReadsOrPairs, 0.1));}
						catch (AssertionError e) {
							if (Constants.verboseLevel> Constants.VERBOSE_NORMAL)
								System.err.println("[ASSERTION] "+getClass().getName()+".outputGFF():\n\tinvariantTestPredSplitFreq= "
										+ invariantTestPredSplitFreq+ ", nrMappingsReadsOrPairs= "+ (pairedEnd?nrMappingsReadsOrPairs*2:nrMappingsReadsOrPairs)
										+ "\n\tlocus: "+ g.trpts[0].getTranscriptID());
						};
						
					} else if (outputBalanced) {
					}
				}
	
				
				// TRANSCRIPTS
				if (outputTranscript|| outputExon|| outputSJunction) {
					float invariantObsAllTx= 0, invariantPredAllTx= 0,
						invariantObsAllEx= 0, invariantPredAllEx= 0;
					for (int i = 0; i < tt.length; i++) {
						++nrTx;
	//					float invariantObsTx= invariantTestObsSplitFreq,
	//					invariantPredTx= invariantTestPredSplitFreq;
						float invariantObsTx= 0, invariantPredTx= 0;
						if (outputTranscript&& !foundTranscripts) {
							if (outputObs|| outputPred) {
								//getGTF(sb, g.trpts[i], solver, g, perM, null, false);	// writer.write
								invariantObsAllTx+= invariantTestObsSplitFreq; //invariantObsTx;
								invariantPredAllTx+= invariantTestPredSplitFreq; // invariantPredTx;
								invariantObsTx= invariantTestObsSplitFreq;
								invariantPredTx= invariantTestPredSplitFreq;
								if (invariantPredTx> 0)
									++nrTxExp;
	
							} else if (outputBalanced) {
								GFFObject obj= GFFObject.createGFFObject(tt[i]);
								sb.append(obj.toString());
								int x= sb.length();
								while(Character.isWhitespace(sb.charAt(--x)))
									sb.delete(x, x+1);
								if (sb.charAt(x)!= ';')
									sb.append("; ");
								else
									sb.append(Constants.SPACE);
								
								sb.append(FluxCapacitorConstants.GTF_ATTRIBUTE_TOKEN_RPKM);
								sb.append(Constants.SPACE);							
								//sb.append(rpkmMap.get(g.trpts[i].getTranscriptID()));
								sb.append(String.format("%1$f", rpkmMap.get(tt[i].getTranscriptID()).floatValue()));	// rgasp parser does not like scientific notation
								sb.append(";\n");
							}
						}
						// EXONS
						float invariantObsEx= 0, invariantPredEx= 0;
						if (outputExon&& !foundExons) {
							Exon[] exons=  tt[i].getExons();
							for (int j = 0; j < exons.length; j++) {
								//getGTF(sb, exons[j], tt[i], g, solver, unsolvedSystem, perM, null, false);
								invariantObsEx+= invariantTestObsSplitFreq;
								invariantPredEx+= invariantTestPredSplitFreq;
							}
						}
						
						// SJ
						if (outputSJunction) {
							Vector<Vector<Edge>> eeV= new Vector<Vector<Edge>>(5,5);
							eeV.add(new Vector<Edge>());
							g.getRPK(tt[i], pairedEnd, Graph.ETYPE_SJ, eeV);
							long[][] sig= new long[][]{g.encodeTset(tt[i])};
							for (int j = 0; j < eeV.elementAt(0).size(); j++) { 
								//getGTF(sb, eeV.elementAt(0).elementAt(j), sig, g, solver, perM);
								invariantObsEx+= invariantTestObsSplitFreq;
								invariantPredEx+= invariantTestPredSplitFreq;
							}
						}
						invariantObsAllEx+= invariantObsEx;
						invariantPredAllEx+= invariantPredEx;
						
						if (outputExon&& outputSJunction&& outputTranscript) {
							try {assert(testInvariant(invariantObsEx, invariantObsTx, 0.05));}	// min: 0.02
							catch (AssertionError e) {
								if (Constants.verboseLevel> Constants.VERBOSE_NORMAL)
									System.err.println("[ASSERTION] "+getClass().getName()+".outputGFF():\n\tinvariantObsEx= "
											+ invariantObsEx+ ", invariantObsTx= "+ invariantObsTx
											+ "\n\tlocus: "+ tt[0].getTranscriptID());
							};
							try {assert(testInvariant(invariantPredEx, invariantPredTx, 0.1));}
							catch (AssertionError e) {
								if (Constants.verboseLevel> Constants.VERBOSE_NORMAL)
									System.err.println("[ASSERTION] "+getClass().getName()+".outputGFF():\n\tinvariantPredEx= "
											+ invariantPredEx+ ", invariantPredTx= "+ invariantPredTx
											+ "\n\tlocus: "+ tt[0].getTranscriptID());
							};
						}
					}
					if (outputTranscript) {
						try {assert(testInvariant(invariantObsAllTx, 
								pairedEnd?nrMappingsReadsOrPairs*2:nrMappingsReadsOrPairs, 0.05));}	// min: 0.01
						catch (AssertionError e) {
							if (Constants.verboseLevel> Constants.VERBOSE_NORMAL)
								System.err.println("[ASSERTION] "+getClass().getName()+".outputGFF():\n\tinvariantObsAllTx= "
										+ invariantObsAllTx+ ", nrMappingsReadsOrPairs= "+ (pairedEnd?nrMappingsReadsOrPairs*2:nrMappingsReadsOrPairs)
										+ "\n\tlocus: "+ tt[0].getTranscriptID());
						};
						try {assert(testInvariant(invariantPredAllTx, 
								pairedEnd?nrMappingsReadsOrPairs*2:nrMappingsReadsOrPairs, 0.1));}
						catch (AssertionError e) {
							if (Constants.verboseLevel> Constants.VERBOSE_NORMAL)
								System.err.println("[ASSERTION] "+getClass().getName()+".outputGFF():\n\tinvariantPredAllTx= "
										+ invariantPredAllTx+ ", nrMappingsReadsOrPairs= "+ (pairedEnd?nrMappingsReadsOrPairs*2:nrMappingsReadsOrPairs)
										+ "\n\tlocus: "+ tt[0].getTranscriptID());
						};
					}
					if (outputExon&& outputSJunction) {
						try {assert(testInvariant(invariantObsAllEx, 
								pairedEnd?nrMappingsReadsOrPairs*2:nrMappingsReadsOrPairs, 0.05));}	// min: 0.02
						catch (AssertionError e) {
							if (Constants.verboseLevel> Constants.VERBOSE_NORMAL)
								System.err.println("[ASSERTION] "+getClass().getName()+".outputGFF():\n\tinvariantObsAllEx= "
										+ invariantObsAllEx+ ", nrMappingsReadsOrPairs= "+ (pairedEnd?nrMappingsReadsOrPairs*2:nrMappingsReadsOrPairs)
										+ "\n\tlocus: "+ tt[0].getTranscriptID());
						};
						try {assert(testInvariant(invariantPredAllEx, 
								pairedEnd?nrMappingsReadsOrPairs*2:nrMappingsReadsOrPairs, 0.1));}
						catch (AssertionError e) {
							if (Constants.verboseLevel> Constants.VERBOSE_NORMAL)
								System.err.println("[ASSERTION] "+getClass().getName()+".outputGFF():\n\tinvariantPredAllEx= "
										+ invariantPredAllEx+ ", nrMappingsReadsOrPairs= "+ (pairedEnd?nrMappingsReadsOrPairs*2:nrMappingsReadsOrPairs)
										+ "\n\tlocus: "+ tt[0].getTranscriptID());
						};
					}
				}
				
				// EVENTS
				if (outputEvent) {
					HashMap<Object,Double> tExpMap= null;
					if (solver!= null) {
						tExpMap= solver.getTrptExprHash();
						Object[] keys= tExpMap.keySet().toArray();	
						for (int i = 0; i < keys.length; i++) {
							if (!(keys[i] instanceof String))
								continue;
							if (tExpMap.get(keys[i])<0)
								tExpMap.put((String) keys[i], 0d);	// TODO ugly
						}
					}
					for (int i = 0; events!= null&& i < events.length; i++) {
						if (outputObs|| outputPred)
							; //getGTF(sb, events[i], g, solver, unsolvedSystem, perM, pv, tExpMap);
						else
							++nrEvents;
						if (outputBalanced) {
							sb.append(events[i].toStringGTF());
							sb.append(Constants.SPACE);
							sb.append("\"");
							boolean allPos= true;
							for (int j = 0; j < events[i].getTranscripts().length; j++) {
								float sum= 0;
								for (int k = 0; k < events[i].getTranscripts()[j].length; k++) 
									sum+= rpkmMap.get(events[i].getTranscripts()[j][k].getTranscriptID());
								sb.append(sum);
								sb.append(",");
								allPos&= (sum> 0);
							}
							if (allPos&& !(outputObs|| outputPred))
								++nrEventsExp;
							sb.replace(sb.length()- 1, sb.length(), "\";\n");
						}
	
					}
				}
	
				// FRAGMENTS and XJUNCTIONS
				if (false&& solver!= null) {
					ArrayList<Edge> cc= new ArrayList<Edge>();
					if (solver!= null) {
						Iterator<Object> iter= solver.getConstraintHash().keySet().iterator();
						while (iter.hasNext()) {
							Object o= iter.next();
							if (o instanceof Edge)
								cc.add((Edge) o);
						}
					}
					Collections.sort(cc, Edge.getDefaultPositionComparator());
	
					Iterator<Edge> iter= cc.iterator();
					while (iter.hasNext()) {
						Edge e= iter.next();
						// no INTRONS
						if ((!(e instanceof SuperEdge))&& (!e.isExonic()))
							continue;
						//getGTF(sb, e, new long[][]{e.getTranscripts()}, g, solver, perM);
					}
				}
					
				try {
					write(sb);
				} catch (Exception e) {
					e.printStackTrace();
				}
				
			}
	
			private boolean testInvariant(double invariant, double reference, double stringency) {
				double delta= Math.abs(reference== 0?invariant: (invariant- reference)/ reference);
				if (delta> stringency) {
					if (invariant<= 0)
						return true; 	// catch 0-predictions
					return false;
				}
				
				return true;
			}
	
			public void setThreadBefore(Thread threadBefore) {
				this.threadBefore = threadBefore;
			}
	
			private void learn_old(Transcript tx, BEDobject2[] beds) {
					
					if (beds== null|| beds.length== 0)
						return;
					
					int elen= tx.getExonicLength();	// TODO this is the effective length
					if (elen< readLenMin)
						return;	// TODO discards reads
					
					int plen= elen- (readLenMin- 1);
					TProfile profile= null;
					if (!uniform) {
						// not any more, now expr x length unique, and IDs
						// profile= func.getProfile(plen, (strand== STRAND_ENABLED), pairedEnd);
						profile= new TProfile(tx.getGene().getGeneID(), plen, (strand== FluxCapacitorConstants.STRAND_ENABLED), pairedEnd);
						func.profiles.add(profile);
					}
					HashMap<CharSequence, BEDobject2[][]> mapPends= new HashMap<CharSequence, BEDobject2[][]>();
					nrReadsSingleLoci+= beds.length;
					// possible read pairs
					// slow
	/*				for (int i = 0; i < beds.length; i++) {
						for (int j = i+1; j < beds.length; j++) {
							if (((beds[i].getStrand()== 0&& beds[j].getStrand()== 0)
									|| (beds[i].getStrand()== -beds[j].getStrand()))
									&& FMRD.matesByName(beds[i].getName(), beds[j].getName()))
								++nrReadsSingleLociPotentialPairs;
						}
					}
	*/				
					for (int i = 0; i < beds.length; i++) {
						
						if (beds[i]== null)
							continue;
						
						int len= beds[i].getLength();
						if (readLenMin< 0|| len< readLenMin)
							readLenMin= len;
						if (len> readLenMax)	// readLenMax< 0|| 
							readLenMax= len;
						int[] ePos= getExonicPos(tx, beds[i], elen);
						
						if (ePos== null) {
							//ePos= getExonicPos(tx, beds[i]);
							++nrReadsSingleLociNoAnnotation;
							continue; 	// doesnt align to the transcript
						}
	
						if (strand== FluxCapacitorConstants.STRAND_SPECIFIC&& beds[i].getStrand()!= tx.getStrand()) {
							++nrMappingsWrongStrand;
							continue;
						}
	
						assert(ePos[1]> ePos[0]);
						if(ePos[0]< 0|| ePos[1]>= elen) {
							//ePos= getExonicPos(tx, beds[i]);
							System.err.println("[ahuuuaaa] learn() does not report an exonic position for read\n\t"+ beds[i]);
							continue;
						}
						
						++nrReadsSingleLociMapped;
						if (pairedEnd) {
	
							byte flag=  (byte) (descriptor.getPairedEndInformation(beds[i].getName())- 1);	// (Fasta.getReadDir(beds[i].getName())- 1))
							CharSequence id= descriptor.getUniqueDescriptor(beds[i].getName());	//Fasta.getReadID(beds[i].getName())
							//CharSequence id= beds[i].getName();
							BEDobject2[][] oo= mapPends.get(id);
							if (oo== null) {
								oo= new BEDobject2[2][];
								mapPends.put(id, oo);
							} 
							if (oo[flag]== null) 
								oo[flag]= new BEDobject2[] {beds[i]};
							else {
								BEDobject2[] op= new BEDobject2[oo[flag].length+ 1];
								System.arraycopy(oo[flag], 0, op, 0, oo[flag].length);
								op[op.length- 1]= beds[i];
								oo[flag]= op;
							}
							
							for (int j = 0; j < oo.length; j++) {
								if (j==flag|| oo[j]== null)
									continue;
								for (int k = 0; k < oo[j].length; k++) {
									int[] ePos2= getExonicPos(tx, oo[j][k], elen);
									if (ePos[0]> ePos2[0]) {
										int[] h= ePos;
										ePos= ePos2;
										ePos2= h;
									}
									if (!uniform)
										profile.addReadPair(ePos[0], ePos2[0], readLenMin);
									addInsertSize(ePos2[1]- ePos[0]+ 1);
									++nrReadsSingleLociPairsMapped;	
	//								if (nrReadsSingleLociPairsMapped> nrReadsSingleLoci)
	//									System.currentTimeMillis();
			//						int delta= ePos2[0]- ePos[0];
								}
							}
								
						} else {
							byte dir= Constants.DIR_FORWARD;
							if (strand== FluxCapacitorConstants.STRAND_ENABLED&& beds[i].getStrand()!= tx.getStrand())
								dir= Constants.DIR_BACKWARD;
							
							if (!uniform)
								profile.addRead(ePos[0],readLenMin,dir);
						}
					}
					
					if ((!uniform)&& profile.getReads()== 0) {
						//profile.eliminate0regs();
						func.profiles.remove(profile);
					}
					//assert(func.getNrReadsInProfiles()== nrReadsLoci);
				}
	
			private GraphLPsolver2 getSolver(Graph g, int mappedReads) {
			
				GraphLPsolver2 solver= new GraphLPsolver2(g, readLenMin, 
						pairedEnd?insertMinMax:null, mappedReads, 
						(strand== FluxCapacitorConstants.STRAND_ENABLED), 
						pairedEnd);
				if (outputLP)
					solver.setFileLPdir(getFileLP());
				solver.costModel= costModel;	// COSTS_LINEAR
				solver.setCostSplit(costSplit);
				solver.setProfile(profile);
				solver.setReadLen(readLenMin);
				solver.setFlow(true);
				solver.costBounds= costBounds;
			
				return solver;
			}
	
			private String getGTF(StringBuilder sb, ASEvent event, Graph g, GraphLPsolver solver, boolean unsolvedSystem, 
						double perM, String pv, HashMap<Object,Double> tExpMap) {
					
			//		for (int i = 0; i < eeV.size(); i++) 
			//			eeV.elementAt(i).removeAllElements();
					Vector<Vector<Edge>> eeV= new Vector<Vector<Edge>>(5,5);
					while (eeV.size()< event.getDimension()) 
						eeV.add(new Vector<Edge>());
					
					g.getRPK(event, pairedEnd, Graph.ETYPE_AL, eeV);
					sb.append(event.toStringGTF());
			
					long[][] sig= new long[event.getDimension()][];
					for (int i = 0; i < sig.length; i++) { 
						sig[i]= g.createAllArray();		
						for (int j = 0; j < sig.length; j++) 
							sig[i]= j== i? sig[i]: 
								Graph.unite(sig[i], g.encodeTset(event.getTranscripts()[j]));
					}
					
					double splitReads= getGTFappend(sb, g, solver, eeV, perM, sig);
					++nrEvents;
					if (splitReads> 0)
						++nrEventsExp;
			
					for (int i = 0; i < eeV.size(); i++) 
						eeV.elementAt(i).removeAllElements();
					
					return sb.toString();
				}
	
			private String getGTF(StringBuilder sb, Edge e, long[][] sig, Graph g, GraphLPsolver solver, 
						double perM) {
					
					Vector<Vector<Edge>> eeV= new Vector<Vector<Edge>>(5,5);
					eeV.add(new Vector<Edge>(1));
					eeV.elementAt(0).add(e);
					
					sb.append(g.trpts[0].getChromosome());
					sb.append("\t");
					Transcript[] tt= g.decodeTset(e.getTranscripts());
			//		if (tt[0].getTranscriptID().equals("ENST00000407980"))
			//			System.currentTimeMillis();
					sb.append(Transcript.getSource(tt));		
					sb.append("\t");
					if (e instanceof SuperEdge) {
						if (((SuperEdge) e).isPend())
							sb.append(FluxCapacitorConstants.GFF_FEATURE_PAIRED);
						else 
							sb.append(FluxCapacitorConstants.GFF_FEATURE_JUNCTION);
					} else
						sb.append(FluxCapacitorConstants.GFF_FEATURE_FRAGMENT);
					sb.append("\t");
					
					int[] frac= e.getFrac(tt[0], readLenMin);
					int start= Math.abs(tt[0].getGenomicPosition(frac[0]));
					int end= Math.abs(tt[0].getGenomicPosition(frac[1]));
					if (start>end) {
						int h= start;
						start= end;
						end= h;
					}
					
					sb.append(Integer.toString(start));
					sb.append("\t");
					sb.append(Integer.toString(end));
					sb.append("\t.\t");
					sb.append(GFFObject.getStrandSymbol(tt[0].getStrand()));
					sb.append("\t.\t");
					
					sb.append(GFFObject.TRANSCRIPT_ID_TAG+" \"");
					for (int j = 0; j < tt.length; j++) {
						sb.append(tt[j].getTranscriptID());
						if (j< tt.length-1)
							sb.append(",");
					}
					sb.append("\";");
					
					sb.append(" ");
					sb.append(GFFObject.GENE_ID_TAG);
					sb.append(" \"");
					sb.append(tt[0].getGene().getGeneID());
					sb.append("\";");
			
					sb.append(" ");
					sb.append("edges");
					sb.append(" \"");
					sb.append(e.toString());
					sb.append("\";");			
					
	//				sb.append(" ");
	//				sb.append(GTF_ATTRIBUTE_LENGTH);
	//				sb.append(" ");	// \"
	//				//sb.append(frac[1]- frac[0]+ 1);
	//				sb.append(getLength(g, eeV.elementAt(0), e.getTranscripts(), false));
	//				sb.append(";");	// \"			
					
					// TODO why not?
					getGTFappend(sb, g, solver, eeV, perM, sig);
		
					// here was the alternative:
	/*				sb.append(" ");
					sb.append(GTF_ATTRIBUTE_TOKEN_OBSV);
					sb.append(" \"");
					sb.append(e.getReadNr());
					sb.append("\";");			
							
					sb.append(" ");
					sb.append(GTF_ATTRIBUTE_TOKEN_PRED);
					sb.append(" \"");
					for (int j = 0; j < tt.length; j++) {
						sb.append(solver.getNFactor()* solver.getTrptExprHash().get(tt[j].getTranscriptID())* 
								solver.getSuperProfileMap().get(tt[j].getTranscriptID()).getAreaFrac(
								GraphLPsolver.bounds2rel(frac, tt[j].getExonicLength()- readLenMin), readLenMin, 
								strandSpecific?TProfile.DIR_BOTH:TProfile.DIR_FORWARD));
						if (j< tt.length-1)
							sb.append(",");
					}
					sb.append("\";");			
			
					sb.append(" ");
					sb.append(GTF_ATTRIBUTE_TOKEN_PRED+"_total");
					sb.append(" \"");
					sb.append(solver.getReads(eeV.elementAt(0), BYTE_0, e.getTranscripts()));
					sb.append("\";");			
			
					
					// expectations
					sb.append(" ");
					sb.append(GTF_ATTRIBUTE_EXPECT);
					sb.append(" \"");
					for (int j = 0; j < tt.length; j++) {
						sb.append(solver.getSuperProfileMap().get(tt[j].getTranscriptID()).getAreaFrac(
								GraphLPsolver.bounds2rel(frac, tt[j].getExonicLength()- readLenMin), readLenMin, 
								strandSpecific?TProfile.DIR_BOTH:TProfile.DIR_FORWARD));
						if (j< tt.length-1)
							sb.append(",");
					}
					sb.append("\";");
					
					// profiles
					sb.append(" ");
					sb.append(GTF_ATTRIBUTE_PROFILE);
					sb.append(" \"");
					for (int j = 0; j < tt.length; j++) {
						Vector<TProfile> v= solver.getSuperProfileMap().get(tt[j].getTranscriptID()).getProfiles();
						for (int i = 0; i < v.size(); i++) {
							sb.append(v.elementAt(i).length());
							if (i< v.size()-1)
								sb.append(":");
						}
						if (j< tt.length-1)
							sb.append(",");
					}
					sb.append("\";");
					
					sb.append("\n");
	*/				
					return sb.toString();
				}
	
			private String getGTF(StringBuilder sb, Exon exon, Transcript t, Graph g, GraphLPsolver solver, boolean unsolvedSystem, 
						double perM, String pv, boolean attributesOnly) {
	
					if (!attributesOnly) {
						GFFObject obj= GFFObject.createGTFObjects(exon, t)[0];
						sb.append(obj.toString());
					}
				
			//		if (eeV.size()< 1)
			//			eeV.add(new Vector<Edge>());
			//		else
			//			eeV.elementAt(0).removeAllElements();
					Vector<Vector<Edge>> eeV= new Vector<Vector<Edge>>(5,5);
					eeV.add(new Vector<Edge>());
					
					//if (g.readCount> 0) // get lengths
					g.getRPK(exon, t, pairedEnd, Graph.ETYPE_AL, eeV);
			
					//containerIntA1A1[0][0]= g.readCount> 0? getLength(eeV.elementAt(0), null, true, false): (exon.getLength()- readLen);
					long[][] sig= new long[][]{g.encodeTset(t)};
					getGTFappend(sb, g, solver, eeV, perM, sig);
					eeV.elementAt(0).removeAllElements(); 
					
					return sb.toString();
				}
	
			private String getGTF(StringBuilder sb, Gene gene, Graph g, GraphLPsolver solver, double perM, String pv) {
				
				//clearEdgeContainer(1);
				Vector<Vector<Edge>> eeV= new Vector<Vector<Edge>>(5,5);
				eeV.add(new Vector<Edge>());
				
				GFFObject obj= GFFObject.createGFFObject(gene);
				sb.append(obj.toString());
				//if (g.readCount> 0) // for getting lengths 
				g.getRPK(gene, pairedEnd, Graph.ETYPE_AL, eeV);
				
				
				//lenExon[0][0]= g.readCount> 0? getLength(eeV.elementAt(0), null): (t.getExonicLength()- readLen);
				//containerIntA1A1[0][0]= getLength(eeV.elementAt(0), null, true, false);
				//debug= true;
				long[][] sig= new long[][]{g.createAllArray()};
				getGTFappend(sb, g, solver, eeV, perM, sig);
				//debug= false;
				return sb.toString();
			}	
			
			private String getGTF(StringBuilder sb, Transcript t, GraphLPsolver solver, Graph g, double perM, String pv, boolean attributesOnly) {
					
					GFFObject obj= GFFObject.createGFFObject(t);
					sb.append(obj.toString());
	
					Vector<Vector<Edge>> eeV= new Vector<Vector<Edge>>(5,5);
					eeV.add(new Vector<Edge>());
					//if (g.readCount> 0) // get lengths 
					g.getRPK(t, pairedEnd, Graph.ETYPE_AL, eeV);
			
					//lenExon[0][0]= g.readCount> 0? getLength(eeV.elementAt(0), null): (t.getExonicLength()- readLen);
					//containerIntA1A1[0][0]= getLength(eeV.elementAt(0), null, true, false);
					long[][] others= new long[1][];
					others[0]= Graph.without(g.createAllArray(), g.encodeTset(new Transcript[] {t}));
					
					long[][] sig= new long[][]{g.encodeTset(t)};	// containerLongA1A[0]
					getGTFappend(sb, g, solver, eeV, perM, sig);
					
					eeV.elementAt(0).removeAllElements();
					
					return sb.toString();
				}
	
			private double getGTFappend(StringBuilder sb, Graph g, GraphLPsolver solver, Vector<Vector<Edge>> eeV, double perM, long[][] tid) {
					
					invariantTestObsSplitFreq= 0; 
					invariantTestPredSplitFreq= 0;
					sb.append(" ");
					Vector<double[]> containerVecLen= new Vector<double[]>();
					for (int x = 0; x < 3; x++) {	// virtual length sum, split, uniq
						boolean output= true;
						if ((x== 0&& !outputAll)|| (x== 1&& !outputSplit)|| (x==2&& !outputUnique))
							output= false;
						if (output) {
							sb.append(FluxCapacitorConstants.GTF_ATTRIBUTE_LENGTH);
							sb.append(FluxCapacitorConstants.GTF_ATTRIBUTE_TOKEN_SEP);
							sb.append(x== 0? FluxCapacitorConstants.GTF_ATTRIBUTE_TOKEN_ALL: (x==1? FluxCapacitorConstants.GTF_ATTRIBUTE_TOKEN_TID: FluxCapacitorConstants.GTF_ATTRIBUTE_TOKEN_EXC));
							sb.append(" ");	// \"
						}
						boolean excl= x==2? true: false;
						for (int i = 0; i < tid.length; i++) {
							long[] sig= x== 0? sigall: tid[i];
							if (i>= containerVecLen.size())
								containerVecLen.add(new double[3]);
							containerVecLen.elementAt(i)[x]= getLength(g, eeV.elementAt(i), sig, excl);
							if (output) {
								sb.append(Double.toString(containerVecLen.elementAt(i)[x]));
								sb.append(",");
							}
						}
						if (output) {
							sb.deleteCharAt(sb.length()- 1);
							sb.append("; ");	// \"
						}
					}
			
					double ms= miss?factor():1d;
					for (int i = 0; i < 3; i++) { // obs, pred, norm
						boolean output1= true;
						if ((i== 0&& !outputObs)|| (i== 1&& !outputPred)|| (i== 2&& !outputBalanced))
							output1= false;
						ms= (i== 0)?Math.round(ms):ms;
						
						ReadStatCalculator calc= (i== 0? FluxCapacitor.this: solver);
						for (int j = 0; j < 3; j++) {	// sum, split, uniq
							boolean output2= true;
							if ((!output1)|| (j== 0&& !outputAll)|| (j== 1&& !outputSplit)|| (j==2&& !outputUnique))
								output2= false;
							boolean excl= j==2? true: false;
							for (int x = 0; x < 2; x++) {	// reads, coverage
								boolean output= true;
								if ((!output2)|| (x== 0&& !outputFreq)
										|| (x== 1&& !outputRfreq)
										|| (x== 2&& !outputRcov))
									output= false;
								if (output) {
									sb.append(i==0? FluxCapacitorConstants.GTF_ATTRIBUTE_TOKEN_OBSV: (i==1? FluxCapacitorConstants.GTF_ATTRIBUTE_TOKEN_PRED:FluxCapacitorConstants.GTF_ATTRIBUTE_TOKEN_BALANCED));
									sb.append(FluxCapacitorConstants.GTF_ATTRIBUTE_TOKEN_SEP);
									sb.append(j== 0?FluxCapacitorConstants.GTF_ATTRIBUTE_TOKEN_ALL:(j==1? FluxCapacitorConstants.GTF_ATTRIBUTE_TOKEN_TID: FluxCapacitorConstants.GTF_ATTRIBUTE_TOKEN_EXC));
									sb.append(FluxCapacitorConstants.GTF_ATTRIBUTE_TOKEN_SEP);				
									sb.append(x== 0?FluxCapacitorConstants.GTF_ATTRIBUTE_TOKEN_READS: FluxCapacitorConstants.GTF_ATTRIBUTE_TOKEN_COV);	// pairedEnd?GTF_ATTRIBUTE_TOKEN_COV:GTF_ATTRIBUTE_TOKEN_RPKM
									sb.append(" ");	// no \"
								}
								for (int m = 0; m < tid.length; m++) {
									long[] sig= j== 0? sigall: tid[m];
			
									if (calc== null) {
										if (output)
											sb.append(FluxCapacitorConstants.VALUE_NA);
									} else {
										
										boolean normalized= i== 2? true: false;
										double val= ms* (j== 0? calc.getReads(eeV.elementAt(m), FluxCapacitorConstants.BYTE_0, sig, normalized):  // ALL								
											calc.getReadsAvg(eeV.elementAt(m), FluxCapacitorConstants.BYTE_0, g, sig, excl, normalized));	// TID, excl 							
										
	//									if (val!= 0&& i>= 1&& j== 2&& containerVecLen.elementAt(m)[j]== 0)
	//										System.currentTimeMillis();
	//									
	//									val= ms* (j== 0? calc.getReads(eeV.elementAt(m), BYTE_0, sig, normalized):  // ALL								
	//										calc.getReadsAvg(eeV.elementAt(m), BYTE_0, g, sig, excl, normalized));	// TID, excl
	//									System.currentTimeMillis();
										
										// TODO only negatives?
										try{assert(true|| val== 0|| val>= 0.0000000001);}catch(AssertionError err){
											if (Constants.verboseLevel> Constants.VERBOSE_NORMAL)
												System.err.println("Encountered value < e-10: "+ val);
											val= 0;	
										}
											
	
										if (x== 0) { //READS
											if (j== 1) {
												if (i== 0)
													invariantTestObsSplitFreq+= val;
												else if (i== 1)
													invariantTestPredSplitFreq+= val;
											}
											// why they should be int
	//										if (i== 0)
	//											sb.append(Integer.toString((int) Math.round(val)));
	//										else
											if (output)
												sb.append(Float.toString((float) val));
										} else {	// RFREQ, COVERAGE
											double length= containerVecLen.elementAt(m)[j];
			//								if (length== 0&& val!= 0) {
			//									System.currentTimeMillis();
			//									getLength(g, eeV.elementAt(0), sig, excl);
			//									val= j== 0? calc.getReads(eeV.elementAt(m), BYTE_0, sig):  // ALL								
			//										calc.getReadsAvg(eeV.elementAt(m), BYTE_0, g, sig, excl);	// TID, excl
			//								}
											if (length== 0) {
												try {assert(val== 0);} catch (AssertionError e){
													if (Constants.verboseLevel> Constants.VERBOSE_NORMAL)
														System.err.println("Value found at 0 length: "+ val);
													val= 0;	
												};
											} else {
												val/= length;
											}
											if (output)
												sb.append(length== 0? FluxCapacitorConstants.FLOAT_STRING_0: Float.toString((float) val));
										} 
									}
									if (output)
										sb.append(",");
								}
								if (output) {
									sb.deleteCharAt(sb.length()- 1);
									sb.append("; ");	// no \"
								}
							}
						}
					}
					
					sb.append("\n"); // here?
					
					return invariantTestPredSplitFreq;
				}
	
			private void outputGFF_save(Graph g, ASEvent[] events, GraphLPsolver solver) {
						++nrLoci;
						if (solver!= null) 
							++nrLociExp;
						double perM= nrReadsAll/ 1000000d;
						// deprecated
							boolean unsolvedSystem= false;	
							double valOF= solver== null?0: solver.getValObjFunc();
							if (valOF> FluxCapacitorConstants.BIG) { 
								++nrUnsolved;
								unsolvedSystem= true;
							}
							String pv= getAttributeOF(valOF, solver, getMappedReadcount());
			
						StringBuilder sb= new StringBuilder();
						// LOCUS TODO genes 
						if (outputGene) {
							getGTF(sb, g.trpts[0].getGene(), g, solver, perM, pv);	
							try {assert(testInvariant(invariantTestObsSplitFreq, 
									pairedEnd?nrMappingsReadsOrPairs*2:nrMappingsReadsOrPairs, 0.05));}	// min: 0.01
							catch (AssertionError e) {
								if (Constants.verboseLevel> Constants.VERBOSE_NORMAL)
									System.err.println("[ASSERTION] "+getClass().getName()+".outputGFF():\n\tinvariantTestObsSplitFreq= "
											+ invariantTestObsSplitFreq+ ", nrMappingsReadsOrPairs= "+ (pairedEnd?nrMappingsReadsOrPairs*2:nrMappingsReadsOrPairs)
											+ "\n\tlocus: "+ g.trpts[0].getTranscriptID());
							};
							try {assert(testInvariant(invariantTestPredSplitFreq, 
									pairedEnd?nrMappingsReadsOrPairs*2:nrMappingsReadsOrPairs, 0.1));}
							catch (AssertionError e) {
								if (Constants.verboseLevel> Constants.VERBOSE_NORMAL)
									System.err.println("[ASSERTION] "+getClass().getName()+".outputGFF():\n\tinvariantTestPredSplitFreq= "
											+ invariantTestPredSplitFreq+ ", nrMappingsReadsOrPairs= "+ (pairedEnd?nrMappingsReadsOrPairs*2:nrMappingsReadsOrPairs)
											+ "\n\tlocus: "+ g.trpts[0].getTranscriptID());
							};
						}
			
						
						// TRANSCRIPTS
						if (outputTranscript|| outputExon|| outputSJunction) {
							float invariantObsAllTx= 0, invariantPredAllTx= 0,
								invariantObsAllEx= 0, invariantPredAllEx= 0;
							for (int i = 0; i < g.trpts.length; i++) {
								++nrTx;
			//					float invariantObsTx= invariantTestObsSplitFreq,
			//					invariantPredTx= invariantTestPredSplitFreq;
								float invariantObsTx= 0, invariantPredTx= 0;
								if (outputTranscript) {
									getGTF(sb, g.trpts[i], solver, g, perM, pv, false);	// writer.write
									invariantObsAllTx+= invariantTestObsSplitFreq; //invariantObsTx;
									invariantPredAllTx+= invariantTestPredSplitFreq; // invariantPredTx;
									invariantObsTx= invariantTestObsSplitFreq;
									invariantPredTx= invariantTestPredSplitFreq;
									if (invariantPredTx> 0)
										++nrTxExp;
								}
								// EXONS
								float invariantObsEx= 0, invariantPredEx= 0;
								if (outputExon) {
									Exon[] exons=  g.trpts[i].getExons();
									for (int j = 0; j < exons.length; j++) {
										getGTF(sb, exons[j], g.trpts[i], g, solver, unsolvedSystem, perM, pv, false);
										invariantObsEx+= invariantTestObsSplitFreq;
										invariantPredEx+= invariantTestPredSplitFreq;
									}
								}
								
								// SJ
								if (outputSJunction) {
									Vector<Vector<Edge>> eeV= new Vector<Vector<Edge>>(5,5);
									eeV.add(new Vector<Edge>());
									g.getRPK(g.trpts[i], pairedEnd, Graph.ETYPE_SJ, eeV);
									long[][] sig= new long[][]{g.encodeTset(g.trpts[i])};
									for (int j = 0; j < eeV.elementAt(0).size(); j++) { 
										getGTF(sb, eeV.elementAt(0).elementAt(j), sig, g, solver, perM);
										invariantObsEx+= invariantTestObsSplitFreq;
										invariantPredEx+= invariantTestPredSplitFreq;
									}
								}
								invariantObsAllEx+= invariantObsEx;
								invariantPredAllEx+= invariantPredEx;
								
								if (outputExon&& outputSJunction&& outputTranscript) {
									try {assert(testInvariant(invariantObsEx, invariantObsTx, 0.05));}	// min: 0.02
									catch (AssertionError e) {
										if (Constants.verboseLevel> Constants.VERBOSE_NORMAL)
											System.err.println("[ASSERTION] "+getClass().getName()+".outputGFF():\n\tinvariantObsEx= "
													+ invariantObsEx+ ", invariantObsTx= "+ invariantObsTx
													+ "\n\tlocus: "+ g.trpts[0].getTranscriptID());
									};
									try {assert(testInvariant(invariantPredEx, invariantPredTx, 0.1));}
									catch (AssertionError e) {
										if (Constants.verboseLevel> Constants.VERBOSE_NORMAL)
											System.err.println("[ASSERTION] "+getClass().getName()+".outputGFF():\n\tinvariantPredEx= "
													+ invariantPredEx+ ", invariantPredTx= "+ invariantPredTx
													+ "\n\tlocus: "+ g.trpts[0].getTranscriptID());
									};
								}
							}
							if (outputTranscript) {
								try {assert(testInvariant(invariantObsAllTx, 
										pairedEnd?nrMappingsReadsOrPairs*2:nrMappingsReadsOrPairs, 0.05));}	// min: 0.01
								catch (AssertionError e) {
									if (Constants.verboseLevel> Constants.VERBOSE_NORMAL)
										System.err.println("[ASSERTION] "+getClass().getName()+".outputGFF():\n\tinvariantObsAllTx= "
												+ invariantObsAllTx+ ", nrMappingsReadsOrPairs= "+ (pairedEnd?nrMappingsReadsOrPairs*2:nrMappingsReadsOrPairs)
												+ "\n\tlocus: "+ g.trpts[0].getTranscriptID());
								};
								try {assert(testInvariant(invariantPredAllTx, 
										pairedEnd?nrMappingsReadsOrPairs*2:nrMappingsReadsOrPairs, 0.1));}
								catch (AssertionError e) {
									if (Constants.verboseLevel> Constants.VERBOSE_NORMAL)
										System.err.println("[ASSERTION] "+getClass().getName()+".outputGFF():\n\tinvariantPredAllTx= "
												+ invariantPredAllTx+ ", nrMappingsReadsOrPairs= "+ (pairedEnd?nrMappingsReadsOrPairs*2:nrMappingsReadsOrPairs)
												+ "\n\tlocus: "+ g.trpts[0].getTranscriptID());
								};
							}
							if (outputExon&& outputSJunction) {
								try {assert(testInvariant(invariantObsAllEx, 
										pairedEnd?nrMappingsReadsOrPairs*2:nrMappingsReadsOrPairs, 0.05));}	// min: 0.02
								catch (AssertionError e) {
									if (Constants.verboseLevel> Constants.VERBOSE_NORMAL)
										System.err.println("[ASSERTION] "+getClass().getName()+".outputGFF():\n\tinvariantObsAllEx= "
												+ invariantObsAllEx+ ", nrMappingsReadsOrPairs= "+ (pairedEnd?nrMappingsReadsOrPairs*2:nrMappingsReadsOrPairs)
												+ "\n\tlocus: "+ g.trpts[0].getTranscriptID());
								};
								try {assert(testInvariant(invariantPredAllEx, 
										pairedEnd?nrMappingsReadsOrPairs*2:nrMappingsReadsOrPairs, 0.1));}
								catch (AssertionError e) {
									if (Constants.verboseLevel> Constants.VERBOSE_NORMAL)
										System.err.println("[ASSERTION] "+getClass().getName()+".outputGFF():\n\tinvariantPredAllEx= "
												+ invariantPredAllEx+ ", nrMappingsReadsOrPairs= "+ (pairedEnd?nrMappingsReadsOrPairs*2:nrMappingsReadsOrPairs)
												+ "\n\tlocus: "+ g.trpts[0].getTranscriptID());
								};
							}
						}
						
						// EVENTS
						if (outputEvent) {
							HashMap<Object,Double> tExpMap= null;
							if (solver!= null) {
								tExpMap= solver.getTrptExprHash();
								Object[] keys= tExpMap.keySet().toArray();	
								for (int i = 0; i < keys.length; i++) {
									if (!(keys[i] instanceof String))
										continue;
									if (tExpMap.get(keys[i])<0)
										tExpMap.put((String) keys[i], 0d);	// TODO ugly
								}
							}
							for (int i = 0; events!= null&& i < events.length; i++) {
								getGTF(sb, events[i], g, solver, unsolvedSystem, perM, pv, tExpMap);
							}
						}
			
						// FRAGMENTS and XJUNCTIONS
						if (false&& solver!= null) {
							ArrayList<Edge> cc= new ArrayList<Edge>();
							if (solver!= null) {
								Iterator<Object> iter= solver.getConstraintHash().keySet().iterator();
								while (iter.hasNext()) {
									Object o= iter.next();
									if (o instanceof Edge)
										cc.add((Edge) o);
								}
							}
							Collections.sort(cc, Edge.getDefaultPositionComparator());
			
							Iterator<Edge> iter= cc.iterator();
							while (iter.hasNext()) {
								Edge e= iter.next();
								// no INTRONS
								if ((!(e instanceof SuperEdge))&& (!e.isExonic()))
									continue;
								getGTF(sb, e, new long[][]{e.getTranscripts()}, g, solver, perM);
							}
						}
							
						try {
							write(sb);
						} catch (Exception e) {
							e.printStackTrace();
						}
						
					}
	
			
			private int[] extend(Transcript tx, BEDobject2[] beds, int[] extension) {
				
				// assert beds sorted by gpos !
				int tstart= tx.get5PrimeEdge(), tend= tx.get3PrimeEdge();
				if (tx.getStrand()< 0) {
					int h= Math.abs(tstart); tstart= Math.abs(tend); tend= h;
				}

				HashMap<CharSequence, BEDobject2> mapMates5= new HashMap<CharSequence, BEDobject2>(),
					mapMates3= new HashMap<CharSequence, BEDobject2>();
				if (tx.getStrand()< 0)
					System.currentTimeMillis();

				// 5' extension
				int i= 0;
				for (; i < beds.length; i++) {
					if (beds[i].getStart()+ 1>= tstart)
						break;
					if (beds[i].getStrand()< 0|| beds[i].getBlockCount()> 1)	// only non-split
						continue; // only reads on forward strand can have a mate falling into transcript
					CharSequence id= descriptor.getUniqueDescriptor(beds[i].getName());	//Fasta.getReadID(beds[i].getName())
					mapMates5.put(id, beds[i]);
				}
				int left= i;
				// 3' extension
				for(i= beds.length- 1; i>= left; --i) {
					if (beds[i].getEnd()<= tend)
						break;
					if (beds[i].getStrand()>= 0|| beds[i].getBlockCount()> 1)	// only non-split
						continue; // only reads on reverse strand can have a mate falling into transcript
					CharSequence id= descriptor.getUniqueDescriptor(beds[i].getName());	//Fasta.getReadID(beds[i].getName())
					mapMates3.put(id, beds[i]);
				}
				int right= i;
				// find mates				
				int min= 0, max= 0;
				for(i= left; i<= right; ++i) {
					boolean contained= contains(tx, beds[i]);
					if (!contained)
						continue;
					CharSequence id= descriptor.getUniqueDescriptor(beds[i].getName());
					BEDobject2 bed= null;
					if (beds[i].getStrand()< 0) {
						bed= mapMates5.get(id);
						if (bed== null)
							continue;
						min= Math.min(min, bed.getStart()+ 1- tstart);
					} else {	// strand>= 0
						bed= mapMates3.get(id);
						if (bed== null)
							continue;
						if (bed.getEnd()- tend> 75)	// TODO
							System.currentTimeMillis();
						max= Math.max(max, bed.getEnd()- tend);
					}
				}
				// TODO additional rescue for pairs that extend on both sides
				
				extension[0]= Math.abs(min);
				extension[1]= max;
				if (tx.getStrand()< 0)
					System.currentTimeMillis();
				if (extension[0]> 0|| extension[1]> 0)	// TODO
					System.currentTimeMillis();		
				if (tx.getStrand()< 0) {
					int h= extension[0]; extension[0]= extension[1]; extension[1]= h;
				}
				return extension;
			}
			
			/**
			 * @param tx
			 * @param bed
			 * @return
			 */
			private boolean contains(Transcript tx, BEDobject2 bed) {
				
//				if (bed.toString().contains("HWUSI-EAS1692:3:30:7649:19668"))
//					System.currentTimeMillis();
				
				Exon[] exons= tx.getExons();
				boolean tsense= tx.getStrand()>= 0;
				int idx= (tsense?0:exons.length- 1), gstart= bed.getStart()+ 1;
				
				// find first exon
				int dstart= tsense?gstart:-gstart;
				for (;idx>= 0&& idx< exons.length;idx+= tsense? 1: -1) {
					if (exons[idx].contains(dstart))
						break;
				}
				if (idx< 0|| idx>= exons.length)
					return false;
				if (bed.getBlockCount()< 2) {
					int dend= tsense?bed.getEnd():-bed.getEnd();
					return exons[idx].contains(dend);
				}
				
				for (int i = 0; i < bed.getBlockCount(); ++i, idx+= tsense? 1: -1) {
					if (idx< 0|| idx>= exons.length)
						return false;
					int bstart= bed.getNextBlockStart()+ gstart,
						bend= bstart+ bed.getNextBlockSize(),
						dend= tsense? bend: -bend;
					dstart= tsense? bstart: -bstart;
						
					if (!(exons[idx].contains(dstart)&& exons[idx].contains(dend)))
						return false;
				}
				return true;
			}
			
			
			private void learn(Transcript tx, BufferedBEDiterator beds) {
							
				if (beds== null)
					return;
				
				// pre-filter reads for non-valid split-maps (null)
				Graph myGraph= getGraph(this.gene);
				map(myGraph, this.gene, beds);
				// DEBUG
				HashMap<CharSequence, CharSequence> map= new HashMap<CharSequence, CharSequence>();
				
				int elen= tx.getExonicLength();	// TODO this is the effective length
//				if (elen< readLenMin)
//					return;	// discards reads
				
				HashMap<CharSequence, BEDobject2[][]> mapPends= null;
				int[] extension= new int[2];	// 5' extension, 3' extension
				extension[0]= 0; extension[1]= 0;
				int extLen= elen;
				if (pairedEnd) {
					mapPends= new HashMap<CharSequence, BEDobject2[][]>();
					//extension= extend(tx, beds, extension);
					extLen+= extension[0]+ extension[1];
				}
				//float rpk= beds.length* 1000f/ elen; 
				UniversalMatrix m= profile.getMatrix(extLen);
				
				nrReadsSingleLoci+= beds.length;
				for (int i = 0; i < beds.length; i++) {
					
					BEDobject2 bed1= beds[i];
					if (bed1== null)
						continue;
					// unique only, bad idea
					//--too little and does not eliminate peaks
//					if (bed1.getScore()> 1)
//						continue;
					CharSequence tag= bed1.getName();
					attributes= descriptor2.getAttributes(tag, attributes);						
					if (stranded) {
						if ((tx.getStrand()== bed1.getStrand()&& attributes.strand== 2)
								|| (tx.getStrand()!= bed1.getStrand()&& attributes.strand== 1))
						++nrMappingsWrongStrand;
						continue;
					}
					
					int bpoint1= getBpoint(tx, bed1);					
//					if (m.sense.length== 1250&& bpoint1== 303&& bed1.getScore()<= 1)
//						System.currentTimeMillis();
					if (bpoint1== Integer.MIN_VALUE) {	// was intron
						++nrReadsSingleLociNoAnnotation;
						continue; 	// doesnt align to the transcript
					}
					bpoint1+= extension[0];
					if(bpoint1< 0|| bpoint1>= extLen) {	// outside tolerated area
						continue;
					}

					int rlen1= bed1.getLength();
					if (readLenMin< 0|| rlen1< readLenMin)
						readLenMin= rlen1;
					if (rlen1> readLenMax) {	// readLenMax< 0|| 
						readLenMax= rlen1;
					}
					
					++nrReadsSingleLociMapped;
					if (pairedEnd) {

						byte flag= attributes.flag;  
							// (byte) (descriptor.getPairedEndInformation(bed1.getName())- 1);	// (Fasta.getReadDir(beds[i].getName())- 1))
						if (flag < 1) {
							System.err.println("\n\tRead ignored, error in readID:\n"+ tag);
							continue;
						}
						--flag; // for array index
						CharSequence id= attributes.id; 	// descriptor.getUniqueDescriptor(bed1.getName());	//Fasta.getReadID(beds[i].getName())
						//CharSequence id= beds[i].getName();
						BEDobject2[][] oo= mapPends.get(id);
						if (oo== null) {
							oo= new BEDobject2[2][];
							mapPends.put(id, oo);
						} 
						if (oo[flag]== null) 
							oo[flag]= new BEDobject2[] {bed1};
						else {
							BEDobject2[] op= new BEDobject2[oo[flag].length+ 1];
							System.arraycopy(oo[flag], 0, op, 0, oo[flag].length);
							op[op.length- 1]= bed1;
							oo[flag]= op;
						}
						// for profiles paired reads only
						boolean found= false;
						boolean dbg= false;
						for (int j = 0; j < oo.length; j++) {
							if (j==flag|| oo[j]== null)
								continue;
							for (int k = 0; k < oo[j].length; k++) {
								BEDobject2 bed2= oo[j][k];
								if (bed2== null)
									continue;
								// unique only, bad idea
								//--too little and does not eliminate peaks
//								if (bed2.getScore()> 1)
//									continue;
								int bpoint2= getBpoint(tx, bed2);
//								if (m.sense.length== 1250&& bpoint2== 303&& bed2.getScore()<= 1)
//									System.currentTimeMillis();
								bpoint2+= extension[0];
								if (bpoint2>= 0&& bpoint2< extLen) {	// inside tolerated area
									int rlen2= bed2.getLength();
//									if (tx.getStrand()< 0)
//										System.currentTimeMillis();
									m.add(bpoint1, bpoint2, rlen1, rlen2, extLen);
									addInsertSize(Math.abs(bpoint2- bpoint1)+ 1);
									
									// count each pair only once, TODO first hit counted
									map.put(bed1.getName(), null);	
									map.put(bed2.getName(), null);
										
									break;
								} 
							}
							if (found)
								break;
						}
					
					} else {
						
						m.add(bpoint1, rlen1, extLen, bed1.getStrand()== tx.getStrand()?Constants.DIR_FORWARD:Constants.DIR_BACKWARD);
					}
				}
				
				nrReadsSingleLociPairsMapped+= map.size()/ 2;
				
				// DEBUG
				if (map.size()!= mapReadOrPairIDs.size()) {
					System.err.println();
					Object[] ooo= map.keySet().toArray();
					Arrays.sort(ooo);
					for (int i = 0; i < ooo.length; i++) {
						if (!mapReadOrPairIDs.contains(ooo[i]))
							System.err.println(ooo[i]);
					}
					
					System.err.println();
					ooo= mapReadOrPairIDs.toArray();
					Arrays.sort(ooo);
					for (int i = 0; i < ooo.length; i++) {
						if (!map.containsKey(ooo[i]))
							System.err.println(ooo[i]);
					}
					
					System.currentTimeMillis();
				} else
					System.currentTimeMillis();
				

			}


			/**
					 * add a SINGLE read
					 * @param regs
					 */
					int mapRead(Graph g, BEDobject2 dobject, boolean force) { 
						// find the edge(s) where the regions align
						// if these do not form a continous chain, create a new edge
						
			//			GFFObject[] gtfs= GFFObject.fromBed(dobject);	// TODO kill intermediary GTFs !!!
			//			DirectedRegion[] regs= new DirectedRegion[gtfs.length];
			//			for (int i = 0; i < regs.length; i++) 
			//				regs[i]= new DirectedRegion(gtfs[i]);
						
						if (force&& mapReadOrPairIDs.contains(dobject.getName())) {
							return 0;
						}
						
						byte flag= getFlag(dobject);  	
						CharSequence ID= getID(dobject); 	
			
						Edge target= g.getEdge(dobject);
						
						if ((target== null)|| ((!pairedEnd)&&  (!(target instanceof SuperEdge))
								&& target.length()< readLenMin))
							return 0;
						if (force) {
							boolean sense= g.trpts[0].getStrand()== dobject.getStrand();
							if (sense)
								target.incrReadNr();
							else
								target.incrRevReadNr();
							mapReadOrPairIDs.add(dobject.getName());
							return 1;
						}
						
						
						byte refStrand= g.trpts[0].getStrand();
						boolean sense= dobject.getStrand()== refStrand;
						byte antiflag= (byte) ((flag==1)?2:1);
						int mapCtr= 0;
						
						
						// add first/single read
						if (pairedEnd) { /* PAIRED END */
			
							//int mappedIDsBefore= mapReadOrPairIDs.size();
							Vector<BEDobject2>[] vv= mapEndsOfPairs.get(ID);
							Vector<BEDobject2> v= null;
							if (vv!= null)
								v= vv[antiflag- 1];
							for (int i = 0; v!= null
								&&i < v.size(); i++) {
								
								BEDobject2 dobject2= v.elementAt(i);
								Edge target2= g.getEdge(dobject2);
								if (target2== null)
									continue;
			
								// check whether they map within isize constraints
								int j = 0;
								for (; target.getSuperEdges()!= null&& j < target.getSuperEdges().size(); j++) {
									if (!target.getSuperEdges().elementAt(j).isPend())
										continue;
									Edge[] ee= target.getSuperEdges().elementAt(j).getEdges();
									int k = 0;
									//for (; k < ee.length&& ee[k]!= target2; k++);	// TODO binarySearch
									// for the case target== target2, better check that there is no other edge
									for (; k < ee.length; k++) {
										if (ee[k]!= target&& ee[k]!= target2)
											break;
									}
										
									if (k== ee.length)
										break;	// common superedge found@ j
								}
								SuperEdge se= null;
								if (target.getSuperEdges()!= null&& j< target.getSuperEdges().size()) 
									se= target.getSuperEdges().elementAt(j);
								else {
									continue;	// not possible paired-end
								}
			
								se.incrReadNr();
								++mapCtr;
								
								try {
									testWriter.write(dobject.toString()+ "\n");
									testWriter.write(dobject2.toString()+ "\n");
								} catch (Exception e) {
									e.printStackTrace();
								}
								
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
			
							if (sense|| (strand!= FluxCapacitorConstants.STRAND_ENABLED)) {
								target.incrReadNr();
								mapCtr= 1;
							} else if (strand!= FluxCapacitorConstants.STRAND_SPECIFIC) {
								target.incrRevReadNr();
								mapCtr= 1;
							} else {
								++nrMappingsWrongStrand;
								mapCtr= 0;
							}
							
							
							
							if (!mapReadOrPairIDs.add(dobject.getName()))
								++nrLocusMultimaps;
			//				if (outputMapped)
			//					writeMappedRead(dobject);
							return mapCtr;
						}
						
					}
		}

	static void printUsage() {
		if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) {
			
			Iterator<String[]> iter= FluxCapacitorConstants.cliExplMap.keySet().iterator();
			int max= 0;
			String[][] ss= new String[FluxCapacitorConstants.cliExplMap.size()][];
			for (int i = 0; iter.hasNext(); i++) {
				ss[i]= iter.next();
				int x= 0;
				for (int j = 0; j < ss[i].length; j++) 
					x+= ss[i][j].length();
				if (x> max)
					max= x;
			}
			
			int tabDist= 4;
			int maxWidth= 60;
			int sep= max+ 3+ tabDist;
			Arrays.sort(ss, new StringArrayByFirstComparator());

			System.err.println("A summary of currently supported command line flags.");
			// + " For a more detailed explanation, see\nhttp://fluxcapacitor.wikidot.com/capacitor:usage\n");
			for (int i = 0; i < ss.length; i++) {
				StringBuilder sb= new StringBuilder("[");
				for (int j = 0; j < ss[i].length; j++) { 
					sb.append(ss[i][j]);
					if (j< ss[i].length- 1)
						sb.append("|");
				}
				sb.append("]");
				int pos= sep- sb.length();
				for (int j = 0; j < pos; j++) 
					sb.append(" ");
				
				String expl= FluxCapacitorConstants.cliExplMap.get(ss[i]);
				pos= sep;
				for (int j = 0; j < expl.length(); j++, pos++) {
					if (pos>= maxWidth) {
						for (; j < expl.length()&&
								(!Character.isWhitespace(expl.charAt(j))); j++) 
							sb.append(expl.charAt(j));
						
						sb.append("\n");
						for (int m = 0; m < sep; m++) 
							sb.append(" ");
						pos= sep;
						continue;
					}
					sb.append(expl.charAt(j));
					if (expl.charAt(j)== '\n') {
						for (int m = 0; m < sep; m++) 
							sb.append(" ");
						pos= sep;
					}
				}
				System.err.println(sb.toString());
			}
			
		}
	}
	
	public static void main(String[] args) {
		
		Execute.initialize(2);
		
		try {
			readProperties(); 
			
			boolean showGUI= false;
			if (args== null|| args.length== 0) {
				if (showGUI) {
					FluxCapacitor.loadLibraries();
					//FluxCapacitorGUI.createGUI();
				} else
					printUsage();
			}
			
			final FluxCapacitor myCapacitor= new FluxCapacitor();
			if (args!= null&& args.length== 1) {
				if (args[0].equalsIgnoreCase("--install")) {
					install();
					System.exit(0);
				} else if (args[0].equalsIgnoreCase("--help")) {
					printUsage();
					System.exit(0);
				} else
					myCapacitor.init2(args);
			}
			
//			myCapacitor.init2(args);
//			if (myCapacitor.isHelpRequested()) {
//				printUsage();
//				System.exit(0);
//			}	
//			wellcome();
//			if (doInstall) {
//				install();
//				System.exit(0);
//			}
			
//			if (!myCapacitor.checkPreliminaries())
//				System.exit(-1);
			
			if (!cheatDisableCleanup) {
				Runtime.getRuntime().addShutdownHook(new Thread("MrProper") {
				    public void run() { 
				    	FileHelper.cleanup(System.getProperty(Constants.PROPERTY_TMPDIR), 
				    			Constants.globalPfx== null?FluxCapacitorConstants.PFX_CAPACITOR+ "."+ myCapacitor.getRunID():Constants.globalPfx+ "_",
				    			null,
				    			Constants.verboseLevel> Constants.VERBOSE_SHUTUP?System.err:null); 
				    }
				});
			} else
				System.err.println("[NOCLEAN] Cleanup disabled!");
			
			// found 1 temporary files with prefix capacitor,...
/*			if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) {
				File dir= new File(System.getProperty(Constants.PROPERTY_TMPDIR));
				String[] fNames= dir.list();
				Vector<File> v= new Vector<File>(), vSort= new Vector<File>();
				for (int i = 0; i < fNames.length; i++) {
					File f= new File(dir+ File.separator+ fNames[i]);
					if (fNames[i].contains(PFX_CAPACITOR))
						v.add(f);
					else if (fNames[i].contains("sort"))
						vSort.add(f);
				}
				
				if (v.size()> 0&& !myCapacitor.force) 
					removeZombies(v, PFX_CAPACITOR);
				if (vSort.size()> 0&& !myCapacitor.force) 
					removeZombies(v, "sort");
				
			}
*/			
			
			int ok= loadLibraries();
			if (ok< 0) 
				exit(-1);
			
			// check
			ok= myCapacitor.checkParameter();
			if (ok< 0) 
				exit(-1);
	
			
			if (Constants.verboseLevel>= Constants.VERBOSE_NORMAL) 
				myCapacitor.printStats(System.err, args);
			
		    // run
			try {
				myCapacitor.run();
			} catch (Throwable XXX) {
				XXX.printStackTrace();
				System.exit(0);
			}
		} catch (Throwable t) {
			if (t instanceof Exception)
				((Exception) t).printStackTrace();
			else if (t instanceof Error)
				((Error) t).printStackTrace();
			else
				System.err.println(t.getMessage());
			if (cheatDoNotExit) { 
				int in= 0;
				while (in!= '\n')
					try {
						in= System.in.read();
					} catch (IOException e) {
						e.printStackTrace();
					}
			}
			
		} finally {
			Execute.shutdown();
		}
	}

	private static void removeZombies(Vector<File> v, String pfx) {
		System.err.println("[ZOMBIE] found "+ v.size()+ " temporary files with prefix "+pfx+"," +
		"\n\tdo you want to remove them (Yes/No/Don't know):");
		boolean yesNo= true;
		if (yesNo) {
			int cnt= 0, cntFail= 0;
			for (int i = 0; i < v.size(); i++) {
				File f= v.elementAt(i);
				boolean failed= false;
				if (f.isDirectory()) {
					if (!FileHelper.rmDir(f)) 
						failed= true;
				} else if (!f.delete()) {
					failed= true;
				}
				if (failed) 
					++cntFail;
				else 
					++cnt;
			}
			if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
				System.err.println("\tremoved "+cnt+" files, failed to remove "
						+ cntFail+ "files.");
		} else
			System.exit(-1);

	}

	public static String[] DEFAULT_PE_SFX= new String[] {"_1", "_2"};
	
	public File fileBED= null, 
		fileGTF= null, 
		fileOut= null, 
		fileOutDir= null, 
		fileOUToriginal= null, 
		fileBEDoriginal= null, 
		fileGTForiginal= null,
		fileMappedReads= null,
		fileNotmappedReads= null,
		fileProfile= null,
		fileLPdir= null,
		fileISize= null;
	
	int readLenMin= 75, readLenMax= -1;
	
	public FluxCapacitor() {		
	}
	
	public int getReadLength() {
		return readLenMin;
	}
	
	static int[] parseInsertSize(String s) {
		int[] a= new int[2];
		String[] ss= s.split(",");
		try {
			a[0]= Integer.parseInt(ss[0]);
			a[1]= Integer.parseInt(ss[1]);
		} catch (Exception e) {
			return null;
		}
		return a;
	}
	
	public static String version= null;
	
	
	/**
	 * @deprecated 
	 * delegated linecount to BEDwrapper
	 * readLength and max/min insert size to profiling
	 * @return
	 */
	boolean prescan() {
		
		long t0= System.currentTimeMillis();
		
		if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
			System.err.println("[PRESCAN] Checking read length(s)");
		nrReadsAll= 0;
		int[] lenMinMax= new int[] {Integer.MAX_VALUE, Integer.MIN_VALUE};
        Log.progressStart("progress ");
		try {
			BufferedReader buffy= new BufferedReader(new FileReader(fileBED));
			File f= File.createTempFile(FluxCapacitorConstants.PFX_CAPACITOR, "prescanIDs");
			BufferedWriter writer= new BufferedWriter(new FileWriter(f));
			long bRead= 0, bTot= fileBED.length();
			int perc= 0, lines= 0;
			for(String s; (s= buffy.readLine())!= null;++nrReadsAll,bRead+= s.length()+1, ++lines) {
                Log.progress(bRead, bTot);
				if (s.charAt(0)== '#')
					continue;

				int i = 0, cnt= 0;

				// check ID: untested, not here, too much overhead
/*				for(++i; s.charAt(i)!= '\t'&& i< s.length(); ++i);	// 3rd sep
				int nameStart= ++i;
				if (i!= s.length()) {
					++cnt;
					for(; s.charAt(i)!= '\t'&& i< s.length(); ++i);	// 4th sep
					if (i!= s.length())
						++cnt; 
				}
				String name= null;
				if (cnt>= 3&& i!= nameStart)
					name= s.substring(nameStart,i);
				
				if (pairedEnd) {
					if (name!= null) {
						int flag= 0;
						if ((flag= FMRD.getPE(name))!= 0) {
							++nrMappingsValid;
							if (flag== 1)
								++nrMappingsP1;
							else if (flag== 2)
								++nrMappingsP2;
						}
						if (i!= s.length())
							++cnt;
					}
				} else {
					if (i!= s.length())
						++nrMappingsValid;
				}
				if (cnt>= 3&& i!= nameStart) {
					writer.write(name);
					writer.write("\n");
				}
*/				
				// get length
				int nowC= cnt;
				for (; i < s.length()&& cnt< (10- nowC); i++) {	// 11 field bsize
					if (s.charAt(i)== '\t')
							++cnt;
				}
				if (cnt== 10) {
					cnt= 0;	// sum
					while(true) {
						int last= i;
						for (; i < s.length()&& s.charAt(i)!= ','&& s.charAt(i)!= '\t'; ++i);
						cnt+= Integer.parseInt(s.substring(last, i));
						if (s.charAt(i)== '\t')
							break;
						else
							++i;
					}
				} else {
					i= 0; cnt= 0;
					int last= 0, start= 0;
					for (; cnt< 3&& i < s.length(); i++) {
						if (s.charAt(i)== '\t') {
							++cnt;
							if (cnt== 2)
								start= Integer.parseInt(s.substring(last, i));
							else if (cnt== 3) {
								cnt= Integer.parseInt(s.substring(last, i))- start;
								break;
							}
							last= i+1;
						}
					}
				}
				if (cnt< lenMinMax[0])
					lenMinMax[0]= cnt;
				if (cnt> lenMinMax[1])
					lenMinMax[1]= cnt;
			}
			buffy.close();
			writer.flush();
			writer.close();
            Log.progressFinish();
		} catch (Exception e) {
            Log.progressFailed("ERROR");
            Log.error("Error : " + e.getMessage(), e);
			return false;
		}
				

        Log.message("\ttook "+((System.currentTimeMillis()- t0)/ 1000)+" sec");
        Log.message("\tfound "+nrReadsAll+" lines");
        Log.message("\treadlength min "+lenMinMax[0]+", max "+ lenMinMax[1]+"\n");

		readLenMin= lenMinMax[0];
		
		return true;
	}

	int init(Method m, String[] args, int p) {
		
		Object res= null;
		try {
			if (m.getParameterTypes()!= null&& m.getParameterTypes().length> 0) {
				if (p+1>= args.length) {
					if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) 
						System.err.println("[HEY] missing parameter for "+ args[p]);
					return -1;
				}
				res= m.invoke(this, args[++p]);
			} else
				res= m.invoke(this, null);
		} catch (Exception e) {
			System.err.println("[OLALA] Could not set parameter "+ args[p]);
			return -1;
		}

		if (res!= null&& res instanceof Boolean&& !((Boolean) res).booleanValue()) {
			if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) 
				System.err.println("[WAAAA] Initialiation failed for parameter "+ args[p]);
			return -1;
		}
		return p;
	}
	
	int init(String[] args) {
		
		String s;
		for (int i = 0; i < args.length; i++) {
			if (args[i].startsWith(FluxCapacitorConstants.CLI_LONG_PFX)) {
				s= args[i].substring(FluxCapacitorConstants.CLI_LONG_PFX.length());
				if (! FluxCapacitorConstants.cliLongMap.containsKey(s)) {
					if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) 
						System.err.println("[OHNO] Unrecognized long option \'"+args[i]+"\'");
					Iterator<String> iter= FluxCapacitorConstants.cliLongMap.keySet().iterator();
					while (iter.hasNext())
						System.err.println("\'"+iter.next()+"\'");
					System.exit(-1);
				}
				Method m= FluxCapacitorConstants.cliLongMap.get(s);
				if ((i= init(m, args, i))< 0) 
					System.exit(-1);
				
			} else if (args[i].startsWith(FluxCapacitorConstants.CLI_SHORT_PFX.toString())){
				s= args[i].substring(1);
				for (int j = 0; j < s.length(); j++) {
					if (! FluxCapacitorConstants.cliShortMap.containsKey(s.charAt(j))) {
						if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) 
							System.err.println("[OHNO] Unrecognized short option \'"+args[i]+"\'");
						System.exit(-1);
					}
					Method m= FluxCapacitorConstants.cliShortMap.get(s.charAt(j));
					if ((i= init(m, args, i))< 0) 
						System.exit(-1);
				}
			}
		}
		
		return 0;
	}
	
	private int checkParameter() {
		
		if (fileBEDoriginal== null|| !fileBEDoriginal.exists()) {
			if (Constants.verboseLevel!= Constants.VERBOSE_SHUTUP) {
				System.err.println("[AIII] I need a input file with aligned reads in order to work.");
				if (fileBED== null)
					System.err.println("\tyou said nothing and this is bad.");
				else
					System.err.println("\tyou said it is "+fileBED+" but it is bad.");
				System.err.println("\tUse the "+FluxCapacitorConstants.CLI_LONG_PFX+FluxCapacitorConstants.CLI_LONG_SRA+" parameter and give me a correct one please.\n");
			}
			return -1;
		}
		if (fileGTForiginal== null|| !fileGTForiginal.exists()) {
			if (Constants.verboseLevel!= Constants.VERBOSE_SHUTUP) {
				System.err.println("[MANO] I need a reference file with the transcripts you want to analyze.");
				if (fileGTF== null)
					System.err.println("\tyou said nothing and this is not good.");
				else
					System.err.println("\tyou said it is "+fileGTF+" but it is not good.");
				System.err.println("\tTry again, using the "+FluxCapacitorConstants.CLI_LONG_PFX+FluxCapacitorConstants.CLI_LONG_REF+" parameter.\n");
			}
			return -1;
		}
		if (fileOut!= null&& fileOut.exists()) {
			if (Constants.verboseLevel>= Constants.VERBOSE_NORMAL) {
				System.err.println("[OUCH] There is already a file at the output location "+fileOut+".");
				System.err.println("\tDo you want me to touch that? (Yes/No/Don't know)");
//				To confirm, write \'easy\' <CR>");
//				String[] expected= new String[] {"easy", "hooo-hooo", "ho-ho-hooo"},
//					proposed= new String[] {
//						"\tSay \'hooo-hooo\'!",
//						"\tSay \'ho-ho-hooo\'!",
//						"\tDo you always everything they tell you? A simple \'yes\' would have been enough."
//				};
//				for (int j = 0; j < expected.length; j++) {
//					String s= readSystemIn().trim();
//					if (s.equalsIgnoreCase("yes")|| s.equalsIgnoreCase("y")) {
//						break;
//					}
//					if (s.equalsIgnoreCase(expected[j])) {
//						System.err.println(proposed[j]);
//					} else {
//						System.err.println("\tI did not get that. Try again.");
//						--j;
//					}
//				}
				boolean yes= true;
				if (!force)
					yes= waitForYesNo();
				if (yes) {
					fileOut.delete();
					if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
						System.err.println("[BUM] I permanently deleted "+fileOut+".\n");
				} else
					System.exit(-1);
			} else 
				return -1;
		}
		
		return 0;
	}
	
	private void printStats(PrintStream p, String[] args) {
		p.println("\n[HEHO] We are set, so let's go!");
		p.print("\tcmd\t"+FluxCapacitorConstants.CLI_CMD);
		for (int i = 0; i < args.length; i++) 
			p.print(" "+args[i]);
		p.println();
		
		try {
			// INPUT
			p.println("\tINPUT");
			p.println("\t"+FluxCapacitorParameters.PAR_ANNOTATION_FILE+"\t"+fileGTForiginal.getCanonicalPath());
			p.println("\t"+FluxCapacitorParameters.PAR_MAPPING_FILE+"\t"+fileBEDoriginal.getCanonicalPath());
			p.println("\tdescriptor\t"+descriptor2);
			//p.println("\t"+CLI_LONG_VERBOSE+"\t"+Constants.VERBOSE_KEYWORDS[Constants.verboseLevel]);
			if (copyLocal)
				p.println("\t"+FluxCapacitorParameters.PAR_COPY_INPUT);
			
			// OUTPUT
			p.println("\tOUTPUT");
			p.println("\tTemporary Folder\t"+ System.getProperty(Constants.PROPERTY_TMPDIR));			
			if (Constants.globalPfx!= null)
				p.println("\t"+FluxCapacitorConstants.CLI_LONG_TPX+"\t"+ Constants.globalPfx);
			p.print("\tQuantification File\t");
			if (fileOut== null)
				p.println("stdout");
			else {
				p.println(fileOut.getCanonicalPath());
				if (compressionOut!= FileHelper.COMPRESSION_NONE)
					p.println("\t"+ FluxCapacitorConstants.CLI_LONG_COMPRESSION+ "\t"+ FileHelper.COMPRESSION_KEYWORDS[compressionOut]);
			}
/*			p.print("\tfeatures:\t");
			if (outputExon)
				p.print("Exons ");
			if (outputSJunction)
				p.print("splice-Junctions ");
			if (outputTranscript)
				p.print("Transcripts ");
			if (outputGene)
				p.print("Genes ");
			if (outputEvent)
				p.print("eVents ");
			p.println();
*/
/*			p.print("\tbases:\t");
			if (outputObs)
				p.print("Observed ");
			if (outputPred)
				p.print("preDicted ");
			if (outputBalanced)
				p.print("Balanced ");
			p.println();
*/			
/*			p.print("\tscopes:\t");
			if (outputAll)
				p.print("All ");
			if (outputSplit)
				p.print("Split ");
			if (outputUnique)
				p.print("Unique ");
			p.println();
*/			
/*			p.print("\tmeasure:\t");
			if (outputFreq)
				p.print("Freq ");
			if (outputSplit)
				p.print("Split ");
			if (outputUnique)
				p.print("Unique ");
			p.println();
*/
/*			if (outputMapped|| outputNotmapped|| outputProfiles|| outputLP) {
				p.print("\tsave:\t");
				if (outputMapped)
					p.print("Mapped-alignments ");
				if (outputNotmapped)
					p.print("Notmapped-alignments ");
				if (outputProfiles)
					p.print("Profiles ");
				if (outputLP)
					p.print("Linear-programs ");
				p.println();
			}
*/
			// ALGORITHM
			p.println("\tADDITIONAL");
			//p.println("\t"+ CLI_LONG_THREAD+" "+ maxThreads);
			p.print("\tRead Distribution\t");
			if (uniform)
				p.println("uniform");
			else if (fileProfile!= null&& fileProfile.exists())
				p.println("from profiles in "+ fileProfile.getAbsolutePath());
			else {
				p.print("profiling is carried out");
				if (fileProfile!= null)
					p.println(" and stored in "+ fileProfile.getAbsolutePath());
				else
					p.println();
			}
			
			if (stranded)
				p.println("\tstrand information considered.");
			if (pairedEnd)
				p.println("\tmate pairing information considered");
			

/*			p.print("\t"+ CLI_LONG_COST_MODEL+" "+GraphLPsolver.COSTS_NAMES[costModel]);
			if (!Double.isNaN(costModelPar))
				p.print(" "+Double.toString(costModelPar));
			p.println();
			p.println("\t"+ CLI_LONG_COST_SPLIT+ " "+costSplit);
			p.print("\t"+ CLI_LONG_COST_BOUNDS+ " ");
			if (costBounds== null)
				p.println(" none.");
			else
				p.println(costBounds[0]+","+costBounds[1]);
			//p.println("\tread length " + readLen);
			p.println("\t"+ CLI_LONG_STRAND+" "+ strandSpecific);
*/
			
		} catch (IOException e) {
			; // :)
		}
		
//		if (pairedEnd)
//			p.println("\t"+CLI_LONG_PAIR+"\t"+insertMinMax[0]+","+insertMinMax[1]);
		//System.err.println("\t"+CLI_LONG_NOISE+"\t"+Float.toString(1- GraphLPsolver.min_read_rest_frac));
		
		p.println();

	}
	
	private final static Random rndLuser= new Random();
	private String errorMissingArgument(String string) {
		return FluxCapacitorConstants.L_USER_COMMENTS[rndLuser.nextInt(FluxCapacitorConstants.L_USER_COMMENTS.length)]
		       + "You forgot to give me an argument for the parameter "
		       + string+ "!";
		
	}

	private static String readSystemIn() {
		StringBuffer sb= new StringBuffer();
		int in;
		while(true) {
			try {
				while((in= System.in.read())!= '\n') {
					sb.append((char) in);
				}
				return sb.toString();
			} catch (Exception e) {
				; // :)
			}
		}
	}
	
	public boolean setFileReads(String path) {
		File f= new File(path);
		if (f.exists()&& !f.isDirectory()) {
			fileBED= f;
			return true;
		}
		return false;
	}
	
	public boolean setFileReference(String path) {
		File f= new File(path);
		if (f.exists()&& !f.isDirectory()) {
			fileGTF= f;
			return true;
		}
		return false;
	}
	
	public void setForce() {
		force= true;
	}
	
	public boolean setNameOutDir(String path) {
		File f= new File(path);
		if (f.exists()&& f.isDirectory()&& f.canWrite()) {
			fileOutDir= f;
			return true;
		}
		return false;

	}
	
	static boolean doInstall= false;
	public static void setInstall() {
		doInstall= true;
	}
	
	static boolean helpRequested= false;
	public static void setHelp() {
		helpRequested= true;
	}
	
	public static boolean isHelpRequested() {
		return helpRequested;
	}
	
	public static boolean setVerbose(String verboseLevel) {
		return Constants.setVerbose(verboseLevel);
	}
	
	static File fileJVMdir;
	public static boolean setJVM(String path) {
		File f= new File(path);
		if (f.exists()&& f.isDirectory()) {
			fileJVMdir= f;
			return true;
		}
		return false;
	}
	
	static File fileLibDir;
	public static boolean setLib(String path) {
		File f= new File(path);
		if (f.exists()&& f.isDirectory()) {
			fileLibDir= f;
			//System.setProperty(CLI_LONG_LIB, f.getAbsolutePath());
			return true;
		}
		return false;
	}
	
	public void setBatch() {
		Constants.verboseLevel= Constants.VERBOSE_SHUTUP;
	}
	
	public void setCompression(String comp) {
		compressionOut= FileHelper.getCompression(comp);
	}

	public void setLocal() {
		copyLocal= true;
	}
	
	public void setTempDir(String tmpDir) {
		if (!tmpDir.endsWith(File.separator))
			tmpDir+= File.separator;
		System.setProperty(Constants.PROPERTY_TMPDIR, tmpDir);
	}
	
	
	public void setTempPfx(String tmpPfx) {
		Constants.globalPfx= tmpPfx;
	}
	
	byte strand= FluxCapacitorConstants.STRAND_NONE;
	public void setStrandSpecific() {
		strand= FluxCapacitorConstants.STRAND_SPECIFIC;
	}
	
	public boolean setThreads(String nrThreads) {
		try {
			int x= Integer.parseInt(nrThreads);
			maxThreads= x;
			return true;
		} catch (NumberFormatException e) {
			return false;
		}
	}
	public boolean setUniformal() {
	
		// guessed
/*		String[] ss= s.split(",");
		try {
			readLenMin= Integer.parseInt(ss[0]);
		} catch (NumberFormatException e) {
			return false;
		}
		insertMinMax= new int[] {-1, -1};
		try {
			insertMinMax[0]= Integer.parseInt(ss[2]);
		} catch (ArrayIndexOutOfBoundsException e) {
			return true;
		} catch (NumberFormatException e) {
			return false;
		}
		try {
			insertMinMax[1]= Integer.parseInt(ss[3]);
		} catch (ArrayIndexOutOfBoundsException e) {
			return true;
		} catch (NumberFormatException e) {
			return false;
		}
*/		
		uniform= true;
		return true;
	}
	public void setPairedEnd() {
		pairedEnd= true;
	}
	
	public boolean setOutput(String s) {
		boolean 
			outputObs= this.outputObs,
			outputPred= this.outputPred,
			outputBalanced= this.outputBalanced,
			outputEvent= this.outputEvent,
			outputMapped= this.outputMapped,
			outputNotmapped= this.outputNotmapped,
			outputLP= this.outputLP,
			outputSorted= this.outputSorted,
			outputProfiles= this.outputProfiles,
			outputAll= this.outputAll,
			outputSplit= this.outputSplit,
			outputUnique= this.outputUnique,
			outputFreq= this.outputFreq,
			outputRfreq= this.outputRfreq,
			outputRcov= this.outputRcov,
			outputExon= this.outputExon,
			outputSJunction= this.outputSJunction,
			outputTranscript= this.outputTranscript,
			outputGene= this.outputGene,
			outputISize= this.outputISize;
			
		for (int i = 0; i < s.length(); i++) {
			if (s.charAt(i)== FluxCapacitorConstants.CLI_OUT_ALL) 
				outputAll= false;
			else if (s.charAt(i)== FluxCapacitorConstants.CLI_OUT_SPLIT)
				outputSplit= false;
			else if (s.charAt(i)== FluxCapacitorConstants.CLI_OUT_UNIQUE)
				outputUnique= false;
			else if (s.charAt(i)== FluxCapacitorConstants.CLI_OUT_OBS)
				outputObs= false;
			else if (s.charAt(i)== FluxCapacitorConstants.CLI_OUT_PRED)
				outputPred= false;
			else if (s.charAt(i)== FluxCapacitorConstants.CLI_OUT_BALANCED)
				outputBalanced= false;
			else if (s.charAt(i)== FluxCapacitorConstants.CLI_OUT_FREQ)
				outputFreq= false;
			else if (s.charAt(i)== FluxCapacitorConstants.CLI_OUT_RFREQ)
				outputRfreq= false;
			else if (s.charAt(i)== FluxCapacitorConstants.CLI_OUT_COVERAGE)
				outputRcov= false;
			else if (s.charAt(i)== FluxCapacitorConstants.CLI_OUT_EXON)
				outputExon= false;
			else if (s.charAt(i)== FluxCapacitorConstants.CLI_OUT_SJUNCTION)
				outputSJunction= false;
			else if (s.charAt(i)== FluxCapacitorConstants.CLI_OUT_TRANSCRIPT)
				outputTranscript= false;
			else if (s.charAt(i)== FluxCapacitorConstants.CLI_OUT_GENE)
				outputGene= false;
			else if (s.charAt(i)== FluxCapacitorConstants.CLI_OUT_EVENTS)
				outputEvent= false;
			else if (s.charAt(i)== FluxCapacitorConstants.CLI_OUT_LP)
				outputLP= false;
			else if (s.charAt(i)== FluxCapacitorConstants.CLI_OUT_PROFILES)
				outputProfiles= false;
			else if (s.charAt(i)== FluxCapacitorConstants.CLI_OUT_MAPPED)
				outputMapped= false;
			else if (s.charAt(i)== FluxCapacitorConstants.CLI_OUT_NOTMAPPED)
				outputNotmapped= false;
			else if (s.charAt(i)== FluxCapacitorConstants.CLI_OUT_ISIZE)
				outputISize= false;
			else if (s.charAt(i)== FluxCapacitorConstants.CLI_OUT_KEEPSORTED)
				outputSorted= false;
			else 
				return false;
		}
		
		this.outputAll= outputAll;
		this.outputSplit= outputSplit;
		this.outputUnique= outputUnique;
		this.outputObs= outputObs;
		this.outputPred= outputPred;
		this.outputBalanced= outputBalanced;
		this.outputFreq= outputFreq;
		this.outputRfreq= outputRfreq;
		this.outputRcov= outputRcov;
		this.outputExon= outputExon;
		this.outputSJunction= outputSJunction;
		this.outputTranscript= outputTranscript;
		this.outputGene= outputGene;
		this.outputEvent= outputEvent;
		this.outputMapped= outputMapped;
		this.outputNotmapped= outputNotmapped;
		this.outputLP= outputLP;
		this.outputProfiles= outputProfiles;
		this.outputISize= outputISize;
		this.outputSorted= outputSorted;

		return true;
	}
	
	
	static {
        try {
            Method m;
            m = FluxCapacitor.class.getDeclaredMethod("setFileReads", new Class[]{String.class});
            FluxCapacitorConstants.cliShortMap.put(FluxCapacitorConstants.CLI_SHORT_SRA, m);
            FluxCapacitorConstants.cliLongMap.put(FluxCapacitorConstants.CLI_LONG_SRA, m);
            FluxCapacitorConstants.cliExplMap.put(new String[]{FluxCapacitorConstants.CLI_SHORT_PFX + FluxCapacitorConstants.CLI_SHORT_SRA.toString(), FluxCapacitorConstants.CLI_LONG_PFX + FluxCapacitorConstants.CLI_LONG_SRA},
                    "set file containing Short Reads Archive (mandatory!)\n");

            m = FluxCapacitor.class.getDeclaredMethod("setFileReference", new Class[]{String.class});
            FluxCapacitorConstants.cliShortMap.put(FluxCapacitorConstants.CLI_SHORT_REF, m);
            FluxCapacitorConstants.cliLongMap.put(FluxCapacitorConstants.CLI_LONG_REF, m);
            FluxCapacitorConstants.cliExplMap.put(new String[]{FluxCapacitorConstants.CLI_SHORT_PFX + FluxCapacitorConstants.CLI_SHORT_REF.toString(), FluxCapacitorConstants.CLI_LONG_PFX + FluxCapacitorConstants.CLI_LONG_REF},
                    "set file with REFerence annotation (mandatory!)\n");

            m = FluxCapacitor.class.getDeclaredMethod("setNameOutDir", new Class[]{String.class});
            FluxCapacitorConstants.cliShortMap.put(FluxCapacitorConstants.CLI_SHORT_FILENAME, m);
            FluxCapacitorConstants.cliLongMap.put(FluxCapacitorConstants.CLI_LONG_FILENAME, m);
            FluxCapacitorConstants.cliExplMap.put(new String[]{FluxCapacitorConstants.CLI_SHORT_PFX + FluxCapacitorConstants.CLI_SHORT_FILENAME.toString(), FluxCapacitorConstants.CLI_LONG_PFX + FluxCapacitorConstants.CLI_LONG_FILENAME},
                    "set output fileName prefix (default stdout)\n");

            m = FluxCapacitor.class.getDeclaredMethod("setForce", null);
            FluxCapacitorConstants.cliShortMap.put(FluxCapacitorConstants.CLI_SHORT_FORCE, m);
            FluxCapacitorConstants.cliLongMap.put(FluxCapacitorConstants.CLI_LONG_FORCE, m);
            FluxCapacitorConstants.cliExplMap.put(new String[]{FluxCapacitorConstants.CLI_SHORT_PFX + FluxCapacitorConstants.CLI_SHORT_FILENAME.toString(), FluxCapacitorConstants.CLI_LONG_PFX + FluxCapacitorConstants.CLI_LONG_FILENAME},
                    "set force (no overwrite checks)\n");

            m = FluxCapacitor.class.getDeclaredMethod("setInstall", (Class[]) null);
            FluxCapacitorConstants.cliLongMap.put(FluxCapacitorConstants.CLI_LONG_INSTALL, m);
            FluxCapacitorConstants.cliExplMap.put(new String[]{FluxCapacitorConstants.CLI_LONG_PFX + FluxCapacitorConstants.CLI_LONG_INSTALL},
                    "installs the basic wrapper script (no reads are mapped)\n");

            m = FluxCapacitor.class.getDeclaredMethod("setJVM", new Class[]{String.class});
            FluxCapacitorConstants.cliLongMap.put(FluxCapacitorConstants.CLI_LONG_JVM, m);
            FluxCapacitorConstants.cliExplMap.put(new String[]{FluxCapacitorConstants.CLI_LONG_PFX + FluxCapacitorConstants.CLI_LONG_JVM},
                    "set a specific Java Virtual Machine home (installation)\n");

            m = FluxCapacitor.class.getDeclaredMethod("setLib", new Class[]{String.class});
            FluxCapacitorConstants.cliLongMap.put(FluxCapacitorConstants.CLI_LONG_LIB, m);
            FluxCapacitorConstants.cliExplMap.put(new String[]{FluxCapacitorConstants.CLI_LONG_PFX + FluxCapacitorConstants.CLI_LONG_LIB},
                    "set path to native libraries (installation)\n");

            m = FluxCapacitor.class.getDeclaredMethod("setBatch", (Class[]) null);
            FluxCapacitorConstants.cliShortMap.put(FluxCapacitorConstants.CLI_SHORT_BATCH, m);
            FluxCapacitorConstants.cliLongMap.put(FluxCapacitorConstants.CLI_LONG_BATCH, m);
            FluxCapacitorConstants.cliExplMap.put(new String[]{FluxCapacitorConstants.CLI_SHORT_PFX + FluxCapacitorConstants.CLI_SHORT_BATCH.toString(), FluxCapacitorConstants.CLI_LONG_PFX + FluxCapacitorConstants.CLI_LONG_BATCH},
                    "set Batch mode, suppresses file checks and stderr communication\n");

            m = FluxCapacitor.class.getDeclaredMethod("setUniformal", (Class[]) null);
            FluxCapacitorConstants.cliShortMap.put(FluxCapacitorConstants.CLI_SHORT_UNIF, m);
            FluxCapacitorConstants.cliLongMap.put(FluxCapacitorConstants.CLI_LONG_UNIF, m);
            FluxCapacitorConstants.cliExplMap.put(new String[]{FluxCapacitorConstants.CLI_SHORT_PFX + FluxCapacitorConstants.CLI_SHORT_UNIF.toString(), FluxCapacitorConstants.CLI_LONG_PFX + FluxCapacitorConstants.CLI_LONG_UNIF},
                    "set uniformal distribution no profiling step is carried out\n");

            m = FluxCapacitor.class.getDeclaredMethod("setPairedEnd", (Class[]) null);
            FluxCapacitorConstants.cliShortMap.put(FluxCapacitorConstants.CLI_SHORT_PAIR, m);
            FluxCapacitorConstants.cliLongMap.put(FluxCapacitorConstants.CLI_LONG_PAIR, m);
            FluxCapacitorConstants.cliExplMap.put(new String[]{FluxCapacitorConstants.CLI_SHORT_PFX + FluxCapacitorConstants.CLI_SHORT_PAIR.toString(), FluxCapacitorConstants.CLI_LONG_PFX + FluxCapacitorConstants.CLI_LONG_PAIR}, "set input paired ends, " +
                    "read name expected in FMRD format (see http://fluxcapacitor.wikidot.com/formats:fmrd)\n");

            m = FluxCapacitor.class.getDeclaredMethod("setProfile", new Class[]{String.class});
            FluxCapacitorConstants.cliLongMap.put(FluxCapacitorConstants.CLI_LONG_PROFILE, m);
            FluxCapacitorConstants.cliExplMap.put(new String[]{FluxCapacitorConstants.CLI_SHORT_PFX + FluxCapacitorConstants.CLI_SHORT_PAIR.toString(), FluxCapacitorConstants.CLI_LONG_PFX + FluxCapacitorConstants.CLI_LONG_PAIR}, "set profile name");

            m = FluxCapacitor.class.getDeclaredMethod("setOutput", new Class[]{String.class});
            FluxCapacitorConstants.cliLongMap.put(FluxCapacitorConstants.CLI_LONG_OUT, m);
            FluxCapacitorConstants.cliShortMap.put(FluxCapacitorConstants.CLI_SHORT_OUT, m);
            FluxCapacitorConstants.cliExplMap.put(new String[]{FluxCapacitorConstants.CLI_SHORT_PFX + FluxCapacitorConstants.CLI_SHORT_OUT.toString(), FluxCapacitorConstants.CLI_LONG_PFX + FluxCapacitorConstants.CLI_LONG_OUT},
                    "select output from [acdefijkgmnoprstuv]\n"
                            + "a All (scope)\n"
                            + "c Coverage (measure)\n"
                            + "d preDiction (base)\n"
                            + "e Exon (feature)\n"
                            + "f Frequency (measure)\n"
                            + "i Insert size (additional, paired-end only)\n"
                            + "j splice Junction (feature)\n"
                            + "k Keepsorted (additional)\n"
                            + "g Gene (feature)\n"
                            + "m Mapped (additional)\n"
                            + "n Notmapped (additional)\n"
                            + "o Observed (base)\n"
                            + "p Profiles (additional)\n"
                            + "r Relative frequency (measure)\n"
                            + "s Split (scope)\n"
                            + "t Transcript (feature)\n"
                            + "u Unique (scope)\n"
                            + "v eVents (feature)\n"
            );
            // g gene, l linear program, m mate-edges, x exon-junctions

            m = FluxCapacitor.class.getDeclaredMethod("setHelp", (Class[]) null);
            FluxCapacitorConstants.cliLongMap.put(FluxCapacitorConstants.CLI_LONG_HELP, m);
            FluxCapacitorConstants.cliShortMap.put(FluxCapacitorConstants.CLI_SHORT_HELP, m);
            FluxCapacitorConstants.cliExplMap.put(new String[]{FluxCapacitorConstants.CLI_SHORT_PFX + FluxCapacitorConstants.CLI_SHORT_HELP.toString(), FluxCapacitorConstants.CLI_LONG_PFX + FluxCapacitorConstants.CLI_LONG_HELP},
                    "print help summary");

            m = FluxCapacitor.class.getDeclaredMethod("setLogLevel", new Class[]{String.class});
            FluxCapacitorConstants.cliLongMap.put(FluxCapacitorConstants.CLI_LONG_VERBOSE, m);
            FluxCapacitorConstants.cliShortMap.put(FluxCapacitorConstants.CLI_SHORT_VERBOSE, m);
            FluxCapacitorConstants.cliExplMap.put(new String[]{FluxCapacitorConstants.CLI_SHORT_PFX + FluxCapacitorConstants.CLI_SHORT_VERBOSE.toString(), FluxCapacitorConstants.CLI_LONG_PFX + FluxCapacitorConstants.CLI_LONG_VERBOSE},
                    "set verbose level (SILENT, VERBOSE, ERRORS, DEBUG)");

            m = FluxCapacitor.class.getDeclaredMethod("setThreads", new Class[]{String.class});
            FluxCapacitorConstants.cliLongMap.put(FluxCapacitorConstants.CLI_LONG_THREAD, m);
            FluxCapacitorConstants.cliShortMap.put(FluxCapacitorConstants.CLI_SHORT_THREAD, m);
            FluxCapacitorConstants.cliExplMap.put(new String[]{FluxCapacitorConstants.CLI_SHORT_PFX + FluxCapacitorConstants.CLI_SHORT_THREAD.toString(), FluxCapacitorConstants.CLI_LONG_PFX + FluxCapacitorConstants.CLI_LONG_THREAD}, "set multi-thread mode, provide number of threads\n" +
                    "(time gain only with complex linear programs, otherwise default=1 recommended)\n");

            m = FluxCapacitor.class.getDeclaredMethod("setLocal", (Class[]) null);
            FluxCapacitorConstants.cliLongMap.put(FluxCapacitorConstants.CLI_LONG_LOCAL, m);
            FluxCapacitorConstants.cliShortMap.put(FluxCapacitorConstants.CLI_SHORT_LOCAL, m);
            FluxCapacitorConstants.cliExplMap.put(new String[]{FluxCapacitorConstants.CLI_SHORT_PFX + FluxCapacitorConstants.CLI_SHORT_LOCAL.toString(), FluxCapacitorConstants.CLI_LONG_PFX + FluxCapacitorConstants.CLI_LONG_LOCAL},
                    "work locally, i.e., copy all files to the temporary directory\n");

            m = FluxCapacitor.class.getDeclaredMethod("setTempDir", new Class[]{String.class});
            FluxCapacitorConstants.cliLongMap.put(FluxCapacitorConstants.CLI_LONG_TMP, m);
            FluxCapacitorConstants.cliExplMap.put(new String[]{FluxCapacitorConstants.CLI_LONG_PFX + FluxCapacitorConstants.CLI_LONG_TMP},
                    "set path to the temporary directory\n");

            m = FluxCapacitor.class.getDeclaredMethod("setTempPfx", new Class[]{String.class});
            FluxCapacitorConstants.cliLongMap.put(FluxCapacitorConstants.CLI_LONG_TPX, m);
            FluxCapacitorConstants.cliExplMap.put(new String[]{FluxCapacitorConstants.CLI_LONG_PFX + FluxCapacitorConstants.CLI_LONG_TPX},
                    "set prefix for temporary files\n");

            m = FluxCapacitor.class.getDeclaredMethod("setCompression", new Class[]{String.class});
            FluxCapacitorConstants.cliShortMap.put(FluxCapacitorConstants.CLI_SHORT_COMPRESSION, m);
            FluxCapacitorConstants.cliLongMap.put(FluxCapacitorConstants.CLI_LONG_COMPRESSION, m);
            FluxCapacitorConstants.cliExplMap.put(new String[]{FluxCapacitorConstants.CLI_SHORT_PFX + FluxCapacitorConstants.CLI_SHORT_COMPRESSION.toString(),
                    FluxCapacitorConstants.CLI_LONG_PFX + FluxCapacitorConstants.CLI_LONG_COMPRESSION},
                    "set compression method for output files (output file, mapped read-mappings, not-mapped read-mappings, insert sizes)\n");

            m = FluxCapacitor.class.getDeclaredMethod("setStrandSpecific", (Class[]) null);
            FluxCapacitorConstants.cliLongMap.put(FluxCapacitorConstants.CLI_LONG_SSPECIFIC, m);
            FluxCapacitorConstants.cliExplMap.put(new String[]{FluxCapacitorConstants.CLI_LONG_PFX + FluxCapacitorConstants.CLI_LONG_SSPECIFIC},
                    "set strand specific reads (default: strand information disregarded/disabled)\n");

        } catch (NoSuchMethodException e) {
            ; // :)
        }
        FluxCapacitorConstants.SHELL_NONE = 0;
    }

	
	private boolean move(File src, File dest) {
		if (FileHelper.move(src, dest)) {
			if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) {
				System.err.println("\t"+ src.getAbsolutePath()+ "\n\t->"+ dest.getAbsolutePath());
			}
			return true;
		} else { 
			if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) {
				System.err.println("\tfailed, output in:\n\t"+ src.getAbsolutePath());
			}
			return false;
		}
	}
	
	private boolean copy(File src, File dest) {
		if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) 
			System.err.println("\t"+ src.getAbsolutePath()+ "\n\t->"+ dest.getAbsolutePath());
		try {
            Log.progressStart("copying");
			FileHelper.fastChannelCopy(src, dest, false);

			Log.progressFinish(StringUtils.OK, true);
			return true;
		} catch (IOException e) {
			if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) 
				System.err.println(e.getMessage());;
			return false;
		}
	}
	void fileFinish() {

		if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
			System.err.println("\n[FINISHING] closing file handles and cleaning up");
		
		// remove temp files: ref and reads
		boolean b= getBedReader().close();
		//b= getFileBED().delete();	// TODO deactivated cleanup
		b= getGTFreader().close();
		//b= getFileGTF().delete();
		
//		if (fileOut!= null)
//			appendFreq();

		if (outputMapped) {
			try {
				getWriterMappedReads().flush();
				getWriterMappedReads().close();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}	
		
		if (outputNotmapped) {
			try {
				getWriterNotmappedReads().flush();
				getWriterNotmappedReads().close();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}

		if (outputLP) {
			FileHelper.setSilent(false); 
			String sfx= FileHelper.getCompressionExtension(FileHelper.COMPRESSION_GZIP);
			File dst= this.pars.fileLPzip; // new File(fileOutDir+ File.separator+ getNameLP()+ (sfx== null? "": Constants.DOT+ sfx));
			if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) 
				System.err.println("\tzipping "+ fileLPdir.getAbsolutePath()
						+"\n\t->"+ dst.getAbsolutePath());
			if (!FileHelper.zipFilesInDir(fileLPdir, dst)) {
				if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
					System.err.println("\n[PROBLEMS] encountered error during zipping, check file.");
			} else {
				FileHelper.rmDir(fileLPdir); // otherwise Mr. Proper
			}
		}
		
// 		if (compressionOut== FileHelper.COMPRESSION_NONE&& !copyLocal)
//			return;

		// move 
		moveOrDeflate(fileOut, fileOUToriginal// new File(fileOutDir+ File.separator+ getNameOut())
				, compressionOut);	

		if (outputMapped) {
			String sfx= FileHelper.getCompressionString(compressionBED);
			moveOrDeflate(getFileMappedReads(), 
					new File(fileOutDir+ File.separator+ getNameMappedReads()+ (sfx== null?"":Constants.DOT+sfx)), 
					compressionBED);
		}
		
		if (outputNotmapped) {
			String sfx= FileHelper.getCompressionString(compressionBED);
			moveOrDeflate(getFileNotMappedReads(), 
					new File(fileOutDir+ File.separator+ getNameNotMappedReads()+ (sfx== null?"":Constants.DOT+sfx)), 
					compressionBED);
		}
		if (outputISize) {
			String sfx= FileHelper.getCompressionString(compressionOut);
			moveOrDeflate(getFileISize(), 
					new File(fileOutDir+ File.separator+ getNameISize()+ (sfx== null?"":Constants.DOT+sfx)), 
					compressionOut);
		}
		// profile never moved
/*		if (outputProfiles) 
			moveOrDeflate(getFileProfile(), 
					new File(fileOutDir+ File.separator+ getNameProfile()), 
					FileHelper.COMPRESSION_NONE);	// already compressed
*/					
		
	}
	
	private boolean moveOrDeflate(File src, File dst,
			byte compression) {
		
		if (dst== null) {
			if (compression!= FileHelper.COMPRESSION_NONE)
				return moveDeflate(src, src, compression, false);
		} else {
			if (compression== FileHelper.COMPRESSION_NONE)
				return move(src, dst);
			else
				return moveDeflate(src, dst, compression, false);
		}
		return true;
	}

	private boolean copyOrDeflate(File src, File dest,
			byte compression) {
		
		if (dest== null) {
			if (compression!= FileHelper.COMPRESSION_NONE)
				return moveDeflate(src, src, compression, true);
		} else {
			if (compression== FileHelper.COMPRESSION_NONE)
				return copy(src, dest);
			else
				return moveDeflate(src, dest, compression, true);
		}
		return true;
	}

	private boolean moveDeflate(File src, File dest,
			byte compression, boolean copy) {
		try {
			if (src.getAbsolutePath().equals(dest.getAbsolutePath())) 
				dest= new File(dest.getAbsolutePath()+ Constants.DOT+ FileHelper.getCompressionExtension(compression));
			if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) 
				System.err.println("\t"+ src.getAbsolutePath()+ "\n\t->"+ dest.getAbsolutePath());			
			FileHelper.deflate(src, dest, compression);
			if (!copy) {
				if (!src.delete())
					return false;
			}
			return true;
		} catch (Exception e) {
			if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) {
				System.err.println("\n[AIII] Problems during deflate: "+ e.getMessage());
				e.printStackTrace();
			}
			return false;
		}
	}

	private boolean copyInflate(File src, File dest,
			byte compression) {
		try {
			if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) {
				System.err.println("\tinflate "+ src.getAbsolutePath()+ "\n\t->"+ dest.getAbsolutePath());
			}
			FileHelper.inflate(src, dest, compression);
			return true;
		} catch (Exception e) {
			if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
				System.err.println("\n[AIII] Problems during inflate");
			return false;
		}
	}
	private static String createID() {
		SimpleDateFormat format= new SimpleDateFormat("yyMMddHHmmssSSSS");
		return format.format(new Date());
	}
	
	public String getRunID() {
		if (runID == null) {
			runID = createID();
		}

		return runID;
	}
	
	private File createTempFile(String id, String ext) {
		String s= System.getProperty(Constants.PROPERTY_TMPDIR)
				+ File.separator
				+ (Constants.globalPfx== null?"":Constants.globalPfx+ "_")
				+ FluxCapacitorConstants.PFX_CAPACITOR
				+ "."
				+ getRunID()
				+ "."+ id
				+ (ext== null?"": "."+ ext);
		File f= new File(s);
		f.deleteOnExit();
		return f;
	}
	
	static String createFileName(String base, byte compression) {
		if (compression== FileHelper.COMPRESSION_NONE) 
			return base;
		if (compression== FileHelper.COMPRESSION_ZIP)
			base+= '.'+ MyFile.SFX_ZIP;
		else if (compression== FileHelper.COMPRESSION_GZIP)
			base+= '.'+ MyFile.SFX_GZ;
		return base;
	}
	
	public boolean fileInitReference() {
		
		isReadyGTF= false;
		gtfReader= null;
		
		// copy or deflate, if necessary: fGTForig -> fGTF
		compressionGTF= FileHelper.getCompression(fileGTForiginal);
		if (compressionGTF== FileHelper.COMPRESSION_NONE&& !copyLocal)
			fileGTF= fileGTForiginal;
		else {
			if (compressionGTF== FileHelper.COMPRESSION_NONE) 
				fileGTF= new File(System.getProperty(Constants.PROPERTY_TMPDIR)+ File.separator+ Constants.getGlobalPfx()+ "_"+ MyFile.getFileNameOnly(fileGTForiginal.getAbsolutePath()));
			else 
				fileGTF= new File(System.getProperty(Constants.PROPERTY_TMPDIR)+ File.separator+ Constants.getGlobalPfx()+ "_"+ MyFile.getFileNameOnly(fileGTForiginal.getAbsolutePath()));
			if (force|| !fileGTF.exists()) {
				if (fileGTF.exists()&& Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
					System.err.println("\texisting file overwritten: "+ fileGTF.getAbsolutePath());
				if(!copyOrDeflate(fileGTForiginal, fileGTF, compressionGTF))
					return false;
			} else {
				if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
					System.err.println("\tfile seems exists and no overwrite forced: "+ fileGTF.getAbsolutePath());
			}
		}
		
		
		// check sorting
		if (cheatDisableFCheck) {
			System.err.println("\tTEST Sort check disabled !!!");
			isSortedGTF= true;
		} else { 
			
			if (getGTFreader().isApplicable()) 
				isSortedGTF= true;
			else {
				isSortedGTF= false;
				File tmpGTF= getGTFreader().createSortedFile();
				if (tmpGTF== null)
					return false;
				//boolean bb= getFileGTF().delete();	// TODO do sth when false..
				this.fileGTF= tmpGTF;
				if (outputSorted&& fileOutDir!= null) {
					String ext= FileHelper.getCompressionExtension(compressionGTF);
					String destFName= fileOutDir+ File.separator+ Constants.getGlobalPfx()+ "_"+ 
										fileGTF.getName()+
										(ext== null?"":Constants.DOT+ ext);
					File destFile= new File(destFName);
					if ((!force)&& !ensureFileCanWrite(destFile))
						return false;
					else if (!copyOrDeflate(tmpGTF, destFile,compressionGTF))
						return false;
				}
			}
		}
		
		// scan file
		if (cheatDisableFCheck) {
			System.err.println("\tTEST GTF scan disabled !!!");
			getGTFreader();
		} else {
			gtfReader= null;	// reinit
			getGTFreader().scanFile();
			if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) {
				System.err.println(Constants.TAB+ getGTFreader().getNrGenes()+ " loci, "
						+ getGTFreader().getNrTranscripts()+ " transcripts, "
						+ getGTFreader().getNrExons()+ " exons.");
			}
		}
		isReadyGTF= true;

		return true;
	}
	
	
	int nrBEDreads= -1, nrBEDmappings= -1;
	private int checkBEDscanMappings= 0;
	public boolean fileInitBED() {

		isReadyBED= false;
		
		// check compression
		byte compressionBED= FileHelper.getCompression(fileBEDoriginal);
		if (compressionBED== FileHelper.COMPRESSION_NONE&& !copyLocal)
			fileBED= fileBEDoriginal;
		else {	// copy/decompress local
			if (compressionBED== FileHelper.COMPRESSION_NONE) 
				fileBED= new File(System.getProperty(Constants.PROPERTY_TMPDIR)+ File.separator+ Constants.getGlobalPfx()+ "_"+ MyFile.getFileNameOnly(fileBEDoriginal.getAbsolutePath()));
			else 
				fileBED= new File(System.getProperty(Constants.PROPERTY_TMPDIR)+ File.separator+ Constants.getGlobalPfx()+ "_"+ MyFile.getFileNameOnly(fileBEDoriginal.getAbsolutePath()));
			if (force|| !fileBED.exists()) {
				if (fileBED.exists()&& Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
					System.err.println("\texisting file overwritten: "+ fileBED.getAbsolutePath());
				if(!copyOrDeflate(fileBEDoriginal, fileBED, compressionBED))
					return false;
			} else {
				if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
					System.err.println("\tfile seems exists and no overwrite forced: "+ fileBED.getAbsolutePath());
			}
		}

		
		// check sorting
		if (cheatDisableFCheck) {
			System.err.println("\tTEST sort check disabled !!!");
			isSortedBED= true;
		} else {  
			
			if (getBedReader().isApplicable()) 
				isSortedBED= true;
			else {
				isSortedBED= false;
				File tmp= getBedReader().sortBED(fileBED);
				if (tmp== null)
					return false;
				//getFileBED().delete();
				this.fileBED= tmp;
				if (outputSorted&& fileOutDir!= null) {
					String destFName= fileOutDir+ File.separator+ Constants.getGlobalPfx()+ "_"+ 
										MyFile.append(fileBED.getName(), "_sorted", false, 
												FileHelper.getCompressionExtension(this.compressionBED));
					File destFile= new File(destFName);
					compressionBED= FileHelper.COMPRESSION_NONE;
					if ((!force)&& !ensureFileCanWrite(destFile))
						return false;
					else if (!copyOrDeflate(tmp, destFile, compressionBED))
						return false;
				}
			}
		}
		
		
		// scan file
		bedWrapper= null;
/*		descriptor= getBedReader().checkReadDescriptor(pairedEnd, stranded);
		if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) {
			System.err.print("\tread descriptor ");
			if (descriptor== null)
				System.err.println("none");
			else if (descriptor instanceof FMRD)
				System.err.println("flux");
			else if (descriptor instanceof SolexaDescriptor)
				System.err.println("solexa");
		}
*/
		if(!getBedReader().checkReadDescriptor(descriptor2))
			return false;
		
		if (cheatDisableFCheck) {
			System.err.println("\tTEST file scan disabled, reads= mappings= 1M !!!");
			checkBEDscanMappings= 100000000;
			nrBEDreads= 100000000;
			nrBEDmappings= 100000000;
		} else {
			getBedReader().scanFile();	
			checkBEDscanMappings= getBedReader().getCountAll();
			nrBEDreads= getBedReader().getCountReads();
			nrBEDmappings= getBedReader().getCountAll();
			if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) {
				System.err.println("\t"+ nrBEDreads+ " reads, "+ nrBEDmappings
						+ " mappings: R-factor "+(getBedReader().getCountAll()/ (float) getBedReader().getCountReads()));
				System.err.println("\t"+ getBedReader().getCountEntire()+ " entire, "+ getBedReader().getCountSplit()
						+ " split mappings ("+ (getBedReader().getCountSplit()* 10f/ getBedReader().getCountAll())+ "%)");
			}
		}
		
		isReadyBED= true;
		
		return true;
	}
	
	public boolean isInputReady() {
		return (isReadyBED&& isReadyGTF); 
	}
	
	private static boolean waitForYesNo() {
		while(true) {
			String s= readSystemIn().trim().toLowerCase();
			if (s.equals("y")|| s.equals("yes")|| s.equals("yeah")) {
				return true;
			}
			if (s.equals("n")|| s.equals("no")|| s.equals("nope")) {
				return false;
			}
			System.err.println("\nDidn't get that, repeat:");
		}
	}
	
	private static boolean ensureFileCanWrite(File file) {
		if (cheatDisableFCheck)
			return true;
		
		boolean returnVal= true;;
		if (file.exists()) {
			if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) {
				System.err.println("[PLOING] this file already exists:\n"+ file.getAbsolutePath());
				System.err.println("\tConfirm overwrite (Yes/No/Don't know)");
				returnVal= waitForYesNo();
			} else
				returnVal= false;
		}
		
		if (!returnVal) {
			if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) {
				System.err.println("[CANNOT] cannot write to:\n\t"+ file.getAbsolutePath());
				System.err.println("\nCiao.");
			}
			System.exit(-1);
		}
		
		return true;
	}

	public int[] getInsertMinMax() {
		
		if (!pairedEnd)
			return null;
		
		//System.err.println("[ISIZE] calc insert size..");
		int min= isizeV.ax[0];
		int max= isizeV.ax[isizeV.size- 1];
		//System.err.println("\t"+ isizeV.size+ " inserts, min="+ min+", max="+ max); 
		int q25= isizeV.getQuartile(0.25f);
		int q75= isizeV.getQuartile(0.75f);
		int iqr= q75- q25;
		
		//get factor from dist q25<>min
		// boxplot: 1.5, between 1.5 and 3.0 mild outlier
		double factor= (q25- min)/ (double) iqr;
		
		int minIn= Math.max(min, (int) Math.ceil(q25- (factor* iqr)));
		int maxIn= Math.min(max, (int) Math.floor(q75+ (factor* iqr)));
		//System.err.println("\tlo="+minIn+", x25="+q25+", iqr="+iqr+", x75="+q75+", hi="+maxIn);
		int[] mm= new int[] {minIn, maxIn};
//		if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
//			System.err.println("\tdetermined insert min "+ mm[0]+", max "+ mm[1]);
		mm[0]-= readLenMin; mm[1]-= readLenMin;
		return mm;
	}
	
	
	int profileNr= 3; 
	BufferedWriter testWriter;
	
	public void run() {
		
		if (!fileInit())
			System.exit(-1);
	
		long t0= System.currentTimeMillis();
		
		//createBins(profileNr);		
		
		if (pairedEnd)
			isizeV= new BinVector();
		
		
		profileStub= null;
/*		new int[BIN_LEN.length+ 1][];
		for (int i = 0; i < profileStub.length; i++) {
			//profileStub[i]= new int[i< BIN_LEN.length?BIN_LEN[i]:10000];
			profileStub[i]= new int[20];
		}
*/		
/*		if (strand== STRAND_ENABLED) {
			profileStubRev= new int[BIN_LEN.length+ 1][];
			for (int i = 0; i < profileStub.length; i++) {
				//profileStubRev[i]= new int[i< BIN_LEN.length?BIN_LEN[i]:10000];
				profileStubRev[i]= new int[20];
			}
		}
*/		
		
		// profiling
		profile= getProfile();
		if (profile== null) {
			exit(-1);
		}
			
		
//		System.exit(0);
		
//		if (map)
//			func.finish();
		
//		try {
//			testWriter= new BufferedWriter(new FileWriter("P:\\rgasp1.3\\HepG2_new\\test.bed"));
//		} catch (IOException e) {
//			e.printStackTrace();
//		}
		
		
		explore(FluxCapacitorConstants.MODE_RECONSTRUCT);
		
//		try {
//			testWriter.flush();
//			testWriter.close();
//		} catch (IOException e) {
//			e.printStackTrace();
//		}
		
		
		fileFinish();
		
		if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) {
			System.err.println("\n[TICTAC] I finished flux in "
					+((System.currentTimeMillis()- t0)/ 1000)+" sec.\nCheers!");
		}
		
		//System.err.println("over "+ GraphLPsolver.nrOverPredicted+", under "+GraphLPsolver.nrUnderPredicted);
	}

	private Profile readProfiles(File fileProfileOriginal) {
		
		try {
			profile= new Profile(this);
			
			ZipFile zf= new ZipFile(fileProfile);
			Enumeration entries= zf.entries();
			String line;
			Vector<Integer> v= new Vector<Integer>();
			Vector<UniversalMatrix> w= new Vector<UniversalMatrix>();
			System.err.println("[LOAD] getting profiles");
			while (entries.hasMoreElements()) {
				ZipEntry ze= (ZipEntry) entries.nextElement();
				BufferedReader buffy = new BufferedReader(
		                new InputStreamReader(zf.getInputStream(ze)));
				int lcount= 0;
				while ((line= buffy.readLine())!= null) 
					++lcount;
				buffy.close();
				v.add(lcount);
				UniversalMatrix m= new UniversalMatrix(lcount);
				buffy = new BufferedReader(
		                new InputStreamReader(zf.getInputStream(ze)));
				lcount= 0;
				while ((line= buffy.readLine())!= null) {
					String[] ss= line.split("\t");
					assert(ss.length== 2);
					m.sense[lcount]= Integer.parseInt(ss[0]);
					m.sums+= m.sense[lcount];
					m.asense[lcount]= Integer.parseInt(ss[1]);
					m.suma+= m.asense[lcount];
					++lcount;
				}
				buffy.close();
				assert(lcount== m.sense.length);
				w.add(m);
			}
			zf.close();
			
			int[] len= new int[v.size()];
			for (int i = 0; i < len.length; i++) 
				len[i]= v.elementAt(i);
			Arrays.sort(len);
			profile.masters= new UniversalMatrix[w.size()];
			for (int i = 0; i < len.length; i++) {
				for (int j = 0; j < len.length; j++) {
					if (len[i]== v.elementAt(j)) {
						profile.masters[i]= w.elementAt(j);
						// check
						for (int n = 0; n < profile.masters[i].getLength(); n++) {
							if(profile.masters[i].asense[n]== 0|| profile.masters[i].sense[n]== 0) {
								if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
									System.err.println("\tprofile with 0-count positions");
								return null;
							}
						}
					}
				}
			}
			System.err.println("\tfound "+ profile.masters.length+" profiles.");
			
			return profile;
			
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		return null;
	}

	static boolean miss= false;
	void appendFreq() {
		
		if (fileOut== null)
			return;
		
		long t0= System.currentTimeMillis();
		
		double[] appendLengths= new double[3];
		double[][] appendReads= new double[3][];
		for (int i = 0; i < appendReads.length; i++) 
			appendReads[i]= new double[3];
		double[][] appendCovs= new double[3][];
		for (int i = 0; i < appendReads.length; i++) 
			appendReads[i]= new double[3];
		try {
			//File tmpFile= File.createTempFile("fcapacitor", "gtf");
			if (!fileOut.exists())
				return;
			
			BufferedReader buffy= new BufferedReader(new FileReader(fileOut));
			File fileTmp= createTempFile(fileOut.getName()+"__append", MyFile.getExtension(fileOut.getName())); 
//				(fileOUToriginal== null)? 
//					File.createTempFile(PFX_CAPACITOR, "gtf"): fileOUToriginal;
			BufferedWriter writer= new BufferedWriter(new FileWriter(fileTmp));
			

            Log.message("\tappending relative measurements");
			Log.progressStart("progress");
			long totBytes= fileOut.length(), bytes= 0;
			int perc= 0;
			String s= null;
			while((s= buffy.readLine())!= null) {
				bytes+= s.length()+ 1;
                Log.progress(bytes, totBytes);
				writer.write(s);
				String[] ss= s.split("\\s");	// TODO kill regexp
				if (ss.length< 8)
					System.err.println("[incomplete line] "+ s);

				boolean event= false;
				if (ss[2].contains("event")) 
					event= true;
				
				// TODO check: ss[2].equals(GFF_FEATURE_JUNCTION)
				if (ss[2].equals(FluxCapacitorConstants.GFF_FEATURE_FRAGMENT)|| ss[2].equals(FluxCapacitorConstants.GFF_FEATURE_PAIRED)) {
					writer.write("\n");
					continue;
				}
				int dim= 0;
				boolean repeat= true;
				StringBuilder[] builder= new StringBuilder[18];	// max nb of additional attributes: addAttributes(), (i* 6)+ (j* 2)+ x;
				for (int i = 0; i < builder.length; i++) 
					builder[i]= new StringBuilder(); 
				while (repeat) {
					for (int i = 0; i < appendReads.length; i++) 
						for (int j = 0; j < appendReads[i].length; j++) 
							appendReads[i][j]= -1;
					for (int i = 0; i < appendLengths.length; i++) 
						appendLengths[i]= -1;
					for (int i = 8; i < ss.length; i+= 2) {
						String[] sss= null;
						String target= null;
						if (ss[i].length()== 0) {
							--i;
							continue;
						}
						
						try {
							target= ss[i+1].substring(0, ss[i+1].length()-1);	// 1, ss[i+1].length()-2
						} catch (Exception e) {
							System.currentTimeMillis();
						}
						
						if (event) {
							sss= target.split(",");
							if (dim== (sss.length- 1)) 
								repeat= false;
						} else
							repeat= false;

							// collect
						try {
							if (ss[i].equals(FluxCapacitorConstants.obsReadsAllTag)) 
								appendReads[0][0]= Double.parseDouble(sss== null?target: sss[dim]);
							else if (ss[i].equals(FluxCapacitorConstants.obsReadsSplitTag)) 
								appendReads[0][1]= Double.parseDouble(sss== null?target: sss[dim]);
							else if (ss[i].equals(FluxCapacitorConstants.obsReadsUniqTag)) 
								appendReads[0][2]= Double.parseDouble(sss== null?target: sss[dim]);
							else if (ss[i].equals(FluxCapacitorConstants.predReadsAllTag)) 
								appendReads[1][0]= Double.parseDouble(sss== null?target: sss[dim]);
							else if (ss[i].equals(FluxCapacitorConstants.predReadsSplitTag)) 
								appendReads[1][1]= Double.parseDouble(sss== null?target: sss[dim]);
							else if (ss[i].equals(FluxCapacitorConstants.predReadsUniqTag)) 
								appendReads[1][2]= Double.parseDouble(sss== null?target: sss[dim]);
							else if (ss[i].equals(FluxCapacitorConstants.normReadsAllTag)) 
								appendReads[1][0]= Double.parseDouble(sss== null?target: sss[dim]);
							else if (ss[i].equals(FluxCapacitorConstants.normReadsSplitTag)) 
								appendReads[1][1]= Double.parseDouble(sss== null?target: sss[dim]);
							else if (ss[i].equals(FluxCapacitorConstants.normReadsUniqTag)) 
								appendReads[1][2]= Double.parseDouble(sss== null?target: sss[dim]);
							else if (ss[i].equals(FluxCapacitorConstants.lengthAllTag)) 
								appendLengths[0]= Double.parseDouble(sss== null?target: sss[dim]);
							else if (ss[i].equals(FluxCapacitorConstants.lengthSplitTag)) 
								appendLengths[1]= Double.parseDouble(sss== null?target: sss[dim]);
							else if (ss[i].equals(FluxCapacitorConstants.lengthUniqTag)) 
								appendLengths[2]= Double.parseDouble(sss== null?target: sss[dim]);
						} catch (NumberFormatException e) {
							if (Constants.verboseLevel> Constants.VERBOSE_NORMAL)
								System.err.println("Error parsing double "+ (sss== null?target: sss[dim]));; // NA
						}
					}
					
					addAttributes(appendLengths, appendReads, builder);
					++dim;
				}	
				
				StringBuilder sb= new StringBuilder();
				for (int i = 0; i < builder.length; i++) {
					if (builder[i]!= null) {
						sb.append(builder[i]);
						sb.append(";");	// \"
					}
				}
				sb.append("\n");
				String t= sb.toString();
				writer.write(t);
			}
			
			buffy.close();
			writer.flush();
			writer.close();
			fileOut.delete();			

            Log.progressFinish();

			if (fileOUToriginal== null) {
				if (!FileHelper.move(fileTmp, fileOut)) {
                    Log.error("[FAILED] Cannot move file!");
					System.exit(-1);
				}
			} else {
				fileOut= fileTmp;
			}
				
				

            Log.progressFinish();

//			if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) { 
//				Constants.progress.setValue(0);
//				Constants.progress.setString("progress");
//			}
//			boolean ok= FileHelper.move(tmpFile, fileOut, Constants.progress);
//			if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) { 
//				Constants.progress.finish();
//				System.err.println("\ttook "+((System.currentTimeMillis()- t0)/ 1000)+" sec.");
//			}
			
		} catch (Exception e) {
			e.printStackTrace();
			return;
		}
		
	}
	
	private void addAttributes(double[] lengths, double[][] reads, StringBuilder[] sb) {
		
		for (int i = 0; i < 3; i++) { // obs, pred, balanced
			String tag= FluxCapacitorConstants.GTF_ATTRIBUTES_BASE[i];
			tag+= FluxCapacitorConstants.GTF_ATTRIBUTE_TOKEN_SEP;
			for (int j = 0; j < 3; j++) {	// all, split, uniq
				String tag2= tag+ FluxCapacitorConstants.GTF_ATTRIBUTES_RESOLUTION[j];
				tag2+= FluxCapacitorConstants.GTF_ATTRIBUTE_TOKEN_SEP;
				for (int x = 0; x < 2; x++) {	// freq, cov
					int pos= (i* 6)+ (j* 2)+ x;
					if (lengths[j]< 0|| reads[i][j]< 0) {
						sb[pos]= null;
						continue;
					}
					String tag3= tag2+ FluxCapacitorConstants.GTF_ATTRIBUTES_MEASUREMENT[x+1];
					if (sb[pos].length()== 0) {
						sb[pos].append(" ");
						sb[pos].append(tag3);
						sb[pos].append(" "); // \"
					} else
						sb[pos].append(",");

					double rfreq= 0;
					//assert((nrReadsMapped== 0)== (reads[i][j]== 0));
					double base= nrReadsMapped;
					if (pairedEnd)
						base*= 2;
					rfreq= reads[i][j]/ base;	// rfreq
					if (x== 1) {	// cov
						if (lengths[j]== 0) {
							// happens for *pred* in areas that are too small for reads
							// now no more..
							try{assert(rfreq== 0);}catch(AssertionError err){
								if (Constants.verboseLevel> Constants.VERBOSE_NORMAL)
									System.err.println("Found 0-lengt for value "+rfreq);
							} 
							rfreq= 0;
						} else
							rfreq/= lengths[j];
						//if (!pairedEnd)
						rfreq*= 1000000000;	// rpkm
					}
					sb[pos].append(Float.toString((float) rfreq));
				}
			}
		}
		
	}
	
	public float calcRPKM(float reads, int len) {
		float rpkm= (float) ((reads/ (double) len)* (1000000000l/ (double) nrBEDreads));
		return rpkm;
	}
	
	public String getCompositeFName() {
		if (fileGTF== null|| fileBED== null)
			return null;
		File f= fileGTF, g= fileBED;
		if (fileGTForiginal!= null)
			f= fileGTForiginal;
		if (fileBED!= null)
			g= fileBED;		
		return MyFile.stripExtension(f.getName())+ "__"
					+ MyFile.stripExtension(g.getName());
	}
	
	private String getNameProfile() {
		String fName= getCompositeFName();
		if (fName== null)
			return null;
		String sfx= FileHelper.getCompressionExtension(compressionProfiles);
		return fName+ FluxCapacitorConstants.SFX_PROFILES+ (sfx== null? "": Constants.DOT+ sfx);
	}
	public File getFileProfile() {
		if (fileProfile == null) {
			String s= getNameProfile();
			if (s== null)
				return null;
			fileProfile=  
				new File(fileOutDir+ File.separator+ getCompositeFName()+ "_profiles."
						+ FileHelper.getCompressionExtension(FileHelper.COMPRESSION_ZIP));
				//createTempFile(s, FileHelper.getCompressionExtension(FileHelper.COMPRESSION_ZIP));
				
//				new File(System.getProperty(Constants.PROPERTY_TMPDIR)+ File.separator 
//				+ fName+ SFX_PROFILES
//				+ "."+ FileHelper.getCompressionExtension(FileHelper.COMPRESSION_ZIP));
		}

		return fileProfile;
	}
	
	private String getNameOut() {
		if (getCompositeFName()== null)
			return null;
		return getCompositeFName()+ Constants.DOT+ FluxCapacitorConstants.SFX_GTF;
	}
	
	public File getFileOut() {
		if (fileOut == null) {
			String fName= getCompositeFName();
			if (fName== null)
				return null;
			fileOut= createTempFile(getNameOut(), null);
		}

		return fileOut;
	}
	
	private String getNameISize() {
		String s= getCompositeFName();
		if (s== null)
			return null;
		return s+ FluxCapacitorConstants.SFX_INSERTSIZE+ Constants.DOT+ "txt";
	}
	
	public File getFileISize() {
		
		if (fileISize == null) {
			fileISize = createTempFile(getNameISize(), null);
		}

		return fileISize;
		
	}
	
	private String getNameLP() {
		return getCompositeFName()+ FluxCapacitorConstants.SFX_LP;
	}
	
	public File getFileLP() {
		if (fileLPdir == null) {
			fileLPdir= createTempFile(getNameLP(), null);
			if (fileLPdir.exists())
				fileLPdir.delete();
			boolean b= fileLPdir.mkdir();
		}

		return fileLPdir;
	}
	
	private String getNameMappedReads() {
		String s= getCompositeFName();
		if (s==  null)
			return null;
		return s+ FluxCapacitorConstants.SFX_MAPPED+ Constants.DOT+ FluxCapacitorConstants.SFX_BED;
	}
	
	
	private String getNameNotMappedReads() {
		String s= getCompositeFName();
		if (s== null)
			return null;
		return s+ FluxCapacitorConstants.SFX_NOTMAPPED+ Constants.DOT+ FluxCapacitorConstants.SFX_BED;
	}
	public File getFileMappedReads() {
		if (fileMappedReads == null) {
			fileMappedReads= createTempFile(
					getNameMappedReads(),
					FluxCapacitorConstants.SFX_BED);
		}

		return fileMappedReads;
	}
		
	public File getFileNotMappedReads() {
		if (fileNotmappedReads == null) {
			fileNotmappedReads = createTempFile(getNameNotMappedReads(), null);
		}

		return fileNotmappedReads;
	}
	
	private void writeProfiles() {
		try {
			long t0= System.currentTimeMillis();
			final String MSG_WRITING_PROFILES= "writing profiles",
						NT= "nt", RPKM= "rpkm", UNDERSCORE= "_";

            Log.progressStart(MSG_WRITING_PROFILES);

			FileOutputStream fos = new FileOutputStream(getFileProfile());
		    ZipOutputStream zos = new ZipOutputStream(fos);

		    UniversalMatrix[] mm= profile.getMasters();
		    for (int i = 0; i < mm.length; i++) {
				String lenString= Integer.toString(mm[i].getLength());
//		    	for (int j = 0; j < mm[i].length; j++) {
					//int expUp= TProfileFunction.BIN_EXP[j];
//					String expString= (j== 0?Integer.toString(Profile.EXP_LO):(j==1?Integer.toString(Profile.EXP_UP):"max"));
					String name= "profile"+ UNDERSCORE+ lenString+ NT;
						//+ UNDERSCORE+ expString;
					ZipEntry ze= new ZipEntry(name);
					zos.putNextEntry(ze);
					zos.write(mm[i].toString().getBytes());
//					zos.write(mm[i][j].toString().getBytes());
					zos.closeEntry();
//				}
		    }
			zos.flush();
			fos.flush();
			zos.close();
            Log.progressFinish(StringUtils.OK, true );

		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	/**
	 * @deprecated
	 * Writes out all profiles, heavy disk activity.
	 * @return
	 */
	private Exception writeTProfiles() {
		try {
			long t0= System.currentTimeMillis();
			final String MSG_WRITING_PROFILES= "writing profiles",
						NT= "nt", RPKM= "rpkm", UNDERSCORE= "_";
            Log.progressStart(MSG_WRITING_PROFILES);

			//TProfile[] t= func.getTProfiles();
//			TSuperProfile[][] supis=
//				func.getMasterProfiles(strand== STRAND_ENABLED, pairedEnd, insertMinMax, readLenMin);
			
			// fileOut.getAbsolutePath()+"_tprofiles.zip"
			FileOutputStream fos = new FileOutputStream(getFileProfile());
		    ZipOutputStream zos = new ZipOutputStream(fos);
//		    for (int i = 0; i < supis.length; i++) {
//		    	int lenUp= TProfileFunction.BIN_LEN[i];
//		    	int lenLo= (i> 0)? TProfileFunction.BIN_LEN[i- 1]: 0;
//		    	int lenMed= lenLo+ ((lenUp- lenLo)/ 2);
//		    	int[][] x= new int[2][];
//		    	for (int j = 0; j < x.length; j++) { 
//					x[j]= new int[lenMed];
//					for (int m = 0; m < x[j].length; m++) 
//						x[j] [m]= 0;
//				}
//				for (int j = 0; j < supis[i].length; j++) {
		    for (int i = 0; i < profileStub.length; i++) {
					try {
						String lenString= i< FluxCapacitorConstants.BIN_LEN.length? Integer.toString(FluxCapacitorConstants.BIN_LEN[i]): "big";
						//int expUp= TProfileFunction.BIN_EXP[j];
						String name= "master"+ UNDERSCORE+ 
							lenString+ NT;
//						+ UNDERSCORE+
//							Integer.toString(expUp)+ RPKM;
//						supis[i][j].project(x);
						ZipEntry ze= new ZipEntry(name);
						zos.putNextEntry(ze);
						StringBuilder buf= new StringBuilder();
						for (int k = 0; k < profileStub[i].length; k++) {
							buf.append(profileStub[i][k]);
							if (strand== FluxCapacitorConstants.STRAND_ENABLED) {
								buf.append('\t');
								buf.append(profileStubRev[i][k]);
							}
							buf.append('\n');
						}
						zos.write(buf.toString().getBytes());
						zos.closeEntry();
					} catch (Exception e) {
						//System.err.println(e);;
						e.printStackTrace();
					}					
			}
		    
/*			int perc= 0;
			for (int i = 0; i < t.length; i++) {
				
				if (i* 10/ t.length> perc) {
					++perc;
					if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) {
						System.err.print("*");
					} else {
						if (Constants.progress!= null) 
							Constants.progress.progress();
					}
				}
				
				try {
					int len= t[i].length();
					String name= t[i].getID()+ UNDERSCORE+ 
						Integer.toString(len)+ NT+ UNDERSCORE+
						Integer.toString((int) calcRPKM(t[i].getReads(), len))+ RPKM;
					ZipEntry ze= new ZipEntry(name);
					zos.putNextEntry(ze);
					zos.write(t[i].m.toByteArray());
					zos.closeEntry();
				} catch (Exception e) {
					if (Constants.verboseLevel> Constants.VERBOSE_ERRORS)
						e.printStackTrace();
					if (Constants.progress!= null)
						Constants.progress.finish(Constants.ERROR, System.currentTimeMillis()- t0);
					return e;
				}
			}
*/				
			
			
			zos.flush();
			fos.flush();
			zos.close();
            Log.progressFinish(StringUtils.OK, true);
			return null;
			
		} catch (Exception e) {
			if (Constants.verboseLevel> Constants.VERBOSE_ERRORS)
				e.printStackTrace();
			return e;
		}
		
	}
	
	public static byte mapFileType= FluxCapacitorConstants.FORMAT_SAM;
	private void writeMapFileSam(Graph g, Edge e, DirectedRegion[] regs, DirectedRegion[][] contRegs) {
		
		return;
		
/*		SAMRecord rec= new SAMRecord(getSammy().getHeader());
		rec.setReadName(regs[0].getID());
		int start= regs[0].get5PrimeEdge(), end= regs[regs.length- 1].get3PrimeEdge();
		if (regs[0].getStrand()< 0) {
			rec.setReadNegativeStrandFlag(true);
			int h= start;
			start= -end;
			end= -h;
		} else
			rec.setReadNegativeStrandFlag(false);
		assert(start<= end&& start> 0&& end> 0);
		rec.setAlignmentStart(start);
		rec.setCigar(Sammy.encodeCigar(regs));
		
		// rec.setAlignmentEnd(end);	// unsupported: derive from cigar
		if (pairedEnd) {
			assert(contRegs[0]!= null);
			start= contRegs[0][0].get5PrimeEdge(); 
			end= contRegs[0][contRegs[0].length- 1].get3PrimeEdge();
			SAMRecord rec2= new SAMRecord(getSammy().getHeader());
			rec2.setReadName(contRegs[0][0].getID());
			if (contRegs[0][0].getStrand()< 0) {
				rec2.setMateNegativeStrandFlag(true);
				int h= start;
				start= -end;
				end= -h;
			} else
				rec2.setMateNegativeStrandFlag(false);
			assert(start<= end&& start> 0&& end> 0);
			rec.setMateAlignmentStart(start);
			rec.setProperPairFlag(true);
			
			rec2.setAlignmentStart(start);
			rec2.setCigar(Sammy.encodeCigar(contRegs[0]));
			
			Transcript[] tt= g.decodeTset(e.getTranscripts());
			HashMap<Integer, String> variantMap= new HashMap<Integer, String>(2);
			for (int i = 0; i < tt.length; i++) {
				int min= Math.min(regs[0].get5PrimeEdge(), contRegs[0][contRegs[0].length- 1].get5PrimeEdge()),
					max= Math.max(regs[0].get3PrimeEdge(), contRegs[0][contRegs[0].length- 1].get3PrimeEdge());
				int estart= tt[i].getExonicPosition(min), eend= tt[i].getExonicPosition(max);
				int isize= eend- estart+ 1;
//				if (regs[0].get5PrimeEdge()< contRegs[0][0].get5PrimeEdge()) {
//					rec.setInferredInsertSize(isize);
//					rec2.setInferredInsertSize(-isize);
//				} else {
//					rec.setInferredInsertSize(-isize);
//					rec2.setInferredInsertSize(isize);
//				}
				if (variantMap.containsKey(isize))
					variantMap.put(isize, variantMap.get(isize)+ "/"+ tt[i].getTranscriptID());
				else
					variantMap.put(isize, tt[i].getTranscriptID());
			}
			StringBuilder sbIsize= new StringBuilder(), sbVariants= new StringBuilder();
			Object[] oo= variantMap.keySet().toArray();
			for (int j = 0; j < oo.length; j++) {
				sbIsize.append((j>0?",":"")+ oo[j]);
				sbVariants.append((j>0?",":"")+ variantMap.get(oo[j]));
			}
			variantMap= null;
			rec.setAttribute(Sammy.OPTION_RNA_ANNOTATION, fileBED.getName());
			rec.setAttribute(Sammy.OPTION_RNA_ISIZE, sbIsize.toString());
			rec.setAttribute(Sammy.OPTION_RNA_IDS, sbVariants.toString());
			synchronized(lock) {
				sammy.addRecord(rec);
				sammy.addRecord(rec2);
			}
		} else
			synchronized(lock) {
				sammy.addRecord(rec);
			}
*/
	}
	
	
	void createBins(int nr, TProfileFunction func) {

		func.getTProfiles();
		
		
	}
	

	
	TProfileFunction func= new TProfileFunction(this); 
	private Vector<Thread> threadPool= new Vector<Thread>();
	int maxThreads= 1;
	Vector<String> origLines= null;
	int checkGTFscanExons= 0;
	boolean copyLocal= false;
	
	boolean 
		outputObs= false, 
		outputPred= false, 
		outputBalanced= true, 
		outputFreq= false, 
		outputRfreq= false, 
		outputRcov= true, 
		outputAll= false, 
		outputSplit= true, 
		outputUnique= false,
		outputExon= false,
		outputUnknown= false,
		outputSJunction= false, 
		outputGene= false, 
		outputTranscript= true, 
		outputEvent= false, 
		outputProfiles= true,	// true
		outputMapped= false,	// true
		outputNotmapped= false,
		outputISize= false,
		outputSorted= true,
		outputLP= false; // true 

	boolean isSortedBED, isSortedGTF, isReadyBED= false, isReadyGTF= false;
	byte compressionBED= FileHelper.COMPRESSION_GZIP, 
		compressionGTF= FileHelper.COMPRESSION_NONE, 
		compressionProfiles= FileHelper.COMPRESSION_GZIP,
		compressionOut= FileHelper.COMPRESSION_NONE;
	
	int eventDim= 2;
	long dbgTimeEmptyGraphs= 0;
	private void solve(Gene gene, BufferedBEDiterator beds, boolean decompose) {
		
		// create LP and solve
		LocusSolver2 lsolver= new LocusSolver2(gene, beds, decompose); 
		if (maxThreads> 1) {
			//Thread outThread= new Thread(lsolver);
			Thread lastThread= getLastThread(); 
			int retry= 0;
			while (threadPool.size()>= maxThreads)
	//			|| (retry< maxThreads&& Runtime.getRuntime().freeMemory()< (0.10* Runtime.getRuntime().maxMemory())))
				try {
					//System.err.println(Runtime.getRuntime().freeMemory()+"<"+ (0.25* Runtime.getRuntime().maxMemory()));
					++retry;
					//System.gc();
					
					//Thread.currentThread().sleep(10); // polling bad
					lastThread.join();
					
	//				if (threadPool.size()< maxThreads&& retry> maxThreads)
	//					break;
				} catch (InterruptedException e) {
					; // :)
				}
			
			lastThread= getLastThread();
			lsolver.setThreadBefore(lastThread);
	//		synchronized(FluxCapacitor.this.threadPool) {
				threadPool.add(lsolver);
	//		}
				lsolver.start(); //;
			
		} else
			lsolver.run();
	}
	
	private Thread getLastThread() {
		synchronized(FluxCapacitor.this.threadPool) {
			if (this.threadPool.size()> 0)
				return this.threadPool.get(this.threadPool.size()- 1);
			else 
				return null;
		}
	}

	
	byte costModel= GraphLPsolver.COSTS_LINEAR, costSplit= 1;
	
	String runID= null;
	
	/*
	 *  The value of this variable will never be cached thread-locally: 
	 *  all reads and writes will go straight to "main memory";
     * Access to the variable acts as though it is enclosed in a 
     * synchronized block, synchronized on itself. 
	 */
	volatile int nrLoci= 0, 
		nrLociExp= 0, 
		nrTx= 0, 
		nrEvents= 0, 
		nrEventsExp= 0, 
		nrTxExp= 0, 
		nrReadsMapped= 0, 
		nrReadsLoci= 0, 
		nrReadsSingleLociMapped= 0, 
		nrReadsSingleLociPairsMapped= 0, 
		nrUnsolved= 0, 
		nrReadsSingleLoci= 0,
		nrSingleTranscriptLoci= 0,	// counted in learn AND decompose redundantly
		nrSingleTranscriptLearn= 0,	
		nrReadsSingleLociNoAnnotation= 0,
		nrReadsSingleLociPotentialPairs= 0,
		nrReadsWrongLength= 0,
		nrMappingsWrongStrand= 0,
		nrMultiMaps= 0,
		nrPairsNoTxEvidence= 0, 
		nrPairsWrongOrientation= 0,
		nrMappingsForced= 0;
	
	int nrMappingsValid= 0,
		nrMappingsP1= 0,
		nrMappingsP2= 0; 	

	
	boolean 
		map= false, 
		decompose= false, 
		uniform= true;
	long nrReadsAll= 0;
	double costModelPar= Double.NaN;
	float[] costBounds= new float[] {0.95f, Float.NaN};	// how much of the original observation can be subs/add
	int[] profileBoundaries;
	private int getBinIdx(int len) {
		int p= Arrays.binarySearch(profileBoundaries, len);
		p= (p<0)?-p-1:p;
		return p;
	}
	
	Profile getProfile() {
		if (uniform) {
			profile= new Profile(this);
			int nr= profile.fill();
		} else {
			if (fileProfile!= null&& fileProfile.exists()) {
				profile= readProfiles(fileProfile);
				if (profile!= null) {
					System.err.println("\tsmoothing..");
					for (int i = 0; i < profile.masters.length; i++) {
						int w= profile.masters[i].sense.length/ 5;
						profile.masters[i].sums=
							Kernel.smoothen(Kernel.KERNEL_EPANECHNIKOV, 
								w, profile.masters[i].sense);
						profile.masters[i].suma=
							Kernel.smoothen(Kernel.KERNEL_EPANECHNIKOV, 
								w, profile.masters[i].asense);
					}
				}
			} 
			if (profile== null) {
				profile= new Profile(this);
				try {
					explore(FluxCapacitorConstants.MODE_LEARN);			
				} catch (Throwable e) {
					e.printStackTrace();
					if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
						System.err.println("[FATAL] Error occured during scanning\n\t"+ e.getMessage());
				}
				writeProfiles();
			}
		}
		
		if (profile== null)
			return null;
		// check
		for (int i = 0; i < profile.masters.length; i++) {
			if (profile.masters[i].hasEmptyPositions())
				profile.masters[i].fill();
		}
		

		return profile;
	}
	
	BinVector isizeV; 
	synchronized void addInsertSize(int isize) {
		isizeV.incrTuple(isize);
	}
	
	private String getAttributeOF(double val, GraphLPsolver solver, int readCount) {

		StringBuilder sb= new StringBuilder(FluxCapacitorConstants.GTF_ATTRIBUTE_PVAL);
		sb.append(" \"");
		if (val > FluxCapacitorConstants.BIG|| solver== null) {
			if (val> FluxCapacitorConstants.BIG)
				System.currentTimeMillis();
			sb.append(" \""+FluxCapacitorConstants.VALUE_NA+ "\";" );
		} else {
			val= val/ (val+ readCount);
			sb.append(StringUtils.fprint(val, 2));
			sb.append("\";");
		}
		
		return sb.toString();
	}
	
	private BufferedWriter writer= null;
	private BufferedWriter getWriter() {
		if (writer == null) {
			try {
				if (getFileOut()== null)
					writer= new BufferedWriter(new OutputStreamWriter(System.out));
				else 
					writer= new BufferedWriter(new FileWriter(getFileOut(), true), bufferSize);
				
			} catch (Exception e) {
				return null;
			}
		}

		return writer;
	}
	
	private ThreadedQWriter qwriter= null;
	private ThreadedQWriter getQWriter() {
		if (qwriter == null) {
			qwriter = new ThreadedQWriter(getWriter());
			qwriter.setLimitBytes(100000000);
			qwriter.start();
		}

		return qwriter;
	}

/*	private Sammy sammy= null;
	private Sammy getSammy() {
		if (sammy == null) {
			sammy = new Sammy(fileBED, fileMappings, false);
			sammy.getHeader();
		}

		return sammy;
	}
*/	
	double getControl(Graph g, Transcript t) {
		Node[] nn= g.getNodesInGenomicOrder();
		long[] part= g.encodeTset(new Transcript[] {t}); // TODO method that takes single transcript
		double sum= 0d;
		for (int i = 0; i < nn.length; i++) {
			for (int j = 0; j < nn[i].getOutEdges().size(); j++) {
				Edge e= nn[i].getOutEdges().elementAt(j);
				if (Graph.isNull(Graph.intersect(e.getTranscripts(), part)))
					continue;
				Transcript[] tt= g.decodeTset(e.getTranscripts());	// TODO method that returns int nr
				sum+= (e.getReadNr()/ (double) tt.length);
				
				for (int k = 0; e.getSuperEdges()!= null&& k < e.getSuperEdges().size(); k++) {
					SuperEdge se= e.getSuperEdges().elementAt(k);
					if (Graph.isNull(Graph.intersect(se.getTranscripts(), part)))
						continue;
					tt= g.decodeTset(se.getTranscripts());
					sum+= (se.getReadNr()/ (double) tt.length);
					
					for (int m = 0; se.getSuperEdges()!= null&& m < se.getSuperEdges().size(); m++) {
						SuperEdge sse= se.getSuperEdges().elementAt(m);
						if (Graph.isNull(Graph.intersect(sse.getTranscripts(), part)))
							continue;
						tt= g.decodeTset(sse.getTranscripts());
						sum+= (sse.getReadNr()/ (double) tt.length);
					}
				}
				
			}
		}
		return sum;
	}
	
	
	// synchronize writes block of output in file !
	void write(StringBuilder sb) throws Exception {
		
		// sync: 680 sec
		// threaded: 680 sec
		// single files: 720 sec
		synchronized (getWriter()) {
			BufferedWriter writer= getWriter();	
			writer.write(sb.toString());
			//writer.flush();
		}
//		getQWriter().add(sb);
//		getQWriter().interrupt();
	}

	private Vector<Vector<Edge>> eeV= new Vector<Vector<Edge>>();
	
	static AtomicLong along= new AtomicLong((0L ^ 0x5DEECE66DL) & ((1L << 48) - 1));
	private static int c= 0;
	private double factor() {
		
		long l = 0;
		for (int i = 26; i <= 27; i++) {
	        long a, b;
	        do {
		    a = along.get();
		    b = (a * 0x5DEECE66DL + 0xBL) & ((1L << 48) - 1);
	        } while (!along.compareAndSet(a, b));
	        int next= (int)(b >>> (48 - i));
	        l+= (i== 26)?((long) next) << 27: next;
		}
        double dbl= (l / (double)(1L << 53))*((c++% 2== 0)?100:1); 
        
		return Math.max(0.001, dbl);
	}
		

	private static final float[] rpkm_1= new float[3], rpkm_2= new float[3];
	Vector<Edge> edgeColl1= new Vector<Edge>(), edgeColl2= new Vector<Edge>();
	int[][] containerIntA1A1= new int[1][];
	{ containerIntA1A1[0]= new int[1]; }
	long[][] containerLongA1A= new long[1][];
	boolean keepTmpSorted= false;
	
	public double getLength(Graph g, Vector<Edge> v, long[] sig, boolean exclusive) {
		double len= 0; 
		for (int i = 0; i < v.size(); i++) {
			Edge e= v.elementAt(i);
			long[] trpts= e.getTranscripts();
			long[] inter= Graph.intersect(trpts, sig);
			if (Graph.isNull(inter)|| (exclusive&& !Graph.equalSet(sig, trpts)))
				continue;
			//len+= v.elementAt(i).length();	// NO, we want possible read pos
			int[] frac= e.getFrac(g.getAnyTranscript(v.elementAt(i).getTranscripts()), readLenMin);
			double len1= frac[1]- frac[0]+ 1;	// not: Math.max(0,...)
			if (len1< 0)
				len1= (frac[1]+ (readLenMin- 1)- frac[0]+ 1)/ (double) readLenMin;
			
			// happens with single reads when:
			// e.g. chr1:101,313,508-101,313,599
			// for the area in uc001dua.2
			if (!(e.getTail().getSite().isLeftFlank()!= e.getTail().getSite().isLeftFlank()
					&& e.getTail().getSite().getPos()+1== e.getHead().getSite().getPos()))
				try {assert(len1>=0);} catch (AssertionError err) {
					if (Constants.verboseLevel> Constants.VERBOSE_NORMAL)
						System.err.println("Found strange length: "+ len1);
					len1= 0;
				}
			
			if (frac.length== 2)
				len+= len1;
			
			// deprecated: paired end, will not happen now.. 
			else {
				for (int j = frac[0]; j <= frac[1]; j++) {
					for (int m = j+ insertMinMax[0]+ readLenMin; 
							m <= j+ insertMinMax[1]+ readLenMin; m++) {
						if (m>= frac[2]&& m <= frac[3])
							++len;
					}
				}
			}
		}
		return len; 
	}
	
	public double getReads(Vector<Edge> v, byte dir, long[] sig, boolean normalized) {
		int sum= 0;
		for (int i = 0; i < v.size(); i++) {
			Edge e= v.elementAt(i);
			long[] inter= Graph.intersect(e.getTranscripts(), sig);
			if (Graph.isNull(inter)|| !e.isExonic())
				continue;
			
			if (pairedEnd) {
				for (int j = 0; e.getSuperEdges()!= null&& j < v.elementAt(i).getSuperEdges().size(); j++) {
					SuperEdge se= e.getSuperEdges().elementAt(j);
					if (!se.isPend())
						continue;
					int cnt= 0;
					for (int k = 0; k < se.getEdges().length; k++) 
						if (se.getEdges()[k]== v.elementAt(i))
							++cnt;
					if (dir>= 0)
						sum+= cnt* se.getReadNr();
					if (dir<= 0)
						sum+= cnt* se.getRevReadNr();
				}
			} else {
				if (dir>= 0)
					sum+= e.getReadNr();
				if (dir<= 0)
					sum+= e.getRevReadNr();
			}
		}
		return sum;
	}
	
	public double getReadsAvg(Vector<Edge> v, byte dir, Graph g, long[] sig, boolean excl, boolean normalized) {
		double sum= 0;
		for (int i = 0; i < v.size(); i++) {
			Edge e= v.elementAt(i);
			long[] trpts= v.elementAt(i).getTranscripts();
			long[] inter= Graph.intersect(trpts, sig);
			if (Graph.isNull(inter)|| (excl&& !Graph.equalSet(sig, trpts))|| !e.isExonic())
				continue;
			double sf= (double) g.decodeCount(v.elementAt(i).getTranscripts());
			int mult=  g.decodeCount(inter);
			
			if (pairedEnd) {
				for (int j = 0; e.getSuperEdges()!= null&& 
						j < e.getSuperEdges().size(); j++) {
					SuperEdge se= e.getSuperEdges().elementAt(j);
					if (!se.isPend())
						continue;
					int cnt= 0;
					for (int k = 0; k < se.getEdges().length; k++) 
						if (se.getEdges()[k]== e)
							++cnt;
					if (dir>= 0)
						sum+= (se.getReadNr()* mult* cnt)/ sf;
					if (dir<= 0)
						sum+= (se.getRevReadNr()* mult* cnt)/ sf;
				}
			} else {
				if (dir>= 0)
					sum+= (e.getReadNr()* mult)/ sf;
				if (dir<= 0)
					sum+= (e.getRevReadNr()* mult)/ sf;
			}
			
			System.currentTimeMillis();
		}
		
		return sum;
	}

	private void append(StringBuilder sb, String s1,
			String s2, String s3, String s4, String s5, String s6, String s7) {
		
		sb.append(s1);
		sb.append(s2);
		sb.append(s3);
		sb.append(s4);
		sb.append(s5);
		sb.append(s6);
		sb.append(s7);
	}

	private Vector<Vector<Edge>> clearEdgeContainer(int nr) {
		for (int i = 0; i< eeV.size()&& i < nr; i++) 
			eeV.elementAt(i).removeAllElements();
		for (int i= eeV.size(); i< nr; ++i)
			eeV.add(new Vector<Edge>());
//		if (containerIntA1A1== null) 
//			containerIntA1A1= new int[1][];
//		if (containerIntA1A1[0]== null) {
//			System.out.println("clear container");
//			containerIntA1A1[0]= new int[1];
//		}
		return eeV;
	}
	
	int[] insertMinMax= null;
	private int[] getExonicPos(Transcript tx, BEDobject2 bed, int tlen) {
		int gstart= bed.getStart();	// getAbsoluteStart();	// fuck 0-base in bed
		int gend= bed.getEnd(); // getAbsoluteEnd();		// last pos not incl in bed, +1-1
		++gstart;	// to normal coordinates
		if (tx.getStrand()< 0) {
			int h= gstart;
			gstart= -gend;
			gend= -h;
		}
		int epos= tx.getExonicPosition(gstart), eposX= tx.getExonicPosition(gend);
		if (epos< 0|| eposX< 0|| epos> tlen|| eposX> tlen)
			return null;	// out of bounds
		
		// TODO check
		int chklen1= bed.getLength();
		int chklen2= eposX- epos+ 1;
		if (chklen1!= chklen2) {
//			if (tx.getStrand()< 0)
//				System.currentTimeMillis();
			return null;	// incompatible split-read 
		}
//		if (bed.getBlockCount()> 1)
//			System.currentTimeMillis();
		
		return new int[] {epos, eposX};
	}

	/**
	 * WARNING: no check whether complete read is contained in transcript
	 * @param tx
	 * @param bed
	 * @return
	 */
	private int getBpoint(Transcript tx, BEDobject2 bed) {
		
		// just depends on genomic position, not on sense/antisense!
		int gpos= bed.getStrand()>= 0? bed.getStart()+ 1: bed.getEnd();	
		int epos= tx.getExonicPosition(gpos);
		
		return epos;
	}
	
	private GFFReader gtfReader;
	public GFFReader getGTFreader() {
		if (gtfReader == null) {
			gtfReader= new GFFReader(fileGTF.getAbsolutePath());
//			if (gtfFirstTime&& (!cheatDisableFCheck)) {
////				System.err.println("[MICHA] Reactivate the sorting filecheck for the gtf.");
//				if (gtfReader.isApplicable()) 
//					isSortedGTF= true;
//				else {
//					isSortedGTF= false;
//					File tmpGTF= gtfReader.createSortedFile();
//					fileGTF= tmpGTF;
//					gtfReader= new GFFReader(fileGTF.getAbsolutePath());					
//				}
//				gtfFirstTime= false;
//			}
				
			gtfReader.setNoIDs(null);
			gtfReader.setReadGene(true);
			gtfReader.setReadFeatures(new String[] {"exon","CDS"});
			gtfReader.setReadAheadTranscripts(1);	// only one locus a time
//			gtfReader.setReadAheadTranscripts(-1);
//			gtfReader.setReadAll(true);
			gtfReader.setGeneWise(true);
			gtfReader.setPrintStatistics(false);
			gtfReader.setReuse(true);
			Transcript.removeGaps= false;
			
			//gtfReader.setReuse(true);
			// chr filter set later on, when bedfile is read
		}
		
		return gtfReader;
	}


	private BEDwrapper bedWrapper; 
	public BEDwrapper getBedReader() {
		if (bedWrapper == null) {
			bedWrapper= new BEDwrapper(fileBED.getAbsolutePath());
		}

		return bedWrapper;
	}
	
	private BufferedWriter writerISize;
	private boolean writeISizes() {
		
		try {
			FileOutputStream fos = new FileOutputStream(getFileISize());
		    ZipOutputStream zos = new ZipOutputStream(fos);
			zos.putNextEntry(new ZipEntry(MyFile.getFileNameWithoutExtension(
					fileISize.getAbsolutePath())));
			BufferedWriter buffy= new BufferedWriter(new OutputStreamWriter(zos));
			buffy.write(isizeV.toString());
			buffy.flush();
			zos.closeEntry();
			zos.flush();
			fos.flush();
			buffy.close();
			return true;
		} catch (Exception e) {
			if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
				e.printStackTrace();
			return false;
		}
	}
	
	private File getFileProfiles(int binIdx) {
		return new File(fileOut.getAbsolutePath()+"_bin_"+profileBoundaries[binIdx]);
	}
	
	

	private int bufferSize= 50000000;
	private BufferedWriter writerMappedReads, writerUnmappedReads;
	private BufferedWriter getWriterMappedReads() {
		if (writerMappedReads == null&& getFileMappedReads()!= null) {
			try {
				writerMappedReads = new BufferedWriter(new FileWriter(fileMappedReads), bufferSize);
			} catch (IOException e) {
				e.printStackTrace();
			}
		}

		return writerMappedReads;
	}
	private BufferedWriter getWriterNotmappedReads() {
		if (writerUnmappedReads == null&& getFileNotMappedReads()!= null) {
			try {
				writerUnmappedReads = new BufferedWriter(new FileWriter(fileNotmappedReads), bufferSize);
			} catch (IOException e) {				
				e.printStackTrace();
			}
		}

		return writerUnmappedReads;
	}
	
	private BufferedWriter writerMappings;
	

	
	private String readLenGuessedFrom;

	int[][] profileStub, profileStubRev;
	/**
	 * profile managing the matrices
	 */
	Profile profile;
	private BufferedBEDiterator readBedFile(Gene gene, int from, int to, byte mode) {
		
		//ByteArrayCharSequence chr= new ByteArrayCharSequence(gene.getChromosome());
		
		if (from> to) {
			System.err.println("reading range error: "+from+","+to);
		}
//		if (gene.getGeneID().equals("chr12:58213712-58240747C")) {
//			System.err.println("\t"+ gene.getGeneID()+" from "+from+" to "+to);
//		}
		
		assert(from>= 0&&to>= 0&&from<= to);
//		for (int i = 0; i < gene.getTranscriptCount(); i++) {
//			if (gene.getTranscripts()[i].getTranscriptID().equals("ENST00000373548")) {
//				System.currentTimeMillis();
//				break;
//			}
//		}

		
		//BEDobject[] beds= getBedReader().read_old(gene.getChromosome(), start, end);
		BEDobject2[] beds= getBedReader().read(gene.getChromosome(), from, to);
		BufferedBEDiterator iter= new BEDiteratorArray(beds);
//		if (gene.getGeneID().equals("chr19:1609293-1652326C"))
//			System.currentTimeMillis();
		if (beds== null)
			return null;
		
		return iter;
		
	}
	
	private int splitBedFile(Gene gene, SyncIOHandler2 handler, OutputStream ostream) {
		
		int start= gene.getStart();
		int end= gene.getEnd();
		if (gene.getStrand()< 0) {
			start= -start;
			end= -end;
		}
		//ByteArrayCharSequence chr= new ByteArrayCharSequence(gene.getChromosome());
		
		assert(start>= 0&&end>= 0&&start<= end);
		
		//BEDobject[] beds= getBedReader().read_old(gene.getChromosome(), start, end);
		return getBedReader().get(gene.getChromosome(), start, end, handler, ostream);
		
	}


	public int addPE(Graph g, int[] insertMinMax, int readLen) {

			// HashMap<String, TProfile> supaMap, 
			
			Edge[] edges= g.getExonicEdgesInGenomicOrder();
			
			int ctr= 0;
			for (int i = 0; i < edges.length; i++) {	// e1
				int p0= edges[i].getHead().getSite().getPos();
				long[] t0= edges[i].getTranscripts();
				for (int j = i; j< edges.length; ++j) {	// e2
					long[] inter_E_E= Graph.intersect(t0, edges[j].getTranscripts());
					if (Graph.isNull(inter_E_E))
						continue;
					int p1= edges[j].getHead().getSite().getPos();
					long[] supp= g.getSupport(edges[i], edges[j], readLen, insertMinMax, inter_E_E);
					// connect edge x edge
					if (!Graph.isNull(supp)&& edges[i].length()>= readLen&& edges[j].length()>= readLen) {
						Edge[] ee= new Edge[] {edges[i], edges[j]};
//						Arrays.sort(ee, g.defaultEdgeCoordComparator);
						g.createPairedEnd(ee, supp);
						++ctr;
						/*double exp= getAllExpectedFracs(g, supaMap, supp, ee, readLen);
						if (exp> 0) {
							g.createPairedEnd(ee, supp);
							++ctr;
						}*/
					}
	
					for (int k = 0; edges[i].getSuperEdges()!= null&& k < edges[i].getSuperEdges().size(); ++k) {	// se1
						SuperEdge se= edges[i].getSuperEdges().elementAt(k);
						if (se.isPend()|| se.getEdges()[0]!= edges[i])
							continue;
						int pSE0= se.getLastEJ();
						long[] t1= se.getTranscripts();
						
						// connect ej x edge
						long[] inter_SE_E= Graph.intersect(se.getTranscripts(), edges[j].getTranscripts());
						if (Graph.isNull(inter_SE_E))
							continue;
						supp= g.getSupport(se, edges[j], readLen, insertMinMax, inter_SE_E);	// pSE0, p1
						if (!Graph.isNull(supp)&& edges[j].length()>= readLen) {
							Edge[] ee= new Edge[] {se, edges[j]};
//							Arrays.sort(ee, g.defaultEdgeCoordComparator);
							g.createPairedEnd(ee, supp);
							++ctr;
							/*double exp= getAllExpectedFracs(g, supaMap, supp, ee, readLen); 
							if (exp> 0) {
								g.createPairedEnd(ee, supp);
								++ctr;
							}*/
						}
	
						for (int m = 0; edges[j].getSuperEdges()!= null&& m < edges[j].getSuperEdges().size(); m++) { // se2
							
							SuperEdge se2= edges[j].getSuperEdges().elementAt(m);
							if (se2.isPend()|| se2.getEdges()[se2.getEdges().length- 1]!= edges[j])
								continue;
							
							long[] inter_E_SE= Graph.intersect(t0, se2.getTranscripts());
							if (Graph.isNull(inter_E_SE))
								continue;
							// connect edge x ej
							int pSE1= se2.getFirstEJ();	//TODOapprox 
							supp= g.getSupport(edges[i], se2, readLen, insertMinMax, inter_E_SE);	// p0, pSE1
							if (!Graph.isNull(supp)&& edges[i].length()>= readLen) {
								Edge[] ee= new Edge[] {edges[i], se2};
//								Arrays.sort(ee, g.defaultEdgeCoordComparator);
								g.createPairedEnd(ee, supp);
								++ctr;
								/*double exp= getAllExpectedFracs(g, supaMap, supp, ee, readLen);
								if (exp> 0) {
									g.createPairedEnd(ee, supp);
									++ctr;
								}*/
							}
							
							// connect ej X ej
							long[] inter_SE_SE= Graph.intersect(t1,se2.getTranscripts());
							supp= g.getSupport(se, se2, readLen, insertMinMax, inter_SE_SE);	// pSE0, pSE1
							if (!Graph.isNull(supp)) {
								Edge[] ee= new Edge[] {se, se2};
//								Arrays.sort(ee, g.defaultEdgeCoordComparator);
								g.createPairedEnd(ee,supp);
								++ctr;
								/*double exp= getAllExpectedFracs(g, supaMap, supp, ee, readLen);							
								if (exp> 0) {
									g.createPairedEnd(ee,supp);
									getAllExpectedFracs(g, supaMap, supp, ee, readLen); // TODO delme
									++ctr;
								}*/
							}
						}
					}
					
						
				}
			}

			return ctr;
		}

	static void readProperties() {
		String wrapper= System.getProperty(Constants.PROPERTY_KEY_WRAPPER_BASE);
		if (wrapper== null)
			wrapper= FluxCapacitorConstants.FNAME_PROPERTIES;
		else
			wrapper+= File.separator+ FluxCapacitorConstants.FNAME_PROPERTIES;
		File pFile= new File(wrapper);		
		Properties props= new Properties();
		try {
			props.load(new FileInputStream(pFile));
		} catch (Exception e) {
			; // :)
		}
	
		if (props.containsKey(FluxCapacitorConstants.PROPERTY_BUILD))
			version= (String) props.get(FluxCapacitorConstants.PROPERTY_BUILD);
		if (props.containsKey(FluxCapacitorConstants.PROPERTY_JDK)) {
			try {
				float v= Float.parseFloat((String) props.get(FluxCapacitorConstants.PROPERTY_JDK));
				String ver= System.getProperty("java.version");
				int p= ver.indexOf('.', 0);
				p= ver.indexOf('.', p+1);
				float v2= Float.parseFloat(ver.substring(0, p));
				if (v2< v) {
					if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
						System.err.println("[BEHIND] Wrong java version, I need "+v2+" but I found "+v);
					System.exit(-1);
				}
			} catch (Exception e) {
				; // :)
			}
			
		}
		
	}

    /**
     * Check operating system and load the native libraries. Exceptions are catched and logged here.
     * Use the return value to check whether loading was successfull
     *
     *
     * @return
     */
	public static int loadLibraries() {
		Log.info("PRE-CHECK","I am checking availability of the required lpsolve JNI libs.");
		VersionInfo lpVer= null;
		try {

            // first check the operating system




			System.loadLibrary(FluxCapacitorConstants.LPSOLVE_LIB_NAME);
			System.loadLibrary(FluxCapacitorConstants.LPSOLVE_JNI_NAME);


			lpVer= LpSolve.lpSolveVersion();
			if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) {
				System.err.println("\t* JNI in java library path");
				System.err.println("\t* successfully loaded lpsolve JNI (version "+lpVer.getMajorversion()+"."+lpVer.getMinorversion()
						+",release "+lpVer.getRelease()+",build "+lpVer.getBuild()+(miss?";":"")+")\n");
			}
			return 0;
		} catch (UnsatisfiedLinkError e) {
			e.printStackTrace();
			if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) 
				System.err.println("\t* there are no lpsolve libraries in the java library path");
		}
		
		return -1;
	}
	
	static int install() {
		
		if (fileJVMdir!= null) {
			File home= fileJVMdir;
			if (!home.exists()) {
				if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
					System.err.println("[NOCOFFEE] Java Home does not exist "+ home);
				System.exit(-1);
			}
			File cmd= new File(home.getAbsolutePath()+ File.separator+ "bin"+ File.separator+ "java");
			if (!cmd.exists()) {
				if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
					System.err.println("[NOJAVA] Java command does not exist "+ cmd);
				System.exit(-1);
			}
			String dir= home.getAbsolutePath();
			try {
				dir= home.getCanonicalPath();
			} catch (Exception e) {
				; // :)
			}
			return install1(dir);
		} else
			return install1(null);
	
	}
	
	static int install1(String jvm) {
		
		if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) 
			System.err.println("[INSTALL] I am installing the program on your system"+ (SystemInspector.checkRuntime()?";":"."));
	
		// guess library path
		String wrapper= System.getProperty(Constants.PROPERTY_KEY_WRAPPER_BASE);
		String basePath= wrapper+ File.separator+ "..";
		if (Constants.verboseLevel>= Constants.VERBOSE_SHUTUP)
			System.err.println("\t* base path: "+ basePath);
		try {
			basePath= new File(basePath).getCanonicalPath();
		} catch (IOException e) {
			if (Constants.verboseLevel>= Constants.VERBOSE_SHUTUP)
				System.err.println("\t* could not create canonical path: "+e.getMessage());
		}	// throws out the /../
		if (Constants.verboseLevel>= Constants.VERBOSE_SHUTUP)
			System.err.println();
		
		String os= SystemInspector.getOSGroupName(); 
		String arch= SystemInspector.getArchGroupName();
		int bits= SystemInspector.getJvmWidth();
		if (Constants.verboseLevel>= Constants.VERBOSE_SHUTUP)
			System.err.println("\t* found OS="+SystemInspector.getOSname()+", ARCH="+SystemInspector.getArchName()+", "+bits+" bit");
		int bits2= getJVMbits(jvm); 
		if (bits2== -1) {
			if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
				System.err.println("\n[AIAIAI] there is a problem with your virtual machine, aborting.");
			System.exit(-1);
		}
		if (bits2< bits) {
			if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
				System.err.println("\t* found JVM with "+bits2+ "bits, switching down.");
			bits= bits2;
		}
		String bit= Integer.toString(bits);
		String nativePath= System.getProperty(FluxCapacitorConstants.CLI_LONG_LIB);
		if (nativePath== null) {
			nativePath= basePath+ File.separator+ FluxCapacitorConstants.SUBDIR_NATIVELIBS+ File.separator+ FluxCapacitorConstants.SUBDIR_LPSOLVE+ File.separator+ 
				os+ File.separator+ arch+ File.separator+ bit;
			try {
				nativePath= new File(nativePath).getCanonicalPath();
			} catch (IOException e) {
				if (Constants.verboseLevel>= Constants.VERBOSE_SHUTUP)
					System.err.println("\t* could not create canonical path: "+e.getMessage());
			}	// throws out the /../
			if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
				System.err.println("\t* guessed library path "+nativePath);
		} else {
			if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
				System.err.println("\t* you provided the path "+nativePath);
		}
		if (Constants.verboseLevel>= Constants.VERBOSE_SHUTUP)
			System.err.println();
	
		// write wrapper scripts
		installWrapper(jvm, basePath, nativePath);
		
		return 0;
	}
	
	private static int getJVMbits(String jvm) {
		
		String cmd= "java -version";
		if (jvm!= null) 
			cmd= jvm+ File.separator+ "bin"+ File.separator+ "java -version";
		
		try {
			Process p= Runtime.getRuntime().exec(cmd);

			BufferedReader buffy= new BufferedReader(new InputStreamReader(p.getErrorStream()));
			int bits= 32;
			boolean ok= false;
			for(String s; (s= buffy.readLine())!= null; ) {
				if (s.contains("java version"))
					ok= true;
				if (s.contains("64-Bit")) {
					bits= 64;
					break;
				}
			}
			buffy.close();

			if (!ok) {
				buffy= new BufferedReader(new InputStreamReader(p.getInputStream()));
				for(String s; (s= buffy.readLine())!= null; ) {
					if (s.contains("java version"))
						ok= true;
					if (s.contains("64-Bit")) {
						bits= 64;
						break;
					}
				}
				buffy.close();
			}

			p.destroy();
			
			if (!ok)
				return -1;
			return bits;
			
		} catch (Exception e) {
			return -1;
		}
		
	}

	static void installWrapper(String jvm, String basePath, String nativePath) {
		try {
			String fp= basePath+ File.separator+ "bin"+ File.separator+ "flux.sh";
			BufferedWriter writer= new BufferedWriter(new FileWriter(fp));
			
			if (SystemInspector.getOSgroup()== SystemInspector.OS_GROUP_WINNT 
				 || SystemInspector.getOSgroup()== SystemInspector.OS_GROUP_VISTA) {
				
				writer.write("SET PATH="+ nativePath+ System.getProperty("path.separator")
						+ "%PATH%"+ System.getProperty("line.separator"));
				
			} else {
				
				byte shell= 0;
				try {
					Process p= Runtime.getRuntime().exec("echo $SHELL");
					BufferedReader buffy= new BufferedReader(new InputStreamReader(p.getInputStream()));
					String s= null;
					while((s= buffy.readLine())!= null) {
						s= s.trim();
						if (s.endsWith("bash"))
							shell= FluxCapacitorConstants.SHELL_BASH;
						else if (s.endsWith("csh"))
							shell= FluxCapacitorConstants.SHELL_CSH;
						else if (s.endsWith("ksh"))
							shell= FluxCapacitorConstants.SHELL_KSH;
					}
					p.destroy();
				} catch (Exception e) {
					; // :)
				}
				
				if (shell== FluxCapacitorConstants.SHELL_CSH)
					writer.write("setenv ");
				else
					writer.write("export ");
				
				if (SystemInspector.getOSgroup()== SystemInspector.OS_GROUP_MACOSX)
					writer.write("DY");
				writer.write("LD_LIBRARY_PATH");
				if (shell== FluxCapacitorConstants.SHELL_CSH)
					writer.write(" ");
				else
					writer.write("=");
				writer.write(nativePath+ System.getProperty("line.separator"));
			}
			if (jvm!= null) {
				if (SystemInspector.getOSgroup()== SystemInspector.OS_GROUP_WINNT 
						 || SystemInspector.getOSgroup()== SystemInspector.OS_GROUP_VISTA)
					writer.write("set JAVA_HOME="+ jvm);
				else 
					writer.write("export JAVA_HOME="+ jvm);
			}
			writer.write(System.getProperty("line.separator"));
			
			String cmd= jvm== null? "java": jvm+ File.separator+ "bin"+ File.separator+ "java";
			writer.write(
				cmd
//				+ " -Xmx1500M"
//				+ " -XX:+AggressiveHeap"
				// http://blogs.sun.com/partnertech/entry/a_short_primer_to_java
				// thasso:
				// http://kirk.blog-city.com/advice_on_jvm_heap_tuning_dont_touch_that_dial.htm
				+ " -Xms500m -Xmx"+ (SystemInspector.getJvmWidth()> 32?"12G": "8G") 
				+ " -XX:MaxNewSize=1500m -XX:NewSize=120m" 
				+ " -XX:+UseParNewGC" 
				+ " -XX:+UseConcMarkSweepGC " 
				+ " -XX:+CMSParallelRemarkEnabled"
				+ " -XX:TargetSurvivorRatio=90"
				+ " -Djava.library.path=\""+ nativePath
				+ "\" -jar \""+ basePath+ File.separator+ "lib"+ File.separator+ "FluxCapacitor.jar\"");
			if (SystemInspector.getOSgroup()== SystemInspector.OS_GROUP_WINNT 
					 || SystemInspector.getOSgroup()== SystemInspector.OS_GROUP_VISTA) 
				writer.write(" %*");
			else
				writer.write(" $@");
			writer.write(System.getProperty("line.separator"));
			writer.flush();
			writer.close();
			
			if (!(SystemInspector.getOSgroup()== SystemInspector.OS_GROUP_WINNT 
					 || SystemInspector.getOSgroup()== SystemInspector.OS_GROUP_VISTA)) {
				try {
					Process p= Runtime.getRuntime().exec("chmod a+x "+fp);
					p.waitFor();
				} catch (Exception e) {
					; // :)
				}
			}
			
		} catch (Exception e) {
			if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
				e.printStackTrace();
		}
	}

	public TProfileFunction getFunc() {
		return func;
	}

	public void setFunc(TProfileFunction func) {
		this.func = func;
	}

	public int getNrSingleTranscriptLoci() {
		return nrSingleTranscriptLoci;
	}

	public int getNrReadsSingleLoci() {
		return nrReadsSingleLoci;
	}

	public int getNrReadsSingleLociMapped() {
		return nrReadsSingleLociMapped;
	}

	boolean fileInit() {
			if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
				System.err.println("\n[INITING] preparing input/output files");
			
			// init file names
			File f= new File(System.getProperty(Constants.PROPERTY_TMPDIR));
			if (!f.exists()) {
				if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
					System.err.println("[AHUUU] Temporary directory does not exist:\n\t"
							+ System.getProperty(Constants.PROPERTY_TMPDIR));
				return false;
			}
			
			boolean returnVal= fileInitReference();
			if (!returnVal)
				return false;
			
			returnVal= fileInitBED();
			if (!returnVal)
				return false;
			
			if (copyLocal) {
				returnVal= fileInitLocalOutput();
				if (!returnVal)
					return false;
			}
			
			return true;
		}

	private boolean fileInitLocalOutput() {
		
		if (fileOut!= null) {
			fileOUToriginal= fileOut;
			String name= MyFile.getFileNameOnly(fileOut.getAbsolutePath());
			String ext= MyFile.getExtension(fileOut.getAbsolutePath());
			fileOut= createTempFile(name, ext);
		}
			
		return true;
	}

	public File getFileBED() {
		return fileBED;
	}

	public File getFileGTF() {
		return fileGTF;
	}

	private BEDobject2[] readBedFile(Gene gene, byte mode) {
		
		int start= gene.getStart();
		int end= gene.getEnd();
		if (gene.getStrand()< 0) {
			start= -start;
			end= -end;
		}
		//ByteArrayCharSequence chr= new ByteArrayCharSequence(gene.getChromosome());
		
		assert(start>= 0&&end>= 0&&start<= end);
		
		//BEDobject[] beds= getBedReader().read_old(gene.getChromosome(), start, end);
		BEDobject2[] beds= getBedReader().read(gene.getChromosome(), start, end);
		if (beds== null)
			return null;
		
		return beds;
		
	}

	FluxCapacitorParameters pars= null;
	void init2(String[] args) {

		if (args== null|| args.length< 1) {
			if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) {
				System.err.println("[MISSING] Please specify a parameter file.");
			}
			System.exit(-1);
		}
		
		File f= new File(args[0]);
		if (!f.exists()) {
			if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) {
				System.err.println("[UPS] Parameter file does not appear to exist "+ f.getAbsolutePath());
			}
		}
		FluxCapacitorParameters pars= FluxCapacitorParameters.create(f);
		if (pars== null|| !pars.check()) 
			System.exit(-1);
		this.fileGTForiginal= pars.fileAnnotation;
		this.fileBEDoriginal= pars.fileMappings;
		this.fileProfile= pars.fileProfile;
		if (fileProfile!= null)
			uniform= false;
		this.pairedEnd= pars.pairedEnd;
		this.stranded= pars.stranded;
		this.descriptor2= pars.descriptor;
		this.fileOut= pars.fileStdOut;
		this.copyLocal= pars.ioInTemp;
		if (pars.fileLPzip!= null) {
			//this.fileLPdir= pars.fileLP;
			this.outputLP= true;
		}
		
		this.outputSorted= pars.fileMappingsSorted!= null; 
		this.pars= pars;
	}

	public boolean explore(byte mode) {
	
			nrSingleTranscriptLoci= 0;
			nrReadsLoci= 0;
			nrReadsMapped= 0; 
			nrReadsWrongLength= 0;
			nrMappingsWrongStrand= 0;
			
			BEDobject2[] leftover= null;
			
			SyncIOHandler2 handler= new SyncIOHandler2(10* 1024* 1024);
			
			if (mode== FluxCapacitorConstants.MODE_LEARN) {
				nrReadsSingleLoci= 0;
				nrReadsSingleLociMapped= 0;
			} 
			
			//System.out.println(System.getProperty("java.library.path"));
			long t0= System.currentTimeMillis();
			try {
				
				Transcript.setEdgeConfidenceLevel(Transcript.ID_SRC_MOST_INCONFIDENT);
				//this.gtfReader= null;
				//GFFReader gtfReader= getGTFreader();
				gtfReader.reset();
				bedWrapper.reset();
				
				if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) {
					if (mode== FluxCapacitorConstants.MODE_LEARN) 
						System.err.println("\n[LEARN] Scanning the input and getting the attributes.");
					else if (mode== FluxCapacitorConstants.MODE_RECONSTRUCT)
						System.err.println("\n[SOLVE] Deconvolving reads of overlapping transcripts.");
				}
				final String profiling= "profiling ", decomposing= "decomposing "; 
	
	            if (mode== FluxCapacitorConstants.MODE_LEARN)
	                Log.progressStart(profiling);
	            else if (mode== FluxCapacitorConstants.MODE_RECONSTRUCT)
	                Log.progressStart(decomposing);
	
	
				
				if (mode== FluxCapacitorConstants.MODE_LEARN) 
					gtfReader.setKeepOriginalLines(false);
				else if (mode== FluxCapacitorConstants.MODE_RECONSTRUCT)
					gtfReader.setKeepOriginalLines(true);
				
				getGTFreader().read();
				Gene[] gene= null, geneNext= getGTFreader().getGenes();
				
				long tlast= System.currentTimeMillis();
				boolean output= false;
		
				String lastChr= null; 
				byte lastStr= 0;
				int lastEnd= -1;
				int tol= this.tolerance; // 1000
	
				if (geneNext!= null) {
					lastChr= geneNext[0].getChromosome();
					lastStr= geneNext[0].getStrand();
				}
				
				Thread readerThread= null;
				int readObjects= 0;
				while (lastChr!= null) {	// MAIN LOOP
					
					
					if ((gene= geneNext)== null)
						break;
					if (mode== FluxCapacitorConstants.MODE_RECONSTRUCT)
						origLines= (Vector<String>) getGTFreader().getVLines().clone();	// TODO make array, trim..
					
					// http://forums.sun.com/thread.jspa?threadID=5171135&tstart=1095
					if (readerThread== null)
						readerThread= new GTFreaderThread();
					//readerThread.start();
					readerThread.run();
	//				while (readerThread!= null&& readerThread.isAlive())
	//					try {
	//						readerThread.join();
	//					} catch (InterruptedException e) {
	//						; // :)
	//					}
					geneNext= getGTFreader().getGenes();
	
					for (int i = 0; (gene!= null)&& i < gene.length; i++) {
						
						//System.gc();
						//Thread.yield();
	//					if (i>= 1500) { 
	//						int c= 0;
	//						while (c!= '\n') {
	//							System.err.println("start?");
	//							c= System.in.read();
	//						}
	//					}
							
						
						// flop strand
						if (lastChr.equals(gene[i].getChromosome())) {
							if (lastStr!= gene[i].getStrand()) {
								//System.err.println(lastChr+" "+lastStr+ " "+ readObjects+ " wrote "+ dbgCntWriteMap +" not "+ dbgCntWriteNonmap);
								readObjects= 0;	
								leftover= null;
								// jump back
								getBedReader().reset(gene[i].getChromosome());
								lastStr= gene[i].getStrand();
								lastEnd= -1;
							}
						} else {						// flop chr
							//System.err.println(lastChr+" "+lastStr+ " "+ readObjects+ " wrote "+ dbgCntWriteMap +" not "+ dbgCntWriteNonmap);
							readObjects= 0;
							leftover= null;
							lastChr= gene[i].getChromosome();
							lastStr= gene[i].getStrand();
							lastEnd= -1;
						}
					
	//					for (int j = 0; j < gene[i].getTranscripts().length; j++) {
	//						if (gene[i].getTranscripts()[j].getTranscriptID().equals("ENST00000391372"))
	//							System.currentTimeMillis();
	//					}
						
						if (gene[i].getTranscriptCount()== 1)
							++nrSingleTranscriptLoci;
						else if (mode== FluxCapacitorConstants.MODE_LEARN)
							continue;	// performance for not reading beds
						
						BufferedBEDiterator beds= null;
	
	/*					File f= File.createTempFile("fluxpfx", ".bed");
						FileOutputStream fos= new FileOutputStream(f);
						handler.addStream(fos);
						fileBED= f;
						bedWrapper= null;
	*/					
						// boundaries
						int start= gene[i].getStart();
						int end= gene[i].getEnd();
						assert(geneNext== null|| geneNext.length== 1);
						
						if (gene[i].getStrand()< 0) {
							start= -start;
							end= -end;
						}
						tol= 0;
						start= Math.max(1, start- tol);
						end= end+ tol;					
	/*					if (lastEnd< 0)
							start= Math.max(1, start- tol);
						else {
							start= Math.max(lastEnd+ 1, start- tol);
						}
						if (geneNext== null|| (!geneNext[0].getChromosome().equals(gene[i].getChromosome()))
								|| (geneNext[0].getStrand()!= gene[i].getStrand()))
							end+= tol;
						else {
							int next= Math.abs(geneNext[0].getStart());
							end= Math.min(end+ tol, end+ ((next- end)/ 2));
						}
						lastEnd= end;
	*/			
						
	//					if (false&& geneNext[0].getGeneID().equals("chr19:1609293-1652326C"))
	//						System.currentTimeMillis();
	
						beds= readBedFile(gene[i], start, end, mode);
						
	/*					if (false&& leftover!= null) {
							BEDobject2[] nuBeds= 
								new BEDobject2[leftover.length+ (beds== null? 0: beds.length)];
							System.arraycopy(leftover, 0, nuBeds, 0, leftover.length);
							if (beds!= null) 
								System.arraycopy(beds, 0, nuBeds, leftover.length, beds.length);
							beds= nuBeds;
							leftover= null;
						}
	*/					
					
	//					if (geneNext[0].getGeneID().equals("chr12:58213712-58240747C"))
	//						System.currentTimeMillis();
						
						if (beds!= null) {
							
	/*						if (false&& geneNext!= null&& geneNext[0].getChromosome().equals(gene[i].getChromosome())
									&& geneNext[0].getStrand()== gene[i].getStrand()) {
	
								int bp= Math.abs(geneNext[0].getStart())- tol;
								if (bp< end) {
									int p= beds.length- 1;
									for(;p>= 0;--p) {
										if (beds[p].getStart()< bp)
											break;
									}
									if (p< 0)
										p= 0;	// take all
									leftover= new BEDobject2[beds.length- p];
									for (int j = p; j < beds.length; j++) 
										leftover[j- p]= beds[j];
									readObjects+= beds.length- p;
								}
							} else
	*/						 	//TODO beds.size() no longer available
								//readObjects+= beds.size();
	//						if (beds.length> 0&& mode== MODE_RECONSTRUCT)
	//							System.err.println(gene[i].toUCSCString()+ " "+ beds.length);
						}
						
	//					if (i>= 1500)
	//						System.err.println("read "+beds.length+" objects");
						
						if (mode== FluxCapacitorConstants.MODE_LEARN&& beds!= null) {
	//						if (Constants.progress!= null) 
	//							Constants.progress.setString(profiling+ gene[i].getGeneID());
							solve(gene[i], beds, false);
						}
						else if (mode== FluxCapacitorConstants.MODE_RECONSTRUCT) {
	
							// check length
	//						for (int k = 0; readLenMin> 0&& beds!= null&& k < beds.length; k++) {
	//							int tmpLen= beds[k].getLength();
	//							
	//							if (tmpLen!= readLenMin) {
	//								
	//								++nrReadsWrongLength;								
	//								//int diff= tmpLen- readLenMin;	// was set to min
	//								beds[k]= null;
	//								/*boolean b= beds[k].trim(beds[k].getStrand()< 0, diff); // (-) always from start, (+) always from end, regardless gene.strand
	//								if (!b) {
	//									if (Constants.verboseLevel> Constants.VERBOSE_NORMAL)
	//										System.err.println("[HEY] mapping length "+tmpLen+" < minReadLen "+readLenMin+"!");
	//									beds[k]= null;
	//								} else try {assert(beds[k].length()== readLenMin);} catch (AssertionError err) {
	//									if (Constants.verboseLevel> Constants.VERBOSE_NORMAL)
	//										System.err.println("[OOPS] failed to trim bed from "+tmpLen+" to "+readLenMin+":\n\t"+ beds[k]);
	//									beds[k]= null;
	//								}*/
	//								
	//							}
	//								
	//						}
							
	//						if (Constants.progress!= null) {
	//							Constants.progress.setString(decomposing);	// + gene[i].getGeneID()
	//						}
							
							solve(gene[i], beds, true); 
						}
							
						if (output) {
							System.out.println(gene[i].getChromosome()+ " "+
									gene[i].getStrand()+
									" cluster "+ gene[i].getGeneID());
									// TODO beds.size() no longer available
									//", "+beds.size()+" reads.");
							if ((lastStr!= gene[i].getStrand()
									||!(lastChr.equals(gene[i].getChromosome())))) {
								long t= System.currentTimeMillis();
								if (lastStr!= 0&& (!lastChr.equals("")))
									System.out.println(lastChr+" "+lastStr+
											" "+((t-tlast)/1000)+" sec.");
								tlast= t;
								lastStr= gene[i].getStrand();
								lastChr= gene[i].getChromosome();
							}		
						}
					}
					//getWriter().flush();
					
					
				}	// end iterate GTF
				
				getBedReader().finish();
				
				while (threadPool.size()> 0&& threadPool.elementAt(0).isAlive())
					try {
						threadPool.elementAt(0).join();
					} catch (Exception e) {
						; //:)
					}
	            Log.progressFinish(StringUtils.OK, true);
				if (checkGTFscanExons> 0&& checkGTFscanExons!= getGTFreader().getNrExons())
					System.err.println("[ERROR] consistency check failed in GTF reader: "+ checkGTFscanExons+ "<>"+ getGTFreader().getNrExons());
				checkGTFscanExons= getGTFreader().getNrExons(); 
				if (checkBEDscanMappings> 0&& checkBEDscanMappings!= getBedReader().getNrLines())
					System.err.println("[ERROR] consistency check failed in BED reader "+ checkBEDscanMappings+ "<>"+ getBedReader().getNrLines());
				//checkBEDscanMappings= getBedReader().getNrLines();
				if (mode== FluxCapacitorConstants.MODE_LEARN) {
	
					if (pairedEnd&& outputISize)
						writeISizes();
					if (pairedEnd&& func.getTProfiles()!= null) {
						insertMinMax= getInsertMinMax();
					}
	
					if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) {
						
						//System.err.println(" OK.");
						System.err.println("\tfirst round finished .. took "+((System.currentTimeMillis()- t0)/ 1000)+ " sec.\n\n\t"
								+ nrSingleTranscriptLoci+" single transcript loci\n\t"							
								+ getBedReader().getNrLines()+ " mappings in file\n\t"
								+ nrReadsSingleLoci+" mappings fall in single transcript loci\n\t"	// these loci(+/-"+tolerance+"nt)\n\t"
								+ nrReadsSingleLociMapped+" mappings map to annotation\n\t"
								+ ((strand== FluxCapacitorConstants.STRAND_SPECIFIC)?nrMappingsWrongStrand+" mappings map to annotation in antisense direction,\n\t":"")
								//+ (pairedEnd?(nrReadsSingleLociPotentialPairs+ " mappings form potential pairs,\n\t"):"")
								+ (pairedEnd?(nrReadsSingleLociPairsMapped* 2)+" mappings in annotation-mapped pairs\n\t":"")
								//+ nrReadsSingleLociNoAnnotation+ " mappings do NOT match annotation,\n\t"
								//+ (uniform?"":func.profiles.size()+" profiles collected\n\t")
								+ readLenMin+ ","+ readLenMax+ " min/max read length\n\t"							
								+ (pairedEnd&& insertMinMax!= null?insertMinMax[0]+","+insertMinMax[1]+" min/max insert size\n\t":""));
						//nrUniqueReads= getBedReader().getNrUniqueLinesRead();
						//System.err.println("\ttotal lines in file "+nrUniqueReads);
						System.err.println();
					}
					
				} else if (mode== FluxCapacitorConstants.MODE_RECONSTRUCT) {
					while (threadPool.size()> 0&& threadPool.elementAt(0).isAlive())
						try {
							threadPool.elementAt(0).join();
						} catch (Exception e) {
							; //:)
						}
					getWriter().flush();
					getWriter().close();
	
	//					if (fileMappings!= null)
	//						getSammy().close();
					
					//assert(nrUniqueReads==getBedReader().getNrUniqueLinesRead());	// take out for cheat
					if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) {
						System.err.println();
						System.err.println("\treconstruction finished .. took "+((System.currentTimeMillis()- t0)/ 1000)+ " sec.\n\n\t"
								+ getBedReader().getNrLines()+" mappings read from file\n\t"
								// no info, reads in redundantly many reads
								//+ nrReadsLoci+" mappings in annotated loci regions\n\t"
								+ nrReadsMapped+ " mapping"+ (pairedEnd?" pairs":"s") +" map to annotation\n"
								+ (pairedEnd?
									"\t"+ nrPairsNoTxEvidence+ " pairs without tx evidence\n"
									+ "\t"+ nrPairsWrongOrientation+ " pairs in wrong orientation\n"
									+ "\t"+ nrMappingsForced+ " single mappings forced\n":"")
								+ ((strand== FluxCapacitorConstants.STRAND_SPECIFIC)?nrMappingsWrongStrand+" mappings map to annotation in antisense direction\n\t":"")
								//+ nrMultiMaps+" mapped multiply.\n\n\t"
								+ (outputGene?"\n\t"+ nrLoci+ " loci, "+ nrLociExp+ " detected":"")
								+ (outputTranscript?"\n\t"+nrTx+" transcripts, "+nrTxExp+" detected":"")
								+ (outputEvent?"\n\t"+ nrEvents+" ASevents of dimension "+eventDim+", "+nrEventsExp+" detected":"")
								+ "\n"
								//+ nrUnsolved+" unsolved systems."
								);
					}					
				}
				
			} catch (Exception e1) {
				e1.printStackTrace();
				return false;
			}
	
			return true;
		}

}
 