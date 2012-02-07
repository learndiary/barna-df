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

package barna.flux.capacitor.reconstruction;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.PipedInputStream;
import java.io.PipedOutputStream;
import java.io.PrintStream;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.Date;
import java.util.Enumeration;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Properties;
import java.util.Vector;
import java.util.concurrent.atomic.AtomicLong;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;
import java.util.zip.ZipOutputStream;

import lpsolve.LpSolve;
import lpsolve.VersionInfo;

import org.cyclopsgroup.jcli.ArgumentProcessor;
import org.cyclopsgroup.jcli.annotation.Cli;
import org.cyclopsgroup.jcli.annotation.Option;

import barna.commons.Execute;
import barna.commons.launcher.CommandLine;
import barna.commons.launcher.FluxTool;
import barna.commons.launcher.HelpPrinter;
import barna.commons.log.Log;
import barna.commons.system.SystemInspector;
import barna.commons.thread.SyncIOHandler2;
import barna.commons.utils.StringUtils;
import barna.flux.capacitor.graph.AnnotationMapper;
import barna.flux.capacitor.graph.MappingsInterface;
import barna.flux.capacitor.reconstruction.FluxCapacitorSettings.AnnotationMapping;
import barna.genome.lpsolver.LPSolverLoader;
import barna.io.AbstractFileIOWrapper;
import barna.io.AnnotationWrapper;
import barna.io.BufferedIterator;
import barna.io.BufferedIteratorDisk;
import barna.io.BufferedIteratorRAM;
import barna.io.FileHelper;
import barna.io.MappingWrapper;
import barna.io.bed.BEDDescriptorComparator;
import barna.io.bed.BEDwrapper;
import barna.io.gtf.GTFwrapper;
import barna.io.rna.UniversalReadDescriptor;
import barna.io.rna.UniversalReadDescriptor.Attributes;
import barna.io.state.MappingWrapperState;
import barna.model.ASEvent;
import barna.model.DirectedRegion;
import barna.model.Exon;
import barna.model.Gene;
import barna.model.Transcript;
import barna.model.bed.BEDobject2;
import barna.model.commons.Coverage;
import barna.model.commons.MyFile;
import barna.model.constants.Constants;
import barna.model.gff.GFFObject;
import barna.model.splicegraph.AbstractEdge;
import barna.model.splicegraph.SimpleEdge;
import barna.model.splicegraph.SplicingGraph;
import barna.model.splicegraph.SuperEdge;



/**
 * Flux Tool that implements the simulation pipeline.
 *
 * @author Micha Sammeth (gmicha@gmail.com)
 *
 */
@Cli(name = "capacitor", description = "Flux Capacitor")
public class FluxCapacitor implements FluxTool<Void>, ReadStatCalculator {

	public static enum SupportedFormatExtensions {
		GTF, GFF, BED,
	}
		
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
	
	public static boolean debug= false;
	public static boolean outputPbClusters= false;
	public boolean pairedEnd= false, stranded= false, force= false;
	int tolerance= 1000;	// +/- gene-near region
	public static boolean 
		cheatDoNotExit= false,
		cheatLearn= false, 
		cheatDisableFCheck= false,
		cheatEnableCleanup= false,
		cheatCopyFile= false,
		doUseLocusNormalization= false;
	
	static void exit(int code) {
		String pfx= "[ASTALAVISTA] ";
		if (code< 0)
			pfx= "[CIAO] ";
		if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) 
			System.err.println(pfx+ "I'm exiting.");
		System.exit(code);
	}
	
	private static class StringArrayByFirstComparator implements Comparator<String[]> {
		public int compare(String[] o1, String[] o2) {
			return o1[0].compareTo(o2[0]);
		}
	}
	
	class LocusSolver extends Thread {
			Gene gene= null;
			ASEvent[] events= null;
			BufferedIterator beds= null;
			boolean decompose= false;
			Thread threadBefore= null;
			int nrMappingsReadsOrPairs;
			HashSet<CharSequence> mapReadOrPairIDs;
			HashMap<CharSequence, Vector<BEDobject2>[]> mapEndsOfPairs;
			long[] sigall= null;
			UniversalReadDescriptor.Attributes attributes= null;
	
			private float invariantTestObsSplitFreq= 0, invariantTestPredSplitFreq= 0;
			
			public LocusSolver(Gene newGene, BufferedIterator newBeds, boolean decompose) {
				//super(newGene.getGeneID());
				
				this.gene= newGene; 
				this.beds= newBeds;
				this.decompose= decompose;
				
				nrMappingsReadsOrPairs= 0;
				mapReadOrPairIDs= new HashSet<CharSequence>();
				if (pairedEnd)
					mapEndsOfPairs = new HashMap<CharSequence, Vector<BEDobject2>[]>();
				attributes= settings.get(FluxCapacitorSettings.READ_DESCRIPTOR).createAttributes();
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
							//SpliceGraph myGraph= getGraph(this.gene);
							AnnotationMapper mapper= new AnnotationMapper(this.gene);
							//map(myGraph, this.gene, this.beds); 
							mapper.map(this.beds, settings);
							nrReadsLoci+= mapper.nrMappingsLocus;
							nrReadsMapped+= mapper.getNrMappingsMapped();
							nrMappingsReadsOrPairs+= mapper.getNrMappingsMapped()/ 2;
							nrPairsNoTxEvidence+= mapper.getNrMappingsNotMappedAsPair();
							nrPairsWrongOrientation+= mapper.getNrMappingsWrongPairOrientation();
							
							GraphLPsolver mySolver= null;
							// != mapReadOrPairIDs.size()> 0, does also count singles
//							if (nrMappingsReadsOrPairs> 0&& this.gene.getTranscriptCount()> 1) {	// OPTIMIZE			
//								mySolver= getSolver(myGraph, nrMappingsReadsOrPairs* 2); // not: getMappedReadcount()
//								mySolver.run();
//							}
							if (mapper.nrMappingsMapped> 0&& this.gene.getTranscriptCount()> 1) {	// OPTIMIZE			
								mySolver= getSolver(mapper, (int) (mapper.nrMappingsMapped* 2)); // not: getMappedReadcount()
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
	
			
			SplicingGraph getGraph(Gene gene) {
					boolean output= false;
					
					// construct graph
				long t0= System.currentTimeMillis();
				
				SplicingGraph myGraph= new SplicingGraph(gene);
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
	
			private int nrLocusMultimaps= 0;
			public int getMappedReadcount() {
				if (pairedEnd)
					return mapReadOrPairIDs.size()/ 2;
				return mapReadOrPairIDs.size();
			}
			
	
			private void outputGFF(AnnotationMapper g, ASEvent[] events, GraphLPsolver solver) {

				// check locus
				++nrLoci;
				if (solver!= null|| nrMappingsReadsOrPairs> 0) 
					++nrLociExp;
				double valOF= solver== null?0: solver.getValObjFunc();
				if (valOF> FluxCapacitorConstants.BIG) { 
					++nrUnsolved;
					Log.warn("Unsolved system: "+ gene.getGeneID());
				}

				// pre-build rpkm hash
				HashMap<String, Double> rpkmMap= null;
				double base= (nrBEDreads< 0? 1: nrBEDreads);
				Transcript[] tt= gene.getTranscripts();
				if (outputBalanced) {
					rpkmMap= new HashMap<String, Double>(tt.length, 1f);
					for (int i = 0; i < tt.length; i++) {
						Transcript tx= tt[i];
						String tid= tt[i].getTranscriptID();
						
						double val= 0d;
						if (solver== null) 									
							val= nrMappingsReadsOrPairs* 2;
						else {
							val= solver.getTrptExprHash().get(tid).doubleValue();
							if (val< 1- costBounds[0]) // 1- 0.95
								val= 0;
						}

						if (val> 0&& !(outputObs|| outputPred))
							++nrTxExp;

						double rpkm= (float) ((val/ (double) tx.getExonicLength())* (1000000000l/ base));
						if (Double.isNaN(rpkm))
							Log.warn("NaN RPKM produced: "+ val+ " / "+ base+ " = "+ rpkm);

						rpkmMap.put(tid, rpkm);
					}
				}
				
				
				// reproduce original
				boolean foundExons= true, foundTranscripts= false;
				if (getGTFreader().isKeepOriginalLines()&& origLines!= null) {
					foundTranscripts= outputGFForiginalLines(g, events, solver, rpkmMap);
				}
				
				StringBuilder sb= new StringBuilder();
				// LOCUS TODO genes 
				if (outputGene) {
					if (outputObs|| outputPred) {
						//getGTF(sb, g.trpts[0].getGene(), g, solver, perM, pv);	
						try {assert(testInvariant(invariantTestObsSplitFreq, 
								pairedEnd?nrMappingsReadsOrPairs*2:nrMappingsReadsOrPairs, 0.05));}	// min: 0.01
						catch (AssertionError e) {
							Log.warn(getClass().getName()+".outputGFF():\n\tinvariantTestObsSplitFreq= "
										+ invariantTestObsSplitFreq+ ", nrMappingsReadsOrPairs= "+ (pairedEnd?nrMappingsReadsOrPairs*2:nrMappingsReadsOrPairs)
										+ "\n\tlocus: "+ g.trpts[0].getTranscriptID());
						};
						try {assert(testInvariant(invariantTestPredSplitFreq, 
								pairedEnd?nrMappingsReadsOrPairs*2:nrMappingsReadsOrPairs, 0.1));}
						catch (AssertionError e) {
							Log.warn(getClass().getName()+".outputGFF():\n\tinvariantTestPredSplitFreq= "
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
						String tid= tt[i].getTranscriptID();
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
									sb.append(" ");
								
								sb.append(FluxCapacitorConstants.GTF_ATTRIBUTE_TOKEN_READS);
								sb.append(" ");
								sb.append(String.format("%1$f", 
										(float) (rpkmMap.get(tid)* tt[i].getExonicLength()* (base/ 1000000000l))));
								sb.append("; ");
								
								sb.append(FluxCapacitorConstants.GTF_ATTRIBUTE_TOKEN_RPKM);
								sb.append(" ");							
								//sb.append(rpkmMap.get(g.trpts[i].getTranscriptID()));
								// avoid scientific notation								
								sb.append(String.format("%1$f", rpkmMap.get(tid).floatValue()));
								sb.append("\n");
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
							Vector<Vector<AbstractEdge>> eeV= new Vector<Vector<AbstractEdge>>(5,5);
							eeV.add(new Vector<AbstractEdge>());
							g.getRPK(tt[i], pairedEnd, SplicingGraph.ETYPE_SJ, eeV);
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
								Log.warn(getClass().getName()+".outputGFF():\n\tinvariantObsEx= "
											+ invariantObsEx+ ", invariantObsTx= "+ invariantObsTx
											+ "\n\tlocus: "+ tt[0].getTranscriptID());
							};
							try {assert(testInvariant(invariantPredEx, invariantPredTx, 0.1));}
							catch (AssertionError e) {
								Log.warn(getClass().getName()+".outputGFF():\n\tinvariantPredEx= "
											+ invariantPredEx+ ", invariantPredTx= "+ invariantPredTx
											+ "\n\tlocus: "+ tt[0].getTranscriptID());
							};
						}
					}
					if (outputTranscript) {
						try {assert(testInvariant(invariantObsAllTx, 
								pairedEnd?nrMappingsReadsOrPairs*2:nrMappingsReadsOrPairs, 0.05));}	// min: 0.01
						catch (AssertionError e) {
							Log.warn(getClass().getName()+".outputGFF():\n\tinvariantObsAllTx= "
										+ invariantObsAllTx+ ", nrMappingsReadsOrPairs= "+ (pairedEnd?nrMappingsReadsOrPairs*2:nrMappingsReadsOrPairs)
										+ "\n\tlocus: "+ tt[0].getTranscriptID());
						};
						try {assert(testInvariant(invariantPredAllTx, 
								pairedEnd?nrMappingsReadsOrPairs*2:nrMappingsReadsOrPairs, 0.1));}
						catch (AssertionError e) {
							Log.warn(getClass().getName()+".outputGFF():\n\tinvariantPredAllTx= "
										+ invariantPredAllTx+ ", nrMappingsReadsOrPairs= "+ (pairedEnd?nrMappingsReadsOrPairs*2:nrMappingsReadsOrPairs)
										+ "\n\tlocus: "+ tt[0].getTranscriptID());
						};
					}
					if (outputExon&& outputSJunction) {
						try {assert(testInvariant(invariantObsAllEx, 
								pairedEnd?nrMappingsReadsOrPairs*2:nrMappingsReadsOrPairs, 0.05));}	// min: 0.02
						catch (AssertionError e) {
							Log.warn(getClass().getName()+".outputGFF():\n\tinvariantObsAllEx= "
										+ invariantObsAllEx+ ", nrMappingsReadsOrPairs= "+ (pairedEnd?nrMappingsReadsOrPairs*2:nrMappingsReadsOrPairs)
										+ "\n\tlocus: "+ tt[0].getTranscriptID());
						};
						try {assert(testInvariant(invariantPredAllEx, 
								pairedEnd?nrMappingsReadsOrPairs*2:nrMappingsReadsOrPairs, 0.1));}
						catch (AssertionError e) {
							Log.warn(getClass().getName()+".outputGFF():\n\tinvariantPredAllEx= "
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
							sb.append(" ");
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
					ArrayList<AbstractEdge> cc= new ArrayList<AbstractEdge>();
					if (solver!= null) {
						Iterator<Object> iter= solver.getConstraintHash().keySet().iterator();
						while (iter.hasNext()) {
							Object o= iter.next();
							if (o instanceof SimpleEdge)
								cc.add((SimpleEdge) o);
						}
					}
					Collections.sort(cc, SimpleEdge.getDefaultPositionComparator());
	
					Iterator<AbstractEdge> iter= cc.iterator();
					while (iter.hasNext()) {
						AbstractEdge e= iter.next();
						// no INTRONS
						if ((!(e instanceof SuperEdge))&& (!e.isExonic()))
							continue;
						//getGTF(sb, e, new long[][]{e.getTranscripts()}, g, solver, perM);
					}
				}
					
				Log.print(sb.toString());
				
			}
	
			private boolean outputGFForiginalLines(AnnotationMapper g,
					ASEvent[] events2, GraphLPsolver solver, HashMap<String, Double> rpkmMap) {
				
				Transcript[] tt= gene.getTranscripts();
				boolean foundTranscripts= false;
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
						
						Log.print(sb.toString());
						
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
						

						Log.print(sb.toString());
					} else if (outputUnknown) {
						Log.print(s+ System.getProperty("line.separator"));
					}
				}
				
				return foundTranscripts;
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
	
			private GraphLPsolver getSolver(AnnotationMapper mapper, int mappedReads) {
			
				GraphLPsolver solver= new GraphLPsolver(mapper, readLenMin, 
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
	
			private String getGTF(StringBuilder sb, ASEvent event, AnnotationMapper g, GraphLPsolver solver, boolean unsolvedSystem, 
						double perM, String pv, HashMap<Object,Double> tExpMap) {
					
			//		for (int i = 0; i < eeV.size(); i++) 
			//			eeV.elementAt(i).removeAllElements();
					Vector<Vector<AbstractEdge>> eeV= new Vector<Vector<AbstractEdge>>(5,5);
					while (eeV.size()< event.getDimension()) 
						eeV.add(new Vector<AbstractEdge>());
					
					g.getRPK(event, pairedEnd, SplicingGraph.ETYPE_AL, eeV);
					sb.append(event.toStringGTF());
			
					long[][] sig= new long[event.getDimension()][];
					for (int i = 0; i < sig.length; i++) { 
						sig[i]= g.createAllArray();		
						for (int j = 0; j < sig.length; j++) 
							sig[i]= j== i? sig[i]: 
								SplicingGraph.unite(sig[i], g.encodeTset(event.getTranscripts()[j]));
					}
					
					double splitReads= getGTFappend(sb, g, solver, eeV, perM, sig);
					++nrEvents;
					if (splitReads> 0)
						++nrEventsExp;
			
					for (int i = 0; i < eeV.size(); i++) 
						eeV.elementAt(i).removeAllElements();
					
					return sb.toString();
				}
	
			private String getGTF(StringBuilder sb, AbstractEdge e, long[][] sig, SplicingGraph g, GraphLPsolver solver, 
						double perM) {
					
					Vector<Vector<AbstractEdge>> eeV= new Vector<Vector<AbstractEdge>>(5,5);
					eeV.add(new Vector<AbstractEdge>(1));
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
	
			private String getGTF(StringBuilder sb, Exon exon, Transcript t, AnnotationMapper g, GraphLPsolver solver, boolean unsolvedSystem, 
						double perM, String pv, boolean attributesOnly) {
	
					if (!attributesOnly) {
						GFFObject obj= GFFObject.createGTFObjects(exon, t)[0];
						sb.append(obj.toString());
					}
				
			//		if (eeV.size()< 1)
			//			eeV.add(new Vector<Edge>());
			//		else
			//			eeV.elementAt(0).removeAllElements();
					Vector<Vector<AbstractEdge>> eeV= new Vector<Vector<AbstractEdge>>(5,5);
					eeV.add(new Vector<AbstractEdge>());
					
					//if (g.readCount> 0) // get lengths
					g.getRPK(exon, t, pairedEnd, SplicingGraph.ETYPE_AL, eeV);
			
					//containerIntA1A1[0][0]= g.readCount> 0? getLength(eeV.elementAt(0), null, true, false): (exon.getLength()- readLen);
					long[][] sig= new long[][]{g.encodeTset(t)};
					getGTFappend(sb, g, solver, eeV, perM, sig);
					eeV.elementAt(0).removeAllElements(); 
					
					return sb.toString();
				}
	
			private String getGTF(StringBuilder sb, Gene gene, AnnotationMapper g, GraphLPsolver solver, double perM, String pv) {
				
				//clearEdgeContainer(1);
				Vector<Vector<AbstractEdge>> eeV= new Vector<Vector<AbstractEdge>>(5,5);
				eeV.add(new Vector<AbstractEdge>());
				
				GFFObject obj= GFFObject.createGFFObject(gene);
				sb.append(obj.toString());
				//if (g.readCount> 0) // for getting lengths 
				g.getRPK(gene, pairedEnd, SplicingGraph.ETYPE_AL, eeV);
				
				
				//lenExon[0][0]= g.readCount> 0? getLength(eeV.elementAt(0), null): (t.getExonicLength()- readLen);
				//containerIntA1A1[0][0]= getLength(eeV.elementAt(0), null, true, false);
				//debug= true;
				long[][] sig= new long[][]{g.createAllArray()};
				getGTFappend(sb, g, solver, eeV, perM, sig);
				//debug= false;
				return sb.toString();
			}	
			
			private String getGTF(StringBuilder sb, Transcript t, GraphLPsolver solver, AnnotationMapper g, double perM, String pv, boolean attributesOnly) {
					
					GFFObject obj= GFFObject.createGFFObject(t);
					sb.append(obj.toString());
	
					Vector<Vector<AbstractEdge>> eeV= new Vector<Vector<AbstractEdge>>(5,5);
					eeV.add(new Vector<AbstractEdge>());
					//if (g.readCount> 0) // get lengths 
					g.getRPK(t, pairedEnd, SplicingGraph.ETYPE_AL, eeV);
			
					//lenExon[0][0]= g.readCount> 0? getLength(eeV.elementAt(0), null): (t.getExonicLength()- readLen);
					//containerIntA1A1[0][0]= getLength(eeV.elementAt(0), null, true, false);
					long[][] others= new long[1][];
					others[0]= SplicingGraph.without(g.createAllArray(), g.encodeTset(new Transcript[] {t}));
					
					long[][] sig= new long[][]{g.encodeTset(t)};	// containerLongA1A[0]
					getGTFappend(sb, g, solver, eeV, perM, sig);
					
					eeV.elementAt(0).removeAllElements();
					
					return sb.toString();
				}
	
			private double getGTFappend(StringBuilder sb, SplicingGraph g, GraphLPsolver solver, Vector<Vector<AbstractEdge>> eeV, double perM, long[][] tid) {
					
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
	
			/**
			 * extends a transcripts by the coordiantes covered by the BED objects
			 * @deprecated no longer used
			 * @param tx
			 * @param beds
			 * @param extension
			 * @return
			 */
			private int[] extend(Transcript tx, BEDobject2[] beds, int[] extension) {
				
				Attributes a= null;
				
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
					a= settings.get(FluxCapacitorSettings.READ_DESCRIPTOR).getAttributes(beds[i].getName(), a);	//Fasta.getReadID(beds[i].getName())
					mapMates5.put(a.id, beds[i]);
				}
				int left= i;
				// 3' extension
				for(i= beds.length- 1; i>= left; --i) {
					if (beds[i].getEnd()<= tend)
						break;
					if (beds[i].getStrand()>= 0|| beds[i].getBlockCount()> 1)	// only non-split
						continue; // only reads on reverse strand can have a mate falling into transcript
					a= settings.get(FluxCapacitorSettings.READ_DESCRIPTOR).getAttributes(beds[i].getName(), a);	//Fasta.getReadID(beds[i].getName())
					mapMates3.put(a.id, beds[i]);
				}
				int right= i;
				// find mates				
				int min= 0, max= 0;
				for(i= left; i<= right; ++i) {
					boolean contained= contains(tx, beds[i]);
					if (!contained)
						continue;
					a= settings.get(FluxCapacitorSettings.READ_DESCRIPTOR).getAttributes(beds[i].getName(), a);
					CharSequence id= a.id; 
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
			
			Coverage coverage= null;
			
			private void learn(Transcript tx, BufferedIterator beds) {
							
				if (beds== null)
					return;
				
				BEDobject2 bed1, bed2;
				UniversalReadDescriptor.Attributes 
					attributes= settings.get(FluxCapacitorSettings.READ_DESCRIPTOR).createAttributes(), 
					attributes2= settings.get(FluxCapacitorSettings.READ_DESCRIPTOR).createAttributes();
				int elen= tx.getExonicLength();	// this is the "effective" length, modify by extensions		
//				if (elen< readLenMin)
//					return;	// discards reads
				
				UniversalMatrix m= profile.getMatrix(elen);
				if (settings.get(FluxCapacitorSettings.COVERAGE_STATS)) {
					if (coverage== null)
						coverage= new Coverage(elen);
					else
						coverage.reset(elen);
				}
				
				while (beds.hasNext()) {
					
					++nrReadsSingleLoci;
					bed1= new BEDobject2(beds.next());
					CharSequence tag= bed1.getName();
					attributes= settings.get(FluxCapacitorSettings.READ_DESCRIPTOR).getAttributes(tag, attributes);	
					if (pairedEnd) {
						if (attributes.flag< 1)
							Log.warn("Read ignored, error in readID: "+ tag);
						if (attributes.flag== 2)	// don't iterate second read
							continue;
					}

					if (stranded) {
						if ((tx.getStrand()== bed1.getStrand()&& attributes.strand== 2)
								|| (tx.getStrand()!= bed1.getStrand()&& attributes.strand== 1))
						++nrMappingsWrongStrand;
						continue;
					}
					
					int bpoint1= getBpoint(tx, bed1);					
					if(bpoint1< 0|| bpoint1>= elen) {	// outside tx area, or intron (Int.MIN_VALUE)
						++nrReadsSingleLociNoAnnotation;
						continue;
					}

					++nrReadsSingleLociMapped;	// the (first) read maps
					
					if (pairedEnd) {

						beds.mark();
						while(beds.hasNext()) {
							bed2= new BEDobject2(beds.next());
							attributes2= settings.get(FluxCapacitorSettings.READ_DESCRIPTOR).getAttributes(bed2.getName(), attributes2);
							if (attributes2== null)
								continue;
							if (!attributes.id.equals(attributes2.id))
								break;						
							if (attributes2.flag== 1)	// not before break, inefficient
								continue;

							int bpoint2= getBpoint(tx, bed2);
							if (bpoint2< 0|| bpoint2>= elen) {
								++nrReadsSingleLociNoAnnotation;
								continue;
							}
								
							// check again strand in case one strand-info had been lost
							if (stranded) {
								if ((tx.getStrand()== bed2.getStrand()&& attributes2.strand== 2)
										|| (tx.getStrand()!= bed2.getStrand()&& attributes2.strand== 1))
								++nrMappingsWrongStrand;
								continue;
							}
							
							// check directionality (sequencing-by-synthesis)
							if (bed1.getStrand()== bed2.getStrand()
									|| (bed1.getStart()< bed2.getStart()&& bed1.getStrand()!= Transcript.STRAND_POS)
									|| (bed2.getStart()< bed1.getStart()&& bed2.getStrand()!= Transcript.STRAND_POS)) {
								nrPairsWrongOrientation+= 2;	
								continue;
							}

							m.add(bpoint1, bpoint2, -1, -1, elen);	// TODO rlen currently not used
							// update coverage
							if (settings.get(FluxCapacitorSettings.COVERAGE_STATS)) {
								if (bpoint1< bpoint2) {
									for (int i = bpoint1; i < bpoint1+ bed1.length(); i++) 
										coverage.increment(i);
									for (int i = bpoint2- bed2.length()+ 1; i <= bpoint2; i++) 
										coverage.increment(i);
								} else {
									for (int i = bpoint2; i < bpoint2+ bed2.length(); i++) 
										coverage.increment(i);
									for (int i = bpoint1- bed1.length()+ 1; i <= bpoint1; i++) 
										coverage.increment(i);
								}
							}
							//addInsertSize(Math.abs(bpoint2- bpoint1)+ 1);	// TODO write out insert size distribution
							
							nrReadsSingleLociPairsMapped+= 2;
							
						}
						beds.reset();
						
					} else {	// single reads						
						m.add(bpoint1, -1, elen, 
								bed1.getStrand()== tx.getStrand()?Constants.DIR_FORWARD:Constants.DIR_BACKWARD);
						// update coverage
						if (settings.get(FluxCapacitorSettings.COVERAGE_STATS)) {
							if (bed1.getStrand()== tx.getStrand()) {
								for (int i = bpoint1; i < bpoint1+ bed1.length(); i++) 
									coverage.increment(i);
							} else {
								for (int i = bpoint1- bed1.length()+ 1; i <= bpoint1; i++) 
									coverage.increment(i);
							}
						}
					}

				} // iterate bed objects
				

				// output coverage stats
				if (settings.get(FluxCapacitorSettings.COVERAGE_STATS)) {
					writeCoverageStats(
							tx.getGene().getGeneID(),
							tx.getTranscriptID(),
							tx.isCoding(),
							tx.getExonicLength(),
							pairedEnd? nrReadsSingleLociPairsMapped: nrReadsSingleLociMapped,
							coverage.getFractionCovered(),
							coverage.getChiSquare(true),
							coverage.getCV(true));
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

			int ok= loadLibraries();
			if (ok< 0) 
				exit(-1);
			

			final FluxCapacitor myCapacitor= new FluxCapacitor();
			myCapacitor.setFile(new File(args[0]));
			
		    // run
			try {
				myCapacitor.call();
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
	
	/**
	 * The parameter file.
	 */
	protected File file= null;
	
	public File fileMappedReads= null,
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
	
	private void printStats() {
		Log.info("HEHO", "We are set, so let's go!");
		// TODO
		// settings.write(Log.logStream);		
		StringBuilder sb;
		
		// INPUT
		Log.info(FluxCapacitorSettings.ANNOTATION_FILE.getName(),
				settings.get(FluxCapacitorSettings.ANNOTATION_FILE).getAbsolutePath());
		Log.info(FluxCapacitorSettings.MAPPING_FILE.getName(), 
				settings.get(FluxCapacitorSettings.MAPPING_FILE).getAbsolutePath());
		Log.info(FluxCapacitorSettings.READ_DESCRIPTOR.getName(), 
				settings.get(FluxCapacitorSettings.READ_DESCRIPTOR).toString());
		// TODO
		//p.println("\t"+CLI_LONG_VERBOSE+"\t"+Constants.VERBOSE_KEYWORDS[Constants.verboseLevel]);
		//if (copyLocal)
		//	p.println("\t"+FluxCapacitorParameters.PAR_COPY_INPUT);
		Log.info(settings.SORT_IN_RAM.getName(),
				Boolean.toString(settings.get(FluxCapacitorSettings.SORT_IN_RAM)));
		
		// OUTPUT
		Log.info(FluxCapacitorSettings.TMP_DIR.getName(),
				settings.get(FluxCapacitorSettings.TMP_DIR).getAbsolutePath());
		// TODO
		//if (Constants.globalPfx!= null)
			//p.println("\t"+FluxCapacitorConstants.CLI_LONG_TPX+"\t"+ Constants.globalPfx);
		sb= new StringBuilder();
		if (settings.get(FluxCapacitorSettings.STDOUT_FILE)== null)
			sb.append("stdout");
		else {
			sb.append(settings.get(FluxCapacitorSettings.STDOUT_FILE).getAbsolutePath());
			if (compressionOut!= FileHelper.COMPRESSION_NONE)
				sb.append("\t"+ FluxCapacitorConstants.CLI_LONG_COMPRESSION+ "\t"+ FileHelper.COMPRESSION_KEYWORDS[compressionOut]);
		}
		Log.info(settings.STDOUT_FILE.getName(), sb.toString());
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
			//p.println("\t"+ CLI_LONG_THREAD+" "+ maxThreads);
//		sb= new StringBuilder("Read Distribution\t");
//		if (uniform)
//			sb.append("uniform");
//		else if (fileProfile!= null&& fileProfile.exists())
//			sb.append("from profiles in "+ fileProfile.getAbsolutePath());
//		else {
//			sb.append("profiling is carried out");
//			if (fileProfile!= null)
//				sb.append(" and stored in "+ fileProfile.getAbsolutePath());
//		}
//		Log.info(sb.toString());
		
		if (stranded)
			Log.info("\tstrand information considered.");
		if (pairedEnd)
			Log.info("\tmate pairing information considered");
			

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
		if(settings.get(FluxCapacitorSettings.INSERT_FILE) != null){
            Log.info("\twriting insert sizes to "+settings.get(FluxCapacitorSettings.INSERT_FILE).getAbsolutePath());
        }
//		if (pairedEnd)
//			p.println("\t"+CLI_LONG_PAIR+"\t"+insertMinMax[0]+","+insertMinMax[1]);
		//System.err.println("\t"+CLI_LONG_NOISE+"\t"+Float.toString(1- GraphLPsolver.min_read_rest_frac));
		
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
	
	
	
	
	void fileFinish() {

		// TODO close input should occur by reader or interface method 
		boolean b= bedWrapper.close();
		b= gtfReader.close();

		if (settings.get(FluxCapacitorSettings.COVERAGE_STATS)) {
			if (FileHelper.move(
					fileTmpCovStats, 
					settings.get(FluxCapacitorSettings.COVERAGE_FILE))) {
				fileTmpCovStats.delete();
				Log.info("Coverage statistics in "+ settings.get(FluxCapacitorSettings.COVERAGE_FILE).getAbsolutePath());
			} else
				Log.warn("Failed to move coverage statistics to "+  
						settings.get(FluxCapacitorSettings.COVERAGE_FILE).getAbsolutePath()+ "\n"
						+ "\tinformation in "+ fileTmpCovStats.getAbsolutePath());
		}

		// TODO close files for non-/mapped reads, insert sizes, LPs, profiles  

		// close output
		if (Log.outputStream!= System.out&& Log.outputStream!= System.err) 
			Log.outputStream.close();
		
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
	
	static String createFileName(String base, byte compression) {
		if (compression== FileHelper.COMPRESSION_NONE) 
			return base;
		if (compression== FileHelper.COMPRESSION_ZIP)
			base+= '.'+ MyFile.SFX_ZIP;
		else if (compression== FileHelper.COMPRESSION_GZIP)
			base+= '.'+ MyFile.SFX_GZ;
		return base;
	}
	
	
	private void fileStats(AnnotationWrapper wrapper) {

		// (3) scan
		((AbstractFileIOWrapper) wrapper).scanFile();
		if(((AbstractFileIOWrapper) wrapper).getNrInvalidLines()> 0)
			Log.warn("Skipped "+ ((AbstractFileIOWrapper) wrapper).getNrInvalidLines()+ " lines.");

		Log.info(Constants.TAB+ wrapper.getNrGenes()+ " loci, "
				+ wrapper.getNrTranscripts()+ " transcripts, "
				+ wrapper.getNrExons()+ " exons.");
	}


	
	int nrBEDreads= -1, nrBEDmappings= -1;
	private int checkBEDscanMappings= 0;
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
	
	@Override
	public Void call() throws Exception {

		// TODO not here
		if (loadLibraries()< 0)
			System.exit(-1);
		
		// load parameters
        if (file == null || !file.exists()) {
            throw new RuntimeException("I have no parameter file and I want to scream!");
        }
        try {
            settings = FluxCapacitorSettings.createSettings(file);
        } catch (Exception e) {
            throw new RuntimeException("Unable to load settings from " + file + "\n\n " + e.getMessage(), e);
        }

        FileHelper.tempDirectory = settings.get(FluxCapacitorSettings.TMP_DIR);
		
		// prepare output files
		if (settings.get(FluxCapacitorSettings.STDOUT_FILE)!= null) {
			File f= settings.get(FluxCapacitorSettings.STDOUT_FILE);
			if (f.exists()&& !CommandLine.confirm(
                    "[CAUTION] I overwrite the output file " +
                            settings.get(FluxCapacitorSettings.STDOUT_FILE).getName() +
                            ", please confirm:\n\t(Yes,No,Don't know)")) {
				exit(-1);
			}
			
			try {
				Log.outputStream= new PrintStream(new FileOutputStream(f));
			} catch (FileNotFoundException e) {
				Log.warn("Cannot write log file to "+ f.getAbsolutePath());	// let it on stderr?!
			}
		}

		// prepare input files
		AbstractFileIOWrapper wrapperAnnotation;
		AbstractFileIOWrapper wrapperMappings;
		if (cheatDisableFCheck) {
			Log.warn("Development run, file check disabled !!!");
			wrapperAnnotation= getWrapper(settings.get(FluxCapacitorSettings.ANNOTATION_FILE));
			wrapperMappings= getWrapper(settings.get(FluxCapacitorSettings.MAPPING_FILE));
		} else {
			wrapperAnnotation=
				fileInit(settings.get(FluxCapacitorSettings.ANNOTATION_FILE));
			fileStats((AnnotationWrapper) wrapperAnnotation);
			
			wrapperMappings= 
				fileInit(settings.get(FluxCapacitorSettings.MAPPING_FILE));
			fileStats((MappingWrapper) wrapperMappings);

		}
		

		// TODO parameters
		pairedEnd= settings.get(FluxCapacitorSettings.ANNOTATION_MAPPING).equals(AnnotationMapping.PAIRED)
				|| settings.get(FluxCapacitorSettings.ANNOTATION_MAPPING).equals(AnnotationMapping.COMBINED);
		stranded= settings.get(FluxCapacitorSettings.ANNOTATION_MAPPING).equals(AnnotationMapping.STRANDED)
				|| settings.get(FluxCapacitorSettings.ANNOTATION_MAPPING).equals(AnnotationMapping.COMBINED);
		
		printStats();
		
		// run
		long t0= System.currentTimeMillis();

		profile= getProfile();
		if (profile== null) {
			exit(-1);
		}
		
		explore(FluxCapacitorConstants.MODE_RECONSTRUCT);

		
		fileFinish();
		
		Log.info("\n[TICTAC] I finished flux in "
				+((System.currentTimeMillis()- t0)/ 1000)+" sec.\nCheers!");
		
		//System.err.println("over "+ GraphLPsolver.nrOverPredicted+", under "+GraphLPsolver.nrUnderPredicted);
		
		return null;
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
	/**
	 * @deprecated 
	 */
	void appendFreq() {

		File fileOut= settings.get(FluxCapacitorSettings.STDOUT_FILE);
		if (fileOut== null)
			return;
		
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
			File fileTmp= createTempFile(null,
					fileOut.getName()+"__append", 
					MyFile.getExtension(fileOut.getName()),
					true); 
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
		float rpkm= (float) ((reads/ (double) len)* (1000000000l/ (double) (nrBEDreads< 0? 1: nrBEDreads)));
		return rpkm;
	}
	
	public String getCompositeFName() {
		File f= settings.get(FluxCapacitorSettings.ANNOTATION_FILE), 
			g= settings.get(FluxCapacitorSettings.MAPPING_FILE);
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
	
	/**
	 * Creates a temporary file in the location provided, iff write access is 
	 * available there. Otherwise the file is created in the custom or system
	 * temporary directory. 
	 * 
	 * @param location a file in the target directory or the directory itself,
	 * may be <code>null</code>
	 * @param name prefix of the file to be created, class name is appended
	 * at the beginning
	 * @param extension (optional) suffix of the temporary file that is created
	 * @param deleteOnExit flag for calling the <code>deleteOnExit()</code> 
	 * method for the file
	 * @return a temporary file according to the specifications
	 */
	protected File createTempFile(File location, String name, String extension, boolean deleteOnExit) {
		
		// get location
		if (location== null)
			location= settings.get(FluxCapacitorSettings.TMP_DIR);
		else {
			if (!location.isDirectory())
				location= location.getParentFile();
			if (!location.canWrite())
				location= settings.get(FluxCapacitorSettings.TMP_DIR);
		}

		// get name
		if (name== null)
			name= getClass().getSimpleName();
		else
			name= getClass().getSimpleName()+ "_"+ name;
		
		File f= null;
		try {
			f= FileHelper.createTempFile(name, extension, location);
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
		
		return createFile(f, deleteOnExit);
	}
	
	/**
	 * Control gateway for file creation from the main class, 
	 * adds a hook for delete on exit in case.
	 * 
	 * @param f the file that has been created
	 * @param deleteOnExit flag to mark for deletion on exit
	 * @return
	 */
	protected File createFile(File f, boolean deleteOnExit) {
		if (deleteOnExit)
			f.deleteOnExit();
		
		return f;
	}
	
	public File getFileProfile() {
		if (fileProfile == null) {
			String s= getNameProfile();
			if (s== null)
				return null;
			
			
			fileProfile=  createTempFile(
					settings.get(FluxCapacitorSettings.MAPPING_FILE), 
					FileHelper.stripExtension(settings.get(FluxCapacitorSettings.MAPPING_FILE).getName()), 
					"prf", 
					false);
				//createTempFile(s, FileHelper.getCompressionExtension(FileHelper.COMPRESSION_ZIP));
				
//				new File(System.getProperty(Constants.PROPERTY_TMPDIR)+ File.separator 
//				+ fName+ SFX_PROFILES
//				+ "."+ FileHelper.getCompressionExtension(FileHelper.COMPRESSION_ZIP));
		}

		return fileProfile;
	}
	
	private String getNameISize() {
		String s= getCompositeFName();
		if (s== null)
			return null;
		return s+ FluxCapacitorConstants.SFX_INSERTSIZE+ Constants.DOT+ "txt";
	}
		
	private String getNameLP() {
		return getCompositeFName()+ FluxCapacitorConstants.SFX_LP;
	}
	
	public File getFileLP() {
		if (fileLPdir == null) {
			fileLPdir= createTempFile(null, getNameLP(), null, false);
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
	
	File fileTmpCovStats= null;
	private BufferedWriter writerTmpCovStats= null;
	
	/**
	 * Writes coverage statistics of a transcript to disk
	 * @param geneID locus identifier
	 * @param transcriptID transcript identifier
	 * @param cds flag to indicate whether transcript has an annotated ORF
	 * @param length (processed) length of the transcript 
	 * @param i
	 * @param fractionCovered
	 * @param chiSquare
	 * @param cv
	 */
	private void writeCoverageStats(String geneID, String transcriptID,
			boolean cds, int length, int nrReads,
			float fracCov, long chiSquare, double cv) {
		
		try {
			if (fileTmpCovStats== null) 
				fileTmpCovStats= FileHelper.createTempFile("tmpCovStats", ".pro");
			if (writerTmpCovStats== null)
				writerTmpCovStats= new BufferedWriter(new FileWriter(fileTmpCovStats));
			writerTmpCovStats.write(
				geneID+ "\t"+ transcriptID+ "\t"+ (cds? "CDS": "NC")+ "\t"
				+ Integer.toString(length)+ "\t"+ Integer.toString(nrReads)+ "\t"
				+ Float.toString(fracCov)+ "\t"+ Long.toString(chiSquare)+ "\t"
				+ Float.toString((float) cv)+ "\n"
			);
			
		} catch (Exception e) {
			throw new RuntimeException(e);
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
	
    /**
     * Set the parameter file
     *
     * @param file parameter file
     */
    @Option(name = "p", longName = "parameter", description = "specify parameter file (PAR file)", displayName = "file", required = true)
    public void setFile(File file) {
        this.file = file;
    }
    
    public boolean validateParameters(HelpPrinter printer, ArgumentProcessor toolArguments) {

        if(isPrintParameters()){
            FluxCapacitorSettings settings = new FluxCapacitorSettings();
            settings.write(System.out);
            return false;
        }

        if (getFile() == null) {
            Log.error("");
            Log.error("No parameter file specified !");
            Log.error("\n");
            printer.print(toolArguments);
            return false;
        }
        if (!getFile().canRead()) {
            Log.error("");
            Log.error("Parameter file " + getFile().getAbsolutePath() + " does not exist or I can not read it!");
            Log.error("\n");
            printer.print(toolArguments);
            return false;
        }

        return true;
    }

	
	private boolean isPrintParameters() {
		return printParameters;
	}

	public static byte mapFileType= FluxCapacitorConstants.FORMAT_SAM;
	private void writeMapFileSam(SplicingGraph g, SimpleEdge e, DirectedRegion[] regs, DirectedRegion[][] contRegs) {
		
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
	private void solve(Gene gene, BufferedIterator beds, boolean decompose) {
		
		// create LP and solve
		LocusSolver lsolver= new LocusSolver(gene, beds, decompose); 
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

	
	byte costModel= GraphLPsolver.COSTS_LINEAR;
	byte costSplit= 1;
	
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

	
	boolean uniform= false;
	long nrReadsAll= 0;
	double costModelPar= Double.NaN;
	float[] costBounds= new float[] {0.95f, Float.NaN};	// how much of the original observation can be subs/add
	int[] profileBoundaries;
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
				//writeProfiles(); // TODO
			}
		}
		
		if (profile== null)
			return null;
		// check
		for (int i = 0; i < profile.getMasters().length; i++) {
			if (profile.getMasters()[i].hasEmptyPositions())
				profile.getMasters()[i].fill();
		}
		

		return profile;
	}
	
	BinVector isizeV; 
	synchronized void addInsertSize(int isize) {
		isizeV.incrTuple(isize);
	}
	
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
		

	Vector<SimpleEdge> edgeColl1= new Vector<SimpleEdge>(), edgeColl2= new Vector<SimpleEdge>();
	int[][] containerIntA1A1= new int[1][];
	{ containerIntA1A1[0]= new int[1]; }
	long[][] containerLongA1A= new long[1][];
	boolean keepTmpSorted= false;
	
	public double getLength(SplicingGraph g, Vector<AbstractEdge> v, long[] sig, boolean exclusive) {
		double len= 0; 
		for (int i = 0; i < v.size(); i++) {
			AbstractEdge e= v.elementAt(i);
			long[] trpts= e.getTranscripts();
			long[] inter= SplicingGraph.intersect(trpts, sig);
			if (SplicingGraph.isNull(inter)|| (exclusive&& !SplicingGraph.equalSet(sig, trpts)))
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
	
	public double getReads(Vector<AbstractEdge> v, byte dir, long[] sig, boolean normalized) {
		int sum= 0;
		for (int i = 0; i < v.size(); i++) {
			AbstractEdge e= v.elementAt(i);
			long[] inter= SplicingGraph.intersect(e.getTranscripts(), sig);
			if (SplicingGraph.isNull(inter)|| !e.isExonic())
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
						sum+= cnt* ((MappingsInterface) se).getMappings().getReadNr();
					if (dir<= 0)
						sum+= cnt* ((MappingsInterface) se).getMappings().getRevReadNr();
				}
			} else {
				if (dir>= 0)
					sum+= ((MappingsInterface) e).getMappings().getReadNr();
				if (dir<= 0)
					sum+= ((MappingsInterface) e).getMappings().getRevReadNr();
			}
		}
		return sum;
	}
	
	public double getReadsAvg(Vector<AbstractEdge> v, byte dir, SplicingGraph g, long[] sig, boolean excl, boolean normalized) {
		double sum= 0;
		for (int i = 0; i < v.size(); i++) {
			AbstractEdge e= v.elementAt(i);
			long[] trpts= v.elementAt(i).getTranscripts();
			long[] inter= SplicingGraph.intersect(trpts, sig);
			if (SplicingGraph.isNull(inter)|| (excl&& !SplicingGraph.equalSet(sig, trpts))|| !e.isExonic())
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
						sum+= (((MappingsInterface) se).getMappings().getReadNr()* mult* cnt)/ sf;
					if (dir<= 0)
						sum+= (((MappingsInterface) se).getMappings().getRevReadNr()* mult* cnt)/ sf;
				}
			} else {
				if (dir>= 0)
					sum+= (((MappingsInterface) e).getMappings().getReadNr()* mult)/ sf;
				if (dir<= 0)
					sum+= (((MappingsInterface) e).getMappings().getRevReadNr()* mult)/ sf;
			}
			
			System.currentTimeMillis();
		}
		
		return sum;
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
		
		// security check, get distance between both exonic coordinates
		int epos2= tx.getExonicPosition(bed.getStrand()>= 0? bed.getEnd(): bed.getStart()+ 1);
		int len= bed.getLength();
		if (readLenMin< 0|| len< readLenMin)
			readLenMin= len;
		if (len> readLenMax) 
			readLenMax= len;

		if (len!= Math.abs(epos- epos2)+ 1)
			return Integer.MIN_VALUE;
		return epos;
	}
	
	private GTFwrapper gtfReader;
	private AbstractFileIOWrapper getWrapperGTF(File inputFile) {
		
		gtfReader= new GTFwrapper(inputFile.getAbsolutePath());
		gtfReader.setNoIDs(null);
		gtfReader.setReadGene(true);
		gtfReader.setReadFeatures(new String[] {"exon","CDS"});
		gtfReader.setReadAheadTranscripts(1);	// only one locus a time
//		gtfReader.setReadAheadTranscripts(-1);
//		gtfReader.setReadAll(true);
		gtfReader.setGeneWise(true);
		gtfReader.setPrintStatistics(false);
		gtfReader.setReuse(true);
		Transcript.removeGaps= false;
		
		return gtfReader;
	}


	/**
	 * @deprecated
	 * @return
	 */
	public GTFwrapper getGTFreader() {
		if (gtfReader == null) {
			gtfReader= new GTFwrapper(settings.get(FluxCapacitorSettings.ANNOTATION_FILE).getAbsolutePath());
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
	private AbstractFileIOWrapper getWrapperBED(File inputFile) {
		bedWrapper= new BEDwrapper(inputFile.getAbsolutePath());

		return bedWrapper;
	}


	/**
	 * @deprecated
	 * @return
	 */
	public BEDwrapper getBedReader() {
		if (bedWrapper == null) {
			bedWrapper= new BEDwrapper(settings.get(FluxCapacitorSettings.MAPPING_FILE).getAbsolutePath());
		}

		return bedWrapper;
	}
	
	private File getFileProfiles(int binIdx) {
		return new File(settings.get(FluxCapacitorSettings.STDOUT_FILE).getAbsolutePath()+"_bin_"+profileBoundaries[binIdx]);
	}
	
	

	int[][] profileStub, profileStubRev;
	/**
	 * profile managing the matrices
	 */
	Profile profile;
	private BufferedIterator readBedFile(Gene gene, int from, int to, byte mode) {
        return readBedFile(gene, from, to, mode, 0, 1);
    }
	private BufferedIterator readBedFile(Gene gene, int from, int to, byte mode, int retryCount, long timeInSeconds) {
		
		if (from> to|| from< 0|| to< 0) 
			throw new RuntimeException("BED reading range error: "+from+" -> "+to);

		// init iterator
		BufferedIterator iter= null;

		try {
			if (settings.get(FluxCapacitorSettings.SORT_IN_RAM)) {
				// memory
				MappingWrapperState state= bedWrapper.read(gene.getChromosome(), from, to);
				if (state.result== null)
					return null;
				BEDobject2[] beds= (BEDobject2[]) state.result;
				Arrays.sort(beds, getDescriptorComparator());
				iter= new BufferedIteratorRAM(beds);
				
			} else {
				
				// read, maintain main thread			
				PipedInputStream  pin= new PipedInputStream();
		        PipedOutputStream pout= new PipedOutputStream(pin);
				Comparator<CharSequence> c= new BEDDescriptorComparator(settings.get(FluxCapacitorSettings.READ_DESCRIPTOR));
				File tmpFile= createTempFile(null, gene.getChromosome()+ ":"+ from+ "-"+ to+ ".", "bed", true);
				BufferedIteratorDisk biter= new BufferedIteratorDisk(pin, tmpFile, c);
				biter.init();
				iter= biter;
				MappingWrapperState state= bedWrapper.read(pout, gene.getChromosome(), from, to);
				pout.flush();
				pout.close();
				if (state.count== 0)
					return null;
			}			
		} catch (IOException e) {
            /**
             * "Resource temporarily unavailable"
             * Catch this exception and try again after sleeping for a while
             */
            if(e.getMessage().contains("Resource temporarily unavailable")){
                if(retryCount < 6){
                    Log.warn("Filesystem reports : 'Resource temporarily unavailable', I am retrying ("+(retryCount+1)+")");
                    try {Thread.sleep(1000 * (timeInSeconds));} catch (InterruptedException e1) {}
                    return readBedFile(gene, from, to, mode, retryCount + 1, timeInSeconds*6);
                }
            }
			throw new RuntimeException(
				"Could not get reads for locus "+ gene.getChromosome()+ ":"+ from+ "-"+ to +", retried " + retryCount + " times", e);
		}
        
		return iter;
		
	}
	
	
	BEDDescriptorComparator comp= null;
	private Comparator<? super BEDobject2> getDescriptorComparator() {
		if (comp == null) {
			comp = new BEDDescriptorComparator(
					settings.get(FluxCapacitorSettings.READ_DESCRIPTOR));
		}

		return comp;
	}

	/**
	 * @deprecated
	 * @param gene
	 * @param handler
	 * @param ostream
	 * @return
	 */
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
            LPSolverLoader.load();
			lpVer= LpSolve.lpSolveVersion();
			if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) {
				Log.info("PRE-CHECK", "\t* successfully loaded lpsolve JNI (version "+lpVer.getMajorversion()+"."+lpVer.getMinorversion()
						+",release "+lpVer.getRelease()+",build "+lpVer.getBuild()+(miss?";":"")+")\n");
			}
			return 0;
		} catch (Exception e) {
            Log.error("Error while loading native libraries: " + e.getMessage());
            Log.error("You can try to set the environment variables " + LPSolverLoader.ENV_JNI + " for \n" +
                    "the JNI library and " + LPSolverLoader.ENV_LIB + " for the shared object library to support \n" +
                    "your operating system. The sources for the LPSolver can be found @ http://sourceforge.net/projects/lpsolve");
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

	FluxCapacitorSettings settings= null;
	protected boolean printParameters;
	

	public boolean explore(byte mode) {
	
			nrSingleTranscriptLoci= 0;
			nrReadsLoci= 0;
			nrReadsMapped= 0; 
			nrReadsWrongLength= 0;
			nrMappingsWrongStrand= 0;
			
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
				
				gtfReader.read();
				Gene[] gene= null, geneNext= gtfReader.getGenes();
				
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
						origLines= (Vector<String>) gtfReader.getVLines().clone();	// TODO make array, trim..
					
					// http://forums.sun.com/thread.jspa?threadID=5171135&tstart=1095
					if (readerThread== null)
						readerThread= new GTFreaderThread();
					//readerThread.start();
					readerThread.run();
					geneNext= gtfReader.getGenes();
	
					for (int i = 0; (gene!= null)&& i < gene.length; i++) {
						
						
						// flop strand
						if (lastChr.equals(gene[i].getChromosome())) {
							if (lastStr!= gene[i].getStrand()) {
								//System.err.println(lastChr+" "+lastStr+ " "+ readObjects+ " wrote "+ dbgCntWriteMap +" not "+ dbgCntWriteNonmap);
								readObjects= 0;	
								// jump back
								bedWrapper.reset(gene[i].getChromosome());
								lastStr= gene[i].getStrand();
								lastEnd= -1;
							}
						} else {						// flop chr
							//System.err.println(lastChr+" "+lastStr+ " "+ readObjects+ " wrote "+ dbgCntWriteMap +" not "+ dbgCntWriteNonmap);
							readObjects= 0;
							lastChr= gene[i].getChromosome();
							lastStr= gene[i].getStrand();
							lastEnd= -1;
						}
					
						if (gene[i].getTranscriptCount()== 1)
							++nrSingleTranscriptLoci;
						else if (mode== FluxCapacitorConstants.MODE_LEARN)
							continue;	// performance for not reading beds
						
						BufferedIterator beds= null;
	
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
						
						beds= readBedFile(gene[i], start, end, mode);
						
						if (mode== FluxCapacitorConstants.MODE_LEARN&& beds!= null) {
							solve(gene[i], beds, false);
						} else if (mode== FluxCapacitorConstants.MODE_RECONSTRUCT) {
							solve(gene[i], beds, true); 
						}
						
						if (beds!= null)
							beds.clear();
							
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
				
				bedWrapper.finish();
				
				while (threadPool.size()> 0&& threadPool.elementAt(0).isAlive())
					try {
						threadPool.elementAt(0).join();
					} catch (Exception e) {
						; //:)
					}
	            Log.progressFinish(StringUtils.OK, true);
				if (checkGTFscanExons> 0&& checkGTFscanExons!= gtfReader.getNrExons())
					System.err.println("[ERROR] consistency check failed in GTF reader: "+ checkGTFscanExons+ "<>"+ gtfReader.getNrExons());
				checkGTFscanExons= gtfReader.getNrExons(); 
				if (checkBEDscanMappings> 0&& checkBEDscanMappings!= bedWrapper.getNrLines())
					System.err.println("[ERROR] consistency check failed in BED reader "+ checkBEDscanMappings+ "<>"+ bedWrapper.getNrLines());
				//checkBEDscanMappings= getBedReader().getNrLines();
				if (mode== FluxCapacitorConstants.MODE_LEARN) {
	
					if (pairedEnd&& func.getTProfiles()!= null) {
						insertMinMax= getInsertMinMax();
					}
					// close coverage writer
					if (settings.get(FluxCapacitorSettings.COVERAGE_STATS))
						try {
							writerTmpCovStats.close();
							writerTmpCovStats= null;
						} catch (Exception e) {
							throw new RuntimeException(e);
						}
						
	
					if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) {
						
						//System.err.println(" OK.");
						System.err.println("\tfirst round finished .. took "+((System.currentTimeMillis()- t0)/ 1000)+ " sec.\n\n\t"
								+ nrSingleTranscriptLoci+" single transcript loci\n\t"							
								+ bedWrapper.getNrLines()+ " mappings in file\n\t"
								+ nrReadsSingleLoci+" mappings fall in single transcript loci\n\t"	// these loci(+/-"+tolerance+"nt)\n\t"
								// counter un-reliable, /2 read is skipped in paired-end mode
								// + nrReadsSingleLociMapped+" mappings map to annotation\n\t"
								+ ((strand== FluxCapacitorConstants.STRAND_SPECIFIC)?nrMappingsWrongStrand+" mappings map to annotation in antisense direction,\n\t":"")
								//+ (pairedEnd?(nrReadsSingleLociPotentialPairs+ " mappings form potential pairs,\n\t"):"")
								+ (pairedEnd?(nrReadsSingleLociPairsMapped)+" mappings in annotation-mapped pairs\n\t":"")
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
	
	//					if (fileMappings!= null)
	//						getSammy().close();
					
					//assert(nrUniqueReads==getBedReader().getNrUniqueLinesRead());	// take out for cheat
					if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) {
						System.err.println();
						System.err.println("\treconstruction finished .. took "+((System.currentTimeMillis()- t0)/ 1000)+ " sec.\n\n\t"
								+ bedWrapper.getNrLines()+" mappings read from file\n\t"
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
				Log.error("Error while iterating loci:", e1);
				throw new RuntimeException(e1);
			}
	
			return true;
		}

	/**
	 * Checks whether file is to be uncompressed and/or sorted.
	 * 
	 * @param inputFile
	 * @param wrapper
	 * @return
	 */
	public AbstractFileIOWrapper fileInit(File inputFile) {

		// (1) unpack, if compressed
		byte cb= FileHelper.getCompression(inputFile);
		if (cb!= FileHelper.COMPRESSION_NONE) {
			File f= new File(FileHelper.stripExtension(inputFile.getAbsolutePath()));
			if (f.exists()) {
				Log.println("Assuming file "+ f.getName()+" is a decompressed version of "+ inputFile.getName());
			} else {
				f= createTempFile(null, FileHelper.getFileNameWithoutExtension(f), FileHelper.getExtension(f), true);
				try {
					FileHelper.inflate(inputFile, f, cb);
				} catch (Exception e) {
					throw new RuntimeException(e);
				}
			}
			inputFile= f;
			inputFile.deleteOnExit();	// carefully
		}
		
		// (2) sort, if needed
		AbstractFileIOWrapper wrapper= getWrapper(inputFile);
		if (!wrapper.isApplicable()) {
			File f= FileHelper.getSortedFile(inputFile);
			File lock= FileHelper.getLockFile(f);

			if (f.exists()&& !lock.exists()) {
				
				Log.warn("Assuming file "+ f.getName()+" is a sorted version of "+ inputFile.getName());
				
			} else {	// we have to sort
				
				boolean lockCreated= false;
				if (settings.get(FluxCapacitorSettings.KEEP_SORTED_FILES)) {	// try to store in original
					
					if (lock.exists()) {	// switch to sorting to temp
						Log.warn("Seems that another process is just sorting file "+ inputFile+
								"\nremove lock file "+ lock.getName()+" if dead leftover."+
								"\nContinuing with sorting to temporary file "+ 
								(f= createTempFile(f, 	// access to non-Temp
										FileHelper.getFileNameWithoutExtension(f), 
										FileHelper.getExtension(f),
										false)).getAbsolutePath());
						
					} else if (!f.getParentFile().canWrite()) {	// sort to temp, but do not delete (parameter)
						Log.warn("Cannot write sorted file to "+ f.getAbsolutePath()+
								"\nContinuing with sorting to temporary file "+ 
								(f= createTempFile(f, // access to non-Temp
										FileHelper.getFileNameWithoutExtension(f), 
										FileHelper.getExtension(f),
										false)).getAbsolutePath());
						
					} else {	// sort to default sorted file
						try {
							lock.createNewFile();
						} catch (Exception e) {
							throw new RuntimeException(e);
						}
						lockCreated= true;
					}
					
				} else {	// do not keep sorted files, sort to temp and delete on exit
					f= createTempFile(null, 
							FileHelper.getFileNameWithoutExtension(f), 
							FileHelper.getExtension(f),
							true);
				}

				// doit
				wrapper.sort(f);
				
				// if locked
				if (lockCreated)
					lock.delete();

				// if unzipped before
				if (cb!= FileHelper.COMPRESSION_NONE)
					inputFile.delete();	// carefully
				
				inputFile= f;
			}
			wrapper= getWrapper(inputFile);
		}
		
		return wrapper;
	}
	
	private void fileStats(MappingWrapper wrapper) {

		if (settings.get(FluxCapacitorSettings.NR_READS_MAPPED)<= 0) {
			// (3) scan
			((AbstractFileIOWrapper) wrapper).scanFile();
			if(((AbstractFileIOWrapper) wrapper).getNrInvalidLines()> 0)
				Log.warn("Skipped "+ ((AbstractFileIOWrapper) wrapper).getNrInvalidLines()+ " lines.");

			checkBEDscanMappings= wrapper.getCountMappings();
			nrBEDreads= wrapper.getCountReads();
			nrBEDmappings= wrapper.getCountMappings();
		} else {
			checkBEDscanMappings= -1;
			nrBEDreads= settings.get(FluxCapacitorSettings.NR_READS_MAPPED);
			nrBEDmappings= -1;
		}
		
		Log.info("\t"+ nrBEDreads+ " reads"
				+ (nrBEDmappings> 0? nrBEDmappings+ " mappings: R-factor "+ (wrapper.getCountMappings()/ (float) wrapper.getCountReads()): ""));
		if (nrBEDmappings> 0)
			Log.info("\t" + wrapper.getCountContinuousMappings() + " entire, " + wrapper.getCountSplitMappings()
	                + " split mappings (" + (wrapper.getCountSplitMappings() * 10f / wrapper.getCountMappings()) + "%)");
		
		// (4) check if read descriptor is applicable
		if (wrapper.isApplicable(settings.get(FluxCapacitorSettings.READ_DESCRIPTOR)))
			Log.info("\tRead descriptor seems OK");
		else {
			String msg= "Read Descriptor "+ settings.get(FluxCapacitorSettings.READ_DESCRIPTOR)
					+ " incompatible with read IDs";
			Log.error(msg);
			throw new RuntimeException(msg);
		}
	}

	private AbstractFileIOWrapper getWrapper(File inputFile) {
		
		String ext= FileHelper.getExtension(inputFile).toUpperCase();
		SupportedFormatExtensions sup= null;
		try {
			sup= SupportedFormatExtensions.valueOf(ext);
		} catch (Exception e) {
			throw new RuntimeException("Unsupported file format "+ ext);
		}
		
		
		switch (sup) {
		case GTF:
		case GFF:
			return getWrapperGTF(inputFile);

		case BED:
			return getWrapperBED(inputFile);
			
		}
		
		return null;	// make compiler happy
	}

	public File getFile() {
		return file;
	}

    /**
     * Enable default parameter printing
     *
     * @param printParameters enable disable
     */
    @Option(name = "o", longName = "printParameters", description = "Print default parameters", required = false)
    public void setPrintParameters(final boolean printParameters) {
        this.printParameters = printParameters;
    }

    
}
 