package fbi.genome.sequencing.rnaseq.reconstructionew;


//import de.luschny.math.factorial.FactorialPrimeSchoenhage;
//import de.luschny.math.arithmetic.Xint;

//import io.Sammy;

import fbi.commons.Log;
import fbi.commons.StringUtils;
import fbi.commons.file.FileHelper;
import fbi.commons.system.SystemInspector;
import fbi.commons.thread.SyncIOHandler2;
import fbi.commons.thread.ThreadedQWriter;
import fbi.genome.io.bed.BEDwrapper;
import fbi.genome.io.gff.GFFReader;
import fbi.genome.io.rna.BARNAdescriptor;
import fbi.genome.io.rna.Descriptor;
import fbi.genome.io.rna.RegExpDescriptor;
import fbi.genome.model.*;
import fbi.genome.model.bed.BEDobject;
import fbi.genome.model.bed.BEDobject2;
import fbi.genome.model.commons.Distribution;
import fbi.genome.model.commons.DoubleVector;
import fbi.genome.model.commons.IntVector;
import fbi.genome.model.commons.MyFile;
import fbi.genome.model.constants.Constants;
import fbi.genome.model.gff.GFFObject;
import fbi.genome.model.splicegraph.Edge;
import fbi.genome.model.splicegraph.Graph;
import fbi.genome.model.splicegraph.Node;
import fbi.genome.model.splicegraph.SuperEdge;
import fbi.genome.sequencing.rnaseq.reconstruction.*;
import fbi.genome.sequencing.rnaseq.reconstruction.gui.FluxCapacitorGUI;
import lpsolve.LpSolve;
import lpsolve.LpSolveException;
import lpsolve.VersionInfo;
import org.apache.commons.cli.Options;

import java.io.*;
import java.lang.reflect.Method;
import java.text.DecimalFormat;
import java.text.SimpleDateFormat;
import java.util.*;
import java.util.concurrent.atomic.AtomicLong;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;
import java.util.zip.ZipOutputStream;

//import gphase.solexa.lp.GraphLPsolverRelease;
//import gphase.solexa.lp.TProfile;
//import gphase.solexa.lp.TProfileFunction;


public class FluxCapacitorNew extends FluxCapacitor {

		/**
		 * +2 offset to the actual return values
		NOMEMORY (-2)  	Out of memory
		OPTIMAL (0) 	An optimal solution was obtained
		SUBOPTIMAL (1) 	The model is sub-optimal. Only happens if there are integer variables and there is already an integer solution found. The solution is not guaranteed the most optimal one.
	
		    * A timeout occured (set via set_timeout or with the -timeout option in lp_solve)
		    * set_break_at_first was called so that the first found integer solution is found (-f option in lp_solve)
		    * set_break_at_value was called so that when integer solution is found that is better than the specified value that it stops (-o option in lp_solve)
		    * set_mip_gap was called (-g/-ga/-gr options in lp_solve) to specify a MIP gap
		    * An abort function is installed (put_abortfunc) and this function returned TRUE
		    * At some point not enough memory could not be allocated 
	
		INFEASIBLE (2) 	The model is infeasible
		UNBOUNDED (3) 	The model is unbounded
		DEGENERATE (4) 	The model is degenerative
		NUMFAILURE (5) 	Numerical failure encountered
		USERABORT (6) 	The abort routine returned TRUE. See put_abortfunc
		TIMEOUT (7) 	A timeout occurred. A timeout was set via set_timeout
		PRESOLVED (9) 	The model could be solved by presolve. This can only happen if presolve is active via set_presolve
		PROCFAIL (10) 	The B&B routine failed
		PROCBREAK (11) 	The B&B was stopped because of a break-at-first (see set_break_at_first) or a break-at-value (see set_break_at_value)
		FEASFOUND (12) 	A feasible B&B solution was found
		NOFEASFOUND (13) 	No feasible B&B solution found
	 */
	public static final String[] RETURN_VERBOSE= new String[] {
		"NOMEMORY (-2)\tOut of memory",
		null,	// -1
		"OPTIMAL (0)\tAn optimal solution was obtained",
		"SUBOPTIMAL (1)\tThe model is sub-optimal.\n" +
		"\tOnly happens if there are integer variables\n" +
		"\tand there is already an integer solution found.\n" +
		"\tThe solution is not guaranteed the most optimal one.",
		"INFEASIBLE (2)\tThe model is infeasible",
		"UNBOUNDED (3)\tThe model is unbounded",
		"DEGENERATE (4)\tThe model is degenerative",
		"NUMFAILURE (5)\tNumerical failure encountered",
		"USERABORT (6)\tThe abort routine returned TRUE. See put_abortfunc",
		"TIMEOUT (7)\tA timeout occurred. A timeout was set via set_timeout",
		null,	// 8
		"PRESOLVED (9)\tThe model could be solved by presolve.\n" +
		"\tThis can only happen if presolve is active via set_presolve",
		"PROCFAIL (10)\tThe B&B routine failed",
		"PROCBREAK (11)\tThe B&B was stopped because of a break-at-first\n" +
		"\t(see set_break_at_first) or a break-at-value (see set_break_at_value)",
		"FEASFOUND (12)\tA feasible B&B solution was found",
		"NOFEASFOUND (13)\tNo feasible B&B solution found"
	};

	private static final int MAX_MATES= 100;
	
	public static byte SHELL_NONE= 0, SHELL_BASH= 1, SHELL_CSH= 2, SHELL_KSH= 3;
	
	public static final String GTF_ATTRIBUTE_LENGTH= "slots",
	GTF_ATTRIBUTE_TOKEN_OBSV= "obsv",
	GTF_ATTRIBUTE_TOKEN_PRED= "pred",
	GTF_ATTRIBUTE_TOKEN_BALANCED= "bal",
	GTF_ATTRIBUTE_TOKEN_ALL= "all",
	GTF_ATTRIBUTE_TOKEN_TID= "split",
	GTF_ATTRIBUTE_TOKEN_EXC= "uniq",
	GTF_ATTRIBUTE_TOKEN_READS= "freq",
	GTF_ATTRIBUTE_TOKEN_RFREQ= "rfreq",
	GTF_ATTRIBUTE_TOKEN_RPKM= "RPKM", // rpkm
	GTF_ATTRIBUTE_TOKEN_COV= "cov", 
	GTF_ATTRIBUTE_TOKEN_FWD= "fwd",
	GTF_ATTRIBUTE_TOKEN_REV= "rev",
	GTF_ATTRIBUTE_TOKEN_BID= "bid",
	GTF_ATTRIBUTE_TOKEN_SEP= "_";
	
	public static final String GTF_ATTRIBUTE_PROFILE= "profile",
		GTF_ATTRIBUTE_EXPECT= "expect";

	static final String[] GTF_ATTRIBUTES_BASE= new String[] {GTF_ATTRIBUTE_TOKEN_OBSV, GTF_ATTRIBUTE_TOKEN_PRED, GTF_ATTRIBUTE_TOKEN_BALANCED};
	static final String[] GTF_ATTRIBUTES_RESOLUTION= new String[] {GTF_ATTRIBUTE_TOKEN_ALL, GTF_ATTRIBUTE_TOKEN_TID, GTF_ATTRIBUTE_TOKEN_EXC};
	static final String[] GTF_ATTRIBUTES_MEASUREMENT= new String[] {GTF_ATTRIBUTE_TOKEN_READS, GTF_ATTRIBUTE_TOKEN_RFREQ, GTF_ATTRIBUTE_TOKEN_RPKM};

	static final String PFX_CAPACITOR= "capacitor", 
		PFX_MAPPED_READS= "mapped",
		PFX_UNMAPPED_READS= "notmapped";
	
	static final char UNDERSCORE= '_';
	
	public static final String PROPERTY_BUILD= "capacitor.build";
	public static final String PROPERTY_JDK= "capacitor.jdk";
	
	Object lock= new Object();

	BARNAdescriptor descriptor= new BARNAdescriptor();
	
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
	
	public static final String SFX_GTF= "gtf", SFX_BED= "bed";
	public static final String VERSION_ID= "1.2";
	
	public static final String SUBDIR_NATIVELIBS= "lib"+File.separator+"native", SUBDIR_LPSOLVE= "lpsolve55", LPSOLVE_LIB_NAME= "lpsolve55", LPSOLVE_JNI_NAME= "lpsolve55j";
	
	public static final String DEBUG_LP_OUT= "C:\\lp_out";
	
	public static final String CLI_CMD= "flux";
	
	public static final String
		CLI_ABBREV_COST_BOUNDS= "cb", 
		CLI_ABBREV_COST_MODEL= "cm", 
		CLI_ABBREV_COST_SPLIT= "cs";

	public static final Character  
	 CLI_SHORT_PFX= '-',
	 CLI_SHORT_ATTRIBUTES= 'a',
	 CLI_SHORT_BATCH= 'b',
	 CLI_SHORT_COMPRESSION= 'c',
	 CLI_SHORT_FORCE= 'f',
	 CLI_SHORT_LOCAL= 'l',
	 CLI_SHORT_FILENAME= 'n',
	 CLI_SHORT_HELP= 'h',
	 CLI_SHORT_OUT= 'o',
	 CLI_SHORT_PAIR= 'p',
	 CLI_SHORT_THREAD= 't',
	 CLI_SHORT_REF= 'r', 
	 CLI_SHORT_SRA= 's',
	 CLI_SHORT_UNIF= 'u',
	 CLI_SHORT_VERBOSE= 'v',
	 CLI_SHORT_EXPERIMENT= 'x'
	 ;
	
	public static final String 
		CLI_LONG_PFX= "--",
		CLI_LONG_ATTRIBUTES= "attributes", 
		CLI_LONG_BATCH= "batch", 
		CLI_LONG_COMPRESSION= "compress",
		CLI_LONG_FORCE= "force",
		CLI_LONG_HELP= "help", 
		CLI_LONG_INSTALL= "install", 
		CLI_LONG_JVM= "jvm", 
		CLI_LONG_LIB= "lib", 
		CLI_LONG_LOCAL= "local", 
		CLI_LONG_FILENAME= "name", 
		CLI_LONG_OUT= "output", 
		CLI_LONG_PAIR= "pair", 
		CLI_LONG_PROFILE= "pro", 
		CLI_LONG_REF= "ref",
		CLI_LONG_SRA= "sra", 
		CLI_LONG_SSPECIFIC= "sp", 
		CLI_LONG_THREAD= "thread",
		CLI_LONG_TMP= "tmp",
		CLI_LONG_TPX= "tpx",
		CLI_LONG_UNIF= "uniform",
		CLI_LONG_VERBOSE= "verbose";
	
	public static boolean debug= false, outputPbClusters= false;
	public boolean pairedEnd= false, stranded= false, force= false;
	int tolerance= 1000;	// +/- gene-near region
	public static boolean 
		cheatDoNotExit= false,
		cheatLearn= false, 
		cheatDisableFCheck= false,
		cheatDisableCleanup= true,
		cheatCopyFile= false;
	
	static void exit(int code) {
		String pfx= "[ASTALAVISTA] ";
		if (code< 0)
			pfx= "[CIAO] ";
		if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) 
			System.err.println(pfx+ "I'm exiting.");
		System.exit(code);
	}
	
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
						+ CLI_SHORT_PFX+ CLI_SHORT_OUT+ "|"+ CLI_LONG_PFX+ CLI_LONG_OUT+ "] [ejtgv]");
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
		
		if (strand== STRAND_SPECIFIC&& pairedEnd) {
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
			BEDobject2[] beds= null;
			boolean decompose= false;
			int nrMappingsReadsOrPairs;
			long[] sigall= null;
			double[] result= null;
	
			private float invariantTestObsSplitFreq= 0, invariantTestPredSplitFreq= 0;
			
			public LocusSolver2() {
				//super(newGene.getGeneID());
				
			}
			
			void init(Gene newGene, BEDobject2[] newBeds, boolean decompose) {
				this.gene= newGene; 
				this.beds= newBeds;
				this.decompose= decompose;
				
				nrMappingsReadsOrPairs= 0;
				if (pairedEnd) {
					mateObjects= new Object[MAX_MATES];
					mateOrient= new byte[MAX_MATES];
					mateInt= new int[MAX_MATES];
					mateP= 0;
				}

			}
			
			
			/**
			 * @deprecated
			 */
			public void run_old() {
					
				try {
					//int mapped= mapReadOrPairIDs.size();
					if (decompose) {
						if (this.gene.getTranscriptCount()== 1) {
							mapTrivial(gene.getTranscripts()[0], beds);
							outputGFF(null, null /*, null*/);
						} else {
							Graph myGraph= getGraph(this.gene);
							map(myGraph, this.gene, this.beds);
							
							GraphLPsolver2 mySolver= null;
							// != mapReadOrPairIDs.size()> 0, does also count singles
							if (nrMappingsReadsOrPairs> 0&& this.gene.getTranscriptCount()> 1) {	// OPTIMIZE			
								mySolver= getSolver(myGraph, nrMappingsReadsOrPairs* 2); // not: getMappedReadcount()
								mySolver.run();
							}
			//				if (this.gene.getTranscriptCount()== 1)
			//					System.currentTimeMillis();
							outputGFF(myGraph, events/*, mySolver*/);
							
						}
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
					FluxCapacitorNew.this.threadPool.remove(this);
	//			}
			}
	
			/**
			 * @deprecated
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
					byte flag=  (byte) 1;// (descriptor.getPairedEndInformation(beds[i].getName())- 1);	// (Fasta.getReadDir(beds[i].getName())- 1))
					if (flag < 0) {
						System.err.println("Error in readID:\n"+ beds[i].getName());
						continue;
					}
					CharSequence id= ""; // descriptor.getUniqueDescriptor(beds[i].getName());	//Fasta.getReadID(beds[i].getName())
					int sep= i, end= i+1;
					for (; end < beds.length; ++end) {
						if (!beds[end].getName().startsWith(id))
							break;
						
						if (sep== i) {
							byte flag2= (byte) 1; // (descriptor.getPairedEndInformation(beds[end].getName())- 1);
							if (flag2< 0) {
								System.err.println("Error in readID:\n"+ beds[i].getName());
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
					if (!(first|| second)) {
						for (int j = i; j < end; j++) 
							beds[j]= null;
					} else {
						//nrMappingsReadsOrPairs+= 2;
					}
					i= end- 1;
				}
				
				// build individual matrix
				int tlen= tx.getExonicLength();
				byte tstrand= tx.getStrand();
				UniversalMatrix m= new UniversalMatrix(Math.max(tlen/ 100, 10));
				for (int i = 0; i < beds.length; i++) {
					if (beds[i]== null)
						continue;
					++nrMappingsReadsOrPairs;
					int p= getBpoint(tx, beds[i]);
					if (p< 0|| p>= tlen) 
						continue;		
					boolean sense= beds[i].getStrand()== tstrand;
					m.add(p, beds[i].getLength(), tlen, sense?Constants.DIR_FORWARD:Constants.DIR_BACKWARD); 
				}
				
				// normalize biases out
				//UniversalMatrix m= profile.getMatrix(tlen);	// no good idea
				if (nrMappingsReadsOrPairs> 0) {
					//better individual matrix, also for saturation biases
					//UniversalMatrix m= profile.getMatrix(tlen);
					
					/*if (tx.getTranscriptID().contains("ENST00000262241")||
							tx.getTranscriptID().contains("ENST00000252818"))
						System.currentTimeMillis();
					*/
					//System.err.println("\n"+ tx.getTranscriptID());
					int w= m.sense.length/ 5;
					m.sums= Kernel.smoothen(Kernel.KERNEL_EPANECHNIKOV, 
							w, m.sense);
					m.suma= Kernel.smoothen(Kernel.KERNEL_EPANECHNIKOV, 
							w, m.asense);
					
					double f= m.getNfactor(0.2d);
					nrMappingsReadsOrPairs*= f;
				}
			}


			Graph getGraph(Gene gene) {
					boolean output= false;
					
					// construct graph
				long t0= System.currentTimeMillis();
				
				Graph myGraph= new Graph(gene);
				myGraph.createDefaultCoordComparator(mapLenMin);
				myGraph.constructGraph();
				//myGraph.collapseFuzzyFlanks();
				
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
	
			/** 
			 * maps reads
			 * @param myGraph
			 * @param g
			 * @param beds
			 * @return
			 */
			int map(Graph myGraph, Gene g, BEDobject2[] beds) {
				
				
				if (beds== null|| beds.length== 0)
					return 0;
				
				Transcript tx= null;
				UniversalMatrix m= null;
				if (myGraph== null) {
					tx= gene.getTranscripts()[0];
					if (decompose)
						m= new UniversalMatrix(Math.max(tx.getExonicLength()/ 100, 10));
					else
						m= profile.getMatrix(tx.getExonicLength()); 
				}
				CharSequence lastID= null;
				int mapCtr= 0;
				for (int x = 0; x< beds.length; ++x) {
					BEDobject2 dobject= beds[x];
					int mode= -1;
					if (stranded|| pairedEnd) {
						CharSequence name= dobject.getName();
						mode= descriptor.getMode(name, fromTo);
						if (pairedEnd) {
							CharSequence ID= name.subSequence(fromTo[0], fromTo[1]);
							if (pairedEnd&& !ID.equals(lastID)) {
								mateP= 0;
								lastID= ID;
							}
						}
					}

					// determine a/sense
					byte rStrand= dobject.getStrand();
					byte orient= (byte) Math.abs(gene.getStrand()- rStrand);
					assert(orient== 0|| orient== 2);
					if (stranded) {	// filter if strand known
						byte sense= descriptor.getStrand(dobject);
						if (sense!= Descriptor.ORIENTATION_UNKNOWN&& sense!= orient)
							continue;
					} 
					
					// map to annotation
					int bpoint1= -1;
					Edge target= null;
					if (myGraph== null) {
						bpoint1= getBpoint(tx, beds[x]);
						if (bpoint1< 0|| bpoint1> tx.getExonicLength())
							continue;
					} else {
						target= myGraph.getEdge2(dobject);
						if (target== null)
							continue;
					}
					
					// increase counters
					boolean mated= false;
					int mate= descriptor.getPairedEndInformation(dobject); //pairedEnd? Mate(mode): Descriptor.MATE_UNKNOWN;
					if (pairedEnd&& mate== Descriptor.MATE_2) {
						assert(mate!= Descriptor.MATE_UNKNOWN);
						
						// try pairing
						// TODO redundancy check (mateP> 1)
						for (int i = 0; i < mateP; i++) {
							
							if (mateOrient[i]== orient) {
								if (gene.getStrand()< 0) {
									if (myGraph== null)
										nrMappingsPairsWrongOrientationCt+= 2;
									else
										nrMappingsPairsWrongOrientationCg+= 2;
								} else {
									if (myGraph== null)
										nrMappingsPairsWrongOrientationWt+= 2;
									else
										nrMappingsPairsWrongOrientationWg+= 2;
								}
								continue;
							}

							// map pair
							if (decompose) {
								if (myGraph== null) { // learn or decompose trivial
									BEDobject2 bed2= (BEDobject2) mateObjects[i];
									int bpoint2= getBpoint(tx, bed2);	// TODO save bpoint before
									//if (bpoint2>= 0&& bpoint2< tx.getExonicLength())	// has to be 
									m.add(mateInt[i], bpoint2, 1, 1, tx.getExonicLength());
								} else {
									Edge target2= (Edge) mateObjects[i];
									SuperEdge se= myGraph.getSuperEdge(target, target2, true, null);
									if (se== null) {
										if (gene.getStrand()< 0) {
											if (myGraph== null)
												nrMappingsPairsWoTxEvidenceCt+= 2;
											else
												nrMappingsPairsWoTxEvidenceCg+= 2;
										} else {
											if (myGraph== null)
												nrMappingsPairsWoTxEvidenceWt+= 2;
											else
												nrMappingsPairsWoTxEvidenceWg+= 2;
										}
										continue;	// no tx evidence
									}
									se.incrReadNr();
									
									// remove anterior single mapping
									if (mateInt[i]== 0) {
										if (mateOrient[i]== Descriptor.ORIENTATION_SENSE)
											target2.decrReadNr();
										else
											target2.decrRevReadNr();
									}
								}
								
							} 
							
							// adjust counters
							mated= true;
							mapCtr+= 2;
							if (mateInt[i]== 0) {
								if (gene.getStrand()< 0) {
									if (mateOrient[i]== Descriptor.ORIENTATION_SENSE) {
										if (myGraph== null)
											--nrMappingsSingleCSt;
										else
											--nrMappingsSingleCSg;
									} else {
										if (myGraph== null)
											--nrMappingsSingleCAt;
										else
											--nrMappingsSingleCAg;
									}
								} else {
									if (mateOrient[i]== Descriptor.ORIENTATION_SENSE) {
										if (myGraph== null)
											--nrMappingsSingleWSt;
										else
											--nrMappingsSingleWSg;
									} else {
										if (myGraph== null)
											--nrMappingsSingleWAt;
										else
											--nrMappingsSingleWAg;
									}
								}
							}
							++mateInt[i];
							
						}
					} 
					
					// map single reads
					if (!mated) {
						if (myGraph== null) {
							if (!pairedEnd)
								m.add(bpoint1, 1, tx.getExonicLength(), 
										(orient== Descriptor.ORIENTATION_SENSE?Constants.DIR_FORWARD:Constants.DIR_BACKWARD));
						} else {
							if (orient== Descriptor.ORIENTATION_SENSE)
								target.incrReadNr();
							else
								target.incrRevReadNr();
						}
						
						if (pairedEnd) {
							if (mate== Descriptor.MATE_1) {
								if (myGraph== null)
									mateObjects[mateP]= beds[x];
								else
									mateObjects[mateP]= target;
								mateOrient[mateP]= orient;
								mateInt[mateP]= 0;
								++mateP;
							}
						} else
							++mapCtr;

						if (gene.getStrand()< 0) {
							if (orient== Descriptor.ORIENTATION_SENSE) {
								if (myGraph== null)
									++nrMappingsSingleCSt;
								else
									++nrMappingsSingleCSg;
							} else {
								if (myGraph== null)
									++nrMappingsSingleCAt;
								else
									++nrMappingsSingleCAg;
							}
						} else {
							if (orient== Descriptor.ORIENTATION_SENSE) {
								if (myGraph== null)
									++nrMappingsSingleWSt;
								else
									++nrMappingsSingleWSg;
							} else {
								if (myGraph== null)
									++nrMappingsSingleWAt;
								else
									++nrMappingsSingleWAg;
							}
						}
					}

				} // iterate beds
				
				return mapCtr;
			}
					
			private int nrLocusMultimaps= 0;
			/**
			 * @deprecated
			 * @param regs
			 */
			int mapRead2(Graph g, BEDobject2 dobject, boolean force) {
				
				HashSet<CharSequence> mapReadOrPairIDs= new HashSet<CharSequence>();
				HashMap<CharSequence, Vector<BEDobject2>[]> mapEndsOfPairs= new HashMap<CharSequence, Vector<BEDobject2>[]>();

				// find the edge(s) where the regions align
				// if these do not form a continous chain, create a new edge
				
	//			GFFObject[] gtfs= GFFObject.fromBed(dobject);	// TODO kill intermediary GTFs !!!
	//			DirectedRegion[] regs= new DirectedRegion[gtfs.length];
	//			for (int i = 0; i < regs.length; i++) 
	//				regs[i]= new DirectedRegion(gtfs[i]);
				
				if (force&& mapReadOrPairIDs.contains(dobject.getName())) {
					return 0;
				}
				
				byte flag= (byte) 1; // getFlag(dobject); 
				if (flag <= 0) {
					System.err.println("Error in readID:\n"+ dobject.getName());
					return 0;
				}
				CharSequence ID= ""; // getID(dobject); 	
	
				if (ID.equals("mouse_7_112_1503_1238")|| ID.equals("mouse_7_9_1185_579"))
					System.currentTimeMillis();
				
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
						if (dobject.getStrand()== dobject2.getStrand()) {
							if (gene.getStrand()< 0)
								;//++nrMappingsPairsWrongOrientationC;
							else
								;//++nrMappingsPairsWrongOrientationW;
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
							if (gene.getStrand()< 0)
								;//++nrMappingsPairsWoTxEvidenceC;
							else
								;//++nrMappingsPairsWoTxEvidenceW;
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
	
					if (sense|| (strand!= STRAND_ENABLED)) {
						target.incrReadNr();
						mapCtr= 1;
					} else if (strand!= STRAND_SPECIFIC) {
						target.incrRevReadNr();
						mapCtr= 1;
					} else {
						if (gene.getStrand()< 0) {
							//++nrMappingsWrongStrandC;
						} else {
							//++nrMappingsWrongStrandW;
						}
						mapCtr= 0;
					}
					
					
					
					if (!mapReadOrPairIDs.add(dobject.getName()))
						++nrLocusMultimaps;
	//				if (outputMapped)
	//					writeMappedRead(dobject);
					return mapCtr;
				}
				
			}
			
			public int getMappedReadcount() {
				HashSet<CharSequence> mapReadOrPairIDs= new HashSet<CharSequence>();

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
			
			
			private void outputSimple() {
				StringBuffer sb= new StringBuffer();
				double sum= nrMappingsLocusMapped;
				if (trptExprHash!= null) {
					sum= 0;
					for (int i = 0; i < gene.getTranscripts().length; i++) 
						sum+= trptExprHash.get(gene.getTranscripts()[i].getTranscriptID());
				}
				for (int i = 0; i < gene.getTranscripts().length; i++) {
					sb.append(gene.getGeneID());
					sb.append("\t");
					String tid= gene.getTranscripts()[i].getTranscriptID();
					sb.append(tid);
					sb.append("\t");
					double val= nrMappingsLocusMapped;
					if (trptExprHash!= null)
						val= trptExprHash.get(gene.getTranscripts()[i].getTranscriptID());
					val/= (gene.getTranscripts()[i].getExonicLength()/ 1000d);	// RPK normalize
					sb.append(val);
					sb.append("\t");
					if (sum== 0)
						sb.append("0.0");
					else
						sb.append(val/ sum);
					sb.append("\n");
				}
				try {
					getWriter().write(sb.toString());
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
			
			/**
			 * @deprecated
			 * @param g
			 * @param events
			 * @param solver
			 */
			private void outputGFF(Graph g, ASEvent[] events) {
				
				++nrLoci;
				if (result!= null) 
					++nrLociExp;
				
				double perM= nrReadsAll/ 1000000d;
				// deprecated
				boolean unsolvedSystem= false;	
				double valOF= result== null?0: result[0];
				if (valOF> BIG) { 
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
						if (result== null) {
							// no reads
							val= nrMappingsLocusMapped;
						} else {
	//						if (solver.getTrptExprHash().get(g.trpts[i].getTranscriptID())!= 0)
	//							System.currentTimeMillis();
							val= (float) (trptExprHash.get(tx.getTranscriptID()).doubleValue());
							if (val< 1- costBounds[0]) // 1- 0.95
								val= 0;
							try {
								assert(tt.length> 1|| val== nrMappingsReadsOrPairs);
							} catch (AssertionError e) {
								System.err.println(val+ " x "+ nrMappingsReadsOrPairs);
								getNFactor();
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
								sb.append(GTF_ATTRIBUTE_TOKEN_RPKM);
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
								sb.append(GTF_ATTRIBUTE_TOKEN_RPKM);
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
								
								sb.append(GTF_ATTRIBUTE_TOKEN_RPKM);
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
/*					if (result!= null) {
						Object[] keys= tExpMap.keySet().toArray();	
						for (int i = 0; i < keys.length; i++) {
							if (!(keys[i] instanceof String))
								continue;
							if (tExpMap.get(keys[i])<0)
								tExpMap.put((String) keys[i], 0d);	// TODO ugly
						}
					}
*/					
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
				if (false&& result!= null) {
					ArrayList<Edge> cc= new ArrayList<Edge>();
					if (result!= null) {
						Iterator<Object> iter= null; //solver.getConstraintHash().keySet().iterator();
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
	
			private GraphLPsolver2 getSolver(Graph g, int mappedReads) {
			
				GraphLPsolver2 solver= new GraphLPsolver2(g, readLenMin, 
						pairedEnd?insertMinMax:null, mappedReads, 
						(strand== STRAND_ENABLED), 
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
							sb.append(GFF_FEATURE_PAIRED);
						else 
							sb.append(GFF_FEATURE_JUNCTION);
					} else
						sb.append(GFF_FEATURE_FRAGMENT);
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
							sb.append(GTF_ATTRIBUTE_LENGTH);
							sb.append(GTF_ATTRIBUTE_TOKEN_SEP);
							sb.append(x== 0? GTF_ATTRIBUTE_TOKEN_ALL: (x==1? GTF_ATTRIBUTE_TOKEN_TID: GTF_ATTRIBUTE_TOKEN_EXC));
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
						
						ReadStatCalculator calc= (i== 0? FluxCapacitorNew.this: solver);
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
									sb.append(i==0? GTF_ATTRIBUTE_TOKEN_OBSV: (i==1? GTF_ATTRIBUTE_TOKEN_PRED:GTF_ATTRIBUTE_TOKEN_BALANCED));
									sb.append(GTF_ATTRIBUTE_TOKEN_SEP);
									sb.append(j== 0?GTF_ATTRIBUTE_TOKEN_ALL:(j==1? GTF_ATTRIBUTE_TOKEN_TID: GTF_ATTRIBUTE_TOKEN_EXC));
									sb.append(GTF_ATTRIBUTE_TOKEN_SEP);				
									sb.append(x== 0?GTF_ATTRIBUTE_TOKEN_READS: GTF_ATTRIBUTE_TOKEN_COV);	// pairedEnd?GTF_ATTRIBUTE_TOKEN_COV:GTF_ATTRIBUTE_TOKEN_RPKM
									sb.append(" ");	// no \"
								}
								for (int m = 0; m < tid.length; m++) {
									long[] sig= j== 0? sigall: tid[m];
			
									if (calc== null) {
										if (output)
											sb.append(VALUE_NA);
									} else {
										
										boolean normalized= i== 2? true: false;
										double val= ms* (j== 0? calc.getReads(eeV.elementAt(m), BYTE_0, sig, normalized):  // ALL								
											calc.getReadsAvg(eeV.elementAt(m), BYTE_0, g, sig, excl, normalized));	// TID, excl 							
										
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
												sb.append(length== 0? FLOAT_STRING_0: Float.toString((float) val));
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
			 * @deprecated
			 * @param g
			 * @param events
			 * @param solver
			 */
			private void outputGFF_save(Graph g, ASEvent[] events, GraphLPsolver solver) {
						++nrLoci;
						if (solver!= null) 
							++nrLociExp;
						double perM= nrReadsAll/ 1000000d;
						// deprecated
							boolean unsolvedSystem= false;	
							double valOF= solver== null?0: solver.getValObjFunc();
							if (valOF> BIG) { 
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
	
			
			/**
			 * @param tx
			 * @param bed
			 * @return
			 */
			private boolean contains(Transcript tx, BEDobject2 bed) {
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
			
			
			/**
			 * @deprecated
			 * @param tx
			 * @param beds
			 */
			private void learn(Transcript tx, BEDobject2[] beds) {
							
				if (beds== null|| beds.length== 0)
					return;
				
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
					if (strand== STRAND_SPECIFIC&& bed1.getStrand()!= tx.getStrand()) {
						//++nrMappingsWrongStrand;
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
					if (rlen1> readLenMax)	// readLenMax< 0|| 
						readLenMax= rlen1;
					
					++nrReadsSingleLociMapped;
					if (pairedEnd) {	// && flag== Descriptor.MATE_2

						CharSequence name= bed1.getName();
						int mode= descriptor.getMode(name, fromTo);
						if (mode <0) {
							System.err.println("Error in readID:\n"+ bed1.getName());
							continue;
						}
						byte flag= descriptor.getPairedEndInformation(bed1); // getMate(mode);
						CharSequence id= name.subSequence(fromTo[0], fromTo[1]);
						
						
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
						for (int j = 0; j < oo.length; j++) {
							if (j==flag|| oo[j]== null)
								continue;
							for (int k = 0; k < oo[j].length; k++) {
								BEDobject2 bed2= oo[j][k];
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
									++nrReadsSingleLociPairsMapped;	
								}
							}
						}
							
					} else {
						
						System.err.println("[DEACTIVATED] no single reads for now");
						if (1== 1)
							System.exit(0);
						int[] ePos= null;
						int binLen= 0;
						TProfile tprofile= null;
						
						byte dir= Constants.DIR_FORWARD;
						if (strand== STRAND_ENABLED) {
							//assert(ePos[0]> 0);	// can be, out of range
							if (!tx.getTranscriptID().startsWith("AT1G04050")) {
								if (ePos[0]> 0&& ePos[1]<= elen) {
									if (beds[i].getStrand()!= tx.getStrand()) {
										dir= Constants.DIR_BACKWARD;
										if (ePos[1]<= elen) {
											//++profileStubRev[binLen][p];
										}
									} else {
										if (ePos[0]> 0) {
											//++profileStub[binLen][p];
										}
									}
								}
							}
						} else if (!uniform)
							tprofile.addRead(ePos[0],readLenMin,dir);
					}
				}
				
			}


			/**
			 * @deprecated
			 * @param regs
			 */
			int mapRead(Graph g, BEDobject2 dobject, boolean force) { 
				
				HashSet<CharSequence> mapReadOrPairIDs= new HashSet<CharSequence>();
				HashMap<CharSequence, Vector<BEDobject2>[]> mapEndsOfPairs= new HashMap<CharSequence, Vector<BEDobject2>[]>();

				// find the edge(s) where the regions align
				// if these do not form a continous chain, create a new edge
				
	//			GFFObject[] gtfs= GFFObject.fromBed(dobject);	// TODO kill intermediary GTFs !!!
	//			DirectedRegion[] regs= new DirectedRegion[gtfs.length];
	//			for (int i = 0; i < regs.length; i++) 
	//				regs[i]= new DirectedRegion(gtfs[i]);
				
				if (force&& mapReadOrPairIDs.contains(dobject.getName())) {
					return 0;
				}
				
				byte flag= 1; // getFlag(dobject);  	
				CharSequence ID= ""; // getID(dobject); 	
	
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
	
					if (sense|| (strand!= STRAND_ENABLED)) {
						target.incrReadNr();
						mapCtr= 1;
					} else if (strand!= STRAND_SPECIFIC) {
						target.incrRevReadNr();
						mapCtr= 1;
					} else {
						//++nrMappingsWrongStrand;
						mapCtr= 0;
					}
					
					
					
					if (!mapReadOrPairIDs.add(dobject.getName()))
						++nrLocusMultimaps;
	//				if (outputMapped)
	//					writeMappedRead(dobject);
					return mapCtr;
				}
				
			}

			public void run() {
							
				try {
					
//					if (!gene.getGeneID().equals("chr14:50256232-50367589C")) 
//						return;
					
					// init
					trptExprHash= null;
					
					// data structure
					Graph myGraph= null;
					if (decompose&& this.gene.getTranscriptCount()> 1) {
						myGraph= getGraph(this.gene);
						if (outputEvent)
							events= myGraph.getEvents(2);
						int nrSJ= 0;
						if (!pairedEnd)
							nrSJ= myGraph.addEJ(27);						
					}
					
					if (beds!= null&& beds.length> 0) {
						
						// mapping
						if (pairedEnd) 
							Arrays.sort(beds, BEDobject2.DEFAULT_ID_COMPARATOR);	// TODO file
						nrMappingsLocusMapped= 0;
						nrMappingsLocusMapped= map(myGraph, this.gene, this.beds);
						// DEBUG
						if (myGraph== null) {
							Graph gg= getGraph(this.gene);
							int chk= map(myGraph, this.gene, this.beds);
							if (chk!= nrMappingsLocusMapped)
								System.currentTimeMillis();
						}
						if (gene.getStrand()< 0) {
							if (myGraph== null)
								nrMappingsPairsMappedCt+= nrMappingsLocusMapped;
							else
								nrMappingsPairsMappedCg+= nrMappingsLocusMapped;
						} else {
							if (myGraph== null)
								nrMappingsPairsMappedWt+= nrMappingsLocusMapped;
							else
								nrMappingsPairsMappedWg+= nrMappingsLocusMapped;
						}
						
						// solving
						if (nrMappingsLocusMapped> 0&& this.gene.getTranscriptCount()> 1) {  
							solve(myGraph);
						}
					}
					
					// output
					//outputGFF(myGraph, null);
					outputSimple();
					
				} catch (Throwable e) {
					System.err.println("\n[ERROR] in locus "+ gene.getGeneID());
					e.printStackTrace();
					System.err.flush();
				}
				
				beds= null;
				gene= null;
			}

		/**
		 * #(Edge,Transcript) x (int[],int[][])
		 * watch out with hash function of edges
		 * @return
		 */
		private HashMap<String, Double> mapCCheck= null; 
		double[] ccheck;
		boolean stopperBalance= false;
		boolean stopperBalance2= false;
		public HashMap<String, Integer> setConstraints_save(Graph g, boolean count) {
			
			//boolean debug= false;
	//		if (g.trpts[0].getTranscriptID().equals("ENST00000323441")) {
	//			debug= true;
	//		}
			
			// transcript constraint variables
			Transcript[] trpts= g.trpts;
			IntVector v= null, w= null;
			DoubleVector u= null;
			HashMap<String, Integer> tMap= null;
			if (count)
				constraintCtr+= trpts.length;
			else { 
				try {
					getLPsolve().setLpName(trpts[0].getTranscriptID());
				} catch (LpSolveException e1) {
					e1.printStackTrace();
				}			
				tMap= new HashMap<String, Integer>(trpts.length* 2);
				for (int i = 0; i < trpts.length; i++) { 
					tMap.put(trpts[i].getTranscriptID(), ++constraintCtr);
					if (debug)
						System.out.println("C"+constraintCtr+"\t"+trpts[i].getTranscriptID());
				}
				v= new IntVector();	// indices for transcript/part 
				w= new IntVector();	// indices for cost function
				u= new DoubleVector();	// observation, bases for cost function
			}
			
			// iterate edges
			Edge[] edges= g.getExonicEdgesInGenomicOrder();
			if (!count) {
				mapCCheck= new HashMap<String, Double>();
				for (int i = 0; i < trpts.length; i++) 
					mapCCheck.put(trpts[i].getTranscriptID(), 0d);
			}
			for (int i = 0; i < edges.length; i++) {
				Edge e= edges[i];
				if (!e.isExonic())
					continue;
				
				// the base edge
				int basElen= e.length();
				Transcript[] tt= g.decodeTset(e.getTranscripts());
				HashMap<Edge, IntVector> mapE= new HashMap<Edge, IntVector>();			
				// sense/anti
				for (int sa= 0; sa< 2; ++sa) {
					for (int x = 0; x < tt.length; ++x) {
						if (!count)
							v.removeAll();
						long[] sig= g.encodeTset(tt[x]);
//						if (e.toString().equals("67937779-67937952^")&& !count)
//							System.currentTimeMillis();
						getConstraints(e, sig, v, mapE, sa== 0, count);
						
						// add transcript constraint
						if (!count) {
							int[] idx= new int[v.length+ 1]; // obs parts+ tx frac
							System.arraycopy(v.vector, 0, idx, 0, v.length);
							idx[idx.length- 1]= tMap.get(tt[x].getTranscriptID());
							double[] val= new double[idx.length];
							Arrays.fill(val, 1d);
							int tlen= tt[x].getExonicLength();
							UniversalMatrix m= profile.getMatrix(tlen);
							double f= m.getFrac(
										tt[x].getExonicPosition(e.getFrac(true)),
										tt[x].getExonicPosition(e.getFrac(false)),
										tlen,
										sa== 0?Constants.DIR_FORWARD:Constants.DIR_BACKWARD);
							//System.err.println(f);
							if (Double.isInfinite(f)|| Double.isNaN(f)) {
								System.err.println("infinite value");
								f= m.getFrac(
										tt[x].getExonicPosition(e.getFrac(true)),
										tt[x].getExonicPosition(e.getFrac(false)),
										tlen,
										sa== 0?Constants.DIR_FORWARD:Constants.DIR_BACKWARD);
							}
							mapCCheck.put(tt[x].getTranscriptID(),
									mapCCheck.get(tt[x].getTranscriptID())+ f);
							val[val.length- 1]= -f;
							if (debug&& !count) {
								StringBuilder sb= new StringBuilder(e.toString());
								sb.append(": ");
								for (int k = 0; k < idx.length; k++) {
									sb.append(val[k]>0?"+":"");
									sb.append(val[k]%1==0?((int) val[k]):val[k]);
									sb.append("C");
									sb.append(idx[k]+" ");
								}
								sb.append("= 0");
								System.out.println(sb);
							}
							addConstraintToLp(idx, val, LpSolve.EQ, 0);
							++restrNr;
							
							// force junction coverage
							if (val.length== 3) {
								assert(val[0]<0&&val[1]>0&&val[2]>0);
								int fixedreadlen= 27;
								double frac= 0.9;
								val[0]= 0d;
								val[1]= frac/ (basElen- fixedreadlen+ 1);
								val[2]= -1d/ (fixedreadlen- 1);
								addConstraintToLp(idx, val, LpSolve.LE, 0);
								++restrNr;
							}							
						}
					}
					// add edge constraints
					Edge[] ee= new Edge[mapE.size()];
					mapE.keySet().toArray(ee);
					
					// total obs
	/*				int sumObs= 0;
					for (int j = 0; j < ee.length; j++) {
						Edge f= ee[j];
						boolean paird= (f instanceof SuperEdge)&& ((SuperEdge) f).isPend();
						int nr= ((paird|| sa== 0)? f.getReadNr(): f.getRevReadNr());
						sumObs+= nr;
					}
	*/				
					for (int j = 0; j < ee.length; j++) {
						Edge f= ee[j];
						
						// paired-end super-edge
						boolean paird= (f instanceof SuperEdge)&& ((SuperEdge) f).isPend();
						int nr= ((paird|| sa== 0)? f.getReadNr(): f.getRevReadNr());
						// coverage or cost calc
/*						int len= (f instanceof SuperEdge)? 27- 1: f.length()- 27+ 1;
						double cov= 0;
						if (len<= 0) {
							assert(nr== 0);
						} else
							cov= nr/ (double) len;
						double cost= 1;
						if (len> 0) {
							cost= 1d+ (10* Math.exp(-nr/ 10));
							//cost= 1d+ (10* Math.pow((nr+ 1)/10d, -4));
						}						
						if (f instanceof SuperEdge) {
							cost= 1d+ (10* Math.exp(-nr/10d));
						} else {
							if (len> 0)
								cost= 1d+ (Math.pow(cov+ 1, -1));
						}
						if (Double.isNaN(cost)|| Double.isInfinite(cost))
							System.currentTimeMillis();
*/							
						v= mapE.remove(f);
						if (count)
							constraintCtr+= 2;
						else {
							int[] idx= new int[v.length+ 2];	// +/-
							System.arraycopy(v.vector, 0, idx, 0, v.length);
							int c= ++constraintCtr;
							// plus on not paired edges at 0-cost, it substracts from obs
							// usual edges cost when in single read mode
							if ((!pairedEnd)|| (pairedEnd&& paird)) {
								w.add(c);
//								if (nr== 0)
//									u.add(Double.NEGATIVE_INFINITY);
//								else
									u.add(-nr);
							}
							idx[idx.length- 2]= c;
							// plus has to be limited, it substracts
							int lim= (paird||(!pairedEnd))? Math.max(nr- 1, 0): nr; // Math.max(nr- 1, 0): nr; 
							try {
								getLPsolve().setUpbo(constraintCtr, lim);
							} catch (LpSolveException e1) {
								e1.printStackTrace();
							}
								
							c= ++constraintCtr;
							// do not limit adding, even with f= 100 unsolvable systems
							// adding reads always costs, also on single edges
							w.add(c);
							u.add(nr);
							idx[idx.length- 1]= c;
							double[] val= new double[idx.length];
							Arrays.fill(val, 1d);
							val[val.length- 1]= -1d;
							if (debug&& !count) {
								StringBuilder sb= new StringBuilder(f.toString());
								sb.append(": ");
								for (int k = 0; k < idx.length; k++) {
									sb.append(val[k]==1?"+C":"-C");
									sb.append(idx[k]+" ");
								}
								sb.append("= "+nr);
								System.out.println(sb);
							}
							addConstraintToLp(idx, val, LpSolve.EQ, nr);
							++restrNr;
						}
					}
				}
			} // end all edges
			
			if (count) {
				if (stopperBalance|| stopperBalance2)
					constraintCtr+= 2;
				return null;
			}
			
			if (stopperBalance2) {
				int[] idx= new int[tMap.size()+ 2];
				Iterator<Integer> iter= tMap.values().iterator();
				int ctr= 0;
				while (iter.hasNext()) {
					idx[ctr++]= iter.next();
				}
				idx[ctr++]= ++constraintCtr;
				idx[ctr++]= ++constraintCtr;
				double[] vals= new double[idx.length];
				Arrays.fill(vals, 1d);
				vals[vals.length- 1]= -1d;
				addConstraintToLp(idx, vals, LpSolve.EQ, nrMappingsLocusMapped);
				++restrNr;
			}
			
			// set objective function/costs
	//		double[] a= createArray(w.toIntArray());	// linear costs
			//double min= Math.exp(-nrMappingsLocusMapped);
			int ext= 0;
			if (stopperBalance)
				ext= 2;
			double[] costs= new double[w.size()+ ext];
			double[] pmPfx= new double[w.size()+ ext];
			double[] dd= u.toDoubleArray();
			for (int i = 0; i < dd.length; i++) {
				dd[i]= Math.abs(dd[i]);
			}
			Distribution dist= new Distribution(dd);
			double med= dist.getMedian();
			if (med== 0)
				med= 1;
			StringBuilder sb= null;
			if (debug) {
				sb= new StringBuilder();
			}
			for (int i = 0; i < costs.length- ext; i++) {
				double x= u.get(i);	// <0 : read substracting
				pmPfx[i]= ((x< 0|| x== Double.NEGATIVE_INFINITY)?-1d:1d);
				//costs[i]= 1+ Math.exp(-Math.log(1+ x));	// neutralizes
				//costs[i]= 1d+ Math.pow(x+1d, -1d/2d);
				//costs[i]= 0;	// 
				//costs[i]= 1d+ Math.log(x+ 1);
				//costs[i]= 1d+ Math.exp(-x);
				//costs[i]= 100d/ (x+ 1);	// last
				
				double c= 1d;
				if (!pairedEnd) {
					double y= (x== Double.NEGATIVE_INFINITY)?0:x;
					c= (10d* Math.exp(-(Math.abs(y)/ med)));
					if (x< 0/*|| x== Double.NEGATIVE_INFINITY*/)
						c+= 1;
				}
				//c= med/ (x+ 1);
				//c= 10d/ (1+ Math.pow(x, 1));
				//c= 10d/ (1+ x);
				//c= x;
				if (Double.isInfinite(c)|| Double.isNaN(c))
					System.currentTimeMillis();
				costs[i]= c;
				
				
				if (debug) { 
					sb.append("C");
					sb.append(w.get(i));
					sb.append("\t");
					sb.append(StringUtils.fprint(c, 3));
					sb.append("\n");
				}
				//costs[i]= 1d+ Math.log(x+ 1);
				//costs[i]= 1d+ Math.sqrt(x+ 1);	// best
				//costs[i]= 1d+ Math.pow(x+ 1, 1d/ 3d);
				
				//costs[i]= (Math.log(x+ 1)/ (x+ 1));	// logdiv
				//costs[i]= 1d/ (1d+ Math.log(x+ 1d));	// divlog
			}
			if (debug) {
				System.out.println("COSTS:");
				System.out.println(sb);
			}
			int[] idx= null;
			if (stopperBalance) {
				idx= new int[w.length+ ext];
				pmPfx[pmPfx.length- 2]= 1;
				pmPfx[pmPfx.length- 1]= -1;
				System.arraycopy(w.vector, 0, idx, 0, w.length);
				idx[pmPfx.length- 2]= ++constraintCtr;
				idx[pmPfx.length- 1]= ++constraintCtr;
				addConstraintToLp(idx, pmPfx, LpSolve.EQ, 0);
				++restrNr;
				
				costs[costs.length- 2]= 1;
				costs[costs.length- 1]= 1;
			} else {
				if (stopperBalance2) {
					idx= new int[w.length+ 2];
					double[] costs2= new double[costs.length+ 2];
					System.arraycopy(w.vector, 0, idx, 0, w.length);
					idx[idx.length- 2]= constraintCtr- 2;
					idx[idx.length- 1]= constraintCtr- 1;
					System.arraycopy(costs, 0, costs2, 0, costs.length);
					costs2[costs2.length- 2]= 1d;
					costs2[costs2.length- 1]= 1d;
					costs= costs2;
				} else
					idx= w.toIntArray();
			}
			
			double[] a= createArray(idx, costs);
			try {
				getLPsolve().setObjFn(a);
				getLPsolve().setMinim();
			} catch (LpSolveException e) {
				e.printStackTrace();
			}
			
			// TODO consistency check
			Object[] oo= mapCCheck.keySet().toArray();
			for (int i = 0; i < oo.length; i++) {
				double val= mapCCheck.get(oo[i]);
				//System.err.println("check "+ val);
				if (Math.abs(2d- val)> 0.2d)
					System.err.println("Fraction inconsistency "+ oo[i]+"\t"+val);
			}
			
			return tMap;
		}

		private HashMap<Object,Double> trptExprHash;
			
		public strictfp int solve(Graph g) {
					
			//		if (costModel== COSTS_LINEAR&& costBounds== null)
			//			costBounds= new int[] {100,100};
					
					long t0= System.currentTimeMillis();
					debug= false;	
					// chr12:67919663-67951290W (CPSF6)
					// chr1:1205831-1217267W (5 NUMERIC) ok
					// chr1:154318993-154376504W (2 INFEASABLE) ok
					// chr1:8843648-8978713C (5 NUMERIC) ok
					// chr1:1374932-1459930W (2 INFEASABLE) ok
					// chr1:44978079-45006026W (KIF2C)
					// chr6:30775563-30793645C (MDC1)
					// chr2:47863725-47887588W (MSH6)
					if (gene.getGeneID().equals("chr2:47863725-47887588W")) {
						debug= true;
						System.err.println();
					} else 
						return -1;
					
					collectStats(g);
					
					int c= 0;
					try {
						setCAll(MODE_COUNT, g);
						c= constraintCtr;
					} catch (LpSolveException e2) {
						e2.printStackTrace();
						return -1;
					}
					//int resSave= restrNr;	// doesnt work properly
					
//					decVal= new double[g.trpts.length+ 2]; 
//					Arrays.fill(decVal, 0, decVal.length- 1, 1d);
//					decVal[decVal.length- 1]= -1d;
					decVal= new double[g.trpts.length* 3];
					for (int i = 0; i < decVal.length; i++) 
						decVal[i]= ((i-1)% 3== 0? -1d: 1d);	// +-+ +-+ +-+
					decIdx= new int[decVal.length];
					boundIdx= new int[2];
					boundVal= new double[2];
					boundVal[0]= 1d;
					boundVal[1]= -1d;
					
					getLPsolve(c);
					constraintCtr= 0;
					restrNr= 0;
						
			/*		if (g.trpts[0].getTranscriptID().startsWith("NM_001168507")||
							g.trpts[0].getTranscriptID().startsWith("NM_173453"))
						debug= true;
			*/			
					
			
					// set up program
					ccheck= new double[g.trpts.length];
					Arrays.fill(ccheck, 0d);
					try {
						int cc= setCAll(MODE_CONSTRAINT, g);
						//int resNow= restrNr;
						//if (resNow!= resSave)
						//	System.currentTimeMillis();
						if (c!= cc)
							System.currentTimeMillis();
						assert(c== cc);
						//assert(resNow== resSave);
					} catch (LpSolveException e1) {
						e1.printStackTrace();
						return -1;
					}
					
					// write out
					String tmpOutFName= null;
					int ret= 0;
					if (fileLPdir!= null) {
						tmpOutFName= fileLPdir+ File.separator
							+ gene.getGeneID().replace(":", "_")
							+ SFX_LPOUT;
			
						try {
							getLPsolve().setOutputfile(tmpOutFName);
						} catch (LpSolveException e) {
							//e.printStackTrace();
							if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
								System.err.println("[FATAL] failed to set lp output to:\n\t"+ tmpOutFName);
						}
					}
					
				
			//		for (int i = 0; i < g.trpts.length; i++) {
			//			checkPercent(g.trpts[i]);
			//		}
					
					try {
						//getLPsolve().printLp();
						if (tmpOutFName== null)
							getLPsolve().setVerbose(LpSolve.CRITICAL);	//shut up ! IMPORTANT, SEVERE, CRITICAL
						t0= System.currentTimeMillis();
						//getLPsolve().setPresolve(arg0, arg1);
						//System.err.println("scaling "+ gene.getGeneID()+ " ("+constraintCtr+","+restrNr+")");
						//getLPsolve().setScaling(LpSolve.SCALE_GEOMETRIC| LpSolve.SCALE_POWER2);
						
						getLPsolve().setScaling(LpSolve.SCALE_GEOMETRIC+ LpSolve.SCALE_EQUILIBRATE + LpSolve.SCALE_INTEGERS);
						getLPsolve().setImprove(LpSolve.IMPROVE_BBSIMPLEX);
						
						ret= getLPsolve().solve();
						//System.err.println("solved "+ gene.getGeneID());
			/*			 NOMEMORY (-2)  	Out of memory
						   OPTIMAL (0) 	An optimal solution was obtained
						SUBOPTIMAL (1) 	The model is sub-optimal. Only happens if there are integer variables and there is already an integer solution found. The solution is not guaranteed the most optimal one.
			
						 							* A timeout occured (set via set_timeout or with the -timeout option in lp_solve)
						 							* set_break_at_first was called so that the first found integer solution is found (-f option in lp_solve)
						 							* set_break_at_value was called so that when integer solution is found that is better than the specified value that it stops (-o option in lp_solve)
						 							* set_mip_gap was called (-g/-ga/-gr options in lp_solve) to specify a MIP gap
						 							* An abort function is installed (put_abortfunc) and this function returned TRUE
						 							* At some point not enough memory could not be allocated 
			
						INFEASIBLE (2) 	The model is infeasible
						UNBOUNDED (3) 	The model is unbounded
						DEGENERATE (4) 	The model is degenerative
						NUMFAILURE (5) 	Numerical failure encountered
						USERABORT (6) 	The abort routine returned TRUE. See put_abortfunc
						TIMEOUT (7) 	A timeout occurred. A timeout was set via set_timeout
						PRESOLVED (9) 	The model could be solved by presolve. This can only happen if presolve is active via set_presolve
						PROCFAIL (10) 	The B&B routine failed
						PROCBREAK (11) 	The B&B was stopped because of a break-at-first (see set_break_at_first) or a break-at-value (see set_break_at_value)
						FEASFOUND (12) 	A feasible B&B solution was found
						NOFEASFOUND (13) 	No feasible B&B solution found
			*/
						if (ret!= 0) {
							System.err.println("insolvable system: "+ gene.getGeneID());
							System.err.println(RETURN_VERBOSE[ret+ 2]);
							debug= true;
						}
						
					} catch (LpSolveException e) {
						e.printStackTrace();
					}
			
					// get transcription expression levels		
					trptExprHash= getResult();	// tMap
					//normalizeBack2LocusExpr(trptExprHash);
					
					if (debug|| tmpOutFName!= null) {
						//getLPsolve().printLp();
						//restrNr= 0;	// Noo, needed for value retrieval
						constraintCtr= 0;
						try {
							int cc= setCAll(MODE_COMPLETE, g);
						} catch (LpSolveException e) {
							e.printStackTrace();
						}
						getLPsolve().printObjective();
						//getLPsolve().printSolution(1);
						printSolution(System.out, g, result, costIdx, costVal);
						System.exit(-1);
					}
						
					getLPsolve().deleteLp();	// closes file outFName
					
					PrintStream p= null;
					if (debug) {
						System.err.flush(); // doesnt work
						p= System.out;
						for (int i = 0; i < g.trpts.length; i++) {
							; //getAllPercent(g.trpts[i]);
						}
					}
					if (tmpOutFName!= null) {
						try {
							p= new PrintStream(new FileOutputStream(tmpOutFName, true));				
						} catch (Exception e) {
							e.printStackTrace();
						}
					}
					
					// output additionally
					if (p!= null) {
						Iterator<Object> idIter= trptExprHash.keySet().iterator();
						while(idIter.hasNext()) {
							Object o= idIter.next();
							if (!(o instanceof String))
								continue;
							String id= (String) o;
							p.println(id+" "+trptExprHash.get(id));
						}
						p.println();
						p.println("Settings:");
						p.println("paired-end\t"+pairedEnd);
						if (costBounds!= null) {
							if (!Double.isNaN(costBounds[0]))
								p.print("cost boundaries:\tlower /"+costBounds[0]);
							if (!Double.isNaN(costBounds[1]))
								p.print(" upper *"+costBounds[1]);
							p.println();
						}
						p.print("costfunc\t");
						//printConstraintHash(p);
						
						p.println();
						p.println("Transcripts:");
						for (int i = 0; i < g.trpts.length; i++) {
							p.print(g.trpts[i].getTranscriptID()+"\t");
							SpliceSite[] ss= g.trpts[i].getSpliceSitesAll();
							for (int j = 0; j < ss.length; j++) 
								p.print(ss[j].toString());
							p.println();
						}
					}
					
					// close file
					if (p!= null)
						p.flush();
					if (tmpOutFName!= null) 
						p.close();
					
			//		System.err.println("solved "+g.trpts[0].getTranscriptID()+": "+g.trpts.length+" trpts, "+constraintCtr+" constr, "+restrNr+" restr"
			//				+ ((System.currentTimeMillis()- t0)/ 1000)+ " sec.");
					if (debug)
						System.currentTimeMillis();
					
					return ret;
				}

		float medExon= 0, medJunc= 0;
		private void collectStats(Graph g) {
			Edge[] edges= g.getExonicEdgesInGenomicOrder();
			IntVector vecExons= new IntVector(), vecJunc= new IntVector();
			for (int i = 0; i < edges.length; i++) {
				Edge e= edges[i];
				if (!e.isExonic())
					continue;
				
				int exSense= e.getReadNr(), exAsense= e.getRevReadNr();
				for (int j = 0; e.getSuperEdges()!= null&& j < e.getSuperEdges().size(); j++) {
					SuperEdge se= (SuperEdge) e.getSuperEdges().elementAt(j);
					int juncSense= 0, juncAsense= 0;
					if (se.getEdges()[0]== e) {
						if (se.isPend())
							exSense+= se.getReadNr();
						else
							juncSense+= se.getReadNr();
					} else if (se.getEdges()[se.getEdges().length- 1]== e){
						if (se.isPend())
							exAsense+= se.getReadNr();
						else
							juncAsense+= se.getRevReadNr();
					}
					
					if (!se.isPend()) {
						for (int k = 0; se.getSuperEdges()!= null&& k < se.getSuperEdges().size(); k++) {
							SuperEdge sse= (SuperEdge) se.getSuperEdges().elementAt(k);
							assert(sse.isPend());
							if (sse.getEdges()[0]== se)
								juncSense+= sse.getReadNr();
							else	// must be last edge
								juncAsense+= sse.getReadNr();
						}
					}
					vecJunc.add(juncSense);
					vecJunc.add(juncAsense);
				}
				vecExons.add(exSense);
				vecExons.add(exAsense);
			}
			Distribution d= new Distribution(vecExons.toIntArray());
			medExon= (float) d.getMedian();
			if (medExon== 0)
				++medExon;
			d= new Distribution(vecJunc.toIntArray());
			medJunc= (float) d.getMedian();
			if (medJunc== 0)
				medJunc= 1;
		}

		private void printSolution(PrintStream p, Graph g, double[] result, int[] costIdx, double[] costVal) {
			
			StringBuilder sb= new StringBuilder("Variable\tValue\tFactor\tCosts:\n");
			double sumCosts= 0d, sumChanges= 0d;
			for (int i = restrNr+ 1; i < result.length; i++) {
				sb.append("C");
				sb.append(Integer.toString(i- restrNr));
				sb.append("\t");
				sb.append(StringUtils.fprint(result[i], 2));
				sumChanges+= result[i];
				sb.append("\t");
				int pos= Arrays.binarySearch(costIdx, (i- restrNr));
				if (pos< 0) {
					sb.append(".\t.\n");
					continue;
				}
				sb.append(StringUtils.fprint(costVal[pos], 2));
				sb.append("\t");
				double costs= costVal[pos]* result[i];
				sumCosts+= costs;
				sb.append(StringUtils.fprint(costVal[pos] * result[i], 2));
				sb.append("\n");
			}
			sb.append("OF\t");
			sb.append(StringUtils.fprint(sumChanges, 2));
			sb.append("\t.\t");
			sb.append(StringUtils.fprint(sumCosts, 2));
			sb.append("\n");
			p.println(sb.toString());
		}

		LpSolve lpSolve;
		LpSolve getLPsolve() {
				if (lpSolve == null) {
/*					if (writeFile)
						try {
							getLPWriter().flush();
							getLPWriter().close();
							lpSolve= LpSolve.readLp(fileLPinput.getAbsolutePath(), LpSolve.IMPORTANT, null);
						} catch (Exception e) {
							e.printStackTrace();
						}
					else
*/						
						try {
		//				Iterator iter= getConstraintHash().keySet().iterator();
		//				int edgeCtr= 0, varCtr= 0;
		//				while (iter.hasNext()) {
		//					Object o= iter.next();
		//					if (o instanceof Edge)
		//						++edgeCtr;
		//					else if (o instanceof Variation)
		//						++varCtr;
		//				}
		//				
							lpSolve = LpSolve.makeLp(0, constraintCtr);	// no 0-column here
							
						} catch (LpSolveException e) {				
							e.printStackTrace();
						}
					
				}
		
				return lpSolve;
			}
		
		LpSolve getLPsolve(int nrC) {
			
			if (lpSolve == null) {
					try {
						lpSolve = LpSolve.makeLp(0, nrC);	// no 0-column here
					} catch (LpSolveException e) {				
						e.printStackTrace();
					}
				
			}
	
			return lpSolve;
		}
		
		LpSolve getLPsolve(int nrR, int nrC) {
			
			if (lpSolve == null) {
					try {
						lpSolve = LpSolve.makeLp(nrR, nrC);	// no 0-column here
					} catch (LpSolveException e) {				
						e.printStackTrace();
					}
				
			}
	
			return lpSolve;
		}

		/**
		 * create an array where certain indices are set to certain values
		 * 
		 * @param constIdx indices of constraints
		 * @param constVal values of constraints
		 * @return
		 */
		private double[] createArray(int[] constIdx, double[] constVal) {
			
			double[] a= new double[constraintCtr+ 1];
			for (int i = 0; i < a.length; i++) 
				a[i]= 0d;
			for (int i = 0; i < constIdx.length; i++) 
				a[constIdx[i]]= constVal[i];		
			return a;
		}

		private HashMap<Object, Double> getResult() {
			
//			int rows= getLPsolve().getNrows();
//			int cols= getLPsolve().getNcolumns();
			result= new double[1+ restrNr+ constraintCtr];
			try {
				getLPsolve().getPrimalSolution(result);
			} catch (LpSolveException e1) {
				e1.printStackTrace();
			}
			//valObjFunc= result[0];	
			HashMap<Object,Double> trptExprHash= new HashMap<Object,Double>();
			
			// normalize fraction rounding errors
			Transcript[] trpts= gene.getTranscripts();
			double sum= 0;
			for (int i = 0; i < trpts.length; i++) { 
				
				// (sense+ asense)/ 2
				double tot= (result[1+ restrNr+ (i* 2)]+ result[1+ restrNr+ (i* 2)+ 1])/ 2d;
				// normalize by summing inaccuracy
				tot/= ccheck[i]/ 2d;
				sum+= tot;
				if (Double.isNaN(tot))
					System.currentTimeMillis();
				trptExprHash.put(trpts[i].getTranscriptID(), tot);
			}
			if (sum== 0)
				return trptExprHash;
			
			// normalize flow network over-/underprediction
			double nfac= nrMappingsLocusMapped/ (double) sum;
			for (int i = 0; sum!= 0&& i < trpts.length; i++) {
				double x= trptExprHash.get(trpts[i].getTranscriptID());
				x*= nfac;
				if (Double.isNaN(x))
					System.currentTimeMillis();
				trptExprHash.put(trpts[i].getTranscriptID(), x);
			}
			
			
			// normalize transcript profile
			for (int i = 0; i < trpts.length; i++) {
				int tlen= trpts[i].getExonicLength();
				UniversalMatrix m= profile.getMatrix(tlen);
				double f= m.getNfactor(0.2d);
				double x= trptExprHash.get(trpts[i].getTranscriptID());
				x*= f;
				if (Double.isNaN(x))
					System.currentTimeMillis();
				trptExprHash.put(trpts[i].getTranscriptID(), x);
			}
			
			// normalize locus??
			
			return trptExprHash;
		}

		private void addConstraintToLp(int[] idx, double[] val, int eq, double cap) {
				
		//		if (writeFile)
		//			writeRowLP(idx, val, eq, cap);
		//		else if (columnWise) {
		//			;
		//		} else 
					//double[] constraint= createArray(idx, val);
					try {
						getLPsolve().addConstraintex(idx.length, val, idx, eq, cap);					
						//getLPsolve().addConstraint(constraint, eq, cap);
					} catch (LpSolveException e) {
						e.printStackTrace();
					}
		
			}

		private void getConstraintsNew(Edge e, long[] sig, IntVector v,
				HashMap<Edge, IntVector> mapE, boolean sense, boolean count) {
			
			// for the edge itself
			IntVector w= mapE.get(e);
			if (count) 
				++constraintCtr;
			else {
				v.add(++constraintCtr);	// for transcript fraction
				if (w== null)
					w= new IntVector();
				w.add(constraintCtr);	// edge consistency
			}
			mapE.put(e, w);
		
			// iterate super-edges
			for (int j = 0; e.getSuperEdges()!= null&& j < e.getSuperEdges().size(); j++) {
				SuperEdge se= e.getSuperEdges().elementAt(j);
				if (Graph.isNull(Graph.intersect(se.getTranscripts(), sig)))
					continue;
				// sense/anti-sense.. e must be first/last in super-edge
				if ((sense&& se.getEdges()[0]!= e)|| ((!sense)&& se.getEdges()[se.getEdges().length- 1]!= e))
					continue;
				if (count)
					++constraintCtr;
				else
					v.add(++constraintCtr);	// for transcript fraction
				w= mapE.get(se);
				if (!count) {
					if (w== null)
						w= new IntVector();
					w.add(constraintCtr); // for edge consistency
				}
				mapE.put(se, w);
				
				if (se.isPend())
					continue;	// no super-edges
				
				for (int k = 0; se.getSuperEdges()!= null&& k < se.getSuperEdges().size(); k++) {
					SuperEdge se2= se.getSuperEdges().elementAt(k);
					assert(se2.isPend());
					if (Graph.isNull(Graph.intersect(se2.getTranscripts(), sig)))
						continue;
					// sense/anti-sense.. e must be first/last in super-edge
					if ((sense&& se2.getEdges()[0]!= se)|| ((!sense)&& se2.getEdges()[se2.getEdges().length- 1]!= se))
						continue;
					if (count)
						++constraintCtr;
					else 
						v.add(++constraintCtr);	// tx
					w= mapE.get(se2);
					if (!count) {
						if (w== null)
							w= new IntVector();
						w.add(constraintCtr); // for edge consistency
					}
					mapE.put(se2, w);
				}
			}
			
			
		}

		private double nFactor= Double.NaN;
		public double getNFactor() {
			if (Double.isNaN(nFactor)|| true) {
				double fictReads= 0;
				Iterator iter= trptExprHash.keySet().iterator();
				while (iter.hasNext()) {
					Object o= iter.next();
					if (!(o instanceof String))
						continue;
					fictReads+= trptExprHash.get(o);
				}
		
				if (fictReads< nrMappingsLocusMapped)
					++nrLociUnderPredicted;
				else
					++nrLociOverPredicted;
				
				// can happen, when reads are where none expected
//				if (fictReads== 0^ nrMappingsObs== 0)
//					System.currentTimeMillis();
				
				// avoid large scaling; was 0, avoid NaN; 0.5 too large
				if (fictReads> 0.000001) {	
					nFactor= nrMappingsLocusMapped/ fictReads;
				} else 
					nFactor= 1d;
			}	
			return nFactor;
		}

		public HashMap<String, Integer> setConstraintsNew(Graph g, boolean count) {
					
					//boolean debug= false;
			//		if (g.trpts[0].getTranscriptID().equals("ENST00000323441")) {
			//			debug= true;
			//		}
					
					// transcript constraint variables
					Transcript[] trpts= g.trpts;
					IntVector v= null, w= null;
					DoubleVector u= null;
					HashMap<String, Integer> tMap= null;
					HashMap<String, IntVector> tDeltas= null;
					if (count)
						constraintCtr+= trpts.length;
					else { 
						try {
							getLPsolve().setLpName(trpts[0].getTranscriptID());
						} catch (LpSolveException e1) {
							e1.printStackTrace();
						}			
						tMap= new HashMap<String, Integer>(trpts.length* 2);
						tDeltas= new HashMap<String, IntVector>(trpts.length* 2);
						for (int i = 0; i < trpts.length; i++) { 
							tMap.put(trpts[i].getTranscriptID(), ++constraintCtr);
							if (debug)
								System.out.println("C"+constraintCtr+"\t"+trpts[i].getTranscriptID());
						}
						v= new IntVector();	// indices for transcript/part 
						w= new IntVector();	// indices for cost function
						u= new DoubleVector();	// observation, bases for cost function
					}
					
					// iterate edges
					Edge[] edges= g.getExonicEdgesInGenomicOrder();
					if (!count) {
						mapCCheck= new HashMap<String, Double>();
						for (int i = 0; i < trpts.length; i++) 
							mapCCheck.put(trpts[i].getTranscriptID(), 0d);
					}
					for (int i = 0; i < edges.length; i++) {
						Edge e= edges[i];
						if (!e.isExonic())
							continue;
						
						// the base edge
						int basElen= e.length();
						Transcript[] tt= g.decodeTset(e.getTranscripts());
						HashMap<Edge, IntVector> mapE= new HashMap<Edge, IntVector>();			
						// sense/anti
						for (int sa= 0; sa< 2; ++sa) {
							for (int x = 0; x < tt.length; ++x) {
								if (!count)
									v.removeAll();
								long[] sig= g.encodeTset(tt[x]);
		//						if (e.toString().equals("67937779-67937952^")&& !count)
		//							System.currentTimeMillis();
								getConstraints(e, sig, v, mapE, sa== 0, count);
								
								// add transcript constraint
								if (count) {
									constraintCtr+= 2;
								} else {		
									IntVector tDeltaV= tDeltas.get(tt[x].getTranscriptID());
									if (tDeltaV== null) 
										tDeltaV= new IntVector(tt[x].getExons().length* 4);	// s, as, je 1 super-edge
									
									// indices
									int[] idx= new int[v.length+ 3]; // sum(contributions) - delta + delta - frac(tx)
									System.arraycopy(v.vector, 0, idx, 0, v.length);
									idx[idx.length- 3]= ++constraintCtr;
									tDeltaV.add(-constraintCtr);	
									// add to cost function here!
									idx[idx.length- 2]= ++constraintCtr;
									tDeltaV.add(constraintCtr);		// add to cost function here!
									tDeltas.put(tt[x].getTranscriptID(), tDeltaV);
									idx[idx.length- 1]= tMap.get(tt[x].getTranscriptID());
									
									// vals
									double[] val= new double[idx.length];
									Arrays.fill(val, 1d);
									val[val.length- 3]= -1d;
									// frac(tx)
									int tlen= tt[x].getExonicLength();
									UniversalMatrix m= profile.getMatrix(tlen);
									double f= m.getFrac(
												tt[x].getExonicPosition(e.getFrac(true)),
												tt[x].getExonicPosition(e.getFrac(false)),
												tlen,
												sa== 0?Constants.DIR_FORWARD:Constants.DIR_BACKWARD);
									//System.err.println(f);
									if (Double.isInfinite(f)|| Double.isNaN(f)) {
										System.err.println("infinite value");
										f= m.getFrac(
												tt[x].getExonicPosition(e.getFrac(true)),
												tt[x].getExonicPosition(e.getFrac(false)),
												tlen,
												sa== 0?Constants.DIR_FORWARD:Constants.DIR_BACKWARD);
									}
									mapCCheck.put(tt[x].getTranscriptID(),
											mapCCheck.get(tt[x].getTranscriptID())+ f);
									val[val.length- 1]= -f;
									
									// restriction
									if (debug&& !count) {
										StringBuilder sb= new StringBuilder(e.toString());
										sb.append(": ");
										for (int k = 0; k < idx.length; k++) {
											sb.append(val[k]>0?"+":"");
											sb.append(val[k]== 1?"":val[k]);
											sb.append("C");
											sb.append(idx[k]+" ");
										}
										sb.append("= 0");
										System.out.println(sb);
									}
									int[] idx2= new int[idx.length- 2];
									System.arraycopy(idx, 0, idx2, 0, idx2.length);
									double[] val2= new double[val.length- 2];
									System.arraycopy(val, 0, val2, 0, val2.length);
									addConstraintToLp(idx, val, LpSolve.EQ, 0);
									++restrNr;

									// upper boundary delta- (variables)
									if (debug&& !count) {
										StringBuilder sb= new StringBuilder("with");
										sb.append(": ");
										for (int k = 0; k < idx2.length; k++) {
											sb.append(val2[k]>0?"+":"");
											sb.append(val2[k]== 1?"":val2[k]);
											sb.append("C");
											sb.append(idx2[k]+" ");
										}
										sb.append(">= 0");
										System.out.println(sb);
									}
									addConstraintToLp(idx2, val2, LpSolve.GE, 0);
									++restrNr;
									
									// force junction coverage ?
								}
							}
							
							
							// add edge constraints
							Edge[] ee= new Edge[mapE.size()];
							mapE.keySet().toArray(ee);
							for (int j = 0; j < ee.length; j++) {
								Edge f= ee[j];
								
								// paired-end super-edge
								boolean paird= (f instanceof SuperEdge)&& ((SuperEdge) f).isPend();
								int nr= ((paird|| sa== 0)? f.getReadNr(): f.getRevReadNr());
								v= mapE.remove(f);
								if (!count) {
									int[] idx= v.toIntArray();
									double[] val= new double[idx.length];
									Arrays.fill(val, 1d);
									if (debug) {
										StringBuilder sb= new StringBuilder(f.toString());
										sb.append(": ");
										for (int k = 0; k < idx.length; k++) {
											sb.append(val[k]==1?"+C":"-C");
											sb.append(idx[k]+" ");
										}
										sb.append("= "+nr);
										System.out.println(sb);
									}
									addConstraintToLp(idx, val, LpSolve.EQ, nr);
									++restrNr;
								}
							}
						}
					} // end all edges
					
					if (count) {
						constraintCtr+= gene.getTranscripts().length* 2;	// plus/minus 
						return null;
					}
					
					// transcript balance
					for (int i = 0; i < gene.getTranscripts().length; i++) {
						StringBuilder sb= null;
						if (debug) 
							sb= new StringBuilder(gene.getTranscripts()[i].getTranscriptID()+ ": ");
						IntVector tDeltaV= tDeltas.get(gene.getTranscripts()[i].getTranscriptID());
						int[] idx= new int[tDeltaV.size()+ 2];
						double[] val= new double[idx.length];
						for (int j = 0; j < tDeltaV.size(); j++) {
							int x= tDeltaV.get(j);
							if (x< 0) {
								val[j]= -1d;
								idx[j]= -x;
							} else {
								val[j]= 1d;
								idx[j]= x;
							}
							if (debug) 
								sb.append(val[j]== 1? "+C"+idx[j]+" ": "-C"+ idx[j]+" ");
						}

						idx[idx.length- 2]=  ++constraintCtr;
						val[val.length- 2]= 1d;
						w.add(constraintCtr);
						idx[idx.length- 1]=  ++constraintCtr;
						val[val.length- 1]= -1d;
						w.add(constraintCtr);
						if (debug) {
							sb.append("+C"+(constraintCtr- 2)+" -C"+(constraintCtr- 1)+ " = 0");
							System.out.println(sb);
						}
						addConstraintToLp(idx, val, LpSolve.EQ, 0);
						++restrNr;						
					}
					
					// set objective function/costs
					double[] costs= new double[w.size()];
					StringBuilder sb= null;
					if (debug) {
						sb= new StringBuilder();
					}
					for (int i = 0; i < costs.length; i++) {
						//double x= u.get(i);	// <0 : read substracting
						//costs[i]= 1+ Math.exp(-Math.log(1+ x));	// neutralizes
						//costs[i]= 1d+ Math.pow(x+1d, -1d/2d);
						//costs[i]= 0;	// 
						//costs[i]= 1d+ Math.log(x+ 1);
						//costs[i]= 1d+ Math.exp(-x);
						//costs[i]= 100d/ (x+ 1);	// last
						
						double c= 1d;
						//c= med/ (x+ 1);
						//c= 10d/ (1+ Math.pow(x, 1));
						//c= 10d/ (1+ x);
						//c= x;
						costs[i]= c;
						
						if (debug) { 
							sb.append("C");
							sb.append(w.get(i));
							sb.append("\t");
							sb.append(StringUtils.fprint(c, 3));
							sb.append("\n");
						}
						//costs[i]= 1d+ Math.log(x+ 1);
						//costs[i]= 1d+ Math.sqrt(x+ 1);	// best
						//costs[i]= 1d+ Math.pow(x+ 1, 1d/ 3d);
						
						//costs[i]= (Math.log(x+ 1)/ (x+ 1));	// logdiv
						//costs[i]= 1d/ (1d+ Math.log(x+ 1d));	// divlog
					}
					if (debug) {
						System.out.println("COSTS:");
						System.out.println(sb);
					}
					int[] idx= w.toIntArray();
					
					double[] a= createArray(idx, costs);
					try {
						getLPsolve().setObjFn(a);
						getLPsolve().setMinim();
					} catch (LpSolveException e) {
						e.printStackTrace();
					}
					
					// TODO consistency check
					Object[] oo= mapCCheck.keySet().toArray();
					for (int i = 0; i < oo.length; i++) {
						double val= mapCCheck.get(oo[i]);
						//System.err.println("check "+ val);
						if (Math.abs(2d- val)> 0.2d)
							System.err.println("Fraction inconsistency "+ oo[i]+"\t"+val);
					}
					
					return tMap;
				}

		public HashMap<String, Integer> setConstraints(Graph g, boolean count) {
					
					//boolean debug= false;
			//		if (g.trpts[0].getTranscriptID().equals("ENST00000323441")) {
			//			debug= true;
			//		}
					
					// transcript constraint variables
					Transcript[] trpts= g.trpts;
					IntVector v= null, w= null;
					DoubleVector u= null;
					HashMap<String, Integer> tMap= null;
					if (count)
						constraintCtr+= trpts.length;
					else { 
						try {
							getLPsolve().setLpName(trpts[0].getTranscriptID());
						} catch (LpSolveException e1) {
							e1.printStackTrace();
						}			
						tMap= new HashMap<String, Integer>(trpts.length* 2);
						for (int i = 0; i < trpts.length; i++) { 
							tMap.put(trpts[i].getTranscriptID(), ++constraintCtr);
							if (debug)
								System.out.println("C"+constraintCtr+"\t"+trpts[i].getTranscriptID());
						}
						v= new IntVector();	// indices for transcript/part 
						w= new IntVector();	// indices for cost function
						u= new DoubleVector();	// observation, bases for cost function
					}
					
					// iterate edges
					Edge[] edges= g.getExonicEdgesInGenomicOrder();
					if (!count) {
						mapCCheck= new HashMap<String, Double>();
						for (int i = 0; i < trpts.length; i++) 
							mapCCheck.put(trpts[i].getTranscriptID(), 0d);
					}
					for (int i = 0; i < edges.length; i++) {
						Edge e= edges[i];
						if (!e.isExonic())
							continue;
						
						// the base edge
						int basElen= e.length();
						Transcript[] tt= g.decodeTset(e.getTranscripts());
						HashMap<Edge, IntVector> mapE= new HashMap<Edge, IntVector>();			
						// sense/anti
						for (int sa= 0; sa< 2; ++sa) {
							for (int x = 0; x < tt.length; ++x) {
								if (!count)
									v.removeAll();
								long[] sig= g.encodeTset(tt[x]);
		//						if (e.toString().equals("67937779-67937952^")&& !count)
		//							System.currentTimeMillis();
								getConstraints(e, sig, v, mapE, sa== 0, count);
								
								// add transcript constraint
								if (!count) {
									int[] idx= new int[v.length+ 1]; // obs parts+ tx frac
									System.arraycopy(v.vector, 0, idx, 0, v.length);
									idx[idx.length- 1]= tMap.get(tt[x].getTranscriptID());
									double[] val= new double[idx.length];
									Arrays.fill(val, 1d);
									int tlen= tt[x].getExonicLength();
									UniversalMatrix m= profile.getMatrix(tlen);
									double f= m.getFrac(
												tt[x].getExonicPosition(e.getFrac(true)),
												tt[x].getExonicPosition(e.getFrac(false)),
												tlen,
												sa== 0?Constants.DIR_FORWARD:Constants.DIR_BACKWARD);
									//System.err.println(f);
									if (Double.isInfinite(f)|| Double.isNaN(f)) {
										System.err.println("infinite value");
										f= m.getFrac(
												tt[x].getExonicPosition(e.getFrac(true)),
												tt[x].getExonicPosition(e.getFrac(false)),
												tlen,
												sa== 0?Constants.DIR_FORWARD:Constants.DIR_BACKWARD);
									}
									mapCCheck.put(tt[x].getTranscriptID(),
											mapCCheck.get(tt[x].getTranscriptID())+ f);
									val[val.length- 1]= -f;
									if (debug&& !count) {
										StringBuilder sb= new StringBuilder(e.toString());
										sb.append(": ");
										for (int k = 0; k < idx.length; k++) {
											sb.append(val[k]>0?"+":"");
											sb.append(val[k]%1==0?((int) val[k]):val[k]);
											sb.append("C");
											sb.append(idx[k]+" ");
										}
										sb.append("= 0");
										System.out.println(sb);
									}
									addConstraintToLp(idx, val, LpSolve.EQ, 0);
									++restrNr;
									
									// force junction coverage
									if (val.length== 3) {
										assert(val[0]<0&&val[1]>0&&val[2]>0);
										int fixedreadlen= 27;
										double frac= 0.9;
										val[0]= 0d;
										val[1]= frac/ (basElen- fixedreadlen+ 1);
										val[2]= -1d/ (fixedreadlen- 1);
										addConstraintToLp(idx, val, LpSolve.LE, 0);
										++restrNr;
									}							
								}
							}
							// add edge constraints
							Edge[] ee= new Edge[mapE.size()];
							mapE.keySet().toArray(ee);
							
							// total obs
			/*				int sumObs= 0;
							for (int j = 0; j < ee.length; j++) {
								Edge f= ee[j];
								boolean paird= (f instanceof SuperEdge)&& ((SuperEdge) f).isPend();
								int nr= ((paird|| sa== 0)? f.getReadNr(): f.getRevReadNr());
								sumObs+= nr;
							}
			*/				
							for (int j = 0; j < ee.length; j++) {
								Edge f= ee[j];
								
								// paired-end super-edge
								boolean paird= (f instanceof SuperEdge)&& ((SuperEdge) f).isPend();
								int nr= ((paird|| sa== 0)? f.getReadNr(): f.getRevReadNr());
								// coverage or cost calc
		/*						int len= (f instanceof SuperEdge)? 27- 1: f.length()- 27+ 1;
								double cov= 0;
								if (len<= 0) {
									assert(nr== 0);
								} else
									cov= nr/ (double) len;
								double cost= 1;
								if (len> 0) {
									cost= 1d+ (10* Math.exp(-nr/ 10));
									//cost= 1d+ (10* Math.pow((nr+ 1)/10d, -4));
								}						
								if (f instanceof SuperEdge) {
									cost= 1d+ (10* Math.exp(-nr/10d));
								} else {
									if (len> 0)
										cost= 1d+ (Math.pow(cov+ 1, -1));
								}
								if (Double.isNaN(cost)|| Double.isInfinite(cost))
									System.currentTimeMillis();
		*/							
								v= mapE.remove(f);
								if (count)
									constraintCtr+= 2;
								else {
									int[] idx= new int[v.length+ 2];	// +/-
									System.arraycopy(v.vector, 0, idx, 0, v.length);
									int c= ++constraintCtr;
									// plus on not paired edges at 0-cost, it substracts from obs
									// usual edges cost when in single read mode
									if ((!pairedEnd)|| (pairedEnd&& paird)) {
										w.add(c);
		//								if (nr== 0)
		//									u.add(Double.NEGATIVE_INFINITY);
		//								else
											u.add(-nr);
									}
									idx[idx.length- 2]= c;
									// plus has to be limited, it substracts
									int lim= (paird||(!pairedEnd))? Math.max(nr- 1, 0): nr; // Math.max(nr- 1, 0): nr; 
									try {
										getLPsolve().setUpbo(constraintCtr, lim);
									} catch (LpSolveException e1) {
										e1.printStackTrace();
									}
										
									c= ++constraintCtr;
									// do not limit adding, even with f= 100 unsolvable systems
									// adding reads always costs, also on single edges
									w.add(c);
									u.add(nr);
									idx[idx.length- 1]= c;
									double[] val= new double[idx.length];
									Arrays.fill(val, 1d);
									val[val.length- 1]= -1d;
									if (debug&& !count) {
										StringBuilder sb= new StringBuilder(f.toString());
										sb.append(": ");
										for (int k = 0; k < idx.length; k++) {
											sb.append(val[k]==1?"+C":"-C");
											sb.append(idx[k]+" ");
										}
										sb.append("= "+nr);
										System.out.println(sb);
									}
									addConstraintToLp(idx, val, LpSolve.EQ, nr);
									++restrNr;
								}
							}
						}
					} // end all edges
					
					if (count) {
						if (stopperBalance|| stopperBalance2)
							constraintCtr+= 2;
						return null;
					}
					
					if (stopperBalance2) {
						int[] idx= new int[tMap.size()+ 2];
						Iterator<Integer> iter= tMap.values().iterator();
						int ctr= 0;
						while (iter.hasNext()) {
							idx[ctr++]= iter.next();
						}
						idx[ctr++]= ++constraintCtr;
						idx[ctr++]= ++constraintCtr;
						double[] vals= new double[idx.length];
						Arrays.fill(vals, 1d);
						vals[vals.length- 1]= -1d;
						addConstraintToLp(idx, vals, LpSolve.EQ, nrMappingsLocusMapped);
						++restrNr;
					}
					
					// set objective function/costs
			//		double[] a= createArray(w.toIntArray());	// linear costs
					//double min= Math.exp(-nrMappingsLocusMapped);
					int ext= 0;
					if (stopperBalance)
						ext= 2;
					double[] costs= new double[w.size()+ ext];
					double[] pmPfx= new double[w.size()+ ext];
					double[] dd= u.toDoubleArray();
					for (int i = 0; i < dd.length; i++) {
						dd[i]= Math.abs(dd[i]);
					}
					Distribution dist= new Distribution(dd);
					double med= dist.getMedian();
					if (med== 0)
						med= 1;
					StringBuilder sb= null;
					if (debug) {
						sb= new StringBuilder();
					}
					for (int i = 0; i < costs.length- ext; i++) {
						double x= u.get(i);	// <0 : read substracting
						pmPfx[i]= ((x< 0|| x== Double.NEGATIVE_INFINITY)?-1d:1d);
						//costs[i]= 1+ Math.exp(-Math.log(1+ x));	// neutralizes
						//costs[i]= 1d+ Math.pow(x+1d, -1d/2d);
						//costs[i]= 0;	// 
						//costs[i]= 1d+ Math.log(x+ 1);
						//costs[i]= 1d+ Math.exp(-x);
						//costs[i]= 100d/ (x+ 1);	// last
						
						double c= 1d;
						if (!pairedEnd) {
							double y= (x== Double.NEGATIVE_INFINITY)?0:x;
							c= (10d* Math.exp(-(Math.abs(y)/ med)));
							if (x< 0/*|| x== Double.NEGATIVE_INFINITY*/)
								c+= 1;
						}
						//c= med/ (x+ 1);
						//c= 10d/ (1+ Math.pow(x, 1));
						//c= 10d/ (1+ x);
						//c= x;
						if (Double.isInfinite(c)|| Double.isNaN(c))
							System.currentTimeMillis();
						costs[i]= c;
						
						
						if (debug) { 
							sb.append("C");
							sb.append(w.get(i));
							sb.append("\t");
							sb.append(StringUtils.fprint(c, 3));
							sb.append("\n");
						}
						//costs[i]= 1d+ Math.log(x+ 1);
						//costs[i]= 1d+ Math.sqrt(x+ 1);	// best
						//costs[i]= 1d+ Math.pow(x+ 1, 1d/ 3d);
						
						//costs[i]= (Math.log(x+ 1)/ (x+ 1));	// logdiv
						//costs[i]= 1d/ (1d+ Math.log(x+ 1d));	// divlog
					}
					if (debug) {
						System.out.println("COSTS:");
						System.out.println(sb);
					}
					int[] idx= null;
					if (stopperBalance) {
						idx= new int[w.length+ ext];
						pmPfx[pmPfx.length- 2]= 1;
						pmPfx[pmPfx.length- 1]= -1;
						System.arraycopy(w.vector, 0, idx, 0, w.length);
						idx[pmPfx.length- 2]= ++constraintCtr;
						idx[pmPfx.length- 1]= ++constraintCtr;
						addConstraintToLp(idx, pmPfx, LpSolve.EQ, 0);
						++restrNr;
						
						costs[costs.length- 2]= 1;
						costs[costs.length- 1]= 1;
					} else {
						if (stopperBalance2) {
							idx= new int[w.length+ 2];
							double[] costs2= new double[costs.length+ 2];
							System.arraycopy(w.vector, 0, idx, 0, w.length);
							idx[idx.length- 2]= constraintCtr- 2;
							idx[idx.length- 1]= constraintCtr- 1;
							System.arraycopy(costs, 0, costs2, 0, costs.length);
							costs2[costs2.length- 2]= 1d;
							costs2[costs2.length- 1]= 1d;
							costs= costs2;
						} else
							idx= w.toIntArray();
					}
					
					double[] a= createArray(idx, costs);
					try {
						getLPsolve().setObjFn(a);
						getLPsolve().setMinim();
					} catch (LpSolveException e) {
						e.printStackTrace();
					}
					
					// TODO consistency check
					Object[] oo= mapCCheck.keySet().toArray();
					for (int i = 0; i < oo.length; i++) {
						double val= mapCCheck.get(oo[i]);
						//System.err.println("check "+ val);
						if (Math.abs(2d- val)> 0.2d)
							System.err.println("Fraction inconsistency "+ oo[i]+"\t"+val);
					}
					
					return tMap;
				}

		private void getConstraints(Edge e, long[] sig, IntVector v,
				HashMap<Edge, IntVector> mapE, boolean sense, boolean count) {
			
			// for the base-edge itself
			IntVector w= mapE.get(e);
			if (count) 
				++constraintCtr;
			else {
				v.add(++constraintCtr);	// for transcript fraction
				if (w== null)
					w= new IntVector();
				w.add(constraintCtr);	// edge consistency
			}
			mapE.put(e, w);
		
			// iterate super-edges
			for (int j = 0; e.getSuperEdges()!= null&& j < e.getSuperEdges().size(); j++) {
				SuperEdge se= e.getSuperEdges().elementAt(j);
				if (Graph.isNull(Graph.intersect(se.getTranscripts(), sig)))
					continue;
				// sense/anti-sense.. e must be first/last in super-edge
				if ((sense&& se.getEdges()[0]!= e)|| ((!sense)&& se.getEdges()[se.getEdges().length- 1]!= e))
					continue;
				if (count)
					++constraintCtr;
				else
					v.add(++constraintCtr);	// for transcript fraction
				w= mapE.get(se);
				if (!count) {
					if (w== null)
						w= new IntVector();
					w.add(constraintCtr); // for edge consistency
				}
				mapE.put(se, w);
				
				if (se.isPend())
					continue;	// no super-edges
				
				for (int k = 0; se.getSuperEdges()!= null&& k < se.getSuperEdges().size(); k++) {
					SuperEdge se2= se.getSuperEdges().elementAt(k);
					assert(se2.isPend());
					if (Graph.isNull(Graph.intersect(se2.getTranscripts(), sig)))
						continue;
					// sense/anti-sense.. e must be first/last in super-edge
					if ((sense&& se2.getEdges()[0]!= se)|| ((!sense)&& se2.getEdges()[se2.getEdges().length- 1]!= se))
						continue;
					if (count)
						++constraintCtr;
					else 
						v.add(++constraintCtr);	// tx
					w= mapE.get(se2);
					if (!count) {
						if (w== null)
							w= new IntVector();
						w.add(constraintCtr); // for edge consistency
					}
					mapE.put(se2, w);
				}
			}
			
			
		}

		int[] cBalTrpt= null;
		int[] boundIdx, decIdx;
		double[] boundVal, decVal;
		int[] costIdx;
		double[] costVal;
		private void printConstraint(PrintStream p, String id, int length, double[] val, int[] idx,
				int eq, double d) {
			StringBuilder sb= new StringBuilder(id);
			sb.append("\t");
			for (int i = 0; i < length; i++) {
				sb.append(val[i]>= 0?"+":"-");
				sb.append(Math.abs(val[i])== 1?"": Float.toString((float) Math.abs(val[i])));
				sb.append("C");
				sb.append(idx[i]);
				sb.append(" ");
			}
			sb.append(eq== LpSolve.EQ?"=":(eq== LpSolve.GE?">=":"<="));
			sb.append(" ");
			sb.append(d);
			p.println(sb.toString());
		}

		int mapLenMin= 27, mapLenMax= 27; // should be minReadLen
		IntVector idxDecV; 
		DoubleVector valDecV;
		private int setCDeconvolution_deltaMulti(boolean count, Graph g, Edge f, Edge e, boolean sense, IntVector[] cSumV, IntVector[] cBalV,
				IntVector costIdxV, DoubleVector costValV, IntVector[][] covIdxV, DoubleVector[][] covValV) throws LpSolveException {
			
			int constraintSave= constraintCtr;
			
			// attributes of e
			int obs= (sense? e.getReadNr(): e.getRevReadNr());
			boolean sj= (e instanceof SuperEdge)&& (!((SuperEdge) e).isPend());
			double effLenMin= calcEffLen(e, f, true);
			double effLenMax= calcEffLen(e, f, false);
			assert(effLenMin> 0|| obs== 0);
			if (!count) {
				if (idxDecV== null) {
					idxDecV= new IntVector(100);
					valDecV= new DoubleVector(100);
				} else {
					idxDecV.removeAll();
					valDecV.removeAll();
				}
					
			}
			
			// all transcript contributions
			for (int j = 0; j < g.trpts.length; j++) {

				if (!g.contains(e.getTranscripts(), j))
					continue;
				if (count) {
					++constraintCtr;
					++restrNr;
					continue;
				}
				
				//int x= ctrTT* 3;	// contribution, d-, d+
				idxDecV.add(++constraintCtr);
				valDecV.add(1d);
				cSumV[j].add(constraintCtr);	// sum of (e_j,t_i)
				if ((sj &&effLenMin> 0)|| ((!sj)&& effLenMax> 0)) {
					covIdxV[0][j].add(constraintCtr);
					covValV[0][j].add(1d/ (sj? -effLenMin: effLenMax));	// inverse
				}
				if ((sj &&effLenMax> 0)|| ((!sj)&& effLenMin> 0)) {
					covIdxV[1][j].add(constraintCtr);
					covValV[1][j].add(1d/ (sj? -effLenMax: effLenMin));	// inverse
				}
				
			}
			
			// obs: deconvolution
			int partSize= constraintCtr- constraintSave; 
			setConstraints(count, obs, idxDecV, valDecV, costIdxV, costValV);

			if (count) {
				++restrNr;
				return partSize;
			}
			int[] decIdx= idxDecV.toIntArray();
			double[] decVal= valDecV.toDoubleArray();
			getLPsolve().addConstraintex(decIdx.length, decVal, decIdx, LpSolve.EQ, obs);
			if (debug)
				printConstraint(System.out, e+ " ("+(sense?"s":"a")+", dec)", decIdx.length, decVal, decIdx, LpSolve.EQ, obs); 
			++restrNr;
			
			return partSize;

		}
		
		private double calcEffLen(Edge e, Edge f, boolean min) {
			
			double effLen= -1d;
			if (e instanceof SuperEdge) {
				SuperEdge se= (SuperEdge) e;
				Edge[] see= se.getEdges();
				if (!se.isPend()) {	// sj
					int interNt= 1;
					for (int i = 1; i < se.getEdges().length- 1; i++) 
						interNt+= se.getEdges()[i].length();
					effLen= (min? 
								Math.min(Math.min(mapLenMin- interNt, see[0].length()),
									see[see.length- 1].length()):
								Math.min(Math.min(mapLenMax- interNt, see[0].length())
									, see[see.length- 1].length()));
				} else {
					//effLenMin= effLenMax= f.length(); // not: - fixedReadLen+ 1;
					int base= f.length();
					effLen= (min? base- mapLenMax+ 1: base- mapLenMin+ 1);
				}
			} else {
				// 0effLenMin= effLenMax= e.length(); // not: - fixedReadLen+ 1;
				int base= e.length();
				effLen= (min? base- mapLenMax+ 1: base- mapLenMin+ 1);
			}

			return effLen;
		}

		private int setConstraints(boolean count, int obs, IntVector idxDecV, DoubleVector valDecV, 
				IntVector costIdxVec, DoubleVector costValVec) throws LpSolveException {
			
			int constraintSave= constraintCtr, maxMinus= obs+ Math.max(obs- 1, 0);
//			long fac= FactorialPrimeSchoenhage.factorial(obs).intValue();
			int lastPlus= -1, lastMinus= -1;
			for (int i = obs+ 1;; i++) {
//				double exp= exp(-i);
//				double pow= pow(i, obs);
//				double cost= fac/(pow* exp);
//				if (cost< 0)
//					System.currentTimeMillis();
				int delta= i- obs;
				double obs1= Math.max(1, obs);
				double cost= exp(-delta/ obs);
				if (cost> 1000&& i> obs+ 1)
					break;
				// d-, adding reads to obs				
				if (count)
					++constraintCtr;
				else {
					getLPsolve().setUpbo(++constraintCtr, 1d);
					idxDecV.add(constraintCtr);
					valDecV.add(-1d);
					costIdxVec.add(constraintCtr);
					costValVec.add(cost);
					lastPlus= constraintCtr;
				}
				
				// d+, substracting reads from obs
				if(i<= maxMinus) {
					if (count) 
						++constraintCtr;
					else {
						getLPsolve().setUpbo(++constraintCtr, 1d);
						idxDecV.add(constraintCtr);
						valDecV.add(1d);
						costIdxVec.add(constraintCtr);
						costValVec.add(cost);
						lastMinus= constraintCtr;
					}
				}
			}
			if (!count)
				getLPsolve().setUpbo(lastPlus, getLPsolve().getInfinite());
			
			return (constraintCtr- constraintSave);
		}
		
		private double e= Math.E;
		private double exp(int i) {
			
			if (i== 0)
				return e;
			int j= i>= 0? i: -i;
			double val= 1d;
			for (int x = 0; x < j; x++) 
				val*= e; 
			if (i< 0)
				val= 1d/ val;
			
			return val;
		}

		private long pow(int base, int exp) {
			assert(exp>= 0);
			if (exp== 0)
				return 1;
			if (exp== 1)
				return base;
			long val= 1l;
			for (int i = 0; i < exp; i++) 
				val*= base;
			return val;
		}

		private double getPoissonProbability(int exp, int obs) {
			// (exp^obs* e^(-exp))/ fact(obs)
			int fac= 0; //FactorialPrimeSchoenhage.factorial(obs).intValue();
			double p= Math.pow(exp, obs)* Math.exp(-exp)/ fac;
			
			return p;
		}
		
		


		private double getCostObs(int obs) {
			double c= Math.exp(-obs);
			return c;
		}
		
		private double getCostObs(int obs, double len) {
			if (len<= 0)
				return 1; 
			double cov= obs/ (double) len;
			double c= Math.exp(-cov);
			//c*= 100;
			return c;
		}
		
		private double getCostObs(boolean sj, int obs, double len) {
			if (len<= 0)
				return 2;
			double c= 1d+ Math.exp(-obs/ (sj?medJunc:medExon));
			//double c= Double.MAX_VALUE;	// force infeasable
/*			double c= 0;
			if (sj)
				c= Math.log(obs+ 1d);
			else 
				c= Math.max(0, -Math.log((obs+ 1)/ (double) len));
			if (c< 0.000001&& c> 0)
				System.currentTimeMillis();
*/
			
			return c;
		}

		private HashMap<Object, Double> getResult(HashMap<String, Integer> tMap) {
			
			result= new double[1+ restrNr+ constraintCtr];
			try {
				getLPsolve().getPrimalSolution(result);
			} catch (LpSolveException e1) {
				e1.printStackTrace();
			}
			//valObjFunc= result[0];	
			HashMap<Object,Double> trptExprHash= new HashMap<Object,Double>();
			
			// normalize fraction rounding errors
			Transcript[] trpts= gene.getTranscripts();
			double sum= 0;
			for (int i = 0; i < trpts.length; i++) { 
				int c= tMap.get(trpts[i].getTranscriptID());
				if (Double.isNaN(c))
					System.currentTimeMillis();
				double x= result[restrNr+ c];
				double tot= mapCCheck.get(trpts[i].getTranscriptID());
				tot/= 2d;
				x/= tot;
				sum+= x;
				if (Double.isNaN(x))
					System.currentTimeMillis();
				trptExprHash.put(trpts[i].getTranscriptID(), x);
			}
			if (sum== 0)
				return trptExprHash;
			
			// normalize flow network over-/underprediction
			double nfac= nrMappingsLocusMapped/ (double) sum;
			for (int i = 0; sum!= 0&& i < trpts.length; i++) {
				double x= trptExprHash.get(trpts[i].getTranscriptID());
				x*= nfac;
				if (Double.isNaN(x))
					System.currentTimeMillis();
				trptExprHash.put(trpts[i].getTranscriptID(), x);
			}
			
			
			// normalize transcript profile
			for (int i = 0; i < trpts.length; i++) {
				int tlen= trpts[i].getExonicLength();
				UniversalMatrix m= profile.getMatrix(tlen);
				double f= m.getNfactor(0.2d);
				double x= trptExprHash.get(trpts[i].getTranscriptID());
				x*= f;
				if (Double.isNaN(x))
					System.currentTimeMillis();
				trptExprHash.put(trpts[i].getTranscriptID(), x);
			}
			
			// normalize locus??
			
			return trptExprHash;
		}

		private int setCDeconvolution_pgauss(byte mode, Graph g, Edge f, Edge e, boolean sense, IntVector[] cSumV, IntVector[] cBalV,
						IntVector costIdxV, DoubleVector costValV, IntVector[][] covIdxV, DoubleVector[][] covValV) throws LpSolveException {
					
					// obs: bounds
					int ctrTT= 0;
					IntVector minusIdxV= null;
					if (mode!= MODE_COUNT) 
						minusIdxV= new IntVector(g.trpts.length);	// TODO reuse
					
					// attributes of e
					int obs= (sense? e.getReadNr(): e.getRevReadNr());
					boolean sj= (e instanceof SuperEdge)&& (!((SuperEdge) e).isPend());
					double effLenMin= calcEffLen(e, f, true);
					double effLenMax= calcEffLen(e, f, false);
					assert(effLenMin> 0|| obs== 0);
					double cost= getCostObs(sj, obs, effLenMax);

					// iterate transcripts
					for (int j = 0; j < g.trpts.length; j++) {
						
						if (!g.contains(e.getTranscripts(), j))
							continue;
						if (mode== MODE_COUNT) {
							++ctrTT;
							++restrNr;
							continue;
						}
						
						int x= ctrTT* 3;
						decIdx[x]= ++constraintCtr;		// c
						cSumV[j].add(constraintCtr);	// Contribution sum for expectation
						//if (se!= null) ???
						boundIdx[0]= constraintCtr;		// delta<= c
						//boundVal[0]= 1d;	// TODO init once
						if ((sj &&effLenMin> 0)|| ((!sj)&& effLenMax> 0)) {
							covIdxV[0][j].add(constraintCtr);
							covValV[0][j].add(1d/ (sj? -effLenMin: effLenMax));	// inverse
						}
						if ((sj &&effLenMax> 0)|| ((!sj)&& effLenMin> 0)) {
							covIdxV[1][j].add(constraintCtr);
							covValV[1][j].add(1d/ (sj? -effLenMax: effLenMin));	// inverse
						}
						
						// d-, adding reads to obs
						decIdx[x+ 1]= ++constraintCtr;
						if (effLenMin<= 0&& (!sj)) {
							boolean found= false;
							for (int k = 0; (!found)&& e.getSuperEdges()!= null&& k < e.getSuperEdges().size(); k++) {
								SuperEdge se2 = (SuperEdge) e.getSuperEdges().elementAt(k);
								if (se2.isPend()|| (sense?se2.getEdges()[0]!= e: se2.getEdges()[se2.getEdges().length- 1]!= e))
									continue;
								if (g.contains(se2.getTranscripts(), j)) {
									found= true;
									break;
								}
							}
							if (found) {	// ensure super-edge before constraining
								getLPsolve().setUpbo(constraintCtr, 0);
								if (false&& mode== MODE_COMPLETE)
									System.out.println(e+ " ("+(sense?"s, ":"a, ")+"ub d+)\tC"+constraintCtr+ " = 0");
							}
						}
						cBalV[j].add(-constraintCtr);
						boundIdx[1]= constraintCtr;
						if (effLenMax> 0) {	// costs only iff there could be mappings
//							costIdxV.add(constraintCtr);
//							costValV.add(cost);
						}
						if (mode== MODE_CONSTRAINT) {
							getLPsolve().addConstraintex(2, boundVal, boundIdx, LpSolve.GE, 0);
							++restrNr;
						} if (false&& mode== MODE_COMPLETE) {
							printConstraint(System.out, e+ " ("+(sense?"s, ":"a, ")+"ub d+)", 2, boundVal, boundIdx, LpSolve.GE, 0);
						}
		
						
						// d+, substracting reads from obs 
						// bounded by c of the tx
						decIdx[x+ 2]= ++constraintCtr;	// d+
						if ((mode!= MODE_COUNT)) 
							minusIdxV.add(constraintCtr);	// last man standing				
						cBalV[j].add(constraintCtr);
						if (effLenMax> 0) {
							cost= getCostObs(sj, obs, effLenMax);
//							costIdxV.add(constraintCtr);
//							costValV.add(cost);
						}
						
						++ctrTT;
					}
					
					// set costs
					int saveCC= constraintCtr;
					setCostsPseudoGaussian(mode, sj, obs, effLenMax, ctrTT* 3, decIdx, costIdxV, costValV);
					int costCC= constraintCtr- saveCC;
					
					// obs: deconvolution
					if (mode!= MODE_COUNT) {
						// limit last man standing
/*						double bound= Math.max(obs/ 100d, 0.01d);
						int[] idx= minusIdxV.toIntArray();
						double[] val= new double[idx.length];
						Arrays.fill(val, 1d);
						if (mode== MODE_CONSTRAINT) {
							getLPsolve().addConstraintex(idx.length, val, idx, LpSolve.LE, bound);
							++restrNr;
						} else if (false&& mode== MODE_COMPLETE) {
							printConstraint(System.out, e+ " ("+(sense?"s, ":"a, ")+"ub d-)", idx.length, val, idx, LpSolve.LE, bound);
						}
*/						
						// deconvolution
						if (mode== MODE_CONSTRAINT) {
							getLPsolve().addConstraintex(ctrTT* 3, decVal, decIdx, LpSolve.EQ, obs);
							++restrNr;
						} else if (mode== MODE_COMPLETE)
							printConstraint(System.out, e+ " ("+(sense?"s, ":"a, ")+"dec)", ctrTT* 3, decVal, decIdx, LpSolve.EQ, obs);
							//completeConstraint(System.out, e+ " ("+(sense?"s, ":"a, ")+"dec)", ctrTT* 3, decVal, decIdx, LpSolve.EQ, obs); 
					}
					
					return (ctrTT* 3)+ costCC;
				}

		private void setCostsGaussian(byte mode, boolean sj, int obs, double effLen, int length, int[] decIdx,
				IntVector costIdxV, DoubleVector costValV) throws LpSolveException {
			
			
		}
		
		private double[] gauss_costs= new double[] 
		                                         //{2.5d, 4d, 20d, 100d};
		                                         {1, 10, 100, 1000, 10000, 100000};
		private void setCostsPseudoGaussian(byte mode, boolean sj, int obs, double effLen, int length, int[] decIdx,
				IntVector costIdxV, DoubleVector costValV) throws LpSolveException {
			
//			if (1== 1)
//				return;
			
			int base= Math.max(obs, 1);
			//double baseCov= base/ (double) Math.max(effLen, 1d);
			//double gaussFac= (1d/ Math.sqrt(2* Math.PI* baseCov));
			double fac= //Math.max(effLen, 1d)/ obs;
				1d+ Math.exp(-obs/ (sj?medJunc:medExon));
			
			int indices= (length/ 3);
			int[] idx= null;
			double[] val= null;
			if (mode!= MODE_COUNT) {
				idx= new int[indices+ gauss_costs.length];
				for (int i = 0; i < indices; i++) 
					idx[i]= decIdx[(i* 3)+ 1];			
				val= new double[idx.length];
				Arrays.fill(val, 0, indices, 1d);
				Arrays.fill(val, indices, val.length, -1d);
			}
			
			// costs d-, adding to obs
			//double deltaObs= base;
			double deltaObs= base/ (double) gauss_costs.length;
			for (int i = 0; i < gauss_costs.length; ++i) {
				if (mode== MODE_COUNT) {
					++constraintCtr;
					continue;
				}
				idx[indices+ i]= ++constraintCtr;
				costIdxV.add(constraintCtr);				
				if (i< gauss_costs.length- 1)
					getLPsolve().setUpbo(constraintCtr, deltaObs);
				costValV.add(fac* gauss_costs[i]);
			}
			if (mode== MODE_CONSTRAINT) {
				getLPsolve().addConstraintex(idx.length, val, idx, LpSolve.EQ, 0);
				++restrNr;
			} else if (mode== MODE_COMPLETE&& false) {
				printConstraint(System.out, "cost", idx.length, val, idx, LpSolve.EQ, 0d);
			}
			
			// costs d+, substracting from obs
			for (int i = 0; mode!= MODE_COUNT&& i < indices; i++) 
				idx[i]= decIdx[(i* 3)+ 2];			
			deltaObs= base/ (double) gauss_costs.length;
			for (int i = 0; i < gauss_costs.length; ++i) {
				if (mode== MODE_COUNT) {
					++constraintCtr;
					continue;
				}
				idx[indices+ i]= ++constraintCtr;
				costIdxV.add(constraintCtr);				
				getLPsolve().setUpbo(constraintCtr, deltaObs);
				costValV.add(fac* gauss_costs[i]);
			}			
			if (mode== MODE_CONSTRAINT) {
				getLPsolve().addConstraintex(idx.length, val, idx, LpSolve.EQ, 0);
				++restrNr;
			} else if (mode== MODE_COMPLETE&& false) {
				printConstraint(System.out, "cost", idx.length, val, idx, LpSolve.EQ, 0d);
			}
		}

		private int setCBase_save(boolean count, Graph g, Edge e, boolean sense, IntVector[] cSumV, IntVector[] cBalV,
				IntVector costIdxV, DoubleVector costValV, IntVector[][] covIdxV, DoubleVector[][] covValV) throws LpSolveException {
			
			assert(!(e instanceof SuperEdge));
			for (int i = 0; (!count)&& i < covIdxV.length; i++) {
				for (int j = 0; j < covValV[i].length; j++) {
					covIdxV[i][j]= new IntVector(2);
					covValV[i][j]= new DoubleVector(2);
				}
			}
			for (int i = 0; (!count)&& i < cSumV.length; i++) 
				cSumV[i].removeAll();
			
			int nrC= setCDeconvolution_save(count, g, null, e, sense, cSumV, cBalV, costIdxV, costValV, covIdxV, covValV);
			
			// iterate superedges
			for (int i = 0; e.getSuperEdges()!= null&& i < e.getSuperEdges().size(); ++i) {
				SuperEdge se= e.getSuperEdges().elementAt(i);
				
				// sense/anti-sense.. e must be first/last in super-edge
				if ((sense&& se.getEdges()[0]!= e)
						|| ((!sense)&& se.getEdges()[se.getEdges().length- 1]!= e))
					continue;
				
				nrC+= setCDeconvolution_save(count, g, e, se, sense, cSumV, cBalV, costIdxV, costValV, covIdxV, covValV);
				
				// TODO 2nd level for paired-end
			}
			if (count) {
				// NP: restrNr+= cSumV.length+ 2;	// < cSumV.length
				return nrC;
			}
			
			// expectation restriction
			for (int i = 0; i < cSumV.length; i++) {
				int len= cSumV[i].size();
				if (len== 0)
					continue;
				int[] idx= new int[len+ 1];	// TODO reuse
				idx[0]= (i+1);
				System.arraycopy(cSumV[i].vector, 0, idx, 1, cSumV[i].length);
				double[] val= new double[idx.length];
				Arrays.fill(val, 1d);
				Transcript tx= g.trpts[i];
				int tlen= tx.getExonicLength();
				UniversalMatrix m= profile.getMatrix(tlen);
				double f= m.getFrac(
							tx.getExonicPosition(e.getFrac(true)),
							tx.getExonicPosition(e.getFrac(false)),
							tlen,
							sense?Constants.DIR_FORWARD:Constants.DIR_BACKWARD);
				if (Double.isInfinite(f)|| Double.isNaN(f)) {
					System.err.println("infinite value");
					f= m.getFrac(
							tx.getExonicPosition(e.getFrac(true)),
							tx.getExonicPosition(e.getFrac(false)),
							tlen,
							sense?Constants.DIR_FORWARD:Constants.DIR_BACKWARD);
				}
				ccheck[i]+= f;
				val[0]= -f;
				
				getLPsolve().addConstraintex(idx.length, val, idx, LpSolve.EQ, 0);
				++restrNr;
				if (debug)
					printConstraint(System.out, e+ " ("+ (sense?"s":"a")+ ", "+ g.trpts[i].getTranscriptID()+ ")", 
							idx.length, val, idx, LpSolve.EQ, 0d); 
			}
			
			// flux preservation, 2 for min/max readlength
			for (int i = 0; (!count)&& i < covIdxV.length; i++) {
				for (int j = 0; j < covValV[i].length; j++) {	// per tx
					// assure there is an exon and a sj
					if (covIdxV[i][j].length< 2) 	
						continue;				
					boolean plus= false, minus= false;
					for (int k = 0; k < covValV[i][j].length; k++) {
						double d= covValV[i][j].get(k);
						plus|= (d> 0);
						minus|= (d< 0);
						if (plus&& minus)
							break;
					}
					if (!(plus&& minus))
						continue;
					getLPsolve().addConstraintex(covIdxV[i][j].length, covValV[i][j].toDoubleArray(), 
							covIdxV[i][j].toIntArray(), i== 0? LpSolve.LE: LpSolve.GE, 0);
					++restrNr;
					if (debug)
						printConstraint(System.out, e+ " (flux)", covIdxV[i][j].length, covValV[i][j].toDoubleArray(), 
								covIdxV[i][j].toIntArray(), i== 0? LpSolve.LE: LpSolve.GE, 0); 
				}
			}
			
			
			return nrC;
		}

		private int setCAll_save(boolean count, Graph g) throws LpSolveException {
			IntVector costIdxV= null;
			DoubleVector costValV= null;
			IntVector[] cSumV= null;
			IntVector[][] cBalV= null; // sense/asense, trpts[]
			IntVector[][] covIdxV= null;	// min/maxCov, trpts[]
			DoubleVector[][] covValV= null;
			Edge[] edges= g.getExonicEdgesInGenomicOrder();
		
			if (!count) {
				costIdxV= new IntVector(2* 2* g.trpts.length* edges.length);
				costValV= new DoubleVector(costIdxV.vector.length);
				cSumV= new IntVector[g.trpts.length];
				for (int i = 0; i < cSumV.length; i++) 
					cSumV[i]= new IntVector(5);	// se per tx in e
				cBalV= new IntVector[2][g.trpts.length];
				for (int i = 0; i < cBalV.length; i++) 
					for (int j = 0; j < cBalV[i].length; j++) 
						cBalV[i][j]= new IntVector(2* edges.length);	// 2 (+/-)* |E|
				covIdxV= new IntVector[2][g.trpts.length];
				covValV= new DoubleVector[2][g.trpts.length];
				for (int i = 0; i < covIdxV.length; i++) {
					for (int j = 0; j < g.trpts.length; j++) {
						covIdxV[i][j]= new IntVector(2);
						covValV[i][j]= new DoubleVector(2);
					}
				}
			}
		
			// iterate edges
			int nrC= 0;
			for (int i = 0; i < g.trpts.length; i++) {
				++constraintCtr;
				if (count)
					++nrC;
				else if (debug)
					System.out.println("C"+ constraintCtr+ "\t"+ g.trpts[i].getTranscriptID());
			}
			for (int i = 0; i < edges.length; i++) {
				Edge e= edges[i];
				if (!e.isExonic())
					continue;
				
				for (int sa = 0; sa < 2; sa++) {
					boolean sense= (sa== 0);
					
					nrC+= setCBase_save(count, g, e, sense, cSumV, count?null:cBalV[sa], costIdxV, costValV, covIdxV, covValV);
		
				}
			}
				
			// transcript balances
			if (count) {
				nrC+= 2* 2* g.trpts.length;	// sense/asense, +/-
				restrNr+= g.trpts.length;
			}
			for (int i = 0; (!count)&& i < cBalV.length; i++) {
				for (int j = 0; j < cBalV[i].length; j++) {
					IntVector v= cBalV[i][j];
					int len= v.size();
					int[] idx= new int[len+ 2];				// TODO reuse
					double[] val= new double[idx.length];	// TODO reuse
					for (int k = 0; k < val.length; k++) {
						idx[k]= Math.abs(v.get(k));
						val[k]= (cBalV[i][j].get(k)> 0? 1d: -1d);	// alternate
					}
					idx[len]= ++constraintCtr;
					val[len]= 1d;
					costIdxV.add(constraintCtr);
					costValV.add(1);
					idx[len+ 1]= ++constraintCtr;
					val[len+ 1]= -1d;
					costIdxV.add(constraintCtr);
					costValV.add(1);
					getLPsolve().addConstraintex(idx.length, val, idx, LpSolve.EQ, 0d);
					if (debug)
						printConstraint(System.out, g.trpts[j].getTranscriptID()+ " (bal)", 
								idx.length, val, idx, LpSolve.EQ, 0d); 
					++restrNr;
				}
			}
		
			
			if (count)
				return nrC;
				
			// OF
			double[] a= new double[constraintCtr+ 1];	// 1-based
			Arrays.fill(a, 0d);
			for (int i = 0; i < costIdxV.size(); ++i) 
				a[costIdxV.get(i)]= costValV.get(i);		
			getLPsolve().setObjFn(a);
			getLPsolve().setMinim();
			costIdx= costIdxV.toIntArray();		// TODO
			costVal= costValV.toDoubleArray();	// TODO
			
			return constraintCtr;
		}

		private void setCBase_deltaMulti(boolean count, Graph g, Edge e, boolean sense, IntVector[] cSumV, IntVector[] cBalV,
						IntVector costIdxV, DoubleVector costValV, IntVector[][] covIdxV, DoubleVector[][] covValV) throws LpSolveException {
					
					int constraintCheck= constraintCtr;
					
					assert(!(e instanceof SuperEdge));
					for (int i = 0; (!count)&& i < covIdxV.length; i++) {
						for (int j = 0; j < covValV[i].length; j++) {
							covIdxV[i][j]= new IntVector(2);
							covValV[i][j]= new DoubleVector(2);
						}
					}
					for (int i = 0; (!count)&& i < cSumV.length; i++) 
						cSumV[i].removeAll();
					
					int partSize= setCDeconvolution_deltaMulti(count, g, null, e, sense, cSumV, cBalV, costIdxV, costValV, covIdxV, covValV);
					
					// iterate superedges
					for (int i = 0; e.getSuperEdges()!= null&& i < e.getSuperEdges().size(); ++i) {
						SuperEdge se= e.getSuperEdges().elementAt(i);
						
						// sense/anti-sense.. e must be first/last in super-edge
						if ((sense&& se.getEdges()[0]!= e)
								|| ((!sense)&& se.getEdges()[se.getEdges().length- 1]!= e))
							continue;
						
						setCDeconvolution_deltaMulti(count, g, e, se, sense, cSumV, cBalV, costIdxV, costValV, covIdxV, covValV);
						
						// TODO 2nd level for paired-end
					}
					
					if (count) {
						restrNr+= partSize;	// only for expectation
						//constraintCtr+= (2* partSize);
						return;
					}
					
					// expectation, transcript d
					for (int i = 0; i < cSumV.length; i++) {
						int len= cSumV[i].size();
						if (len== 0)
							continue;	// t[i] not in e_trpts
						
						// values
						double[] val= new double[len+ 1];	// + 3 for d+/d-
						Arrays.fill(val, 2, val.length, 1d);
						//val[1]= -1d;
						Transcript tx= g.trpts[i];
						int tlen= tx.getExonicLength();
						UniversalMatrix m= profile.getMatrix(tlen);
						double f= m.getFrac(
									tx.getExonicPosition(e.getFrac(true)),
									tx.getExonicPosition(e.getFrac(false)),
									tlen,
									sense?Constants.DIR_FORWARD:Constants.DIR_BACKWARD);
						if (Double.isInfinite(f)|| Double.isNaN(f)) {
							System.err.println("infinite value");
							f= m.getFrac(
									tx.getExonicPosition(e.getFrac(true)),
									tx.getExonicPosition(e.getFrac(false)),
									tlen,
									sense?Constants.DIR_FORWARD:Constants.DIR_BACKWARD);
						}
						ccheck[i]+= f;
						val[0]= -f;
						//double cost= Math.exp(-f);
						
						// indices
						int[] idx= new int[len+ 1];	// 3 for d+/d-
						idx[0]= (i+1);
//						idx[1]= ++constraintCtr;	// d_exp-(e,t[i])
//						costIdxV.add(constraintCtr);
//						costValV.add(cost);
//						cBalV[i].add(-constraintCtr);	// tx balance
						// bound idx probably super-fluent as by C_tx>= 0
		/*				getLPsolve().addConstraintex(2, boundVal, boundIdx, LpSolve.GE, 0);
						++restrNr;
						if (debug)
							printConstraint(System.out, e+ " ("+(sense?"s, ":"a, ")+"ub d+)", 2, boundVal, boundIdx, LpSolve.GE, 0);
		*/					
		
//						idx[2]= ++constraintCtr;	// d_exp-(e,t[i])
//						costIdxV.add(constraintCtr);
//						costValV.add(cost);
//						cBalV[i].add(constraintCtr);	// tx balance
						System.arraycopy(cSumV[i].vector, 0, idx, 1, cSumV[i].length);	// 3 for d+/d-
						
						getLPsolve().addConstraintex(idx.length, val, idx, LpSolve.EQ, 0);
						++restrNr;
						if (debug)
							printConstraint(System.out, e+ " ("+ (sense?"s":"a")+ ", "+ g.trpts[i].getTranscriptID()+ ")", 
									idx.length, val, idx, LpSolve.EQ, 0d); 
					}
					
					// flux preservation, 2 for min/max readlength
					for (int i = 0; (!count)&& i < covIdxV.length; i++) {
						for (int j = 0; j < covValV[i].length; j++) {	// per tx
							// assure there is an exon and a sj
							if (covIdxV[i][j].length< 2) 	
								continue;				
							boolean plus= false, minus= false;
							for (int k = 0; k < covValV[i][j].length; k++) {
								double d= covValV[i][j].get(k);
								plus|= (d> 0);
								minus|= (d< 0);
								if (plus&& minus)
									break;
							}
							if (!(plus&& minus))
								continue;
							getLPsolve().addConstraintex(covIdxV[i][j].length, covValV[i][j].toDoubleArray(), 
									covIdxV[i][j].toIntArray(), i== 0? LpSolve.LE: LpSolve.GE, 0);
							++restrNr;
							if (debug)
								printConstraint(System.out, e+ " ("+ (sense?"s":"a")+ ", flux)", covIdxV[i][j].length, covValV[i][j].toDoubleArray(), 
										covIdxV[i][j].toIntArray(), i== 0? LpSolve.LE: LpSolve.GE, 0); 
						}
					}
					
				}

		private int setCAll_deltaMulti(boolean count, Graph g) throws LpSolveException {
					IntVector costIdxV= null;
					DoubleVector costValV= null;
					IntVector[] cSumV= null;
					IntVector[][] cBalV= null; // sense/asense, trpts[]
					IntVector[][] covIdxV= null;	// min/maxCov, trpts[]
					DoubleVector[][] covValV= null;
					Edge[] edges= g.getExonicEdgesInGenomicOrder();
		
					if (!count) {
						costIdxV= new IntVector(2* 2* g.trpts.length* edges.length);
						costValV= new DoubleVector(costIdxV.vector.length);
						cSumV= new IntVector[g.trpts.length];
						for (int i = 0; i < cSumV.length; i++) 
							cSumV[i]= new IntVector(5);	// se per tx in e
						cBalV= new IntVector[2][g.trpts.length];
						for (int i = 0; i < cBalV.length; i++) 
							for (int j = 0; j < cBalV[i].length; j++) 
								cBalV[i][j]= new IntVector(2* edges.length);	// 2 (+/-)* |E|
						covIdxV= new IntVector[2][g.trpts.length];
						covValV= new DoubleVector[2][g.trpts.length];
						for (int i = 0; i < covIdxV.length; i++) {
							for (int j = 0; j < g.trpts.length; j++) {
								covIdxV[i][j]= new IntVector(2);
								covValV[i][j]= new DoubleVector(2);
							}
						}
					}
		
					// iterate edges
					for (int i = 0; i < g.trpts.length; i++) {
						++constraintCtr;
						if (!count) {
							getLPsolve().setLowbo(constraintCtr, 1d);
							getLPsolve().setSemicont(constraintCtr, true);
							if (debug) {
								System.out.println("C"+ constraintCtr+ "\t"+ g.trpts[i].getTranscriptID());
								System.out.println("C"+ constraintCtr+ " >= 1 | 0");
							}
						}
						
					}
					for (int i = 0; i < edges.length; i++) {
						Edge e= edges[i];
						if (!e.isExonic())
							continue;
						
						for (int sa = 0; sa < 2; sa++) {
							boolean sense= (sa== 0);
							
							setCBase_deltaMulti(count, g, e, sense, cSumV, count?null:cBalV[sa], costIdxV, costValV, covIdxV, covValV);
		
						}
					}
						
					// transcript balances
					if (count) {
						restrNr+= 2* g.trpts.length; // sense/ asense
					}
					for (int i = 0; (!count)&& i < cBalV.length; i++) {
						for (int j = 0; j < cBalV[i].length; j++) {
							IntVector v= cBalV[i][j];
							int len= v.size();
							int[] idx= new int[len];				// TODO reuse
							double[] val= new double[idx.length];	// TODO reuse
							for (int k = 0; k < val.length; k++) {
								idx[k]= Math.abs(v.get(k));
								val[k]= (cBalV[i][j].get(k)> 0? 1d: -1d);	// alternate
							}
		//					idx[len]= ++constraintCtr;
		//					val[len]= 1d;
		//					costIdxV.add(constraintCtr);
		//					costValV.add(1);
		//					idx[len+ 1]= ++constraintCtr;
		//					val[len+ 1]= -1d;
		//					costIdxV.add(constraintCtr);
		//					costValV.add(1);
							
							//getLPsolve().addConstraintex(idx.length, val, idx, LpSolve.EQ, 0d);
							getLPsolve().addConstraintex(idx.length, val, idx, LpSolve.LE, 1d);
							getLPsolve().addConstraintex(idx.length, val, idx, LpSolve.GE, -1d);
							if (debug)
								System.out.println(g.trpts[j].getTranscriptID()+ " (bal)\t= [-1,1]");
		//						printConstraint(System.out, g.trpts[j].getTranscriptID()+ " (bal)", 
		//								idx.length, val, idx, LpSolve.LE, 1d); 
							++restrNr;
							++restrNr;
						}
					}
				
			
					
					if (count)
						return constraintCtr;
						
					// OF
					double[] a= new double[constraintCtr+ 1];	// 1-based
					Arrays.fill(a, 0d);
					for (int i = 0; i < costIdxV.size(); ++i) 
						a[costIdxV.get(i)]= costValV.get(i);		
					getLPsolve().setObjFn(a);
					getLPsolve().setMinim();
					costIdx= costIdxV.toIntArray();		// TODO
					costVal= costValV.toDoubleArray();	// TODO
					
					//assert(nrC== constraintCtr);
					return constraintCtr;
				}

		private int setCAll_2dim(boolean count, Graph g) throws LpSolveException {
					IntVector costIdxV= null;
					DoubleVector costValV= null;
					IntVector[] cSumV= null;
					IntVector[][] cBalV= null; // sense/asense, trpts[]
					IntVector[][] covIdxV= null;	// min/maxCov, trpts[]
					DoubleVector[][] covValV= null;
					Edge[] edges= g.getExonicEdgesInGenomicOrder();
		
					if (!count) {
						costIdxV= new IntVector(2* 2* g.trpts.length* edges.length);
						costValV= new DoubleVector(costIdxV.vector.length);
						cSumV= new IntVector[g.trpts.length];
						for (int i = 0; i < cSumV.length; i++) 
							cSumV[i]= new IntVector(5);	// se per tx in e
						cBalV= new IntVector[2][g.trpts.length];
						for (int i = 0; i < cBalV.length; i++) 
							for (int j = 0; j < cBalV[i].length; j++) 
								cBalV[i][j]= new IntVector(2* edges.length);	// 2 (+/-)* |E|
						covIdxV= new IntVector[2][g.trpts.length];
						covValV= new DoubleVector[2][g.trpts.length];
						for (int i = 0; i < covIdxV.length; i++) {
							for (int j = 0; j < g.trpts.length; j++) {
								covIdxV[i][j]= new IntVector(2);
								covValV[i][j]= new DoubleVector(2);
							}
						}
					}
		
					// iterate edges
					int nrC= 0;
					for (int i = 0; i < g.trpts.length; i++) {
						++constraintCtr;
						++nrC;
						if (!count) {
							getLPsolve().setLowbo(constraintCtr, 1d);
							getLPsolve().setSemicont(constraintCtr, true);
							if (debug) {
								System.out.println("C"+ constraintCtr+ "\t"+ g.trpts[i].getTranscriptID());
								System.out.println("C"+ constraintCtr+ " >= 1 | 0");
							}
						}
						
					}
					for (int i = 0; i < edges.length; i++) {
						Edge e= edges[i];
						if (!e.isExonic())
							continue;
						
						for (int sa = 0; sa < 2; sa++) {
							boolean sense= (sa== 0);
							
							nrC+= setCBase_2dim(count, g, e, sense, cSumV, count?null:cBalV[sa], costIdxV, costValV, covIdxV, covValV);
		
						}
					}
						
					// transcript balances
					if (count) {
						restrNr+= 2* g.trpts.length; // sense/ asense
					}
					for (int i = 0; (!count)&& i < cBalV.length; i++) {
						for (int j = 0; j < cBalV[i].length; j++) {
							IntVector v= cBalV[i][j];
							int len= v.size();
							int[] idx= new int[len];				// TODO reuse
							double[] val= new double[idx.length];	// TODO reuse
							for (int k = 0; k < val.length; k++) {
								idx[k]= Math.abs(v.get(k));
								val[k]= (cBalV[i][j].get(k)> 0? 1d: -1d);	// alternate
							}
		//					idx[len]= ++constraintCtr;
		//					val[len]= 1d;
		//					costIdxV.add(constraintCtr);
		//					costValV.add(1);
		//					idx[len+ 1]= ++constraintCtr;
		//					val[len+ 1]= -1d;
		//					costIdxV.add(constraintCtr);
		//					costValV.add(1);
							
							//getLPsolve().addConstraintex(idx.length, val, idx, LpSolve.EQ, 0d);
							getLPsolve().addConstraintex(idx.length, val, idx, LpSolve.LE, 1d);
							getLPsolve().addConstraintex(idx.length, val, idx, LpSolve.GE, -1d);
							if (debug)
								System.out.println(g.trpts[j].getTranscriptID()+ " (bal)\t= [-1,1]");
		//						printConstraint(System.out, g.trpts[j].getTranscriptID()+ " (bal)", 
		//								idx.length, val, idx, LpSolve.LE, 1d); 
							++restrNr;
							++restrNr;
						}
					}
				
			
					
					if (count)
						return nrC;
						
					// OF
					double[] a= new double[constraintCtr+ 1];	// 1-based
					Arrays.fill(a, 0d);
					for (int i = 0; i < costIdxV.size(); ++i) 
						a[costIdxV.get(i)]= costValV.get(i);		
					getLPsolve().setObjFn(a);
					getLPsolve().setMinim();
					costIdx= costIdxV.toIntArray();		// TODO
					costVal= costValV.toDoubleArray();	// TODO
					
					//assert(nrC== constraintCtr);
					return constraintCtr;
				}

		private int setCBase_2dim(boolean count, Graph g, Edge e, boolean sense, IntVector[] cSumV, IntVector[] cBalV,
						IntVector costIdxV, DoubleVector costValV, IntVector[][] covIdxV, DoubleVector[][] covValV) throws LpSolveException {
					
					int constraintCheck= constraintCtr;
					
					assert(!(e instanceof SuperEdge));
					for (int i = 0; (!count)&& i < covIdxV.length; i++) {
						for (int j = 0; j < covValV[i].length; j++) {
							covIdxV[i][j]= new IntVector(2);
							covValV[i][j]= new DoubleVector(2);
						}
					}
					for (int i = 0; (!count)&& i < cSumV.length; i++) 
						cSumV[i].removeAll();
					
					int nrC= setCDeconvolution_2dim(count, g, null, e, sense, cSumV, cBalV, costIdxV, costValV, covIdxV, covValV);
					int partSize= nrC- 2;
					
					// iterate superedges
					for (int i = 0; e.getSuperEdges()!= null&& i < e.getSuperEdges().size(); ++i) {
						SuperEdge se= e.getSuperEdges().elementAt(i);
						
						// sense/anti-sense.. e must be first/last in super-edge
						if ((sense&& se.getEdges()[0]!= e)
								|| ((!sense)&& se.getEdges()[se.getEdges().length- 1]!= e))
							continue;
						
						nrC+= setCDeconvolution_2dim(count, g, e, se, sense, cSumV, cBalV, costIdxV, costValV, covIdxV, covValV);
						
						// TODO 2nd level for paired-end
					}
					
					nrC+= (2* partSize);
					if (count) {
						restrNr+= partSize;	// only for expectation
						constraintCtr+= (2* partSize);
						return nrC;
					}
					
					// expectation, transcript d
					for (int i = 0; i < cSumV.length; i++) {
						int len= cSumV[i].size();
						if (len== 0)
							continue;	// t[i] not in e_trpts
						
						// values
						double[] val= new double[len+ 3];
						Arrays.fill(val, 2, val.length, 1d);
						val[1]= -1d;
						Transcript tx= g.trpts[i];
						int tlen= tx.getExonicLength();
						UniversalMatrix m= profile.getMatrix(tlen);
						double f= m.getFrac(
									tx.getExonicPosition(e.getFrac(true)),
									tx.getExonicPosition(e.getFrac(false)),
									tlen,
									sense?Constants.DIR_FORWARD:Constants.DIR_BACKWARD);
						if (Double.isInfinite(f)|| Double.isNaN(f)) {
							System.err.println("infinite value");
							f= m.getFrac(
									tx.getExonicPosition(e.getFrac(true)),
									tx.getExonicPosition(e.getFrac(false)),
									tlen,
									sense?Constants.DIR_FORWARD:Constants.DIR_BACKWARD);
						}
						ccheck[i]+= f;
						val[0]= -f;
						double cost= Math.exp(-(partSize- 1));
						
						// indices
						int[] idx= new int[len+ 3];	// TODO reuse
						idx[0]= (i+1);
						idx[1]= ++constraintCtr;	// d_exp-(e,t[i])
						costIdxV.add(constraintCtr);
						costValV.add(cost);
						cBalV[i].add(-constraintCtr);	// tx balance
						// bound idx probably super-fluent as by C_tx>= 0
		/*				getLPsolve().addConstraintex(2, boundVal, boundIdx, LpSolve.GE, 0);
						++restrNr;
						if (debug)
							printConstraint(System.out, e+ " ("+(sense?"s, ":"a, ")+"ub d+)", 2, boundVal, boundIdx, LpSolve.GE, 0);
		*/					
		
						idx[2]= ++constraintCtr;	// d_exp-(e,t[i])
						costIdxV.add(constraintCtr);
						costValV.add(cost);
						cBalV[i].add(constraintCtr);	// tx balance
						System.arraycopy(cSumV[i].vector, 0, idx, 3, cSumV[i].length);
						
						getLPsolve().addConstraintex(idx.length, val, idx, LpSolve.EQ, 0);
						++restrNr;
						if (debug)
							printConstraint(System.out, e+ " ("+ (sense?"s":"a")+ ", "+ g.trpts[i].getTranscriptID()+ ")", 
									idx.length, val, idx, LpSolve.EQ, 0d); 
					}
					
					// flux preservation, 2 for min/max readlength
					for (int i = 0; (!count)&& i < covIdxV.length; i++) {
						for (int j = 0; j < covValV[i].length; j++) {	// per tx
							// assure there is an exon and a sj
							if (covIdxV[i][j].length< 2) 	
								continue;				
							boolean plus= false, minus= false;
							for (int k = 0; k < covValV[i][j].length; k++) {
								double d= covValV[i][j].get(k);
								plus|= (d> 0);
								minus|= (d< 0);
								if (plus&& minus)
									break;
							}
							if (!(plus&& minus))
								continue;
							getLPsolve().addConstraintex(covIdxV[i][j].length, covValV[i][j].toDoubleArray(), 
									covIdxV[i][j].toIntArray(), i== 0? LpSolve.LE: LpSolve.GE, 0);
							++restrNr;
							if (debug)
								printConstraint(System.out, e+ " ("+ (sense?"s":"a")+ ", flux)", covIdxV[i][j].length, covValV[i][j].toDoubleArray(), 
										covIdxV[i][j].toIntArray(), i== 0? LpSolve.LE: LpSolve.GE, 0); 
						}
					}
					
		//			if (constraintCheck+ nrC!= constraintCtr)
		//				System.currentTimeMillis();
					assert(constraintCheck+ nrC== constraintCtr);
					return nrC;
				}

		private int setCDeconvolution_2dim(boolean count, Graph g, Edge f, Edge e, boolean sense, IntVector[] cSumV, IntVector[] cBalV,
						IntVector costIdxV, DoubleVector costValV, IntVector[][] covIdxV, DoubleVector[][] covValV) throws LpSolveException {
					
					int constraintCheck= constraintCtr;
					
					// attributes of e
					int ctrTT= 0;
					int obs= (sense? e.getReadNr(): e.getRevReadNr());
					boolean sj= (e instanceof SuperEdge)&& (!((SuperEdge) e).isPend());
					double effLenMin= calcEffLen(e, f, true);
					double effLenMax= calcEffLen(e, f, false);
					assert(effLenMin> 0|| obs== 0);
					double cost= getCostObs(sj, obs, effLenMax);
		
					// all transcript contributions
					for (int j = 0; j < g.trpts.length; j++) {
		
						if (!g.contains(e.getTranscripts(), j))
							continue;
						if (count) {
							++ctrTT;
							++constraintCtr;
							++restrNr;
							continue;
						}
						
						//int x= ctrTT* 3;	// contribution, d-, d+
						decIdx[ctrTT++]= ++constraintCtr;		// c
						cSumV[j].add(constraintCtr);	// sum of (e_j,t_i)
						if ((sj &&effLenMin> 0)|| ((!sj)&& effLenMax> 0)) {
							covIdxV[0][j].add(constraintCtr);
							covValV[0][j].add(1d/ (sj? -effLenMin: effLenMax));	// inverse
						}
						if ((sj &&effLenMax> 0)|| ((!sj)&& effLenMin> 0)) {
							covIdxV[1][j].add(constraintCtr);
							covValV[1][j].add(1d/ (sj? -effLenMax: effLenMin));	// inverse
						}
					}
					
					// obs: deconvolution
					if (count) {
						restrNr+= 2;
						constraintCtr+= 2;
					} else {
						
						Arrays.fill(decVal, 0, ctrTT+ 1, 1d);
						
						// d+, substract reads from obs 
						decIdx[ctrTT]= ++constraintCtr;	// d+
						if (effLenMax> 0) {	// costs only iff there could be mappings
							costIdxV.add(constraintCtr);
							costValV.add(cost);
						}
						double bound= Math.max(obs- 1d, 0d);	// last man standing
						// ERROR in set_upbo: status = -1 (Model has not been optimized)
						// getLPsolve().setUpbo(constraintCtr, bound);
						// if (debug)
						//	System.out.println(e+ " ("+(sense?"s, ":"a, ")+"ub d-)\t<= "+ bound);
						getLPsolve().addConstraintex(1, new double[] {1d}, new int[] {constraintCtr}, 
								LpSolve.LE, bound);
						++restrNr;
						
						// d-, adding reads to obs
						decIdx[ctrTT+ 1]= ++constraintCtr;
						decVal[ctrTT+ 1]= -1d;
						if (effLenMax> 0) {	// costs only iff there could be mappings
							costIdxV.add(constraintCtr);
							costValV.add(cost);
						}
						
						// deconvolution
						getLPsolve().addConstraintex(ctrTT+ 2, decVal, decIdx, LpSolve.EQ, obs);
						if (debug)
							printConstraint(System.out, e+ " ("+(sense?"s":"a")+", dec)", ctrTT+ 2, decVal, decIdx, LpSolve.EQ, obs); 
						++restrNr;
					}
					
		//			if (constraintCheck+ ctrTT+ 2!= constraintCtr)
		//				System.currentTimeMillis();
					assert(constraintCheck+ ctrTT+ 2== constraintCtr);
					return (ctrTT+ 2);
		
				}

		private int setCAll(byte mode, Graph g) throws LpSolveException {
			IntVector costIdxV= null;
			DoubleVector costValV= null;
			IntVector[] cSumV= null;
			IntVector[][] cBalV= null; // sense/asense, trpts[]
			IntVector[][] covIdxV= null;	// min/maxCov, trpts[]
			DoubleVector[][] covValV= null;
			Edge[] edges= g.getExonicEdgesInGenomicOrder();
		
			if (mode!= MODE_COUNT) {
				costIdxV= new IntVector(2* 2* g.trpts.length* edges.length);
				costValV= new DoubleVector(costIdxV.vector.length);
				cSumV= new IntVector[g.trpts.length];
				for (int i = 0; i < cSumV.length; i++) 
					cSumV[i]= new IntVector(5);	// se per tx in e
				cBalV= new IntVector[2][g.trpts.length];
				for (int i = 0; i < cBalV.length; i++) 
					for (int j = 0; j < cBalV[i].length; j++) 
						cBalV[i][j]= new IntVector(2* edges.length);	// 2 (+/-)* |E|
				covIdxV= new IntVector[2][g.trpts.length];
				covValV= new DoubleVector[2][g.trpts.length];
				for (int i = 0; i < covIdxV.length; i++) {
					for (int j = 0; j < g.trpts.length; j++) {
						covIdxV[i][j]= new IntVector(2);
						covValV[i][j]= new DoubleVector(2);
					}
				}
			}
		
			// iterate transcripts
			int nrC= 0;
			for (int i = 0; i < g.trpts.length; i++) {
				constraintCtr+= 2;	// sense, anti-sense
				if (mode== MODE_COUNT)
					nrC+= 2;
				else if (mode== MODE_COMPLETE) {
					System.out.println(g.trpts[i].getTranscriptID()+ "\t"+ "C"+ constraintCtr);
				} else if (mode== MODE_CONSTRAINT) {
					int[] idx= new int[2];
					idx[0]= constraintCtr- 1;
					idx[1]= constraintCtr;
					
					Transcript tx= g.trpts[i];
					int tlen= tx.getExonicLength();
					UniversalMatrix m= profile.getMatrix(tlen);
					double[] val= new double[2];
					val[0]= -1d;
					val[1]= m.sums/ (double) m.suma;
					getLPsolve().addConstraintex(idx.length, val, idx, LpSolve.EQ, 0d);
				}
				if (mode!= MODE_COMPLETE)
					++restrNr;
			}
			// iterate edges
			for (int i = 0; i < edges.length; i++) {
				Edge e= edges[i];
				if (!e.isExonic())
					continue;
				
				for (int sa = 0; sa < 2; sa++) {
					boolean sense= (sa== 0);
					setCBase(mode, g, e, sense, cSumV, (mode== MODE_COUNT)?null:cBalV[sa], costIdxV, costValV, covIdxV, covValV);
				}
			}
				
			// transcript balances
			if (mode== MODE_COUNT) {
				constraintCtr+= 2* 2* g.trpts.length;	// sense/asense, +/-
				restrNr+= g.trpts.length;
			}
			IntVector allLenV= new IntVector(g.trpts.length);
			for (int i = 0; (mode!= MODE_COUNT)&& i < cBalV.length; i++) {
				int tlen= g.trpts[i].getExonicLength();
				allLenV.add(tlen);
			}
			Distribution d= new Distribution(allLenV.toIntArray());
			double medLen= d.getMedian();
			for (int i = 0; (mode!= MODE_COUNT)&& i < cBalV.length; i++) {
				for (int j = 0; j < cBalV[i].length; j++) {
					IntVector v= cBalV[i][j];
					int len= v.size();
					int[] idx= new int[len+ 2];				// TODO reuse
					double[] val= new double[idx.length];	// TODO reuse
					for (int k = 0; k < val.length; k++) {
						idx[k]= Math.abs(v.get(k));
						val[k]= (cBalV[i][j].get(k)> 0? 1d: -1d);	// alternate
					}
					int tlen= g.trpts[i].getExonicLength();

					//
					double cost= 10d* (1d+ (tlen/ medLen));
					idx[len]= ++constraintCtr;
					val[len]= 1d;
					costIdxV.add(constraintCtr);
					costValV.add(cost);
					idx[len+ 1]= ++constraintCtr;
					val[len+ 1]= -1d;
					costIdxV.add(constraintCtr);
					costValV.add(cost);
					if (mode== MODE_CONSTRAINT) {
						getLPsolve().addConstraintex(idx.length, val, idx, LpSolve.EQ, 0d);
						++restrNr;
					} else if (mode== MODE_COMPLETE) {
						printConstraint(System.out, g.trpts[j].getTranscriptID()+ " (bal)", 
								idx.length, val, idx, LpSolve.EQ, 0d);
					}
				}
			}
		
			if (mode== MODE_COUNT)
				return nrC;
				
			// OF
			double[] a= new double[constraintCtr+ 1];	// 1-based
			Arrays.fill(a, 0d);
			for (int i = 0; i < costIdxV.size(); ++i)  {
				double val= costValV.get(i);
				val= val== 0? 1000: val;
				a[costIdxV.get(i)]= val;	
			}
			getLPsolve().setObjFn(a);
			getLPsolve().setMinim();
			costIdx= costIdxV.toIntArray();		// TODO
			costVal= costValV.toDoubleArray();	// TODO
			
			return constraintCtr;
		}

		private int setCBase_pgauss(byte mode, Graph g, Edge e, boolean sense, IntVector[] cSumV, IntVector[] cBalV,
				IntVector costIdxV, DoubleVector costValV, IntVector[][] covIdxV, DoubleVector[][] covValV) throws LpSolveException {
			
			assert(!(e instanceof SuperEdge));
			for (int i = 0; (mode!= MODE_COUNT)&& i < covIdxV.length; i++) {
				for (int j = 0; j < covValV[i].length; j++) {
					covIdxV[i][j]= new IntVector(2);
					covValV[i][j]= new DoubleVector(2);
				}
			}
			for (int i = 0; (mode!= MODE_COUNT)&& i < cSumV.length; i++) 
				cSumV[i].removeAll();
			
			int nrC= setCDeconvolution_pgauss(mode, g, null, e, sense, cSumV, cBalV, costIdxV, costValV, covIdxV, covValV);
			
			// iterate superedges
			for (int i = 0; e.getSuperEdges()!= null&& i < e.getSuperEdges().size(); ++i) {
				SuperEdge se= e.getSuperEdges().elementAt(i);
				
				// sense/anti-sense.. e must be first/last in super-edge
				if ((sense&& se.getEdges()[0]!= e)
						|| ((!sense)&& se.getEdges()[se.getEdges().length- 1]!= e))
					continue;
				
				nrC+= setCDeconvolution_pgauss(mode, g, e, se, sense, cSumV, cBalV, costIdxV, costValV, covIdxV, covValV);
				
				// TODO 2nd level for paired-end
			}
			if (mode== MODE_COUNT) {
				// NP: restrNr+= cSumV.length+ 2;	// < cSumV.length
				return nrC;
			}
			
			// expectation restriction
			for (int i = 0; i < cSumV.length; i++) {
				int len= cSumV[i].size();
				if (len== 0)
					continue;
				int[] idx= new int[len+ 1];	// TODO reuse
				idx[0]= sense? ((i* 2)+1): ((i* 2)+2);
				System.arraycopy(cSumV[i].vector, 0, idx, 1, cSumV[i].length);
				double[] val= new double[idx.length];
				Arrays.fill(val, 1d);
				Transcript tx= g.trpts[i];
				int tlen= tx.getExonicLength();
				UniversalMatrix m= profile.getMatrix(tlen);
				double f= m.getFrac(
							tx.getExonicPosition(e.getFrac(true)),
							tx.getExonicPosition(e.getFrac(false)),
							tlen,
							sense?Constants.DIR_FORWARD:Constants.DIR_BACKWARD);
				if (Double.isInfinite(f)|| Double.isNaN(f)) {
					System.err.println("infinite value");
					f= m.getFrac(
							tx.getExonicPosition(e.getFrac(true)),
							tx.getExonicPosition(e.getFrac(false)),
							tlen,
							sense?Constants.DIR_FORWARD:Constants.DIR_BACKWARD);
				}
				ccheck[i]+= f;
				val[0]= -f;
				
				if (mode== MODE_CONSTRAINT) {
					getLPsolve().addConstraintex(idx.length, val, idx, LpSolve.EQ, 0);
					++restrNr;
				} else if (false&& mode== MODE_COMPLETE) {
					printConstraint(System.out, e+ " ("+ (sense?"s":"a")+ ", "+ g.trpts[i].getTranscriptID()+ ")", 
							idx.length, val, idx, LpSolve.EQ, 0d);
				}
			}
			
			// flux preservation, 2 for min/max readlength
			for (int i = 0; (mode!= MODE_COUNT)&& i < covIdxV.length; i++) {
				for (int j = 0; j < covValV[i].length; j++) {	// per tx
					// assure there is an exon and a sj
					if (covIdxV[i][j].length< 2) 	
						continue;				
					boolean plus= false, minus= false;
					for (int k = 0; k < covValV[i][j].length; k++) {
						double d= covValV[i][j].get(k);
						plus|= (d> 0);
						minus|= (d< 0);
						if (plus&& minus)
							break;
					}
					if (!(plus&& minus))
						continue;
					if (mode== MODE_CONSTRAINT) {
						getLPsolve().addConstraintex(covIdxV[i][j].length, covValV[i][j].toDoubleArray(), 
								covIdxV[i][j].toIntArray(), i== 0? LpSolve.LE: LpSolve.GE, 0);
						++restrNr;
					} else if (false&& mode== MODE_COMPLETE) {
						printConstraint(System.out, e+ " (flux)", covIdxV[i][j].length, covValV[i][j].toDoubleArray(), 
								covIdxV[i][j].toIntArray(), i== 0? LpSolve.LE: LpSolve.GE, 0);
					}
				}
			}
			
			
			return nrC;
		}

		private int setCDeconvolution_save(boolean count, Graph g, Edge f, Edge e, boolean sense, IntVector[] cSumV, IntVector[] cBalV,
						IntVector costIdxV, DoubleVector costValV, IntVector[][] covIdxV, DoubleVector[][] covValV) throws LpSolveException {
					
					// obs: bounds
					int ctrTT= 0;
					int obs= (sense? e.getReadNr(): e.getRevReadNr());
					double cost= getCostObs(obs);
					IntVector minusIdxV= null;
					if (!count) 
						minusIdxV= new IntVector(g.trpts.length);	// TODO reuse
					for (int j = 0; j < g.trpts.length; j++) {
						
						if (!g.contains(e.getTranscripts(), j))
							continue;
						if (count) {
							++ctrTT;
							++restrNr;
							continue;
						}
						
						SuperEdge se= null;
						boolean sj= false;
						double effLenMin= 0, effLenMax= 0;
						if (e instanceof SuperEdge) {
							se= (SuperEdge) e;
							Edge[] see= se.getEdges();
							sj= (se!= null&& !se.isPend());
							if (sj) {
								int interNt= 1;
								for (int i = 1; i < se.getEdges().length- 1; i++) 
									interNt+= se.getEdges()[i].length();
								effLenMin= Math.min(Math.min(mapLenMin- interNt, see[0].length()),
												see[see.length- 1].length());
								effLenMax= Math.min(Math.min(mapLenMax- interNt, see[0].length())
												, see[see.length- 1].length());
							} else {
								//effLenMin= effLenMax= f.length(); // not: - fixedReadLen+ 1;
								int base= f.length();
		/*						int lim= sense? g.trpts[j].getExonicLength()- g.trpts[j].getExonicPosition(f.getHead().getSite().getPos()):
									g.trpts[j].getExonicPosition(f.getTail().getSite().getPos());
								effLenMin= base+ Math.min(lim, mapLenMin);
								effLenMax= base+ Math.min(lim, mapLenMax);
		*/
								effLenMin= base- mapLenMax+ 1;
								effLenMax= base- mapLenMin+ 1;
							}
						} else {
							// 0effLenMin= effLenMax= e.length(); // not: - fixedReadLen+ 1;
							int base= e.length();
		/*					int lim= sense? g.trpts[j].getExonicLength()- g.trpts[j].getExonicPosition(e.getHead().getSite().getPos()):
								g.trpts[j].getExonicPosition(e.getTail().getSite().getPos());
							effLenMin= base+ Math.min(lim, mapLenMin);
							effLenMax= base+ Math.min(lim, mapLenMax);
		*/					
							effLenMin= base- mapLenMax+ 1;
							effLenMax= base- mapLenMin+ 1;
						}
						assert(effLenMin> 0|| obs== 0);
						
						int x= ctrTT* 3;
						
						decIdx[x]= ++constraintCtr;		// c
						cSumV[j].add(constraintCtr);
						//if (se!= null) ???
						boundIdx[0]= constraintCtr;
						//boundVal[0]= 1d;	// TODO init once
						if ((sj &&effLenMin> 0)|| ((!sj)&& effLenMax> 0)) {
							covIdxV[0][j].add(constraintCtr);
							covValV[0][j].add(1d/ (sj? -effLenMin: effLenMax));	// inverse
						}
						if ((sj &&effLenMax> 0)|| ((!sj)&& effLenMin> 0)) {
							covIdxV[1][j].add(constraintCtr);
							covValV[1][j].add(1d/ (sj? -effLenMax: effLenMin));	// inverse
						}
						
						// d-, adding reads to obs
						decIdx[x+ 1]= ++constraintCtr;
						if (effLenMin<= 0&& (!sj)) {
							
							boolean found= false;
							for (int k = 0; (!found)&& e.getSuperEdges()!= null&& k < e.getSuperEdges().size(); k++) {
								SuperEdge se2 = (SuperEdge) e.getSuperEdges().elementAt(k);
								if (se2.isPend()|| (sense?se2.getEdges()[0]!= e: se2.getEdges()[se2.getEdges().length- 1]!= e))
									continue;
								if (g.contains(se2.getTranscripts(), j)) {
									found= true;
									break;
								}
							}
							if (found) {
		//						&& (e.getTail().getInEdges().elementAt(0).getTail()!= g.root)
		//						&& (e.getHead().getOutEdges().elementAt(0).getHead()!= g.leaf)
		//						&& e.toString().equals("154351122-154351137[")) {	// NA if not possible, SJ
		//					if (e.getSuperEdges()== null&& e.getSuperEdges().size()== 0)
		//						System.currentTimeMillis();
		//					Transcript[] ttt= g.decodeTset(e.getTranscripts());
		//					assert(e.getSuperEdges()!= null&& e.getSuperEdges().size()> 0);
								getLPsolve().setUpbo(constraintCtr, 0);
								if (debug)
									System.out.println(e+ " ("+(sense?"s, ":"a, ")+"ub d+)\tC"+constraintCtr+ " = 0");
							}
						}
						cBalV[j].add(-constraintCtr);
						boundIdx[1]= constraintCtr;
						//boundVal[1]= -1d;	// TODO init once
						if (effLenMax> 0) {	// costs only iff there could be mappings
							cost= getCostObs(sj, obs, effLenMax);
							if (cost< 1)
								System.currentTimeMillis();
							costIdxV.add(constraintCtr);
							costValV.add(cost);
						}
						// cov constraint only with assigned Cs
		/*				covIdxV[0].add(constraintCtr);
						covValV[0].add(sj? effLenMin: -effLenMin);	// inverse
						covIdxV[1].add(constraintCtr);
						covValV[1].add(sj? effLenMax: -effLenMax);	// inverse
		*/				
						getLPsolve().addConstraintex(2, boundVal, boundIdx, LpSolve.GE, 0);
						++restrNr;
						if (debug)
							printConstraint(System.out, e+ " ("+(sense?"s, ":"a, ")+"ub d+)", 2, boundVal, boundIdx, LpSolve.GE, 0);
		
						
						// d-, adding reads to obs 
						// bounded by c of the tx
						decIdx[x+ 2]= ++constraintCtr;	// d+
						if ((!count)) 
							minusIdxV.add(constraintCtr);	// last man standing				
						cBalV[j].add(constraintCtr);
						if (effLenMax> 0) {
							cost= getCostObs(sj, obs, effLenMax);
							if (cost< 1)
								System.currentTimeMillis();
							costIdxV.add(constraintCtr);
							costValV.add(cost);
						}
						// cov constraint only with assigned Cs
		/*				covIdxV[0].add(constraintCtr);
						covValV[0].add(sj? -effLenMin: effLenMin);	// inverse
						covIdxV[1].add(constraintCtr);
						covValV[1].add(sj? -effLenMax: effLenMax);	// inverse
		*/				
						
						++ctrTT;
					}
					
					// obs: deconvolution
					if (!count) {
						// limit last man standing
						double bound= Math.max(obs- 1d, 0d);
						int[] idx= minusIdxV.toIntArray();
						double[] val= new double[idx.length];
						Arrays.fill(val, 1d);
						getLPsolve().addConstraintex(idx.length, val, idx, LpSolve.LE, bound);
						++restrNr;
						if (debug)
							printConstraint(System.out, e+ " ("+(sense?"s, ":"a, ")+"ub d-)", idx.length, val, idx, LpSolve.LE, bound);
						
						// deconvolution
						getLPsolve().addConstraintex(ctrTT* 3, decVal, decIdx, LpSolve.EQ, obs);
						if (debug)
							printConstraint(System.out, e+ " ("+(sense?"s, ":"a, ")+"dec)", ctrTT* 3, decVal, decIdx, LpSolve.EQ, obs); 
					}
					++restrNr;
					
					return (ctrTT* 3);
				}

		private void completeConstraint(PrintStream p, String id, int length, double[] val, int[] idx,
				int eq, double d) {
			StringBuilder sb= new StringBuilder(id);
			sb.append("\t");
			for (int i = 0; i < length; i++) {
				double res= result[restrNr+ idx[i]];
				double var= val[i]* res;				
				sb.append(var> 0?"+":(var< 0? "-": (val[i]>= 0? "+": "-")));
				sb.append(Math.abs(var)== 1?"": StringUtils.fprint(Math.abs(var), 2));
//				sb.append("C");
//				sb.append(idx[i]);
				sb.append(" ");
				
				int pos= Arrays.binarySearch(costIdx, idx[i]);
				if (pos>= 0) {
					sb.append("(*");
					sb.append(StringUtils.fprint(costVal[pos], 2));
					sb.append("= ");
					sb.append(StringUtils.fprint(costVal[pos] * res, 2));
					sb.append(") ");
				}
			}
			sb.append(eq== LpSolve.EQ?"=":(eq== LpSolve.GE?">=":"<="));
			sb.append(" ");
			sb.append(d);
			p.println(sb.toString());
			
		}

		private int setCAll_pgauss(byte mode, Graph g) throws LpSolveException {
			IntVector costIdxV= null;
			DoubleVector costValV= null;
			IntVector[] cSumV= null;
			IntVector[][] cBalV= null; // sense/asense, trpts[]
			IntVector[][] covIdxV= null;	// min/maxCov, trpts[]
			DoubleVector[][] covValV= null;
			Edge[] edges= g.getExonicEdgesInGenomicOrder();
		
			if (mode!= MODE_COUNT) {
				costIdxV= new IntVector(2* 2* g.trpts.length* edges.length);
				costValV= new DoubleVector(costIdxV.vector.length);
				cSumV= new IntVector[g.trpts.length];
				for (int i = 0; i < cSumV.length; i++) 
					cSumV[i]= new IntVector(5);	// se per tx in e
				cBalV= new IntVector[2][g.trpts.length];
				for (int i = 0; i < cBalV.length; i++) 
					for (int j = 0; j < cBalV[i].length; j++) 
						cBalV[i][j]= new IntVector(2* edges.length);	// 2 (+/-)* |E|
				covIdxV= new IntVector[2][g.trpts.length];
				covValV= new DoubleVector[2][g.trpts.length];
				for (int i = 0; i < covIdxV.length; i++) {
					for (int j = 0; j < g.trpts.length; j++) {
						covIdxV[i][j]= new IntVector(2);
						covValV[i][j]= new DoubleVector(2);
					}
				}
			}
		
			// iterate transcripts
			int nrC= 0;
			for (int i = 0; i < g.trpts.length; i++) {
				constraintCtr+= 2;	// sense, anti-sense
				if (mode== MODE_COUNT)
					nrC+= 2;
				else if (mode== MODE_COMPLETE) {
					System.out.println(g.trpts[i].getTranscriptID()+ "\t"+ "C"+ constraintCtr);
				} else if (mode== MODE_CONSTRAINT) {
					int[] idx= new int[2];
					idx[0]= constraintCtr- 1;
					idx[1]= constraintCtr;
					
					Transcript tx= g.trpts[i];
					int tlen= tx.getExonicLength();
					UniversalMatrix m= profile.getMatrix(tlen);
					double[] val= new double[2];
					val[0]= -1d;
					val[1]= m.sums/ (double) m.suma;
					getLPsolve().addConstraintex(idx.length, val, idx, LpSolve.EQ, 0d);
				}
				if (mode!= MODE_COMPLETE)
					++restrNr;
			}
			// iterate edges
			for (int i = 0; i < edges.length; i++) {
				Edge e= edges[i];
				if (!e.isExonic())
					continue;
				
				for (int sa = 0; sa < 2; sa++) {
					boolean sense= (sa== 0);
					
					nrC+= setCBase_pgauss(mode, g, e, sense, cSumV, (mode== MODE_COUNT)?null:cBalV[sa], costIdxV, costValV, covIdxV, covValV);
		
				}
			}
				
			// transcript balances
			if (mode== MODE_COUNT) {
				nrC+= 2* 2* g.trpts.length;	// sense/asense, +/-
				restrNr+= g.trpts.length;
			}
			IntVector allLenV= new IntVector(g.trpts.length);
			for (int i = 0; (mode!= MODE_COUNT)&& i < cBalV.length; i++) {
				int tlen= g.trpts[i].getExonicLength();
				allLenV.add(tlen);
			}
			Distribution d= new Distribution(allLenV.toIntArray());
			double medLen= d.getMedian();
			for (int i = 0; (mode!= MODE_COUNT)&& i < cBalV.length; i++) {
				for (int j = 0; j < cBalV[i].length; j++) {
					IntVector v= cBalV[i][j];
					int len= v.size();
					int[] idx= new int[len+ 2];				// TODO reuse
					double[] val= new double[idx.length];	// TODO reuse
					for (int k = 0; k < val.length; k++) {
						idx[k]= Math.abs(v.get(k));
						val[k]= (cBalV[i][j].get(k)> 0? 1d: -1d);	// alternate
					}
					int tlen= g.trpts[i].getExonicLength();
		
					//
					double cost= 10d* (1d+ (tlen/ medLen));
					idx[len]= ++constraintCtr;
					val[len]= 1d;
					costIdxV.add(constraintCtr);
					costValV.add(cost);
					idx[len+ 1]= ++constraintCtr;
					val[len+ 1]= -1d;
					costIdxV.add(constraintCtr);
					costValV.add(cost);
					if (mode== MODE_CONSTRAINT) {
						getLPsolve().addConstraintex(idx.length, val, idx, LpSolve.EQ, 0d);
						++restrNr;
					} else if (mode== MODE_COMPLETE) {
						printConstraint(System.out, g.trpts[j].getTranscriptID()+ " (bal)", 
								idx.length, val, idx, LpSolve.EQ, 0d);
					}
				}
			}
		
			if (mode== MODE_COUNT)
				return nrC;
				
			// OF
			double[] a= new double[constraintCtr+ 1];	// 1-based
			Arrays.fill(a, 0d);
			for (int i = 0; i < costIdxV.size(); ++i) 
				a[costIdxV.get(i)]= costValV.get(i);		
			getLPsolve().setObjFn(a);
			getLPsolve().setMinim();
			costIdx= costIdxV.toIntArray();		// TODO
			costVal= costValV.toDoubleArray();	// TODO
			
			return constraintCtr;
		}

		private void setCBase(byte mode, Graph g, Edge e, boolean sense, IntVector[] cSumV, IntVector[] cBalV,
				IntVector costIdxV, DoubleVector costValV, IntVector[][] covIdxV, DoubleVector[][] covValV) throws LpSolveException {
			
			assert(!(e instanceof SuperEdge));
			for (int i = 0; (mode!= MODE_COUNT)&& i < covIdxV.length; i++) {
				for (int j = 0; j < covValV[i].length; j++) {
					covIdxV[i][j]= new IntVector(2);
					covValV[i][j]= new DoubleVector(2);
				}
			}
			for (int i = 0; (mode!= MODE_COUNT)&& i < cSumV.length; i++) 
				cSumV[i].removeAll();
			
			int obs= (sense? e.getReadNr(): e.getRevReadNr());
			assert(!(e instanceof SuperEdge));
			double effLenMin= calcEffLen(e, null, true);
			double effLenMax= calcEffLen(e, null, false);
			assert(effLenMin> 0|| obs== 0);
			double cov= obs/ effLenMin;
			
			int decLen= setCDeconvolution(mode, g, null, -1d, null, -1,  
					e, sense, cSumV, cBalV, costIdxV, costValV, covIdxV, covValV);
			int[] decBase= null;
			if (mode!= MODE_COUNT)
				decBase= decIdx.clone();
			int decPos= 0;
			
			// iterate superedges
			for (int i = 0; e.getSuperEdges()!= null&& i < e.getSuperEdges().size(); ++i) {
				
				SuperEdge se= e.getSuperEdges().elementAt(i);
				if (se.isPend())
					continue;	// TODO 2nd level for paired-end
				if ((sense&& se.getEdges()[0]!= e)
						|| ((!sense)&& se.getEdges()[se.getEdges().length- 1]!= e))
					continue;
				
				decPos= setCDeconvolution(mode, g, e, cov, decBase, decPos, 
						se, sense, cSumV, cBalV, costIdxV, costValV, covIdxV, covValV);
				
			}
			if (mode== MODE_CONSTRAINT && decPos!= decLen) {
				assert(decPos== 0);
				assert((e.getTail().getSite().isTSS()&& !sense)|| 
						(e.getHead().getSite().isTES()&& sense));
				int lenExon= e.length()- readLenMax;
				assert(lenExon> 0);
				double[] vals= getCostGaussian(obs, lenExon, -1, -1, readLenMax, 0.01, gaussVals);
				double cost= vals[0]* readLenMax/ lenExon;
				double ub= vals[1]* lenExon/ readLenMax;
				while(decPos< decLen) {
					for (int i = 0; i < 2; i++) {
						costIdxV.add(decBase[++decPos]);
						costValV.add(cost);
						getLPsolve().setUpbo(decPos, ub);
					}
					++decPos;
				}
				assert(decPos== decLen);
			}
			if (mode== MODE_COUNT) 
				return;
			
			// expectation restriction
			for (int i = 0; i < cSumV.length; i++) {
				int len= cSumV[i].size();
				if (len== 0)
					continue;
				int[] idx= new int[len+ 1];	// TODO reuse
				idx[0]= sense? ((i* 2)+1): ((i* 2)+2);
				System.arraycopy(cSumV[i].vector, 0, idx, 1, cSumV[i].length);
				double[] val= new double[idx.length];
				Arrays.fill(val, 1d);
				Transcript tx= g.trpts[i];
				int tlen= tx.getExonicLength();
				UniversalMatrix m= profile.getMatrix(tlen);
				double f= m.getFrac(
							tx.getExonicPosition(e.getFrac(true)),
							tx.getExonicPosition(e.getFrac(false)),
							tlen,
							sense?Constants.DIR_FORWARD:Constants.DIR_BACKWARD);
				if (Double.isInfinite(f)|| Double.isNaN(f)) {
					System.err.println("infinite value");
					f= m.getFrac(
							tx.getExonicPosition(e.getFrac(true)),
							tx.getExonicPosition(e.getFrac(false)),
							tlen,
							sense?Constants.DIR_FORWARD:Constants.DIR_BACKWARD);
				}
				ccheck[i]+= f;
				val[0]= -f;
				
				if (mode== MODE_CONSTRAINT) {
					getLPsolve().addConstraintex(idx.length, val, idx, LpSolve.EQ, 0);
					++restrNr;
				} else if (false&& mode== MODE_COMPLETE) {
					printConstraint(System.out, e+ " ("+ (sense?"s":"a")+ ", "+ g.trpts[i].getTranscriptID()+ ")", 
							idx.length, val, idx, LpSolve.EQ, 0d);
				}
			}
			
			// flux preservation, 2 for min/max readlength
			for (int i = 0; (mode!= MODE_COUNT)&& i < covIdxV.length; i++) {
				for (int j = 0; j < covValV[i].length; j++) {	// per tx
					// assure there is an exon and a sj
					if (covIdxV[i][j].length< 2) 	
						continue;				
					boolean plus= false, minus= false;
					for (int k = 0; k < covValV[i][j].length; k++) {
						double d= covValV[i][j].get(k);
						plus|= (d> 0);
						minus|= (d< 0);
						if (plus&& minus)
							break;
					}
					if (!(plus&& minus))
						continue;
					if (mode== MODE_CONSTRAINT) {
						getLPsolve().addConstraintex(covIdxV[i][j].length, covValV[i][j].toDoubleArray(), 
								covIdxV[i][j].toIntArray(), i== 0? LpSolve.LE: LpSolve.GE, 0);
						++restrNr;
					} else if (false&& mode== MODE_COMPLETE) {
						printConstraint(System.out, e+ " (flux)", covIdxV[i][j].length, covValV[i][j].toDoubleArray(), 
								covIdxV[i][j].toIntArray(), i== 0? LpSolve.LE: LpSolve.GE, 0);
					}
				}
			}
			
		}

		private double[] getCostGaussian(int obsExon, int lenExon, int obsSJ, int lenSJ, 
				int readLenHere, double pthold, double[] result) {
			
			if (obsSJ== -1&& lenSJ== -1) {
				double x= Math.sqrt(-2d* Math.log(Math.sqrt(2d* Math.PI)* pthold));
				double y= 1d/ Math.sqrt(2d* Math.PI);
				double slope= (y- pthold)/ x;
				result[0]= slope;
				result[1]= x;
				return result;
			}
			
			// mean and variance
			double covExon= (obsExon* readLenHere)/ (double) lenExon;	// cov_nt, reads homogenized to common length
			double covSJ= (obsSJ* readLenHere)/ (double) lenSJ;
			double u= Math.abs(covExon- covSJ)/ 2d;	// mean
			double sig2= (Math.abs(covExon- u)+ Math.abs(covSJ- u))/ 2d;
			sig2*= sig2;	// variance, square of mean(abs. deviations)
			
			// intercept on threshold(p), value p(u)
			double sqrt2PiSig2= Math.sqrt(2d* Math.PI* sig2);
			double x= Math.sqrt(-2d* sig2* Math.log(Math.sqrt(sqrt2PiSig2* pthold)));
			assert(x>= 0);
			result[1]= x;
			double y= 1d/ sqrt2PiSig2;
			
			// -slope, -delta(y)/ delta(x) 
			double slope= (y- pthold)/ x; 
			assert(slope> 0);
			result[0]= slope;
			return result;
		}
		
		private double[] gaussVals= new double[2];
		private int setCDeconvolution(byte mode, Graph g, Edge f, double cov, 
				int[] decBase, int decPos, Edge e, boolean sense, IntVector[] cSumV, IntVector[] cBalV,
				IntVector costIdxV, DoubleVector costValV, IntVector[][] covIdxV, DoubleVector[][] covValV) 
		throws LpSolveException {
							
					// obs: bounds
					int ctrTT= 0;
					IntVector minusIdxV= null;
					if (mode!= MODE_COUNT) 
						minusIdxV= new IntVector(g.trpts.length);	// TODO reuse
					
					// attributes of e
					int obs= (sense? e.getReadNr(): e.getRevReadNr());
					boolean sj= (e instanceof SuperEdge)&& (!((SuperEdge) e).isPend());
					double effLenMin= calcEffLen(e, f, true);
					double effLenMax= calcEffLen(e, f, false);
					assert(effLenMin> 0|| obs== 0);

					// iterate transcripts
					if (f!= null&& f.toString().equals("47884063-47884268^"))
						System.currentTimeMillis();
					for (int j = 0; j < g.trpts.length; j++) {
						
						if (!g.contains(e.getTranscripts(), j))
							continue;
						if (mode== MODE_COUNT) {
							++ctrTT;
							constraintCtr+= 2;
							++restrNr;
							continue;
						}
						
						int x= ctrTT* 3;
						decIdx[x]= ++constraintCtr;		// c
						cSumV[j].add(constraintCtr);	// Contribution sum for expectation
						//if (se!= null) ???
						boundIdx[0]= constraintCtr;		// delta<= c
						//boundVal[0]= 1d;	// TODO init once
						if ((sj &&effLenMin> 0)|| ((!sj)&& effLenMax> 0)) {
							covIdxV[0][j].add(constraintCtr);
							covValV[0][j].add(1d/ (sj? -effLenMin: effLenMax));	// inverse
						}
						if ((sj &&effLenMax> 0)|| ((!sj)&& effLenMin> 0)) {
							covIdxV[1][j].add(constraintCtr);
							covValV[1][j].add(1d/ (sj? -effLenMax: effLenMin));	// inverse
						}
						
						// costs, symmetrical
						int readLenHere= 0;
						double costSJ= 0d, ubReadsSJ= 0d;
						if (sj&& mode== MODE_CONSTRAINT) {
							int obsExon= sense? f.getReadNr(): f.getRevReadNr();
							readLenHere= mapLenMax;	// ~eff. length of junction
							if (effLenMin< 0)
								readLenHere+= effLenMin;
							assert(readLenHere> 0);
							int lenExon= f.length()- readLenHere;
							assert(lenExon> 0);
							gaussVals= getCostGaussian(
									obsExon, lenExon, obs, readLenHere- 1, readLenHere, 0.01d, gaussVals);
							double cost= gaussVals[0];
							double costUB= gaussVals[1];
							// costs for exon d
							double costEx= (cost* readLenHere)/ lenExon,
									ubReadsEx= (costUB* lenExon/ readLenHere);
//							if (decPos+ 1>= decBase.length)
//								System.currentTimeMillis();
							costIdxV.add(decBase[++decPos]);
							costValV.add(costEx);
							getLPsolve().setUpbo(decPos, ubReadsEx);
							costIdxV.add(decBase[++decPos]);
							costValV.add(costEx);
							getLPsolve().setUpbo(decPos, ubReadsEx);
							++decPos;	// set for next
							// cost for SJ
							costSJ= (cost* readLenHere)/ (readLenHere- 1);
							ubReadsSJ= costUB* (readLenHere- 1)/ readLenHere;
						}
						
						// d-, adding reads to obs
						decIdx[x+ 1]= ++constraintCtr;
						cBalV[j].add(-constraintCtr);
						boundIdx[1]= constraintCtr;
						if (sj) {
							costIdxV.add(constraintCtr);
							costValV.add(costSJ);
							getLPsolve().setUpbo(constraintCtr, ubReadsSJ);
						}
						if (mode== MODE_CONSTRAINT) {	// bound to contribution
							getLPsolve().addConstraintex(2, boundVal, boundIdx, LpSolve.GE, 0);
							++restrNr;
						} if (false&& mode== MODE_COMPLETE) {
							printConstraint(System.out, e+ " ("+(sense?"s, ":"a, ")+"ub d+)", 2, boundVal, boundIdx, LpSolve.GE, 0);
						}
		
						
						// d+, substracting reads from obs 
						decIdx[x+ 2]= ++constraintCtr;	// d+
						if ((mode!= MODE_COUNT)) 
							minusIdxV.add(constraintCtr);	// last man standing				
						cBalV[j].add(constraintCtr);
						if (sj) {
							costIdxV.add(constraintCtr);
							costValV.add(costSJ);
							getLPsolve().setUpbo(constraintCtr, ubReadsSJ);
						}
						
						++ctrTT;
					}	/* end tx */
					
					
					// deconvolution
					if (mode!= MODE_COUNT) {
						if (mode== MODE_CONSTRAINT) {
							getLPsolve().addConstraintex(ctrTT* 3, decVal, decIdx, LpSolve.EQ, obs);
							++restrNr;
						} else if (mode== MODE_COMPLETE)
							printConstraint(System.out, e+ " ("+(sense?"s, ":"a, ")+"dec)", ctrTT* 3, decVal, decIdx, LpSolve.EQ, obs);
							//completeConstraint(System.out, e+ " ("+(sense?"s, ":"a, ")+"dec)", ctrTT* 3, decVal, decIdx, LpSolve.EQ, obs); 
					}
					
					return sj?decPos:(ctrTT* 3);
				}

		private Object[] mateObjects= null;
		private byte[] mateOrient= null;
		private int[] mateInt= null;
		private int mateP;
		private int txStrand= -1;
		private int[] fromTo= new int[2];
		private int nrMappingsLocusMapped= 0;
		private int constraintCtr= 0, restrNr= 0;
	} /* end LocusSolver */

	class Partition implements Comparable {
		long[] partition;
		DoubleVector coV;
		double mean= -1d, sig2= -1d;
		int reads= 0;
		
		public Partition(long[] partition) {
			this.partition= partition;
			coV= new DoubleVector();
		}
		public void add(double coverage) {
			assert(mean< 0&& sig2< 0);
			coV.add(coverage);
		}
		
		public void addReads(int reads) {
			this.reads+= reads;
		}
		public double getMean() {
			if (mean< 0) 
				calcStats();
			return mean;
		}
		public double getSig2() {
			if (sig2< 0) 
				calcStats();
			return sig2;
		}
		
		@Override
		public String toString() {
			if (partition== null)
				return "null";
			StringBuilder sb= new StringBuilder(Long.toString(partition[0]));
			for (int i = 1; i < partition.length; i++) {
				sb.append(",");
				sb.append(Long.toString(partition[i]));
			}
			return sb.toString();
		}
		private void calcStats() {
			
			double[] vals= coV.toDoubleArray();
			// median more robust?
			// + exp edges less weight
			// - lim->0 in lowly covered genes
			// Distribution d= new Distribution(vals);
			// mean= d.getMedian();
			mean = 0d;
			for (int i = 0; i < vals.length; i++) 
				mean+= vals[i];
			mean/= vals.length;
			
			sig2= 0;
			if (vals.length== 1)
				sig2= 1;	// set normal distribution for single observation
			else {	// N> 1
				for (int i = 0; i < vals.length; i++) {
					double x= vals[i]- mean;
					sig2+= x* x;
				}
				sig2/= vals.length- 1;	// TODO wiki says 2* vals.length
			}
		}

		public int compareTo(Object o) {
			if (!(o instanceof Partition)) 
				throw new IllegalArgumentException();
			Partition part= (Partition) o;
			if (partition.length!= part.partition.length)
				return (partition.length- part.partition.length);
			for (int i = 0; i < partition.length; i++) {
				if (partition[i]!= part.partition[i]) {
					long x= (partition[i]- part.partition[i]);
					if (x< 0)
						return -1;	// otherwise int overflow
					else if (x> 0)
						return 1;
				}
			}
			return 0;
		}
		
		@Override
		public boolean equals(Object obj) {
			return (compareTo(obj)== 0);
		}
		
	}
	
	class LocusSolver3 extends Thread {
		
				Gene gene= null;
				ASEvent[] events= null;
				BEDobject2[] beds= null;
				boolean decompose= false;
				int nrMappingsReadsOrPairs;
				long[] sigall= null;
				double[] result= null;
		
				private float invariantTestObsSplitFreq= 0, invariantTestPredSplitFreq= 0;
				
				public LocusSolver3() {
					//super(newGene.getGeneID());
					
				}
				
				void init(Gene newGene, BEDobject2[] newBeds, boolean decompose) {
					this.gene= newGene; 
					this.beds= newBeds;
					this.decompose= decompose;
					
					nrMappingsReadsOrPairs= 0;
					if (pairedEnd) {
						mateObjects= new Object[MAX_MATES];
						mateOrient= new byte[MAX_MATES];
						mateInt= new int[MAX_MATES];
						mateP= 0;
					}
	
				}
				
				
				/**
				 * @deprecated
				 */
				public void run_old() {
						
					try {
						//int mapped= mapReadOrPairIDs.size();
						if (decompose) {
							if (this.gene.getTranscriptCount()== 1) {
								mapTrivial(gene.getTranscripts()[0], beds);
								outputGFF(null, null /*, null*/);
							} else {
								Graph myGraph= getGraph(this.gene);
								map(myGraph, this.gene, this.beds);
								
								GraphLPsolver2 mySolver= null;
								// != mapReadOrPairIDs.size()> 0, does also count singles
								if (nrMappingsReadsOrPairs> 0&& this.gene.getTranscriptCount()> 1) {	// OPTIMIZE			
									mySolver= getSolver(myGraph, nrMappingsReadsOrPairs* 2); // not: getMappedReadcount()
									mySolver.run();
								}
				//				if (this.gene.getTranscriptCount()== 1)
				//					System.currentTimeMillis();
								outputGFF(myGraph, events/*, mySolver*/);
								
							}
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
					FluxCapacitorNew.this.threadPool.remove(this);
		//			}
				}
		
				/**
				 * @deprecated
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
						byte flag=  (byte) 1;// (descriptor.getPairedEndInformation(beds[i].getName())- 1);	// (Fasta.getReadDir(beds[i].getName())- 1))
						if (flag < 0) {
							System.err.println("Error in readID:\n"+ beds[i].getName());
							continue;
						}
						CharSequence id= ""; // descriptor.getUniqueDescriptor(beds[i].getName());	//Fasta.getReadID(beds[i].getName())
						int sep= i, end= i+1;
						for (; end < beds.length; ++end) {
							if (!beds[end].getName().startsWith(id))
								break;
							
							if (sep== i) {
								byte flag2= (byte) 1; // (descriptor.getPairedEndInformation(beds[end].getName())- 1);
								if (flag2< 0) {
									System.err.println("Error in readID:\n"+ beds[i].getName());
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
						if (!(first|| second)) {
							for (int j = i; j < end; j++) 
								beds[j]= null;
						} else {
							//nrMappingsReadsOrPairs+= 2;
						}
						i= end- 1;
					}
					
					// build individual matrix
					int tlen= tx.getExonicLength();
					byte tstrand= tx.getStrand();
					UniversalMatrix m= new UniversalMatrix(Math.max(tlen/ 100, 10));
					for (int i = 0; i < beds.length; i++) {
						if (beds[i]== null)
							continue;
						++nrMappingsReadsOrPairs;
						int p= getBpoint(tx, beds[i]);
						if (p< 0|| p>= tlen) 
							continue;		
						boolean sense= beds[i].getStrand()== tstrand;
						m.add(p, beds[i].getLength(), tlen, sense?Constants.DIR_FORWARD:Constants.DIR_BACKWARD); 
					}
					
					// normalize biases out
					//UniversalMatrix m= profile.getMatrix(tlen);	// no good idea
					if (nrMappingsReadsOrPairs> 0) {
						//better individual matrix, also for saturation biases
						//UniversalMatrix m= profile.getMatrix(tlen);
						
						/*if (tx.getTranscriptID().contains("ENST00000262241")||
								tx.getTranscriptID().contains("ENST00000252818"))
							System.currentTimeMillis();
						*/
						//System.err.println("\n"+ tx.getTranscriptID());
						int w= m.sense.length/ 5;
						m.sums= Kernel.smoothen(Kernel.KERNEL_EPANECHNIKOV, 
								w, m.sense);
						m.suma= Kernel.smoothen(Kernel.KERNEL_EPANECHNIKOV, 
								w, m.asense);
						
						double f= m.getNfactor(0.2d);
						nrMappingsReadsOrPairs*= f;
					}
				}
	
	
				Graph getGraph(Gene gene) {
						boolean output= false;
						
						// construct graph
					long t0= System.currentTimeMillis();
					
					Graph myGraph= new Graph(gene);
					myGraph.createDefaultCoordComparator(mapLenMin);
					myGraph.constructGraph();
					//myGraph.collapseFuzzyFlanks();
					
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
		
				/** 
				 * maps reads
				 * @param myGraph
				 * @param g
				 * @param beds
				 * @return
				 */
				int map(Graph myGraph, Gene g, BEDobject2[] beds) {
					
					
					if (beds== null|| beds.length== 0)
						return 0;
					
					Transcript tx= null;
					UniversalMatrix m= null;
					if (myGraph== null) {
						tx= gene.getTranscripts()[0];
						if (decompose)
							m= new UniversalMatrix(Math.max(tx.getExonicLength()/ 100, 10));
						else
							m= profile.getMatrix(tx.getExonicLength()); 
					}
					CharSequence lastID= null;
					int mapCtr= 0;
					for (int x = 0; x< beds.length; ++x) {
						BEDobject2 dobject= beds[x];
						int mode= -1;
						if (stranded|| pairedEnd) {
							CharSequence name= dobject.getName();
							mode= descriptor.getMode(name, fromTo);
							if (pairedEnd) {
								CharSequence ID= name.subSequence(fromTo[0], fromTo[1]);
								if (pairedEnd&& !ID.equals(lastID)) {
									mateP= 0;
									lastID= ID;
								}
							}
						}
	
						// determine a/sense
						byte rStrand= dobject.getStrand();
						byte orient= (byte) Math.abs(gene.getStrand()- rStrand);
						assert(orient== 0|| orient== 2);
						if (stranded) {	// filter if strand known
							byte sense= descriptor.getStrand(dobject);
							if (sense!= Descriptor.ORIENTATION_UNKNOWN&& sense!= orient) 
								continue;
						} 
						
						// map to annotation
						int bpoint1= -1;
						Edge target= null;
						if (myGraph== null) {
							bpoint1= getBpoint(tx, beds[x]);
							if (bpoint1< 0|| bpoint1> tx.getExonicLength())
								continue;
						} else {
							target= myGraph.getEdge2(dobject);
							if (target== null)
								continue;
						}
						
						// increase counters
						boolean mated= false;
						int mate= descriptor.getPairedEndInformation(dobject); // pairedEnd? descriptor.getMate(mode): Descriptor.MATE_UNKNOWN;
						if (pairedEnd&& mate== Descriptor.MATE_2) {
							assert(mate!= Descriptor.MATE_UNKNOWN);
							
							// try pairing
							// TODO redundancy check (mateP> 1)
							for (int i = 0; i < mateP; i++) {
								
								if (mateOrient[i]== orient) {
									if (gene.getStrand()< 0) {
										if (myGraph== null)
											nrMappingsPairsWrongOrientationCt+= 2;
										else
											nrMappingsPairsWrongOrientationCg+= 2;
									} else {
										if (myGraph== null)
											nrMappingsPairsWrongOrientationWt+= 2;
										else
											nrMappingsPairsWrongOrientationWg+= 2;
									}
									continue;
								}
	
								// map pair
								if (decompose) {
									if (myGraph== null) { // learn or decompose trivial
										BEDobject2 bed2= (BEDobject2) mateObjects[i];
										int bpoint2= getBpoint(tx, bed2);	// TODO save bpoint before
										//if (bpoint2>= 0&& bpoint2< tx.getExonicLength())	// has to be 
										m.add(mateInt[i], bpoint2, 1, 1, tx.getExonicLength());
									} else {
										Edge target2= (Edge) mateObjects[i];
										SuperEdge se= myGraph.getSuperEdge(target, target2, true, null);
										if (se== null) {
											if (gene.getStrand()< 0) {
												if (myGraph== null)
													nrMappingsPairsWoTxEvidenceCt+= 2;
												else
													nrMappingsPairsWoTxEvidenceCg+= 2;
											} else {
												if (myGraph== null)
													nrMappingsPairsWoTxEvidenceWt+= 2;
												else
													nrMappingsPairsWoTxEvidenceWg+= 2;
											}
											continue;	// no tx evidence
										}
										se.incrReadNr();
										
										// remove anterior single mapping
										if (mateInt[i]== 0) {
											if (mateOrient[i]== Descriptor.ORIENTATION_SENSE)
												target2.decrReadNr();
											else
												target2.decrRevReadNr();
										}
									}
									
								} 
								
								// adjust counters
								mated= true;
								mapCtr+= 2;
								if (mateInt[i]== 0) {
									if (gene.getStrand()< 0) {
										if (mateOrient[i]== Descriptor.ORIENTATION_SENSE) {
											if (myGraph== null)
												--nrMappingsSingleCSt;
											else
												--nrMappingsSingleCSg;
										} else {
											if (myGraph== null)
												--nrMappingsSingleCAt;
											else
												--nrMappingsSingleCAg;
										}
									} else {
										if (mateOrient[i]== Descriptor.ORIENTATION_SENSE) {
											if (myGraph== null)
												--nrMappingsSingleWSt;
											else
												--nrMappingsSingleWSg;
										} else {
											if (myGraph== null)
												--nrMappingsSingleWAt;
											else
												--nrMappingsSingleWAg;
										}
									}
								}
								++mateInt[i];
								
							}
						} 
						
						// map single reads
						if (!mated) {
							if (myGraph== null) {
								if (!pairedEnd)
									m.add(bpoint1, 1, tx.getExonicLength(), 
											(orient== Descriptor.ORIENTATION_SENSE?Constants.DIR_FORWARD:Constants.DIR_BACKWARD));
							} else {
								if (orient== Descriptor.ORIENTATION_SENSE)
									target.incrReadNr();
								else
									target.incrRevReadNr();
							}
							
							if (pairedEnd) {
								if (mate== Descriptor.MATE_1) {
									if (myGraph== null)
										mateObjects[mateP]= beds[x];
									else
										mateObjects[mateP]= target;
									mateOrient[mateP]= orient;
									mateInt[mateP]= 0;
									++mateP;
								}
							} else
								++mapCtr;
	
							if (gene.getStrand()< 0) {
								if (orient== Descriptor.ORIENTATION_SENSE) {
									if (myGraph== null)
										++nrMappingsSingleCSt;
									else
										++nrMappingsSingleCSg;
								} else {
									if (myGraph== null)
										++nrMappingsSingleCAt;
									else
										++nrMappingsSingleCAg;
								}
							} else {
								if (orient== Descriptor.ORIENTATION_SENSE) {
									if (myGraph== null)
										++nrMappingsSingleWSt;
									else
										++nrMappingsSingleWSg;
								} else {
									if (myGraph== null)
										++nrMappingsSingleWAt;
									else
										++nrMappingsSingleWAg;
								}
							}
						}
	
					} // iterate beds
					
					return mapCtr;
				}
						
				private int nrLocusMultimaps= 0;
				/**
				 * @deprecated
				 * @param regs
				 */
				int mapRead2(Graph g, BEDobject2 dobject, boolean force) {
					
					HashSet<CharSequence> mapReadOrPairIDs= new HashSet<CharSequence>();
					HashMap<CharSequence, Vector<BEDobject2>[]> mapEndsOfPairs= new HashMap<CharSequence, Vector<BEDobject2>[]>();
	
					// find the edge(s) where the regions align
					// if these do not form a continous chain, create a new edge
					
		//			GFFObject[] gtfs= GFFObject.fromBed(dobject);	// TODO kill intermediary GTFs !!!
		//			DirectedRegion[] regs= new DirectedRegion[gtfs.length];
		//			for (int i = 0; i < regs.length; i++) 
		//				regs[i]= new DirectedRegion(gtfs[i]);
					
					if (force&& mapReadOrPairIDs.contains(dobject.getName())) {
						return 0;
					}
					
					byte flag= (byte) 1; // getFlag(dobject); 
					if (flag <= 0) {
						System.err.println("Error in readID:\n"+ dobject.getName());
						return 0;
					}
					CharSequence ID= ""; // getID(dobject); 	
		
					if (ID.equals("mouse_7_112_1503_1238")|| ID.equals("mouse_7_9_1185_579"))
						System.currentTimeMillis();
					
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
							if (dobject.getStrand()== dobject2.getStrand()) {
								if (gene.getStrand()< 0)
									;//++nrMappingsPairsWrongOrientationC;
								else
									;//++nrMappingsPairsWrongOrientationW;
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
								if (gene.getStrand()< 0)
									;//++nrMappingsPairsWoTxEvidenceC;
								else
									;//++nrMappingsPairsWoTxEvidenceW;
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
		
						if (sense|| (strand!= STRAND_ENABLED)) {
							target.incrReadNr();
							mapCtr= 1;
						} else if (strand!= STRAND_SPECIFIC) {
							target.incrRevReadNr();
							mapCtr= 1;
						} else {
							if (gene.getStrand()< 0) {
								//++nrMappingsWrongStrandC;
							} else {
								//++nrMappingsWrongStrandW;
							}
							mapCtr= 0;
						}
						
						
						
						if (!mapReadOrPairIDs.add(dobject.getName()))
							++nrLocusMultimaps;
		//				if (outputMapped)
		//					writeMappedRead(dobject);
						return mapCtr;
					}
					
				}
				
				public int getMappedReadcount() {
					HashSet<CharSequence> mapReadOrPairIDs= new HashSet<CharSequence>();
	
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
				
				
				private void outputSimple() {
					StringBuffer sb= new StringBuffer();
					double sum= nrMappingsLocusMapped;
					if (trptExprHash!= null) {
						sum= 0;
						for (int i = 0; i < gene.getTranscripts().length; i++) 
							sum+= trptExprHash.get(gene.getTranscripts()[i].getTranscriptID());
					}
					for (int i = 0; i < gene.getTranscripts().length; i++) {
						sb.append(gene.getGeneID());
						sb.append("\t");
						String tid= gene.getTranscripts()[i].getTranscriptID();
						sb.append(tid);
						sb.append("\t");
						double val= nrMappingsLocusMapped;
						if (trptExprHash!= null)
							val= trptExprHash.get(gene.getTranscripts()[i].getTranscriptID());
						val/= (gene.getTranscripts()[i].getExonicLength()/ 1000d);	// RPK normalize
						sb.append(val);
						sb.append("\t");
						if (sum== 0)
							sb.append("0.0");
						else
							sb.append(val/ sum);
						sb.append("\n");
					}
					try {
						getWriter().write(sb.toString());
					} catch (IOException e) {
						e.printStackTrace();
					}
				}
				
				/**
				 * @deprecated
				 * @param g
				 * @param events
				 * @param solver
				 */
				private void outputGFF(Graph g, ASEvent[] events) {
					
					++nrLoci;
					if (result!= null) 
						++nrLociExp;
					
					double perM= nrReadsAll/ 1000000d;
					// deprecated
					boolean unsolvedSystem= false;	
					double valOF= result== null?0: result[0];
					if (valOF> BIG) { 
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
							if (result== null) {
								// no reads
								val= nrMappingsLocusMapped;
							} else {
		//						if (solver.getTrptExprHash().get(g.trpts[i].getTranscriptID())!= 0)
		//							System.currentTimeMillis();
								val= (float) (trptExprHash.get(tx.getTranscriptID()).doubleValue());
								if (val< 1- costBounds[0]) // 1- 0.95
									val= 0;
								try {
									assert(tt.length> 1|| val== nrMappingsReadsOrPairs);
								} catch (AssertionError e) {
									System.err.println(val+ " x "+ nrMappingsReadsOrPairs);
									getNFactor();
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
									sb.append(GTF_ATTRIBUTE_TOKEN_RPKM);
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
									sb.append(GTF_ATTRIBUTE_TOKEN_RPKM);
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
									
									sb.append(GTF_ATTRIBUTE_TOKEN_RPKM);
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
	/*					if (result!= null) {
							Object[] keys= tExpMap.keySet().toArray();	
							for (int i = 0; i < keys.length; i++) {
								if (!(keys[i] instanceof String))
									continue;
								if (tExpMap.get(keys[i])<0)
									tExpMap.put((String) keys[i], 0d);	// TODO ugly
							}
						}
	*/					
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
					if (false&& result!= null) {
						ArrayList<Edge> cc= new ArrayList<Edge>();
						if (result!= null) {
							Iterator<Object> iter= null; //solver.getConstraintHash().keySet().iterator();
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
		
				private GraphLPsolver2 getSolver(Graph g, int mappedReads) {
				
					GraphLPsolver2 solver= new GraphLPsolver2(g, readLenMin, 
							pairedEnd?insertMinMax:null, mappedReads, 
							(strand== STRAND_ENABLED), 
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
								sb.append(GFF_FEATURE_PAIRED);
							else 
								sb.append(GFF_FEATURE_JUNCTION);
						} else
							sb.append(GFF_FEATURE_FRAGMENT);
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
								sb.append(GTF_ATTRIBUTE_LENGTH);
								sb.append(GTF_ATTRIBUTE_TOKEN_SEP);
								sb.append(x== 0? GTF_ATTRIBUTE_TOKEN_ALL: (x==1? GTF_ATTRIBUTE_TOKEN_TID: GTF_ATTRIBUTE_TOKEN_EXC));
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
							
							ReadStatCalculator calc= (i== 0? FluxCapacitorNew.this: solver);
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
										sb.append(i==0? GTF_ATTRIBUTE_TOKEN_OBSV: (i==1? GTF_ATTRIBUTE_TOKEN_PRED:GTF_ATTRIBUTE_TOKEN_BALANCED));
										sb.append(GTF_ATTRIBUTE_TOKEN_SEP);
										sb.append(j== 0?GTF_ATTRIBUTE_TOKEN_ALL:(j==1? GTF_ATTRIBUTE_TOKEN_TID: GTF_ATTRIBUTE_TOKEN_EXC));
										sb.append(GTF_ATTRIBUTE_TOKEN_SEP);				
										sb.append(x== 0?GTF_ATTRIBUTE_TOKEN_READS: GTF_ATTRIBUTE_TOKEN_COV);	// pairedEnd?GTF_ATTRIBUTE_TOKEN_COV:GTF_ATTRIBUTE_TOKEN_RPKM
										sb.append(" ");	// no \"
									}
									for (int m = 0; m < tid.length; m++) {
										long[] sig= j== 0? sigall: tid[m];
				
										if (calc== null) {
											if (output)
												sb.append(VALUE_NA);
										} else {
											
											boolean normalized= i== 2? true: false;
											double val= ms* (j== 0? calc.getReads(eeV.elementAt(m), BYTE_0, sig, normalized):  // ALL								
												calc.getReadsAvg(eeV.elementAt(m), BYTE_0, g, sig, excl, normalized));	// TID, excl 							
											
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
													sb.append(length== 0? FLOAT_STRING_0: Float.toString((float) val));
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
				 * @deprecated
				 * @param g
				 * @param events
				 * @param solver
				 */
				private void outputGFF_save(Graph g, ASEvent[] events, GraphLPsolver solver) {
							++nrLoci;
							if (solver!= null) 
								++nrLociExp;
							double perM= nrReadsAll/ 1000000d;
							// deprecated
								boolean unsolvedSystem= false;	
								double valOF= solver== null?0: solver.getValObjFunc();
								if (valOF> BIG) { 
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
		
				
				/**
				 * @param tx
				 * @param bed
				 * @return
				 */
				private boolean contains(Transcript tx, BEDobject2 bed) {
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
				
				
				/**
				 * @deprecated
				 * @param tx
				 * @param beds
				 */
				private void learn(Transcript tx, BEDobject2[] beds) {
								
					if (beds== null|| beds.length== 0)
						return;
					
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
						if (strand== STRAND_SPECIFIC&& bed1.getStrand()!= tx.getStrand()) {
							//++nrMappingsWrongStrand;
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
						if (rlen1> readLenMax)	// readLenMax< 0|| 
							readLenMax= rlen1;
						
						++nrReadsSingleLociMapped;
						if (pairedEnd) {	// && flag== Descriptor.MATE_2
	
							CharSequence name= bed1.getName();
							int mode= descriptor.getMode(name, fromTo);
							if (mode <0) {
								System.err.println("Error in readID:\n"+ bed1.getName());
								continue;
							}
							byte flag= descriptor.getPairedEndInformation(bed1); // getMate(mode);
							CharSequence id= name.subSequence(fromTo[0], fromTo[1]);
							
							
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
							for (int j = 0; j < oo.length; j++) {
								if (j==flag|| oo[j]== null)
									continue;
								for (int k = 0; k < oo[j].length; k++) {
									BEDobject2 bed2= oo[j][k];
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
										++nrReadsSingleLociPairsMapped;	
									}
								}
							}
								
						} else {
							
							System.err.println("[DEACTIVATED] no single reads for now");
							if (1== 1)
								System.exit(0);
							int[] ePos= null;
							int binLen= 0;
							TProfile tprofile= null;
							
							byte dir= Constants.DIR_FORWARD;
							if (strand== STRAND_ENABLED) {
								//assert(ePos[0]> 0);	// can be, out of range
								if (!tx.getTranscriptID().startsWith("AT1G04050")) {
									if (ePos[0]> 0&& ePos[1]<= elen) {
										if (beds[i].getStrand()!= tx.getStrand()) {
											dir= Constants.DIR_BACKWARD;
											if (ePos[1]<= elen) {
												//++profileStubRev[binLen][p];
											}
										} else {
											if (ePos[0]> 0) {
												//++profileStub[binLen][p];
											}
										}
									}
								}
							} else if (!uniform)
								tprofile.addRead(ePos[0],readLenMin,dir);
						}
					}
					
				}
	
	
				/**
				 * @deprecated
				 * @param regs
				 */
				int mapRead(Graph g, BEDobject2 dobject, boolean force) { 
					
					HashSet<CharSequence> mapReadOrPairIDs= new HashSet<CharSequence>();
					HashMap<CharSequence, Vector<BEDobject2>[]> mapEndsOfPairs= new HashMap<CharSequence, Vector<BEDobject2>[]>();
	
					// find the edge(s) where the regions align
					// if these do not form a continous chain, create a new edge
					
		//			GFFObject[] gtfs= GFFObject.fromBed(dobject);	// TODO kill intermediary GTFs !!!
		//			DirectedRegion[] regs= new DirectedRegion[gtfs.length];
		//			for (int i = 0; i < regs.length; i++) 
		//				regs[i]= new DirectedRegion(gtfs[i]);
					
					if (force&& mapReadOrPairIDs.contains(dobject.getName())) {
						return 0;
					}
					
					byte flag= 1; // getFlag(dobject);  	
					CharSequence ID= ""; // getID(dobject); 	
		
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
		
						if (sense|| (strand!= STRAND_ENABLED)) {
							target.incrReadNr();
							mapCtr= 1;
						} else if (strand!= STRAND_SPECIFIC) {
							target.incrRevReadNr();
							mapCtr= 1;
						} else {
							//++nrMappingsWrongStrand;
							mapCtr= 0;
						}
						
						
						
						if (!mapReadOrPairIDs.add(dobject.getName()))
							++nrLocusMultimaps;
		//				if (outputMapped)
		//					writeMappedRead(dobject);
						return mapCtr;
					}
					
				}
	
				public void run() {
								
					try {
						
//						if (!gene.getGeneID().equals("chr14:23685936-23706451W")) 
//							return;
						
						// init
						trptExprHash= null;
						
						// data structure
						Graph myGraph= null;
						if (decompose&& this.gene.getTranscriptCount()> 1) {
							myGraph= getGraph(this.gene);
							if (outputEvent)
								events= myGraph.getEvents(2);
							int nrSJ= 0;
							if (!pairedEnd)
								nrSJ= myGraph.addEJ(27);						
						}
						
						if (beds!= null&& beds.length> 0) {
							
							// mapping
							if (pairedEnd) 
								Arrays.sort(beds, BEDobject2.DEFAULT_ID_COMPARATOR);	// TODO file
							nrMappingsLocusMapped= 0;
							nrMappingsLocusMapped= map(myGraph, this.gene, this.beds);
							// DEBUG
							if (myGraph== null) {
								Graph gg= getGraph(this.gene);
								int chk= map(myGraph, this.gene, this.beds);
								if (chk!= nrMappingsLocusMapped)
									System.currentTimeMillis();
							}
							if (gene.getStrand()< 0) {
								if (myGraph== null)
									nrMappingsPairsMappedCt+= nrMappingsLocusMapped;
								else
									nrMappingsPairsMappedCg+= nrMappingsLocusMapped;
							} else {
								if (myGraph== null)
									nrMappingsPairsMappedWt+= nrMappingsLocusMapped;
								else
									nrMappingsPairsMappedWg+= nrMappingsLocusMapped;
							}
							
							// solving
							if (nrMappingsLocusMapped> 0&& this.gene.getTranscriptCount()> 1) {  
								solve(myGraph);
							}
						}
						
						// output
						//outputGFF(myGraph, null);
						outputSimple();
						
					} catch (Throwable e) {
						System.err.println("\n[ERROR] in locus "+ gene.getGeneID());
						e.printStackTrace();
						System.err.flush();
					}
					
					beds= null;
					gene= null;
				}
	
			/**
			 * #(Edge,Transcript) x (int[],int[][])
			 * watch out with hash function of edges
			 * @return
			 */
			private HashMap<String, Double> mapCCheck= null; 
			double[] ccheck;
			boolean stopperBalance= false;
			boolean stopperBalance2= false;
			private HashMap<Object,Double> trptExprHash;
				
			public strictfp int solve(Graph g) {
						
				//		if (costModel== COSTS_LINEAR&& costBounds== null)
				//			costBounds= new int[] {100,100};
						
						long t0= System.currentTimeMillis();
						debug= false;	
						// chr12:67919663-67951290W (CPSF6)
						// chr1:1205831-1217267W (5 NUMERIC) ok
						// chr1:154318993-154376504W (2 INFEASABLE) ok
						// chr1:8843648-8978713C (5 NUMERIC) ok
						// chr1:1374932-1459930W (2 INFEASABLE) ok
						// chr1:44978079-45006026W (KIF2C)
						// chr6:30775563-30793645C (MDC1)
						// chr2:47863725-47887588W (MSH6)
						// chr1:46541890-46555035W (funny gene with 100:10 s/a mappings in yaspo)
//						if (gene.getGeneID().equals("chr6:30775563-30793645C")) {
//							debug= true;
//							System.err.println();
//						} else 
//							return -1;
						
						collectStats(g);
						
						int cCtrCount= 0, rCtrCount= 0;
						try {
							setCAll(MODE_COUNT, g);
							cCtrCount= constraintCtr;
							rCtrCount= restrNr;
						} catch (LpSolveException e2) {
							e2.printStackTrace();
							return -1;
						}
						//int resSave= restrNr;	// doesnt work properly
						
	//					decVal= new double[g.trpts.length+ 2]; 
	//					Arrays.fill(decVal, 0, decVal.length- 1, 1d);
	//					decVal[decVal.length- 1]= -1d;
						decVal= new double[g.trpts.length* 3];
						for (int i = 0; i < decVal.length; i++) 
							decVal[i]= ((i-1)% 3== 0? -1d: 1d);	// +-+ +-+ +-+
						decIdx= new int[decVal.length];
						boundIdx= new int[2];
						boundVal= new double[2];
						boundVal[0]= 1d;
						boundVal[1]= -1d;
						
						getLPsolve(cCtrCount);
							
				/*		if (g.trpts[0].getTranscriptID().startsWith("NM_001168507")||
								g.trpts[0].getTranscriptID().startsWith("NM_173453"))
							debug= true;
				*/			
						
				
						// set up program
						ccheck= new double[g.trpts.length];
						Arrays.fill(ccheck, 0d);
						try {
							constraintCtr= 0;
							restrNr= 0;
							setCAll(MODE_CONSTRAINT, g);	// CONSTRAINT
							if (constraintCtr!= cCtrCount|| restrNr!= rCtrCount) {
								System.err.println("\n\tERROR constraint: Cs complete "+ constraintCtr+ " (<>"+ cCtrCount+"),\n\t" +
										"Rs complete "+ restrNr+ " (<>"+ rCtrCount+ ")");
								System.exit(-1);
							}
							//assert(resNow== resSave);
						} catch (LpSolveException e1) {
							e1.printStackTrace();
							return -1;
						}
						
						// write out
						String tmpOutFName= null;
						int ret= 0;
						if (fileLPdir!= null) {
							tmpOutFName= fileLPdir+ File.separator
								+ gene.getGeneID().replace(":", "_")
								+ SFX_LPOUT;
				
							try {
								getLPsolve().setOutputfile(tmpOutFName);
							} catch (LpSolveException e) {
								//e.printStackTrace();
								if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
									System.err.println("[FATAL] failed to set lp output to:\n\t"+ tmpOutFName);
							}
						}
						
					
				//		for (int i = 0; i < g.trpts.length; i++) {
				//			checkPercent(g.trpts[i]);
				//		}
						
						try {
							//getLPsolve().printLp();
							if (tmpOutFName== null)
								getLPsolve().setVerbose(LpSolve.CRITICAL);	//shut up ! IMPORTANT, SEVERE, CRITICAL
							t0= System.currentTimeMillis();
							//getLPsolve().setPresolve(arg0, arg1);
							//System.err.println("scaling "+ gene.getGeneID()+ " ("+constraintCtr+","+restrNr+")");
							//getLPsolve().setScaling(LpSolve.SCALE_GEOMETRIC| LpSolve.SCALE_POWER2);
							
							getLPsolve().setScaling(LpSolve.SCALE_GEOMETRIC+ LpSolve.SCALE_EQUILIBRATE + LpSolve.SCALE_INTEGERS);
							getLPsolve().setImprove(LpSolve.IMPROVE_BBSIMPLEX);
							
							ret= getLPsolve().solve();
							//System.err.println("solved "+ gene.getGeneID());
				/*			 NOMEMORY (-2)  	Out of memory
							   OPTIMAL (0) 	An optimal solution was obtained
							SUBOPTIMAL (1) 	The model is sub-optimal. Only happens if there are integer variables and there is already an integer solution found. The solution is not guaranteed the most optimal one.
				
							 							* A timeout occured (set via set_timeout or with the -timeout option in lp_solve)
							 							* set_break_at_first was called so that the first found integer solution is found (-f option in lp_solve)
							 							* set_break_at_value was called so that when integer solution is found that is better than the specified value that it stops (-o option in lp_solve)
							 							* set_mip_gap was called (-g/-ga/-gr options in lp_solve) to specify a MIP gap
							 							* An abort function is installed (put_abortfunc) and this function returned TRUE
							 							* At some point not enough memory could not be allocated 
				
							INFEASIBLE (2) 	The model is infeasible
							UNBOUNDED (3) 	The model is unbounded
							DEGENERATE (4) 	The model is degenerative
							NUMFAILURE (5) 	Numerical failure encountered
							USERABORT (6) 	The abort routine returned TRUE. See put_abortfunc
							TIMEOUT (7) 	A timeout occurred. A timeout was set via set_timeout
							PRESOLVED (9) 	The model could be solved by presolve. This can only happen if presolve is active via set_presolve
							PROCFAIL (10) 	The B&B routine failed
							PROCBREAK (11) 	The B&B was stopped because of a break-at-first (see set_break_at_first) or a break-at-value (see set_break_at_value)
							FEASFOUND (12) 	A feasible B&B solution was found
							NOFEASFOUND (13) 	No feasible B&B solution found
				*/
							if (ret!= 0) {
								//System.err.println("insolvable system: "+ gene.getGeneID()+ "\t("+ g.trpts.length+ " transcripts)");
								//System.err.println(RETURN_VERBOSE[ret+ 2]);
								return 0;
								//  debug= true;
							} else {
//								System.err.println("solved "+ gene.getGeneID());
								getResult();
								for (int i = 0; i < g.trpts.length; ++i) {
									double sense= result[restrNr+ (i* 4)+ 1];
									double asense= result[restrNr+ (i* 4)+ 2];									
									double realQ= sense/ asense;
									if (sense== 0&& asense== 0)
										realQ= 1;
									int tlen= g.trpts[i].getExonicLength();
									UniversalMatrix m= profile.getMatrix(tlen);
									double expQ= m.sums/ (double) m.suma;
									if (Math.abs(realQ- expQ)> 1&& (sense> 10|| asense> 10))
										System.err.println("\n"+gene.getGeneID() 
												+"\n\tWARNING: S/A imbalance "+ realQ+ " ("+ expQ+") with "+ sense+ ","+ asense);
/*									if (expQ!= 1d)
										System.currentTimeMillis();
									System.err.print(realQ+ ","+ expQ+ "\t");
*/									
								}
//								System.err.println();
							}
							
						} catch (LpSolveException e) {
							e.printStackTrace();
						}
				
						// get transcription expression levels		
						trptExprHash= getResult();	// tMap
						//normalizeBack2LocusExpr(trptExprHash);
						debug= true;
						if (debug|| tmpOutFName!= null) {
							getLPsolve().printLp();
							//restrNr= 0;	// Noo, needed for value retrieval
							try {
								constraintCtr= 0;
								restrNr= 0;
								setCAll(MODE_COMPLETE, g);
								if (constraintCtr!= cCtrCount|| restrNr!= rCtrCount) {
									System.err.println("\n\tERROR complete: Cs complete "+ constraintCtr+ " (<>"+ cCtrCount+"),\n\t" +
											"Rs complete "+ restrNr+ " (<>"+ rCtrCount+ ")");
									System.exit(-1);
								}
							} catch (LpSolveException e) {
								e.printStackTrace();
							}
							getLPsolve().printObjective();
							//getLPsolve().printSolution(1);
							printSolution(System.out, g, result, costIdx, costVal);
							System.exit(-1);
						}
							
						getLPsolve().deleteLp();	// closes file outFName
						
						PrintStream p= null;
						if (debug) {
							System.err.flush(); // doesnt work
							p= System.out;
							for (int i = 0; i < g.trpts.length; i++) {
								; //getAllPercent(g.trpts[i]);
							}
						}
						if (tmpOutFName!= null) {
							try {
								p= new PrintStream(new FileOutputStream(tmpOutFName, true));				
							} catch (Exception e) {
								e.printStackTrace();
							}
						}
						
						// output additionally
						if (p!= null) {
							Iterator<Object> idIter= trptExprHash.keySet().iterator();
							while(idIter.hasNext()) {
								Object o= idIter.next();
								if (!(o instanceof String))
									continue;
								String id= (String) o;
								p.println(id+" "+trptExprHash.get(id));
							}
							p.println();
							p.println("Settings:");
							p.println("paired-end\t"+pairedEnd);
							if (costBounds!= null) {
								if (!Double.isNaN(costBounds[0]))
									p.print("cost boundaries:\tlower /"+costBounds[0]);
								if (!Double.isNaN(costBounds[1]))
									p.print(" upper *"+costBounds[1]);
								p.println();
							}
							p.print("costfunc\t");
							//printConstraintHash(p);
							
							p.println();
							p.println("Transcripts:");
							for (int i = 0; i < g.trpts.length; i++) {
								p.print(g.trpts[i].getTranscriptID()+"\t");
								SpliceSite[] ss= g.trpts[i].getSpliceSitesAll();
								for (int j = 0; j < ss.length; j++) 
									p.print(ss[j].toString());
								p.println();
							}
						}
						
						// close file
						if (p!= null)
							p.flush();
						if (tmpOutFName!= null) 
							p.close();
						
				//		System.err.println("solved "+g.trpts[0].getTranscriptID()+": "+g.trpts.length+" trpts, "+constraintCtr+" constr, "+restrNr+" restr"
				//				+ ((System.currentTimeMillis()- t0)/ 1000)+ " sec.");
						if (debug)
							System.currentTimeMillis();
						
						return ret;
					}
	
			float medExon= 0, medJunc= 0;
			long sumSense= 0, sumAsense= 0;
			private void collectStats(Graph g) {
				sumSense= 0; sumAsense= 0;
				Edge[] edges= g.getExonicEdgesInGenomicOrder();
				IntVector vecExons= new IntVector(), vecJunc= new IntVector();
				for (int i = 0; i < edges.length; i++) {
					Edge e= edges[i];
					if (!e.isExonic())
						continue;
					
					int exSense= e.getReadNr(), exAsense= e.getRevReadNr();
					for (int j = 0; e.getSuperEdges()!= null&& j < e.getSuperEdges().size(); j++) {
						SuperEdge se= (SuperEdge) e.getSuperEdges().elementAt(j);
						int juncSense= 0, juncAsense= 0;
						if (se.getEdges()[0]== e) {
							if (se.isPend())
								exSense+= se.getReadNr();
							else
								juncSense+= se.getReadNr();
						} else if (se.getEdges()[se.getEdges().length- 1]== e){
							if (se.isPend())
								exAsense+= se.getReadNr();
							else
								juncAsense+= se.getRevReadNr();
						}
						
						if (!se.isPend()) {
							for (int k = 0; se.getSuperEdges()!= null&& k < se.getSuperEdges().size(); k++) {
								SuperEdge sse= (SuperEdge) se.getSuperEdges().elementAt(k);
								assert(sse.isPend());
								if (sse.getEdges()[0]== se)
									juncSense+= sse.getReadNr();
								else	// must be last edge
									juncAsense+= sse.getReadNr();
							}
						}
						vecJunc.add(juncSense);
						vecJunc.add(juncAsense);
						sumSense+= juncSense;
						sumAsense+= juncAsense;
					}
					vecExons.add(exSense);
					vecExons.add(exAsense);
					sumSense+= exSense;
					sumAsense+= exAsense;
				}
				Distribution d= new Distribution(vecExons.toIntArray());
				medExon= (float) d.getMedian();
				if (medExon== 0)
					++medExon;
				d= new Distribution(vecJunc.toIntArray());
				medJunc= (float) d.getMedian();
				if (medJunc== 0)
					medJunc= 1;
			}
	
			private void printSolution(PrintStream p, Graph g, double[] result, int[] costIdx, double[] costVal) {
				
				StringBuilder sb= new StringBuilder("Variable\tValue\tFactor\tCosts:\n");
				double sumCosts= 0d, sumChanges= 0d;
				for (int i = restrNr+ 1; i < result.length; i++) {
					sb.append("C");
					sb.append(Integer.toString(i- restrNr));
					sb.append("\t");
					sb.append(StringUtils.fprint(result[i], 2));
					sumChanges+= result[i];
					sb.append("\t");
					int pos= Arrays.binarySearch(costIdx, (i- restrNr));
					if (pos< 0) {
						sb.append(".\t.\n");
						continue;
					}
					sb.append(StringUtils.fprint(costVal[pos], 2));
					sb.append("\t");
					double costs= costVal[pos]* result[i];
					sumCosts+= costs;
					sb.append(StringUtils.fprint(costVal[pos] * result[i], 2));
					sb.append("\n");
				}
				sb.append("OF\t");
				sb.append(StringUtils.fprint(sumChanges, 2));
				sb.append("\t.\t");
				sb.append(StringUtils.fprint(sumCosts, 2));
				sb.append("\n");
				p.println(sb.toString());
			}
	
			LpSolve lpSolve;
			LpSolve getLPsolve() {
					if (lpSolve == null) {
	/*					if (writeFile)
							try {
								getLPWriter().flush();
								getLPWriter().close();
								lpSolve= LpSolve.readLp(fileLPinput.getAbsolutePath(), LpSolve.IMPORTANT, null);
							} catch (Exception e) {
								e.printStackTrace();
							}
						else
	*/						
							try {
			//				Iterator iter= getConstraintHash().keySet().iterator();
			//				int edgeCtr= 0, varCtr= 0;
			//				while (iter.hasNext()) {
			//					Object o= iter.next();
			//					if (o instanceof Edge)
			//						++edgeCtr;
			//					else if (o instanceof Variation)
			//						++varCtr;
			//				}
			//				
								lpSolve = LpSolve.makeLp(0, constraintCtr);	// no 0-column here
								
							} catch (LpSolveException e) {				
								e.printStackTrace();
							}
						
					}
			
					return lpSolve;
				}
			
			LpSolve getLPsolve(int nrC) {
				
				if (lpSolve == null) {
						try {
							lpSolve = LpSolve.makeLp(0, nrC);	// no 0-column here
						} catch (LpSolveException e) {				
							e.printStackTrace();
						}
					
				}
		
				return lpSolve;
			}
			
			LpSolve getLPsolve(int nrR, int nrC) {
				
				if (lpSolve == null) {
						try {
							lpSolve = LpSolve.makeLp(nrR, nrC);	// no 0-column here
						} catch (LpSolveException e) {				
							e.printStackTrace();
						}
					
				}
		
				return lpSolve;
			}
	
			/**
			 * create an array where certain indices are set to certain values
			 * 
			 * @param constIdx indices of constraints
			 * @param constVal values of constraints
			 * @return
			 */
			private double[] createArray(int[] constIdx, double[] constVal) {
				
				double[] a= new double[constraintCtr+ 1];
				for (int i = 0; i < a.length; i++) 
					a[i]= 0d;
				for (int i = 0; i < constIdx.length; i++) 
					a[constIdx[i]]= constVal[i];		
				return a;
			}
	
			private HashMap<Object, Double> getResult() {
				
	//			int rows= getLPsolve().getNrows();
	//			int cols= getLPsolve().getNcolumns();
				result= new double[1+ restrNr+ constraintCtr];
				try {
					getLPsolve().getPrimalSolution(result);
				} catch (LpSolveException e1) {
					e1.printStackTrace();
				}
				//valObjFunc= result[0];	
				HashMap<Object,Double> trptExprHash= new HashMap<Object,Double>();
				
				// normalize fraction rounding errors
				Transcript[] trpts= gene.getTranscripts();
				double sum= 0;
				for (int i = 0; i < trpts.length; i++) { 
					
					// (sense+ asense)/ 2
					double tot= (result[1+ restrNr+ (i* 2)]+ result[1+ restrNr+ (i* 2)+ 1])/ 2d;
					// normalize by summing inaccuracy
					tot/= ccheck[i]/ 2d;
					sum+= tot;
					if (Double.isNaN(tot))
						System.currentTimeMillis();
					trptExprHash.put(trpts[i].getTranscriptID(), tot);
				}
				if (sum== 0)
					return trptExprHash;
				
				// normalize flow network over-/underprediction
				double nfac= nrMappingsLocusMapped/ (double) sum;
				for (int i = 0; sum!= 0&& i < trpts.length; i++) {
					double x= trptExprHash.get(trpts[i].getTranscriptID());
					x*= nfac;
					if (Double.isNaN(x))
						System.currentTimeMillis();
					trptExprHash.put(trpts[i].getTranscriptID(), x);
				}
				
				
				// normalize transcript profile
				for (int i = 0; i < trpts.length; i++) {
					int tlen= trpts[i].getExonicLength();
					UniversalMatrix m= profile.getMatrix(tlen);
					double f= m.getNfactor(0.2d);
					double x= trptExprHash.get(trpts[i].getTranscriptID());
					x*= f;
					if (Double.isNaN(x))
						System.currentTimeMillis();
					trptExprHash.put(trpts[i].getTranscriptID(), x);
				}
				
				// normalize locus??
				
				return trptExprHash;
			}
	
			private void addConstraintToLp(int[] idx, double[] val, int eq, double cap) {
					
			//		if (writeFile)
			//			writeRowLP(idx, val, eq, cap);
			//		else if (columnWise) {
			//			;
			//		} else 
						//double[] constraint= createArray(idx, val);
						try {
							getLPsolve().addConstraintex(idx.length, val, idx, eq, cap);					
							//getLPsolve().addConstraint(constraint, eq, cap);
						} catch (LpSolveException e) {
							e.printStackTrace();
						}
			
				}
	
			private void getConstraintsNew(Edge e, long[] sig, IntVector v,
					HashMap<Edge, IntVector> mapE, boolean sense, boolean count) {
				
				// for the edge itself
				IntVector w= mapE.get(e);
				if (count) 
					++constraintCtr;
				else {
					v.add(++constraintCtr);	// for transcript fraction
					if (w== null)
						w= new IntVector();
					w.add(constraintCtr);	// edge consistency
				}
				mapE.put(e, w);
			
				// iterate super-edges
				for (int j = 0; e.getSuperEdges()!= null&& j < e.getSuperEdges().size(); j++) {
					SuperEdge se= e.getSuperEdges().elementAt(j);
					if (Graph.isNull(Graph.intersect(se.getTranscripts(), sig)))
						continue;
					// sense/anti-sense.. e must be first/last in super-edge
					if ((sense&& se.getEdges()[0]!= e)|| ((!sense)&& se.getEdges()[se.getEdges().length- 1]!= e))
						continue;
					if (count)
						++constraintCtr;
					else
						v.add(++constraintCtr);	// for transcript fraction
					w= mapE.get(se);
					if (!count) {
						if (w== null)
							w= new IntVector();
						w.add(constraintCtr); // for edge consistency
					}
					mapE.put(se, w);
					
					if (se.isPend())
						continue;	// no super-edges
					
					for (int k = 0; se.getSuperEdges()!= null&& k < se.getSuperEdges().size(); k++) {
						SuperEdge se2= se.getSuperEdges().elementAt(k);
						assert(se2.isPend());
						if (Graph.isNull(Graph.intersect(se2.getTranscripts(), sig)))
							continue;
						// sense/anti-sense.. e must be first/last in super-edge
						if ((sense&& se2.getEdges()[0]!= se)|| ((!sense)&& se2.getEdges()[se2.getEdges().length- 1]!= se))
							continue;
						if (count)
							++constraintCtr;
						else 
							v.add(++constraintCtr);	// tx
						w= mapE.get(se2);
						if (!count) {
							if (w== null)
								w= new IntVector();
							w.add(constraintCtr); // for edge consistency
						}
						mapE.put(se2, w);
					}
				}
				
				
			}
	
			private double nFactor= Double.NaN;
			public double getNFactor() {
				if (Double.isNaN(nFactor)|| true) {
					double fictReads= 0;
					Iterator iter= trptExprHash.keySet().iterator();
					while (iter.hasNext()) {
						Object o= iter.next();
						if (!(o instanceof String))
							continue;
						fictReads+= trptExprHash.get(o);
					}
			
					if (fictReads< nrMappingsLocusMapped)
						++nrLociUnderPredicted;
					else
						++nrLociOverPredicted;
					
					// can happen, when reads are where none expected
	//				if (fictReads== 0^ nrMappingsObs== 0)
	//					System.currentTimeMillis();
					
					// avoid large scaling; was 0, avoid NaN; 0.5 too large
					if (fictReads> 0.000001) {	
						nFactor= nrMappingsLocusMapped/ fictReads;
					} else 
						nFactor= 1d;
				}	
				return nFactor;
			}
	
			public HashMap<String, Integer> setConstraints(Graph g, boolean count) {
						
						//boolean debug= false;
				//		if (g.trpts[0].getTranscriptID().equals("ENST00000323441")) {
				//			debug= true;
				//		}
						
						// transcript constraint variables
						Transcript[] trpts= g.trpts;
						IntVector v= null, w= null;
						DoubleVector u= null;
						HashMap<String, Integer> tMap= null;
						if (count)
							constraintCtr+= trpts.length;
						else { 
							try {
								getLPsolve().setLpName(trpts[0].getTranscriptID());
							} catch (LpSolveException e1) {
								e1.printStackTrace();
							}			
							tMap= new HashMap<String, Integer>(trpts.length* 2);
							for (int i = 0; i < trpts.length; i++) { 
								tMap.put(trpts[i].getTranscriptID(), ++constraintCtr);
								if (debug)
									System.out.println("C"+constraintCtr+"\t"+trpts[i].getTranscriptID());
							}
							v= new IntVector();	// indices for transcript/part 
							w= new IntVector();	// indices for cost function
							u= new DoubleVector();	// observation, bases for cost function
						}
						
						// iterate edges
						Edge[] edges= g.getExonicEdgesInGenomicOrder();
						if (!count) {
							mapCCheck= new HashMap<String, Double>();
							for (int i = 0; i < trpts.length; i++) 
								mapCCheck.put(trpts[i].getTranscriptID(), 0d);
						}
						for (int i = 0; i < edges.length; i++) {
							Edge e= edges[i];
							if (!e.isExonic())
								continue;
							
							// the base edge
							int basElen= e.length();
							Transcript[] tt= g.decodeTset(e.getTranscripts());
							HashMap<Edge, IntVector> mapE= new HashMap<Edge, IntVector>();			
							// sense/anti
							for (int sa= 0; sa< 2; ++sa) {
								for (int x = 0; x < tt.length; ++x) {
									if (!count)
										v.removeAll();
									long[] sig= g.encodeTset(tt[x]);
			//						if (e.toString().equals("67937779-67937952^")&& !count)
			//							System.currentTimeMillis();
									getConstraints(e, sig, v, mapE, sa== 0, count);
									
									// add transcript constraint
									if (!count) {
										int[] idx= new int[v.length+ 1]; // obs parts+ tx frac
										System.arraycopy(v.vector, 0, idx, 0, v.length);
										idx[idx.length- 1]= tMap.get(tt[x].getTranscriptID());
										double[] val= new double[idx.length];
										Arrays.fill(val, 1d);
										int tlen= tt[x].getExonicLength();
										UniversalMatrix m= profile.getMatrix(tlen);
										double f= m.getFrac(
													tt[x].getExonicPosition(e.getFrac(true)),
													tt[x].getExonicPosition(e.getFrac(false)),
													tlen,
													sa== 0?Constants.DIR_FORWARD:Constants.DIR_BACKWARD);
										//System.err.println(f);
										if (Double.isInfinite(f)|| Double.isNaN(f)) {
											System.err.println("infinite value");
											f= m.getFrac(
													tt[x].getExonicPosition(e.getFrac(true)),
													tt[x].getExonicPosition(e.getFrac(false)),
													tlen,
													sa== 0?Constants.DIR_FORWARD:Constants.DIR_BACKWARD);
										}
										mapCCheck.put(tt[x].getTranscriptID(),
												mapCCheck.get(tt[x].getTranscriptID())+ f);
										val[val.length- 1]= -f;
										if (debug&& !count) {
											StringBuilder sb= new StringBuilder(e.toString());
											sb.append(": ");
											for (int k = 0; k < idx.length; k++) {
												sb.append(val[k]>0?"+":"");
												sb.append(val[k]%1==0?((int) val[k]):val[k]);
												sb.append("C");
												sb.append(idx[k]+" ");
											}
											sb.append("= 0");
											System.out.println(sb);
										}
										addConstraintToLp(idx, val, LpSolve.EQ, 0);
										++restrNr;
										
										// force junction coverage
										if (val.length== 3) {
											assert(val[0]<0&&val[1]>0&&val[2]>0);
											int fixedreadlen= 27;
											double frac= 0.9;
											val[0]= 0d;
											val[1]= frac/ (basElen- fixedreadlen+ 1);
											val[2]= -1d/ (fixedreadlen- 1);
											addConstraintToLp(idx, val, LpSolve.LE, 0);
											++restrNr;
										}							
									}
								}
								// add edge constraints
								Edge[] ee= new Edge[mapE.size()];
								mapE.keySet().toArray(ee);
								
								// total obs
				/*				int sumObs= 0;
								for (int j = 0; j < ee.length; j++) {
									Edge f= ee[j];
									boolean paird= (f instanceof SuperEdge)&& ((SuperEdge) f).isPend();
									int nr= ((paird|| sa== 0)? f.getReadNr(): f.getRevReadNr());
									sumObs+= nr;
								}
				*/				
								for (int j = 0; j < ee.length; j++) {
									Edge f= ee[j];
									
									// paired-end super-edge
									boolean paird= (f instanceof SuperEdge)&& ((SuperEdge) f).isPend();
									int nr= ((paird|| sa== 0)? f.getReadNr(): f.getRevReadNr());
									// coverage or cost calc
			/*						int len= (f instanceof SuperEdge)? 27- 1: f.length()- 27+ 1;
									double cov= 0;
									if (len<= 0) {
										assert(nr== 0);
									} else
										cov= nr/ (double) len;
									double cost= 1;
									if (len> 0) {
										cost= 1d+ (10* Math.exp(-nr/ 10));
										//cost= 1d+ (10* Math.pow((nr+ 1)/10d, -4));
									}						
									if (f instanceof SuperEdge) {
										cost= 1d+ (10* Math.exp(-nr/10d));
									} else {
										if (len> 0)
											cost= 1d+ (Math.pow(cov+ 1, -1));
									}
									if (Double.isNaN(cost)|| Double.isInfinite(cost))
										System.currentTimeMillis();
			*/							
									v= mapE.remove(f);
									if (count)
										constraintCtr+= 2;
									else {
										int[] idx= new int[v.length+ 2];	// +/-
										System.arraycopy(v.vector, 0, idx, 0, v.length);
										int c= ++constraintCtr;
										// plus on not paired edges at 0-cost, it substracts from obs
										// usual edges cost when in single read mode
										if ((!pairedEnd)|| (pairedEnd&& paird)) {
											w.add(c);
			//								if (nr== 0)
			//									u.add(Double.NEGATIVE_INFINITY);
			//								else
												u.add(-nr);
										}
										idx[idx.length- 2]= c;
										// plus has to be limited, it substracts
										int lim= (paird||(!pairedEnd))? Math.max(nr- 1, 0): nr; // Math.max(nr- 1, 0): nr; 
										try {
											getLPsolve().setUpbo(constraintCtr, lim);
										} catch (LpSolveException e1) {
											e1.printStackTrace();
										}
											
										c= ++constraintCtr;
										// do not limit adding, even with f= 100 unsolvable systems
										// adding reads always costs, also on single edges
										w.add(c);
										u.add(nr);
										idx[idx.length- 1]= c;
										double[] val= new double[idx.length];
										Arrays.fill(val, 1d);
										val[val.length- 1]= -1d;
										if (debug&& !count) {
											StringBuilder sb= new StringBuilder(f.toString());
											sb.append(": ");
											for (int k = 0; k < idx.length; k++) {
												sb.append(val[k]==1?"+C":"-C");
												sb.append(idx[k]+" ");
											}
											sb.append("= "+nr);
											System.out.println(sb);
										}
										addConstraintToLp(idx, val, LpSolve.EQ, nr);
										++restrNr;
									}
								}
							}
						} // end all edges
						
						if (count) {
							if (stopperBalance|| stopperBalance2)
								constraintCtr+= 2;
							return null;
						}
						
						if (stopperBalance2) {
							int[] idx= new int[tMap.size()+ 2];
							Iterator<Integer> iter= tMap.values().iterator();
							int ctr= 0;
							while (iter.hasNext()) {
								idx[ctr++]= iter.next();
							}
							idx[ctr++]= ++constraintCtr;
							idx[ctr++]= ++constraintCtr;
							double[] vals= new double[idx.length];
							Arrays.fill(vals, 1d);
							vals[vals.length- 1]= -1d;
							addConstraintToLp(idx, vals, LpSolve.EQ, nrMappingsLocusMapped);
							++restrNr;
						}
						
						// set objective function/costs
				//		double[] a= createArray(w.toIntArray());	// linear costs
						//double min= Math.exp(-nrMappingsLocusMapped);
						int ext= 0;
						if (stopperBalance)
							ext= 2;
						double[] costs= new double[w.size()+ ext];
						double[] pmPfx= new double[w.size()+ ext];
						double[] dd= u.toDoubleArray();
						for (int i = 0; i < dd.length; i++) {
							dd[i]= Math.abs(dd[i]);
						}
						Distribution dist= new Distribution(dd);
						double med= dist.getMedian();
						if (med== 0)
							med= 1;
						StringBuilder sb= null;
						if (debug) {
							sb= new StringBuilder();
						}
						for (int i = 0; i < costs.length- ext; i++) {
							double x= u.get(i);	// <0 : read substracting
							pmPfx[i]= ((x< 0|| x== Double.NEGATIVE_INFINITY)?-1d:1d);
							//costs[i]= 1+ Math.exp(-Math.log(1+ x));	// neutralizes
							//costs[i]= 1d+ Math.pow(x+1d, -1d/2d);
							//costs[i]= 0;	// 
							//costs[i]= 1d+ Math.log(x+ 1);
							//costs[i]= 1d+ Math.exp(-x);
							//costs[i]= 100d/ (x+ 1);	// last
							
							double c= 1d;
							if (!pairedEnd) {
								double y= (x== Double.NEGATIVE_INFINITY)?0:x;
								c= (10d* Math.exp(-(Math.abs(y)/ med)));
								if (x< 0/*|| x== Double.NEGATIVE_INFINITY*/)
									c+= 1;
							}
							//c= med/ (x+ 1);
							//c= 10d/ (1+ Math.pow(x, 1));
							//c= 10d/ (1+ x);
							//c= x;
							if (Double.isInfinite(c)|| Double.isNaN(c))
								System.currentTimeMillis();
							costs[i]= c;
							
							
							if (debug) { 
								sb.append("C");
								sb.append(w.get(i));
								sb.append("\t");
								sb.append(StringUtils.fprint(c, 3));
								sb.append("\n");
							}
							//costs[i]= 1d+ Math.log(x+ 1);
							//costs[i]= 1d+ Math.sqrt(x+ 1);	// best
							//costs[i]= 1d+ Math.pow(x+ 1, 1d/ 3d);
							
							//costs[i]= (Math.log(x+ 1)/ (x+ 1));	// logdiv
							//costs[i]= 1d/ (1d+ Math.log(x+ 1d));	// divlog
						}
						if (debug) {
							System.out.println("COSTS:");
							System.out.println(sb);
						}
						int[] idx= null;
						if (stopperBalance) {
							idx= new int[w.length+ ext];
							pmPfx[pmPfx.length- 2]= 1;
							pmPfx[pmPfx.length- 1]= -1;
							System.arraycopy(w.vector, 0, idx, 0, w.length);
							idx[pmPfx.length- 2]= ++constraintCtr;
							idx[pmPfx.length- 1]= ++constraintCtr;
							addConstraintToLp(idx, pmPfx, LpSolve.EQ, 0);
							++restrNr;
							
							costs[costs.length- 2]= 1;
							costs[costs.length- 1]= 1;
						} else {
							if (stopperBalance2) {
								idx= new int[w.length+ 2];
								double[] costs2= new double[costs.length+ 2];
								System.arraycopy(w.vector, 0, idx, 0, w.length);
								idx[idx.length- 2]= constraintCtr- 2;
								idx[idx.length- 1]= constraintCtr- 1;
								System.arraycopy(costs, 0, costs2, 0, costs.length);
								costs2[costs2.length- 2]= 1d;
								costs2[costs2.length- 1]= 1d;
								costs= costs2;
							} else
								idx= w.toIntArray();
						}
						
						double[] a= createArray(idx, costs);
						try {
							getLPsolve().setObjFn(a);
							getLPsolve().setMinim();
						} catch (LpSolveException e) {
							e.printStackTrace();
						}
						
						// TODO consistency check
						Object[] oo= mapCCheck.keySet().toArray();
						for (int i = 0; i < oo.length; i++) {
							double val= mapCCheck.get(oo[i]);
							//System.err.println("check "+ val);
							if (Math.abs(2d- val)> 0.2d)
								System.err.println("Fraction inconsistency "+ oo[i]+"\t"+val);
						}
						
						return tMap;
					}
	
			private void getConstraints(Edge e, long[] sig, IntVector v,
					HashMap<Edge, IntVector> mapE, boolean sense, boolean count) {
				
				// for the base-edge itself
				IntVector w= mapE.get(e);
				if (count) 
					++constraintCtr;
				else {
					v.add(++constraintCtr);	// for transcript fraction
					if (w== null)
						w= new IntVector();
					w.add(constraintCtr);	// edge consistency
				}
				mapE.put(e, w);
			
				// iterate super-edges
				for (int j = 0; e.getSuperEdges()!= null&& j < e.getSuperEdges().size(); j++) {
					SuperEdge se= e.getSuperEdges().elementAt(j);
					if (Graph.isNull(Graph.intersect(se.getTranscripts(), sig)))
						continue;
					// sense/anti-sense.. e must be first/last in super-edge
					if ((sense&& se.getEdges()[0]!= e)|| ((!sense)&& se.getEdges()[se.getEdges().length- 1]!= e))
						continue;
					if (count)
						++constraintCtr;
					else
						v.add(++constraintCtr);	// for transcript fraction
					w= mapE.get(se);
					if (!count) {
						if (w== null)
							w= new IntVector();
						w.add(constraintCtr); // for edge consistency
					}
					mapE.put(se, w);
					
					if (se.isPend())
						continue;	// no super-edges
					
					for (int k = 0; se.getSuperEdges()!= null&& k < se.getSuperEdges().size(); k++) {
						SuperEdge se2= se.getSuperEdges().elementAt(k);
						assert(se2.isPend());
						if (Graph.isNull(Graph.intersect(se2.getTranscripts(), sig)))
							continue;
						// sense/anti-sense.. e must be first/last in super-edge
						if ((sense&& se2.getEdges()[0]!= se)|| ((!sense)&& se2.getEdges()[se2.getEdges().length- 1]!= se))
							continue;
						if (count)
							++constraintCtr;
						else 
							v.add(++constraintCtr);	// tx
						w= mapE.get(se2);
						if (!count) {
							if (w== null)
								w= new IntVector();
							w.add(constraintCtr); // for edge consistency
						}
						mapE.put(se2, w);
					}
				}
				
				
			}
	
			int[] cBalTrpt= null;
			int[] boundIdx, decIdx;
			double[] boundVal, decVal;
			int[] costIdx;
			double[] costVal;
			private void printConstraint(PrintStream p, String id, int length, double[] val, int[] idx,
					int eq, double d) {
				StringBuilder sb= new StringBuilder(id);
				sb.append("\t");
				for (int i = 0; i < length; i++) {
					sb.append(val[i]>= 0?"+":"-");
					sb.append(Math.abs(val[i])== 1?"": Float.toString((float) Math.abs(val[i])));
					sb.append("C");
					sb.append(idx[i]);
					sb.append(" ");
				}
				sb.append(eq== LpSolve.EQ?"=":(eq== LpSolve.GE?">=":"<="));
				sb.append(" ");
				sb.append(d);
				p.println(sb.toString());
			}
	
			int mapLenMin= 27, mapLenMax= 27; // should be minReadLen
			IntVector idxDecV; 
			DoubleVector valDecV;
			private double calcEffLen(Edge e, Edge f, boolean min) {
				
				double effLen= -1d;
				if (e instanceof SuperEdge) {
					SuperEdge se= (SuperEdge) e;
					Edge[] see= se.getEdges();
					if (!se.isPend()) {	// sj
						int interNt= 1;
						for (int i = 1; i < se.getEdges().length- 1; i++) 
							interNt+= se.getEdges()[i].length();
						effLen= (min? 
									Math.min(Math.min(mapLenMin- interNt, see[0].length()),
										see[see.length- 1].length()):
									Math.min(Math.min(mapLenMax- interNt, see[0].length())
										, see[see.length- 1].length()));
					} else {
						//effLenMin= effLenMax= f.length(); // not: - fixedReadLen+ 1;
						int base= f.length();
						effLen= (min? base- mapLenMax+ 1: base- mapLenMin+ 1);
					}
				} else {
					// 0effLenMin= effLenMax= e.length(); // not: - fixedReadLen+ 1;
					int base= e.length();
					effLen= (min? base- mapLenMax+ 1: base- mapLenMin+ 1);
				}
	
				return effLen;
			}
	
			private int setConstraints(boolean count, int obs, IntVector idxDecV, DoubleVector valDecV, 
					IntVector costIdxVec, DoubleVector costValVec) throws LpSolveException {
				
				int constraintSave= constraintCtr, maxMinus= obs+ Math.max(obs- 1, 0);
	//			long fac= FactorialPrimeSchoenhage.factorial(obs).intValue();
				int lastPlus= -1, lastMinus= -1;
				for (int i = obs+ 1;; i++) {
	//				double exp= exp(-i);
	//				double pow= pow(i, obs);
	//				double cost= fac/(pow* exp);
	//				if (cost< 0)
	//					System.currentTimeMillis();
					int delta= i- obs;
					double obs1= Math.max(1, obs);
					double cost= exp(-delta/ obs);
					if (cost> 1000&& i> obs+ 1)
						break;
					// d-, adding reads to obs				
					if (count)
						++constraintCtr;
					else {
						getLPsolve().setUpbo(++constraintCtr, 1d);
						idxDecV.add(constraintCtr);
						valDecV.add(-1d);
						costIdxVec.add(constraintCtr);
						costValVec.add(cost);
						lastPlus= constraintCtr;
					}
					
					// d+, substracting reads from obs
					if(i<= maxMinus) {
						if (count) 
							++constraintCtr;
						else {
							getLPsolve().setUpbo(++constraintCtr, 1d);
							idxDecV.add(constraintCtr);
							valDecV.add(1d);
							costIdxVec.add(constraintCtr);
							costValVec.add(cost);
							lastMinus= constraintCtr;
						}
					}
				}
				if (!count)
					getLPsolve().setUpbo(lastPlus, getLPsolve().getInfinite());
				
				return (constraintCtr- constraintSave);
			}
			
			private double e= Math.E;
			private double exp(int i) {
				
				if (i== 0)
					return e;
				int j= i>= 0? i: -i;
				double val= 1d;
				for (int x = 0; x < j; x++) 
					val*= e; 
				if (i< 0)
					val= 1d/ val;
				
				return val;
			}
	
			private long pow(int base, int exp) {
				assert(exp>= 0);
				if (exp== 0)
					return 1;
				if (exp== 1)
					return base;
				long val= 1l;
				for (int i = 0; i < exp; i++) 
					val*= base;
				return val;
			}
	
			private double getPoissonProbability(int exp, int obs) {
				// (exp^obs* e^(-exp))/ fact(obs)
				int fac= 0; //FactorialPrimeSchoenhage.factorial(obs).intValue();
				double p= Math.pow(exp, obs)* Math.exp(-exp)/ fac;
				
				return p;
			}
			
			
	
	
			private double getCostObs(int obs) {
				double c= Math.exp(-obs);
				return c;
			}
			
			private double getCostObs(int obs, double len) {
				if (len<= 0)
					return 1; 
				double cov= obs/ (double) len;
				double c= Math.exp(-cov);
				//c*= 100;
				return c;
			}
			
			private double getCostObs(boolean sj, int obs, double len) {
				if (len<= 0)
					return 2;
				double c= 1d+ Math.exp(-obs/ (sj?medJunc:medExon));
				//double c= Double.MAX_VALUE;	// force infeasable
	/*			double c= 0;
				if (sj)
					c= Math.log(obs+ 1d);
				else 
					c= Math.max(0, -Math.log((obs+ 1)/ (double) len));
				if (c< 0.000001&& c> 0)
					System.currentTimeMillis();
	*/
				
				return c;
			}
	
			private HashMap<Object, Double> getResult(HashMap<String, Integer> tMap) {
				
				result= new double[1+ restrNr+ constraintCtr];
				try {
					getLPsolve().getPrimalSolution(result);
				} catch (LpSolveException e1) {
					e1.printStackTrace();
				}
				//valObjFunc= result[0];	
				HashMap<Object,Double> trptExprHash= new HashMap<Object,Double>();
				
				// normalize fraction rounding errors
				Transcript[] trpts= gene.getTranscripts();
				double sum= 0;
				for (int i = 0; i < trpts.length; i++) { 
					int c= tMap.get(trpts[i].getTranscriptID());
					if (Double.isNaN(c))
						System.currentTimeMillis();
					double x= result[restrNr+ c];
					double tot= mapCCheck.get(trpts[i].getTranscriptID());
					tot/= 2d;
					x/= tot;
					sum+= x;
					if (Double.isNaN(x))
						System.currentTimeMillis();
					trptExprHash.put(trpts[i].getTranscriptID(), x);
				}
				if (sum== 0)
					return trptExprHash;
				
				// normalize flow network over-/underprediction
				double nfac= nrMappingsLocusMapped/ (double) sum;
				for (int i = 0; sum!= 0&& i < trpts.length; i++) {
					double x= trptExprHash.get(trpts[i].getTranscriptID());
					x*= nfac;
					if (Double.isNaN(x))
						System.currentTimeMillis();
					trptExprHash.put(trpts[i].getTranscriptID(), x);
				}
				
				
				// normalize transcript profile
				for (int i = 0; i < trpts.length; i++) {
					int tlen= trpts[i].getExonicLength();
					UniversalMatrix m= profile.getMatrix(tlen);
					double f= m.getNfactor(0.2d);
					double x= trptExprHash.get(trpts[i].getTranscriptID());
					x*= f;
					if (Double.isNaN(x))
						System.currentTimeMillis();
					trptExprHash.put(trpts[i].getTranscriptID(), x);
				}
				
				// normalize locus??
				
				return trptExprHash;
			}
	
			private void setCostsGaussian(byte mode, boolean sj, int obs, double effLen, int length, int[] decIdx,
					IntVector costIdxV, DoubleVector costValV) throws LpSolveException {
				
				
			}
			
			private double[] gauss_costs= new double[] 
			                                         //{2.5d, 4d, 20d, 100d};
			                                         {1, 10, 100, 1000, 10000, 100000};
			private void setCostsPseudoGaussian(byte mode, boolean sj, int obs, double effLen, int length, int[] decIdx,
					IntVector costIdxV, DoubleVector costValV) throws LpSolveException {
				
	//			if (1== 1)
	//				return;
				
				int base= Math.max(obs, 1);
				//double baseCov= base/ (double) Math.max(effLen, 1d);
				//double gaussFac= (1d/ Math.sqrt(2* Math.PI* baseCov));
				double fac= //Math.max(effLen, 1d)/ obs;
					1d+ Math.exp(-obs/ (sj?medJunc:medExon));
				
				int indices= (length/ 3);
				int[] idx= null;
				double[] val= null;
				if (mode!= MODE_COUNT) {
					idx= new int[indices+ gauss_costs.length];
					for (int i = 0; i < indices; i++) 
						idx[i]= decIdx[(i* 3)+ 1];			
					val= new double[idx.length];
					Arrays.fill(val, 0, indices, 1d);
					Arrays.fill(val, indices, val.length, -1d);
				}
				
				// costs d-, adding to obs
				//double deltaObs= base;
				double deltaObs= base/ (double) gauss_costs.length;
				for (int i = 0; i < gauss_costs.length; ++i) {
					if (mode== MODE_COUNT) {
						++constraintCtr;
						continue;
					}
					idx[indices+ i]= ++constraintCtr;
					costIdxV.add(constraintCtr);				
					if (i< gauss_costs.length- 1)
						getLPsolve().setUpbo(constraintCtr, deltaObs);
					costValV.add(fac* gauss_costs[i]);
				}
				if (mode== MODE_CONSTRAINT) {
					getLPsolve().addConstraintex(idx.length, val, idx, LpSolve.EQ, 0);
					++restrNr;
				} else if (mode== MODE_COMPLETE&& false) {
					printConstraint(System.out, "cost", idx.length, val, idx, LpSolve.EQ, 0d);
				}
				
				// costs d+, substracting from obs
				for (int i = 0; mode!= MODE_COUNT&& i < indices; i++) 
					idx[i]= decIdx[(i* 3)+ 2];			
				deltaObs= base/ (double) gauss_costs.length;
				for (int i = 0; i < gauss_costs.length; ++i) {
					if (mode== MODE_COUNT) {
						++constraintCtr;
						continue;
					}
					idx[indices+ i]= ++constraintCtr;
					costIdxV.add(constraintCtr);				
					getLPsolve().setUpbo(constraintCtr, deltaObs);
					costValV.add(fac* gauss_costs[i]);
				}			
				if (mode== MODE_CONSTRAINT) {
					getLPsolve().addConstraintex(idx.length, val, idx, LpSolve.EQ, 0);
					++restrNr;
				} else if (mode== MODE_COMPLETE&& false) {
					printConstraint(System.out, "cost", idx.length, val, idx, LpSolve.EQ, 0d);
				}
			}
	
			final double COSTS_INFINITY= 100000;
			private int setCAll(byte mode, Graph g) throws LpSolveException {
				
				boolean balActive= true, costActive= true;
				
				IntVector costIdxV= null, locIdxV= null;
				DoubleVector costValV= null;
				IntVector[][] cBalV= null; // sense/asense, trpts[]
				// always init, for correct C/R counting
				partPos= 0;
				Partition[] part= getPartitions(g);
				assert(partPos> 0);

				IntVector[] cSumV= new IntVector[g.trpts.length];
				for (int i = 0; i < cSumV.length; i++) 
					cSumV[i]= new IntVector(5);	// se per tx in e
				IntVector[] fluxIdxV= new IntVector[g.trpts.length];	// min/maxCov, trpts[]
				DoubleVector[] fluxValV= new DoubleVector[g.trpts.length];
				for (int i = 0; i < fluxValV.length; i++) {
					fluxIdxV[i]= new IntVector();
					fluxValV[i]= new DoubleVector();
				}
				Edge[] edges= g.getExonicEdgesInGenomicOrder();
				if (mode!= MODE_COUNT) {
					locIdxV= new IntVector();
					costIdxV= new IntVector(2* 2* g.trpts.length* edges.length);
					costValV= new DoubleVector(costIdxV.vector.length);					

					cBalV= new IntVector[2][partPos];
					for (int i = 0; i < cBalV.length; i++) 
						for (int j = 0; j < cBalV[i].length; j++) 
							cBalV[i][j]= new IntVector(2* edges.length);	// 2 (+/-)* |E|
				}
			
				// iterate transcripts
				double obsSratio= sumSense/ (double) (sumSense+ sumAsense);
//				if (gene.getGeneID().equals("chr1:46541890-46555035W")&& mode== MODE_CONSTRAINT)
//					System.currentTimeMillis();
				for (int i = 0; i < g.trpts.length; i++) {					
					if (mode== MODE_COUNT)
						constraintCtr+= 4;	// sense, anti-sense
					else if (mode== MODE_COMPLETE) {
						++constraintCtr;
						System.out.println("(s,"+ g.trpts[i].getTranscriptID()+ ")\tC"+ constraintCtr);
						++constraintCtr;
						System.out.println("(a,"+ g.trpts[i].getTranscriptID()+ ")\tC"+ constraintCtr);
						++constraintCtr;
						System.out.println("(a,"+ g.trpts[i].getTranscriptID()+ ")\tC"+ constraintCtr);
						++constraintCtr;
						System.out.println("(a,"+ g.trpts[i].getTranscriptID()+ ")\tC"+ constraintCtr);
					} else if (mode== MODE_CONSTRAINT)	{
						
						Transcript tx= g.trpts[i];
						int tlen= tx.getExonicLength();
						UniversalMatrix m= profile.getMatrix(tlen);
						double ratio= m.sums/ (double) (m.suma+ m.sums);
						double sig2= obsSratio- ratio;
						sig2*= sig2;	// square
						//sig2/= 2;		// mean, for variance
						double ub= 1, cost= COSTS_INFINITY; 
						if (sig2> 0) {
							ub= getGaussUB(sig2, 0.01);
							cost= getGaussCost(sig2, ub, 0.01); 	// 10d
						}
						cost= 10000;
						
						int[] idx= new int[4];
						idx[0]= ++constraintCtr;
						idx[1]= ++constraintCtr;
						idx[2]= ++constraintCtr;
						idx[3]= ++constraintCtr;
						costIdxV.add(idx[2]);
						costValV.add(cost);
						costIdxV.add(idx[3]);
						costValV.add(cost);
						
						double[] val= new double[idx.length];
						val[0]= -1d;
						val[1]= m.sums/ (double) m.suma;
						val[2]= 1d;
						val[3]= -1d;
						getLPsolve().addConstraintex(idx.length, val, idx, LpSolve.EQ, 0d);
					}
					++restrNr;
				}
				
				// iterate edges
				for (int i = 0; i < edges.length; i++) {
					Edge e= edges[i];
					if (!e.isExonic())
						continue;
					
					for (int sa = 0; sa < 2; sa++) {
						boolean sense= (sa== 0);
						setCBase(mode, g, e, sense, 
								locIdxV, cSumV, (mode== MODE_COUNT)?null:cBalV[sa], costIdxV, costValV, fluxIdxV, fluxValV);
					}
				}
					
				// partition balances
				if (mode== MODE_COUNT) {
					//constraintCtr+= 2* 2* partPos;	// sense/asense, +/-
					if (balActive)
						restrNr+= 2* partPos;
				}
				IntVector allLenV= new IntVector(g.trpts.length);
				for (int i = 0; (mode!= MODE_COUNT)&& i < cBalV.length; i++) {
					int tlen= g.trpts[i].getExonicLength();
					allLenV.add(tlen);
				}
				Distribution d= new Distribution(allLenV.toIntArray());
				double medLen= d.getMedian();
				for (int i = 0; (mode!= MODE_COUNT)&& balActive&& i < cBalV.length; i++) {
					for (int j = 0; j < cBalV[i].length; j++) {
						IntVector v= cBalV[i][j];
						int len= v.size();
						int[] idx= new int[len+ 2];				// TODO reuse
						double[] val= new double[idx.length];	// TODO reuse
						for (int k = 0; k < val.length; k++) {
							idx[k]= Math.abs(v.get(k));
							val[k]= (cBalV[i][j].get(k)> 0? 1d: -1d);	// alternate
						}
	
						//
						//double ub= part[j].reads/ 2d;
						//double cost= getPoisCost(ub, part[i].reads, 0.01);
						
						//double sig2= 1d; // part[i].reads/ 2d;
						//double ub= getGaussUB(sig2, 0.01);
						//double cost= 1d/ getGaussCost(sig2, ub, 0.01);
						
						double cost= 100;
						idx[len]= ++constraintCtr;
						val[len]= 1d;
						costIdxV.add(constraintCtr);
						costValV.add(cost);
						idx[len+ 1]= ++constraintCtr;
						val[len+ 1]= -1d;
						costIdxV.add(constraintCtr);
						costValV.add(cost);
						
						if (mode== MODE_CONSTRAINT) {
							getLPsolve().addConstraintex(idx.length, val, idx, LpSolve.EQ, 0d);
						} else if (mode== MODE_COMPLETE) {
							printConstraint(System.out, "part ["+ part[j]+ "] (bal)", 
									idx.length, val, idx, LpSolve.EQ, 0d);
						}
						++restrNr;
					}
				}
			
				if (mode== MODE_COUNT) {
					constraintCtr+= partPos* 2* 2;
					return constraintCtr;
				}
					
				// OF
				double[] a= new double[constraintCtr+ 1];	// 1-based
				Arrays.fill(a, 0d);
				for (int i = 0; i < costIdxV.size(); ++i)  {
					double val= 0d;
					if (!costActive) {
						costValV.set(i, val);
					} else {
						val= costValV.get(i);
						val= val== 0? 1000: val;
					}
					a[costIdxV.get(i)]= val;	
				}
				getLPsolve().setObjFn(a);
				getLPsolve().setMinim();
				costIdx= costIdxV.toIntArray();		// TODO
				costVal= costValV.toDoubleArray();	// TODO
				
				return constraintCtr;
			}
	
			/**
			 * @deprecated
			 * @param ub
			 * @param reads
			 * @param d
			 * @return
			 */
			private double getPoisCost(double ub, int reads, double d) {

				// TODO from de.luschny factorial package:
/*				Xint xxx= new Xint(new Apint(1, 10));
				Xint fact0= null; // FactorialPrimeSchoenhage.factorial(reads);
				xxx.divide(fact0);
				Xint pow= new Xint(new Apint((long) Math.pow(reads, reads)));
				xxx.multiply(pow);
				double p0= xxx.doubleValue();
				double exp= Math.exp(-reads);
				p0*= exp;
				
				int ubInt= (int) Math.round(reads+ ub);
				xxx= new Xint(new Apint(1));
				Xint fact1= null; // FactorialPrimeSchoenhage.factorial(ubInt);
				xxx.divide(fact1);
				pow= new Xint(new Apint((long) Math.pow(reads, ubInt)));
				double p1= xxx.doubleValue();
				p1*= exp;
				
				// slope= -1/ (dy/dx)
				double dy= p0- p1;
				if (dy<= 0)
					System.currentTimeMillis();
				assert(dy> 0);
				double slope= -1d/ (dy/ ub);
				return slope;
*/
				return 0d;
			}

			/**
			 * @deprecated 
			 * @param obsN
			 * @param d
			 * @return
			 */
			private double getPoisUB(int obsN, double d) {
				
				// difficult, got to invert factorial
				return 0d;
			}

			Partition[] part;
			int partPos;
			// TODO from de.luschny factorial package:
			//FactorialPrimeSchoenhage fps= new FactorialPrimeSchoenhage();
			private Partition[] getPartitions(Graph g) {
				
				double[][] res= new double[g.trpts.length][];
				for (int i = 0; i < res.length; i++) 	
					res[i]= new double[2]; 	// u, sig2
				
				// calc max partition number
//				if (g.trpts[0].getTranscriptID().equals("ENST00000379407"))
//					System.currentTimeMillis();
				int min, max, maxPart, txCount= g.trpts.length;
				maxPart= txCount+ 1;	// all partition
				long multiplier, factorial;
				for (int i = 2; i <= txCount- 1; i++) {
					max= (txCount- i);
					if (i< max) 
						min= i;
					else {
						min= max; max= i;
					}
					factorial= 0; // TODO from de.luschny factorial package: fps.factorial(min).longValue();
					multiplier= txCount;
					for (int j = (max+1); j < txCount; j++) 
						multiplier*= j;
					maxPart+= multiplier/ factorial;
					if (maxPart< 0)
						break;	// TODO overflow, change to observed part
				}
				if (maxPart< 0|| maxPart> 1000)
					maxPart= 1000;
				part= new Partition[maxPart];
				
				// iterate edges
				Edge[] ee= g.getExonicEdgesInGenomicOrder();				
				for (int i = 0; i < ee.length; i++) {
					Edge e= ee[i];
					if (!e.isExonic())
						continue;

					//UniversalMatrix m= profile.getMatrix(tlen);
					Partition px= new Partition(e.getTranscripts());
					int p= Constants.binarySearch(part, 0, partPos, px);
					if (p>= 0)
						px= part[p];
					else {
						p= -(p+ 1);
						if (p< partPos)
							System.arraycopy(part, p, part, p+1, partPos- p);
						++partPos;
						part[p]= px;
					}
					Transcript tx= g.getAnyTranscript(e.getTranscripts());
					
					for (int x = 0; x < 2; x++) {
						boolean sense= (x== 0);
						double ntCov= e.getNtCoverage(tx, 
								sense? Constants.DIR_FORWARD: Constants.DIR_BACKWARD, mapLenMax);
						int reads= sense? e.getReadNr(): e.getRevReadNr();
						if (Double.isNaN(ntCov))
							assert(reads== 0); // System.currentTimeMillis();
						else {
							px.add(ntCov);
							px.addReads(reads);
						}
						for (int j= 0; e.getSuperEdges()!= null&& j< e.getSuperEdges().size(); ++j) {
							
							SuperEdge se= e.getSuperEdges().elementAt(j);
//							if (se.toString().equals("895929-896001^896122-896249^"))
//								System.currentTimeMillis();
							if (se.isPend())
								continue;	// TODO 2nd level for paired-end
							if ((sense&& se.getEdges()[0]!= e)
									|| ((!sense)&& se.getEdges()[se.getEdges().length- 1]!= e))
								continue;
							
							Partition ppx= new Partition(se.getTranscripts());
							p= Constants.binarySearch(part, 0, partPos, ppx);
							if (p>= 0)
								ppx= part[p];
							else {
								p= -(p+ 1);
								if (p< partPos) {
									System.arraycopy(part, p, part, p+1, partPos- p);
								}
								++partPos;
								part[p]= ppx;
							}
							Transcript ttx= g.getAnyTranscript(se.getTranscripts());
							ntCov= se.getNtCoverage(ttx,
									sense? Constants.DIR_FORWARD: Constants.DIR_BACKWARD, mapLenMax);
							if (ntCov< 0) {
								ntCov= se.getNtCoverage(ttx,
										sense? Constants.DIR_FORWARD: Constants.DIR_BACKWARD, mapLenMax);
							}
							reads= sense? se.getReadNr(): se.getRevReadNr();
							if (Double.isNaN(ntCov))
								assert(reads== 0);
							else {
								ppx.add(ntCov);
								ppx.addReads(reads);
							}
						}
					}
				}
				return part;
			}

			private void completeConstraint(PrintStream p, String id, int length, double[] val, int[] idx,
					int eq, double d) {
				StringBuilder sb= new StringBuilder(id);
				sb.append("\t");
				for (int i = 0; i < length; i++) {
					double res= result[restrNr+ idx[i]];
					double var= val[i]* res;				
					sb.append(var> 0?"+":(var< 0? "-": (val[i]>= 0? "+": "-")));
					sb.append(Math.abs(var)== 1?"": StringUtils.fprint(Math.abs(var), 2));
	//				sb.append("C");
	//				sb.append(idx[i]);
					sb.append(" ");
					
					int pos= Arrays.binarySearch(costIdx, idx[i]);
					if (pos>= 0) {
						sb.append("(*");
						sb.append(StringUtils.fprint(costVal[pos], 2));
						sb.append("= ");
						sb.append(StringUtils.fprint(costVal[pos] * res, 2));
						sb.append(") ");
					}
				}
				sb.append(eq== LpSolve.EQ?"=":(eq== LpSolve.GE?">=":"<="));
				sb.append(" ");
				sb.append(d);
				p.println(sb.toString());
				
			}
	
			private void setCBase(byte mode, Graph g, Edge e, boolean sense, 
					IntVector locIdxV, IntVector[] cSumV, IntVector[] cBalV,
					IntVector costIdxV, DoubleVector costValV, 
					IntVector[] fluxIdxV, DoubleVector[] fluxValV) throws LpSolveException {
				
				assert(!(e instanceof SuperEdge));
				boolean fluxActive= true;
				
				for (int i = 0; i < fluxIdxV.length; i++) { 
					fluxIdxV[i].removeAll();
					fluxValV[i].removeAll();
				}
				for (int i = 0; i < cSumV.length; i++) 
					cSumV[i].removeAll();	// sum of contributions
				
				// base edge
				setCDeconPartition(mode, g, null, e, sense, 
						locIdxV, cSumV, cBalV, costIdxV, costValV, fluxIdxV, fluxValV);
				
				int[] decBase= null;
				if (mode!= MODE_COUNT)
					decBase= decIdx.clone();
				int decPos= 0;
				
				// iterate superedges
				for (int i = 0; e.getSuperEdges()!= null&& i < e.getSuperEdges().size(); ++i) {
					
					SuperEdge se= e.getSuperEdges().elementAt(i);
					if (se.isPend())
						continue;	// TODO 2nd level for paired-end
					if ((sense&& se.getEdges()[0]!= e)
							|| ((!sense)&& se.getEdges()[se.getEdges().length- 1]!= e))
						continue;
					
					setCDeconPartition(mode, g, e, se, sense, 
							locIdxV, cSumV, cBalV, costIdxV, costValV, fluxIdxV, fluxValV);
					
				}				
				
				// expectation restriction
				for (int i = 0; i < cSumV.length; i++) {
					int len= cSumV[i].size();
					if (len== 0)
						continue;
					if (mode== MODE_COUNT) {
						++restrNr;
						continue;
					}
					int[] idx= new int[len+ 1];	// TODO reuse
					idx[0]= sense? ((i* 4)+1): ((i* 4)+2);
					System.arraycopy(cSumV[i].vector, 0, idx, 1, cSumV[i].length);
					double[] val= new double[idx.length];
					Arrays.fill(val, 1d);
					Transcript tx= g.trpts[i];
					int tlen= tx.getExonicLength();
					UniversalMatrix m= profile.getMatrix(tlen);
					double f= m.getFrac(
								tx.getExonicPosition(e.getFrac(true)),
								tx.getExonicPosition(e.getFrac(false)),
								tlen,
								sense?Constants.DIR_FORWARD:Constants.DIR_BACKWARD);
					if (Double.isInfinite(f)|| Double.isNaN(f)) {
						System.err.println("infinite value");
						f= m.getFrac(
								tx.getExonicPosition(e.getFrac(true)),
								tx.getExonicPosition(e.getFrac(false)),
								tlen,
								sense?Constants.DIR_FORWARD:Constants.DIR_BACKWARD);
					}
					ccheck[i]+= f;
					val[0]= -f;
					
					if (mode== MODE_CONSTRAINT) {
						getLPsolve().addConstraintex(idx.length, val, idx, LpSolve.EQ, 0);
					} else if (mode== MODE_COMPLETE) {
						printConstraint(System.out, e+ " ("+ (sense?"s":"a")+ ", "+ g.trpts[i].getTranscriptID()+ ")", 
								idx.length, val, idx, LpSolve.EQ, 0d);
					}
					++restrNr;
				}
				
				// flux preservation, 2 for min/max readlength
				for (int j = 0; fluxActive&& j < fluxIdxV.length; j++) {	// per tx
					// assure there is an exon and a sj
					if (fluxValV[j].length< 2) 	
						continue;				
					boolean plus= false, minus= false;
					for (int k = 0; k < fluxValV[j].length; k++) {
						double d= fluxValV[j].get(k);
						plus|= (d> 0);
						minus|= (d< 0);
						if (plus&& minus)
							break;
					}
					if (!(plus&& minus))
						continue;
					if (mode== MODE_CONSTRAINT) {
						getLPsolve().addConstraintex(fluxIdxV[j].length, fluxValV[j].toDoubleArray(), 
								fluxIdxV[j].toIntArray(), LpSolve.EQ, 0);	// i== 0? LpSolve.LE: LpSolve.GE
					} else if (mode== MODE_COMPLETE) {
						printConstraint(System.out, e+ " (flux)", fluxIdxV[j].length, fluxValV[j].toDoubleArray(), 
								fluxIdxV[j].toIntArray(), LpSolve.EQ, 0);	// i== 0? LpSolve.LE: LpSolve.GE
					}
					++restrNr;
				}
				
			}
	
			private double getGaussUB(double sig2, double pthold) {
				
				double sqrt2PiSig2= Math.sqrt(2d* Math.PI* sig2);
				double prod= sqrt2PiSig2* pthold;
				double log= Math.log(prod);
				double prod2= -2d* sig2* log;
				double x= Math.sqrt(Math.abs(prod2));
				if (x< 0|| Double.isInfinite(x)|| Double.isNaN(x))
					System.currentTimeMillis();
				assert(x>= 0);
				return x;
			}
			
			private double[] gaussVals= new double[2];
			private void setCDeconPartition(byte mode, Graph g, Edge f, Edge e, boolean sense, 
							IntVector locIdxV, IntVector[] cSumV, IntVector[] cBalV,
							IntVector costIdxV, DoubleVector costValV, 
							IntVector[] fluxIdxV, DoubleVector[] fluxValV) 
					throws LpSolveException {
				
				boolean ubActive= true, partActive= true;
				
				// attributes of e
				int nbTx= g.getTranscriptNb(e.getTranscripts());
				int ctrTT= -1;
				int obs= (sense? e.getReadNr(): e.getRevReadNr());
				boolean sj= (e instanceof SuperEdge)&& (!((SuperEdge) e).isPend());
				//double effLenMin= calcEffLen(e, f, true);
				//double effLenMax= calcEffLen(e, f, false);
				//assert(effLenMin> 0|| obs== 0);

				// DECONVOLUTION constraints f.a. tx
				int ntLength= -1, mapLength= -1;
				double fac= -1;
				for (int j = 0; j < g.trpts.length; j++) {
					
					if (!g.contains(e.getTranscripts(), j))
						continue;
					
					if (ntLength< 0) {	// init once
						Transcript tx= g.trpts[j];
						mapLength= e.getMapLength(tx, mapLenMax);
						ntLength= e.getNtLength(tx, mapLength);
						fac= mapLength/ (double) ntLength;
					}
					
					++ctrTT;
					cSumV[j].add(++constraintCtr);	// sum(contributions)
					fluxIdxV[j].add(constraintCtr);	// flux preservation
					fluxValV[j].add(1d/ (sj? -fac: fac));	// inverse
					if (mode!= MODE_COUNT) {
						decIdx[ctrTT]= constraintCtr;		// c
						decVal[ctrTT]= 1d;
					}
				}	/* end tx */
//				if (ctrTT+ 1!= nbTx)
//					System.currentTimeMillis();
				assert(ctrTT+ 1== nbTx);
				
				// locus deviation, DECONVOLUTION
				Partition px= new Partition(e.getTranscripts());
				int p= Constants.binarySearch(part, 0, partPos, px);				
				assert(p>= 0); 
				if (mode== MODE_COUNT) {
					constraintCtr+= 2;
				} else {
					decIdx[nbTx]= ++constraintCtr; 	// +dL..
					decVal[nbTx]= 1d;
					cBalV[p].add(constraintCtr);
					locIdxV.add(-constraintCtr);		// ..substracting from obs
					decIdx[nbTx+ 1]= ++constraintCtr;	// -dL
					decVal[nbTx+ 1]= -1d;					
					cBalV[p].add(-constraintCtr);
					locIdxV.add(constraintCtr);			// ..adding
					if (mode== MODE_CONSTRAINT) 
						getLPsolve().addConstraintex(nbTx+ 2, decVal, decIdx, LpSolve.EQ, obs);
					else
						printConstraint(System.out, e.toString()+ " ("+(sense?"s, ":"a, ")+"dec)", nbTx+ 2, decVal, decIdx, LpSolve.EQ, obs);
				}
				++restrNr;	// DECONVOLUTION
				
				
				// costs according to observed PARTITION coverages
				for (int i = p; partActive&& i >= 0; --i) {
					
					if (i!= p&& !Graph.intersect(part[i].partition, px.partition).equals(
							part[i].partition))
						continue;	// not completely contained
					
					if (mode== MODE_COUNT) {
						constraintCtr+= 2; // +/- dP, reordering
						++restrNr;
						continue;
					}
					
					// sum deconvoluted coverages in partition
					nbTx= g.getTranscriptNb(part[i].partition);
					int tIdx= g.getNextTxIdx(part[i].partition, -1);
					for (int j = 0; tIdx>= 0&& j < nbTx; j++) {
						decIdx[j]= cSumV[tIdx].get(cSumV[tIdx].size()- 1);
						decVal[j]= fac;
						tIdx= g.getNextTxIdx(part[i].partition, tIdx);
					}
					double u= part[i].getMean();
					double sig= part[i].getSig2();
					double cost= COSTS_INFINITY, ub= 1;
					if (sig> 0) {
						ub= getGaussUB(sig, 0.01);
						cost= getGaussCost(sig, ub, 0.01);
					}
					decIdx[nbTx]= ++constraintCtr;
					decVal[nbTx]= 1d;	// +delta(part)
					if (ubActive) {
						if (mode== MODE_CONSTRAINT)
							getLPsolve().setUpbo(constraintCtr, ub);
						else if (mode== MODE_COMPLETE)
							System.out.println("C"+ constraintCtr+ " <= "+ ub);
					} else {
						if (mode== MODE_COMPLETE)
							System.out.println("NO ub C"+constraintCtr);
					}
					costIdxV.add(constraintCtr);
					costValV.add(cost);
					decIdx[nbTx+ 1]= ++constraintCtr;
					decVal[nbTx+ 1]= -1d;	// -delta(part)
					if (ubActive) {
						if (mode== MODE_CONSTRAINT)
							getLPsolve().setUpbo(constraintCtr, Math.min(ub, u));	
						else if (mode== MODE_COMPLETE)
							System.out.println("C"+ constraintCtr+ " <= "+ ub);
					} else {
						if (mode== MODE_COMPLETE)
							System.out.println("NO ub C"+constraintCtr);
					}
					costIdxV.add(constraintCtr);
					costValV.add(cost);
					
					if (mode== MODE_CONSTRAINT)
						getLPsolve().addConstraintex(nbTx+ 2, decVal, decIdx, LpSolve.EQ, u);
					else if (mode== MODE_COMPLETE)
						printConstraint(System.out, e.toString()+ " ("+(sense?"s, ":"a, ")+"cost)", nbTx+ 2, decVal, decIdx, LpSolve.EQ, u);
					++restrNr;

				}	/* end part */
				
			}

			private double[] getCostGaussian(int obsExon, int lenExon, int obsSJ, int lenSJ, 
					int readLenHere, double pthold, double[] result) {
				
				if (obsSJ== -1&& lenSJ== -1) {
					double x= Math.sqrt(-2d* Math.log(Math.sqrt(2d* Math.PI)* pthold));
					double y= 1d/ Math.sqrt(2d* Math.PI);
					double slope= (y- pthold)/ x;
					result[0]= slope;
					result[1]= x;
					return result;
				}
				
				// mean and variance
				double covExon= (obsExon* readLenHere)/ (double) lenExon;	// cov_nt, reads homogenized to common length
				double covSJ= (obsSJ* readLenHere)/ (double) lenSJ;
				double u= Math.abs(covExon- covSJ)/ 2d;	// mean
				double sig2= (Math.abs(covExon- u)+ Math.abs(covSJ- u))/ 2d;
				sig2*= sig2;	// variance, square of mean(abs. deviations)
				
				// intercept on threshold(p), value p(u)
				double sqrt2PiSig2= Math.sqrt(2d* Math.PI* sig2);
				double x= Math.sqrt(-2d* sig2* Math.log(Math.sqrt(sqrt2PiSig2* pthold)));
				assert(x>= 0);
				result[1]= x;
				double y= 1d/ sqrt2PiSig2;
				
				// -slope, -delta(y)/ delta(x) 
				double slope= (y- pthold)/ x; 
				assert(slope> 0);
				result[0]= slope;
				return result;
			}

			private double getGaussCost(double sig2, double x, double pthold) {
				
				double sqrt2PiSig2= Math.sqrt(2d* Math.PI* sig2);
				double pthold2= 1d/ sqrt2PiSig2* Math.exp(-(x* x)/ (2* sig2)); 
				double y= 1d/ sqrt2PiSig2;
				// slope= delta(y)/ delta(x) 
				double slope= (pthold2- y)/ x; 
				if (slope> 0|| Double.isNaN(slope)|| Double.isInfinite(slope))
					System.currentTimeMillis();
				assert(slope<= 0);
				// for costs, perpendicular
				slope= -1d/ slope;
				return slope;
			}

			private Object[] mateObjects= null;
			private byte[] mateOrient= null;
			private int[] mateInt= null;
			private int mateP;
			private int txStrand= -1;
			private int[] fromTo= new int[2];
			private int nrMappingsLocusMapped= 0;
			private int constraintCtr= 0, restrNr= 0;
			
		} /* end LocusSolver */

	static void printUsage() {
		if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) {
			
			Iterator<String[]> iter= cliExplMap.keySet().iterator();
			int max= 0;
			String[][] ss= new String[cliExplMap.size()][];
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
				
				String expl= cliExplMap.get(ss[i]);
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
	
	private static final String HELP_HEADER= "", HELP_FOOTER= ""; 
	private static void wellcome() {
		if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) 
			System.err.println("\n[HELLO] I am the FLUX CAPACITOR (build "+version+ "), nice to meet you!\n");
	}
	
	public static void main(String[] args) {
		
		try {
			readProperties(); 
			
			boolean showGUI= false;
			if (args== null|| args.length== 0) {
				if (showGUI) {
					FluxCapacitor.loadLibraries();
					FluxCapacitorGUI.createGUI();
				} else
					printUsage();
			}
			
			final FluxCapacitorNew myCapacitor= new FluxCapacitorNew();
			myCapacitor.init(args);
			if (myCapacitor.isHelpRequested()) {
				printUsage();
				System.exit(0);
			}
	
			wellcome();
			if (doInstall) {
				install();
				System.exit(0);
			}
			
			if (!myCapacitor.checkPreliminaries())
				System.exit(-1);
			
			if (!cheatDisableCleanup) {
				Runtime.getRuntime().addShutdownHook(new Thread("MrProper") {
				    public void run() { 
				    	FileHelper.cleanup(System.getProperty(Constants.PROPERTY_TMPDIR), 
				    			Constants.globalPfx== null?PFX_CAPACITOR+ "."+ myCapacitor.getRunID():Constants.globalPfx+ "_",
				    			null,
				    			Constants.verboseLevel> Constants.VERBOSE_SHUTUP?System.err:null); 
				    }
				});
			} else
				System.err.println("[NOCLEAN] Cleanup disabled!");
			
			
			if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) {
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
		fileProfileOriginal= null,
		fileLPdir= null,
		fileISize= null;
	
	int readLenMin= 75, readLenMax= -1;
	
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
	
	static private Options options= null;
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
			File f= File.createTempFile(PFX_CAPACITOR, "prescanIDs");
			BufferedWriter writer= new BufferedWriter(new FileWriter(f));
			long bRead= 0, bTot= fileBED.length();
			int perc= 0, lines= 0;
			for(String s; (s= buffy.readLine())!= null;++nrReadsAll,bRead+= s.length()+1, ++lines) {
                Log.progress(bRead,  bTot);

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
			e.printStackTrace();
			return false;
		}
				
		if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) {
			System.err.println("\ttook "+((System.currentTimeMillis()- t0)/ 1000)+" sec");
			System.err.println("\tfound "+nrReadsAll+" lines");
			System.err.println("\treadlength min "+lenMinMax[0]+", max "+ lenMinMax[1]+"\n");
		}
		readLenMin= lenMinMax[0];
		
		return true;
	}

	int init(Method m, String[] args, int p) {
		
		Object res= null;

		try {
			if (m.getParameterTypes()!= null&& m.getParameterTypes().length> 0) {
				int i = p+ 1; 
				while (i < args.length&& i- p<= m.getParameterTypes().length) {
					if (args[i].charAt(0)== CLI_SHORT_PFX|| args[i].startsWith(CLI_LONG_PFX))
						break;
					++i;
				}
				--i;
				int nrPars= i- p, expPars= m.getParameterTypes().length;				
				int diff= expPars- nrPars;
				if (diff< 0)
					throw new IllegalArgumentException("Too many arguments for method "+ m.getName());
				Object[] mArgs= new Object[expPars];
				for (int j = 0; j < nrPars; j++) 
					mArgs[j]= args[p+ 1+ j];
				for (int j = nrPars; j < expPars; j++) 
					mArgs[j]= null;
				res= m.invoke(this, mArgs);
				 
			} else	// nrPars== 0
				res= m.invoke(this, null);
			
		} catch (Exception e) {
			System.err.println("[OLALA] Could not set parameter "+ args[p]+
					"\n\t"+ e.getMessage());
			e.printStackTrace();
			System.exit(-1);
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
			if (args[i].startsWith(CLI_LONG_PFX)) {
				s= args[i].substring(CLI_LONG_PFX.length());
				if (! cliLongMap.containsKey(s)) {
					if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) 
						System.err.println("[OHNO] Unrecognized long option \'"+args[i]+"\'");
					Iterator<String> iter= cliLongMap.keySet().iterator();
					while (iter.hasNext())
						System.err.println("\'"+iter.next()+"\'");
					System.exit(-1);
				}
				Method m= cliLongMap.get(s);
				if ((i= init(m, args, i))< 0) 
					System.exit(-1);
				
			} else if (args[i].startsWith(CLI_SHORT_PFX.toString())){
				s= args[i].substring(1);
				for (int j = 0; j < s.length(); j++) {
					if (! cliShortMap.containsKey(s.charAt(j))) {
						if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) 
							System.err.println("[OHNO] Unrecognized short option \'"+args[i]+"\'");
						System.exit(-1);
					}
					Method m= cliShortMap.get(s.charAt(j));
					if ((i= init(m, args, i))< 0) 
						System.exit(-1);
				}
			}
		}
		
		return 0;
	}
	
	private int checkParameter() {
		
		if (fileBED== null|| !fileBED.exists()) {
			if (Constants.verboseLevel!= Constants.VERBOSE_SHUTUP) {
				System.err.println("[AIII] I need a input file with aligned reads in order to work, ok?!");
				if (fileBED== null)
					System.err.println("\tyou said nothing and this is bad.");
				else
					System.err.println("\tyou said it is "+fileBED+" but it is bad.");
				System.err.println("\tUse the "+CLI_LONG_PFX+CLI_LONG_SRA+" parameter and give me a correct one please.\n");
			}
			return -1;
		}
		if (fileGTF== null|| !fileGTF.exists()) {
			if (Constants.verboseLevel!= Constants.VERBOSE_SHUTUP) {
				System.err.println("[MANO] I need a reference file with the transcripts you want to analyze.");
				if (fileGTF== null)
					System.err.println("\tyou said nothing and this is not good.");
				else
					System.err.println("\tyou said it is "+fileGTF+" but it is not good.");
				System.err.println("\tTry again, using the "+CLI_LONG_PFX+CLI_LONG_REF+" parameter.\n");
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
		p.print("\tcmd\t"+CLI_CMD);
		for (int i = 0; i < args.length; i++) 
			p.print(" "+args[i]);
		p.println();
		
		try {
			// INPUT
			p.println("\t== INPUT ==");
			p.println("\t"+CLI_LONG_REF+"\t"+fileGTF.getCanonicalPath());
			p.println("\t"+CLI_LONG_SRA+"\t"+fileBED.getCanonicalPath());
			p.println("\tDescriptor\t"+descriptor.toString());
			if (pairedEnd|| stranded) {
				p.print("\tAttributes:\t");
				if (pairedEnd)
					p.print("paired-end ");
				if (stranded)
					p.print("stranded ");
				p.println();
			}
			if (fileProfileOriginal!= null&& fileProfileOriginal.exists())
				p.println("\t"+CLI_LONG_PROFILE+"\t"+fileProfileOriginal.getCanonicalPath());				
			p.println("\t"+CLI_LONG_VERBOSE+"\t"+Constants.VERBOSE_KEYWORDS[Constants.verboseLevel]);
			if (copyLocal)
				p.println("\t"+CLI_LONG_LOCAL+" copies");
			// OUTPUT
			p.println("\t== OUTPUT ==");
			p.println("\t"+CLI_LONG_TMP+"\t"+ System.getProperty("java.io.tmpdir"));
			if (Constants.globalPfx!= null)
				p.println("\t"+CLI_LONG_TPX+"\t"+ Constants.globalPfx);
			p.print("\t"+CLI_LONG_OUT+"\t");
			if (fileOutDir== null)
				p.println("stdout");
			else {
				p.println(fileOutDir.getCanonicalPath());
				if (compressionOut!= FileHelper.COMPRESSION_NONE)
					p.println("\t"+ CLI_LONG_COMPRESSION+ "\t"+ FileHelper.COMPRESSION_KEYWORDS[compressionOut]);
			}
//			if (fileProfileOriginal== null|| !fileProfileOriginal.exists())
//				p.println("\t"+CLI_LONG_PROFILE+"\t"+fileProfileOriginal.getCanonicalPath());				

			p.print("\tfeatures:\t");
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

			p.print("\tbases:\t");
			if (outputObs)
				p.print("observed ");
			if (outputPred)
				p.print("predicted ");
			if (outputBalanced)
				p.print("balanced ");
			p.println();
			
			p.print("\tscopes:\t");
			if (outputAll)
				p.print("all ");
			if (outputSplit)
				p.print("split ");
			if (outputUnique)
				p.print("unique ");
			p.println();
			
			p.print("\tmeasure:\t");
			if (outputFreq)
				p.print("frequency ");
			if (outputRfreq)
				p.print("rel.frequency ");
			if (outputRcov)
				p.print("rel.coverage ");
			p.println();

			if (outputMapped|| outputNotmapped|| outputProfiles|| outputLP) {
				p.print("\tsave:\t");
				if (outputMapped)
					p.print("mapped mappings ");
				if (outputNotmapped)
					p.print("not-mapped mappings ");
				if (outputProfiles)
					p.print("profiles ");
				if (outputLP)
					p.print("linear-programs ");
				p.println();
			}

			// ALGORITHM
			p.println("\tSETTINGS");
			p.println("\t"+ CLI_LONG_THREAD+" "+ maxThreads);
			p.print("\t"+CLI_LONG_UNIF+"\t");
			if (uniform)
				p.println("yes, uniform read distribution");
			else
				p.println("no, read distribution profiles");
			
		} catch (IOException e) {
			; // :)
		}
		
		p.println();

	}
	
	private final static String[] L_USER_COMMENTS= new String[] {"[OOOPS] ", "[HEOO] ", "[PLONG] ", "[BAOOO] "};
	private final static Random rndLuser= new Random();
	private String errorMissingArgument(String string) {
		return L_USER_COMMENTS[rndLuser.nextInt(L_USER_COMMENTS.length)]
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
	
	public boolean setAttributes(String attrib, String regExp) {
		if (attrib== null|| attrib.length()== 0)
			return false;
		attrib= attrib.toLowerCase();
		for (int i = 0; i < attrib.length(); i++) {
			if (attrib.charAt(i)== Descriptor.CHAR_ID_PAIRED)
				pairedEnd= true;
			else if (attrib.charAt(i)== Descriptor.CHAR_ID_STRANDED)
				stranded= true;
			else
				return false;
		}
		
		if (regExp== null) {
			descriptor= new BARNAdescriptor();
			return true;
		} else {
			RegExpDescriptor desc= new RegExpDescriptor();
			desc.init(attrib, regExp);
			//descriptor= desc;
		}
		return true; // descriptor.isValid();
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
	
	public void setProfile(String profileName) {
		this.fileProfileOriginal= new File(profileName);
	}
	
	public void setTempPfx(String tmpPfx) {
		Constants.globalPfx= tmpPfx;
	}
	
	public static final byte STRAND_NONE= 0, STRAND_ENABLED= 1, STRAND_SPECIFIC= 2;
	byte strand= STRAND_NONE;
	public void setStrandSpecific() {
		strand= STRAND_SPECIFIC;
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
	
	public static final char 
		CLI_OUT_ALL= 'a', 
		CLI_OUT_BALANCED= 'b',
		CLI_OUT_COVERAGE= 'c',
		CLI_OUT_PRED= 'd',
		// e
		CLI_OUT_RFREQ= 'f', // fraction
		CLI_OUT_GENE= 'g',
		// h
		CLI_OUT_ISIZE= 'i', 
		CLI_OUT_SJUNCTION= 'j',	// junctions in general.. separate later 
		CLI_OUT_KEEPSORTED= 'k', 
		CLI_OUT_LP= 'l', 
		CLI_OUT_MAPPED= 'm', 
		CLI_OUT_NOTMAPPED= 'n', 
		CLI_OUT_OBS= 'o', 
		CLI_OUT_PROFILES= 'p',
		CLI_OUT_UNIQUE= 'q',	// WAS: u
		CLI_OUT_FREQ= 'r',	// reads 
		CLI_OUT_SPLIT= 's', 
		CLI_OUT_TRANSCRIPT= 't',
		CLI_OUT_LOCUS= 'u',	// locus, cluster
		CLI_OUT_EVENTS= 'v',
		// w
		CLI_OUT_EXON= 'x',	// WAS: e 
		CLI_OUT_INTRON= 'y';
		// z

	
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
			if (s.charAt(i)== CLI_OUT_ALL) 
				outputAll= false;
			else if (s.charAt(i)== CLI_OUT_SPLIT)
				outputSplit= false;
			else if (s.charAt(i)== CLI_OUT_UNIQUE)
				outputUnique= false;
			else if (s.charAt(i)== CLI_OUT_OBS)
				outputObs= false;
			else if (s.charAt(i)== CLI_OUT_PRED)
				outputPred= false;
			else if (s.charAt(i)== CLI_OUT_BALANCED)
				outputBalanced= false;
			else if (s.charAt(i)== CLI_OUT_FREQ)
				outputFreq= false;
			else if (s.charAt(i)== CLI_OUT_RFREQ)
				outputRfreq= false;
			else if (s.charAt(i)== CLI_OUT_COVERAGE)
				outputRcov= false;
			else if (s.charAt(i)== CLI_OUT_EXON)
				outputExon= false;
			else if (s.charAt(i)== CLI_OUT_SJUNCTION)
				outputSJunction= false;
			else if (s.charAt(i)== CLI_OUT_TRANSCRIPT)
				outputTranscript= false;
			else if (s.charAt(i)== CLI_OUT_GENE)
				outputGene= false;
			else if (s.charAt(i)== CLI_OUT_EVENTS)
				outputEvent= false;
			else if (s.charAt(i)== CLI_OUT_LP)
				outputLP= false;
			else if (s.charAt(i)== CLI_OUT_PROFILES)
				outputProfiles= false;
			else if (s.charAt(i)== CLI_OUT_MAPPED)
				outputMapped= false;
			else if (s.charAt(i)== CLI_OUT_NOTMAPPED)
				outputNotmapped= false;
			else if (s.charAt(i)== CLI_OUT_ISIZE)
				outputISize= false;
			else if (s.charAt(i)== CLI_OUT_KEEPSORTED)
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
	
	
//	protected static HashMap<String, Method> cliLongMap= new HashMap<String, Method>(); 
//	protected static HashMap<Character, Method> cliShortMap= new HashMap<Character, Method>();
//	protected static HashMap<String[], String> cliExplMap= new HashMap<String[], String>();
	static {
		try {
			Method m;
			
			m= FluxCapacitor.class.getDeclaredMethod("setAttributes", new Class[] {String.class, String.class});
			cliShortMap.put(CLI_SHORT_ATTRIBUTES, m);
			cliLongMap.put(CLI_LONG_ATTRIBUTES, m);
			cliExplMap.put(new String[] {CLI_SHORT_PFX+ CLI_SHORT_ATTRIBUTES.toString(), CLI_LONG_PFX+ CLI_LONG_ATTRIBUTES}, 
					"set the attributes of experiment, using\n"
					+ "p for paired-end\n"
					+ "s for stranded\n"
					+ "and (optionally) an regular expression\n"
					+ "to tokenize the read identifier, for instance the expression\n"
					+ CLI_SHORT_PFX+ CLI_SHORT_ATTRIBUTES.toString()+ " ps /([12])_strand([12])$\n"
					+ "assigns pair mates by the symbols \'1\',\'2\' between '/' and \'_strand\' "
					+ "and the strandedness according to the symbols \'1\' (sense) respectively "
					+ "\'2\' (antisense) between \'_strand\' and the end of the read identifier. "
					+ "For more information, have a look at a regular expression documentation, e.g.,\n"
					+ "http://java.sun.com/j2se/1.4.2/docs/api/java/util/regex/Pattern.html");
			
			m= FluxCapacitor.class.getDeclaredMethod("setFileReads", new Class[] {String.class});
			cliShortMap.put(CLI_SHORT_SRA, m);
			cliLongMap.put(CLI_LONG_SRA, m);
			cliExplMap.put(new String[] {CLI_SHORT_PFX+ CLI_SHORT_SRA.toString(), CLI_LONG_PFX+ CLI_LONG_SRA}, 
					"set file containing Short Reads Archive (mandatory!)\n");
			
			m= FluxCapacitor.class.getDeclaredMethod("setFileReference", new Class[] {String.class});
			cliShortMap.put(CLI_SHORT_REF, m);
			cliLongMap.put(CLI_LONG_REF, m);
			cliExplMap.put(new String[] {CLI_SHORT_PFX+ CLI_SHORT_REF.toString(), CLI_LONG_PFX+ CLI_LONG_REF}, 
					"set file with REFerence annotation (mandatory!)\n");
			
			m= FluxCapacitor.class.getDeclaredMethod("setNameOutDir", new Class[] {String.class});
			cliShortMap.put(CLI_SHORT_FILENAME, m);
			cliLongMap.put(CLI_LONG_FILENAME, m);
			cliExplMap.put(new String[] {CLI_SHORT_PFX+ CLI_SHORT_FILENAME.toString(), CLI_LONG_PFX+ CLI_LONG_FILENAME}, 
					"set output fileName prefix (default stdout)\n");
			
			m= FluxCapacitor.class.getDeclaredMethod("setForce", null);
			cliShortMap.put(CLI_SHORT_FORCE, m);
			cliLongMap.put(CLI_LONG_FORCE, m);
			cliExplMap.put(new String[] {CLI_SHORT_PFX+ CLI_SHORT_FILENAME.toString(), CLI_LONG_PFX+ CLI_LONG_FILENAME}, 
					"set force (no overwrite checks)\n");
			
			m= FluxCapacitor.class.getDeclaredMethod("setInstall", (Class[]) null);
			cliLongMap.put(CLI_LONG_INSTALL, m);
			cliExplMap.put(new String[] {CLI_LONG_PFX+ CLI_LONG_INSTALL}, 
					"installs the basic wrapper script (no reads are mapped)\n");
			
			m= FluxCapacitor.class.getDeclaredMethod("setJVM", new Class[] {String.class});
			cliLongMap.put(CLI_LONG_JVM, m);
			cliExplMap.put(new String[] {CLI_LONG_PFX+ CLI_LONG_JVM}, 
					"set a specific Java Virtual Machine home (installation)\n");

			m= FluxCapacitor.class.getDeclaredMethod("setLib", new Class[] {String.class});
			cliLongMap.put(CLI_LONG_LIB, m);
			cliExplMap.put(new String[] {CLI_LONG_PFX+ CLI_LONG_LIB}, 
					"set path to native libraries (installation)\n");

			m= FluxCapacitor.class.getDeclaredMethod("setBatch", (Class[]) null);
			cliShortMap.put(CLI_SHORT_BATCH, m);
			cliLongMap.put(CLI_LONG_BATCH, m);
			cliExplMap.put(new String[] {CLI_SHORT_PFX+ CLI_SHORT_BATCH.toString(), CLI_LONG_PFX+ CLI_LONG_BATCH}, 
					"set Batch mode, suppresses file checks and stderr communication\n");
			
			m= FluxCapacitor.class.getDeclaredMethod("setUniformal", (Class[]) null);
			cliShortMap.put(CLI_SHORT_UNIF, m);
			cliLongMap.put(CLI_LONG_UNIF, m);
			cliExplMap.put(new String[] {CLI_SHORT_PFX+ CLI_SHORT_UNIF.toString(), CLI_LONG_PFX+ CLI_LONG_UNIF}, 
					"set uniformal distribution no profiling step is carried out\n");
			
			m= FluxCapacitor.class.getDeclaredMethod("setPairedEnd", (Class[]) null);
			cliShortMap.put(CLI_SHORT_PAIR, m);
			cliLongMap.put(CLI_LONG_PAIR, m);
			cliExplMap.put(new String[] {CLI_SHORT_PFX+ CLI_SHORT_PAIR.toString(), CLI_LONG_PFX+ CLI_LONG_PAIR}, "set input paired ends, " +
					"read name expected in FMRD format (see http://fluxcapacitor.wikidot.com/formats:fmrd)\n");

			m= FluxCapacitor.class.getDeclaredMethod("setProfile", new Class[] {String.class});
			cliLongMap.put(CLI_LONG_PROFILE, m);
			cliExplMap.put(new String[] {CLI_SHORT_PFX+ CLI_SHORT_PAIR.toString(), CLI_LONG_PFX+ CLI_LONG_PAIR}, "set profile name");
			
			m= FluxCapacitor.class.getDeclaredMethod("setOutput", new Class[] {String.class});
			cliLongMap.put(CLI_LONG_OUT, m);
			cliShortMap.put(CLI_SHORT_OUT, m);
			cliExplMap.put(new String[] {CLI_SHORT_PFX+ CLI_SHORT_OUT.toString(), CLI_LONG_PFX+ CLI_LONG_OUT}, 
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

			m= FluxCapacitor.class.getDeclaredMethod("setHelp", (Class[]) null);
			cliLongMap.put(CLI_LONG_HELP, m);
			cliShortMap.put(CLI_SHORT_HELP, m);
			cliExplMap.put(new String[] {CLI_SHORT_PFX+ CLI_SHORT_HELP.toString(), CLI_LONG_PFX+ CLI_LONG_HELP},
					"print help summary");
			
			m= FluxCapacitor.class.getDeclaredMethod("setLogLevel", new Class[] {String.class});
			cliLongMap.put(CLI_LONG_VERBOSE, m);
			cliShortMap.put(CLI_SHORT_VERBOSE, m);
			cliExplMap.put(new String[] {CLI_SHORT_PFX+ CLI_SHORT_VERBOSE.toString(), CLI_LONG_PFX+ CLI_LONG_VERBOSE},
					"set verbose level (SILENT, VERBOSE, ERRORS, DEBUG)");

			m= FluxCapacitor.class.getDeclaredMethod("setThreads", new Class[] {String.class});
			cliLongMap.put(CLI_LONG_THREAD, m);
			cliShortMap.put(CLI_SHORT_THREAD, m);
			cliExplMap.put(new String[] {CLI_SHORT_PFX+ CLI_SHORT_THREAD.toString(), CLI_LONG_PFX+ CLI_LONG_THREAD}, "set multi-thread mode, provide number of threads\n" +
					"(time gain only with complex linear programs, otherwise default=1 recommended)\n");
			
			m= FluxCapacitor.class.getDeclaredMethod("setLocal", (Class[]) null);
			cliLongMap.put(CLI_LONG_LOCAL, m);
			cliShortMap.put(CLI_SHORT_LOCAL, m);
			cliExplMap.put(new String[] {CLI_SHORT_PFX+ CLI_SHORT_LOCAL.toString(), CLI_LONG_PFX+ CLI_LONG_LOCAL}, 
					"work locally, i.e., copy all files to the temporary directory\n");
			
			m= FluxCapacitor.class.getDeclaredMethod("setTempDir", new Class[] {String.class});
			cliLongMap.put(CLI_LONG_TMP, m);
			cliExplMap.put(new String[] {CLI_LONG_PFX+ CLI_LONG_TMP}, 
					"set path to the temporary directory\n");
			
			m= FluxCapacitor.class.getDeclaredMethod("setTempPfx", new Class[] {String.class});
			cliLongMap.put(CLI_LONG_TPX, m);
			cliExplMap.put(new String[] {CLI_LONG_PFX+ CLI_LONG_TPX}, 
					"set prefix for temporary files\n");
			
			m= FluxCapacitor.class.getDeclaredMethod("setCompression", new Class[] {String.class});
			cliShortMap.put(CLI_SHORT_COMPRESSION, m);
			cliLongMap.put(CLI_LONG_COMPRESSION, m);
			cliExplMap.put(new String[] {CLI_SHORT_PFX+ CLI_SHORT_COMPRESSION.toString(), 
						CLI_LONG_PFX+ CLI_LONG_COMPRESSION}, 
					"set compression method for output files (output file, mapped read-mappings, not-mapped read-mappings, insert sizes)\n");
			
			m= FluxCapacitor.class.getDeclaredMethod("setStrandSpecific", (Class[]) null);
			cliLongMap.put(CLI_LONG_SSPECIFIC, m);
			cliExplMap.put(new String[] {CLI_LONG_PFX+ CLI_LONG_SSPECIFIC}, 
					"set strand specific reads (default: strand information disregarded/disabled)\n");
			
		} catch (NoSuchMethodException e) {
			; // :)
		}
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
			File dst= new File(fileOutDir+ File.separator+ getNameLP()+ (sfx== null? "": Constants.DOT+ sfx));
			if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) 
				System.err.println("\tzipping "+ fileLPdir.getAbsolutePath()
						+"\n\t->"+ dst.getAbsolutePath());
			if (!FileHelper.zip(fileLPdir, dst)) {
				if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
					System.err.println("\n[PROBLEMS] encountered error during zipping, check file.");
			} else {
				FileHelper.rmDir(fileLPdir); // otherwise Mr. Proper
			}
		}
		
// 		if (compressionOut== FileHelper.COMPRESSION_NONE&& !copyLocal)
//			return;

		// move 
		moveOrDeflate(fileOut, 
				new File(fileOutDir+ File.separator+ getNameOut()), compressionOut);

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
		if (outputProfiles) 
			moveOrDeflate(getFileProfile(), 
					new File(fileOutDir+ File.separator+ getNameProfile()), 
					FileHelper.COMPRESSION_NONE);	// already compressed
		
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
				+ PFX_CAPACITOR
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
		fileGTForiginal= fileGTF;
		
		// check compression
		compressionGTF= FileHelper.getCompression(fileGTF);
		String deflatedFName= System.getProperty(Constants.PROPERTY_TMPDIR);
		if (compressionGTF== FileHelper.COMPRESSION_NONE) 
			deflatedFName+= fileGTF.getName();
		else
			deflatedFName+= MyFile.stripExtension(fileGTF.getName());
		File deflatedFile= new File(deflatedFName);
		if (!deflatedFile.exists())	// TODO specific runID
			if (!copyOrDeflate(fileGTF, deflatedFile, compressionGTF))
				return false;
		this.fileGTF= deflatedFile;
		
		// check sorting
		if (1== 1) {
			System.err.println("\tTEST Sort check disabled !!!");
			isSortedGTF= true;
		} else if (getGTFreader().isApplicable()) 
			isSortedGTF= true;
		else {
			isSortedGTF= false;
			File tmpGTF= getGTFreader().createSortedFile();
			boolean bb= getFileGTF().delete();	// TODO do sth when false..
			this.fileGTF= tmpGTF;
			if (outputSorted&& fileOutDir!= null) {
				String ext= FileHelper.getCompressionExtension(compressionGTF);
				String destFName= fileOutDir+ File.separator+ 
									fileGTF.getName()+
									(ext== null?"":Constants.DOT+ ext);
				File destFile= new File(destFName);
				if ((!force)&& !ensureFileCanWrite(destFile))
					return false;
				else if (!copyOrDeflate(tmpGTF, destFile,compressionGTF))
					return false;
			}
		}
		
		// scan file
		if (1== 1) {
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
		if (getFileBED()== null)
			return false;
		
		fileBEDoriginal= getFileBED();
		
		// check compression
		byte compressionBED= FileHelper.getCompression(fileBED); // global var for output
		String deflatedFName= System.getProperty(Constants.PROPERTY_TMPDIR);
		File deflatedFile= null;
		if (compressionBED== FileHelper.COMPRESSION_NONE) { 
			deflatedFName+= fileBED.getName();
			deflatedFile= new File(deflatedFName);
			if (!deflatedFile.exists())	// TODO make specific runID
				copy(fileBED, deflatedFile);
		} else {
			deflatedFName+= MyFile.stripExtension(fileBED.getName());
			deflatedFile= new File(deflatedFName);
			if (!deflatedFile.exists())	// TODO make specific runID
				if (!copyInflate(fileBED, deflatedFile, compressionBED))
					return false;		
		}
		this.fileBED= deflatedFile;
		
		// check sorting
		if (1== 1) {
			System.err.println("\tTEST sort check disabled !!!");
			isSortedBED= true;
		} else if (getBedReader().isApplicable()) 
			isSortedBED= true;
		else {
			isSortedBED= false;
			File tmp= getBedReader().sortBED(fileBED);
			getFileBED().delete();
			this.fileBED= tmp;
			if (outputSorted&& fileOutDir!= null) {
				String destFName= fileOutDir+ File.separator+ 
									MyFile.append(fileBED.getName(), "_sorted", false, 
											FileHelper.getCompressionExtension(this.compressionBED));
				File destFile= new File(destFName);
				if ((!force)&& !ensureFileCanWrite(destFile))
					return false;
				else if (!copyOrDeflate(tmp, destFile,this.compressionBED))
					return false;
			}
		}
		
		
		// scan file
		bedWrapper= null;
/*		descriptor= getBedReader().checkReadDescriptor(pairedEnd);
		if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) {
			System.err.print("\tread descriptor ");
			if (descriptor== null)
				System.err.println("none");
			else if (descriptor instanceof FMRD)
				System.err.println("flux");
			else if (descriptor instanceof SolexaDescriptor)
				System.err.println("solexa");
		}
		if (pairedEnd&& descriptor== null) { 
			if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
				System.err.println("[FMRD] No paired end information in FMRD format, see\n\thttp://fluxcapacitor.wikidot.com/formats:fmrd-c");
			exit(-1);
		} 
*/		

		if (1== 1) {
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
	
	public static final String SFX_INFLATE= "__inflated";
	public static final String SFX_MAPPED= "_mapped", SFX_NOTMAPPED= "_notmapped", SFX_PROFILES= "_profiles", SFX_LP= "_lp", SFX_INSERTSIZE= "_insertsize";
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
		
		Constants.verboseLevel= Constants.VERBOSE_NORMAL+ 1;
		
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
		
		profile= new Profile(this);
		if (uniform) {
			int nr= profile.fill();
		} else {
			if (fileProfileOriginal!= null&& fileProfileOriginal.exists()) {
				System.err.println();
				profile= readProfiles(fileProfileOriginal);
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
				// dont save, keep original data
			} else {
				try {
					explore(MODE_LEARN);			
				} catch (Throwable e) {
					e.printStackTrace();
					if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
						System.err.println("[FATAL] Error occured during scanning\n\t"+ e.getMessage());
				}
				if (outputProfiles)
					writeProfiles();
			}

		}
		
//		System.exit(0);
		
//		if (map)
//			func.finish();
		
//		try {
//			testWriter= new BufferedWriter(new FileWriter("P:\\rgasp1.3\\HepG2_new\\test.bed"));
//		} catch (IOException e) {
//			e.printStackTrace();
//		}
		
		
		explore(MODE_RECONSTRUCT);
		
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
			
			ZipFile zf= new ZipFile(fileProfileOriginal);
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

	private static final String obsReadsAllTag= GTF_ATTRIBUTE_TOKEN_OBSV+ GTF_ATTRIBUTE_TOKEN_SEP+ GTF_ATTRIBUTE_TOKEN_ALL+ GTF_ATTRIBUTE_TOKEN_SEP+ GTF_ATTRIBUTE_TOKEN_READS, 
	obsReadsSplitTag= GTF_ATTRIBUTE_TOKEN_OBSV+ GTF_ATTRIBUTE_TOKEN_SEP+ GTF_ATTRIBUTE_TOKEN_TID+ GTF_ATTRIBUTE_TOKEN_SEP+ GTF_ATTRIBUTE_TOKEN_READS,				
	obsReadsUniqTag= GTF_ATTRIBUTE_TOKEN_OBSV+ GTF_ATTRIBUTE_TOKEN_SEP+ GTF_ATTRIBUTE_TOKEN_EXC+ GTF_ATTRIBUTE_TOKEN_SEP+ GTF_ATTRIBUTE_TOKEN_READS, 
	predReadsAllTag= GTF_ATTRIBUTE_TOKEN_PRED+ GTF_ATTRIBUTE_TOKEN_SEP+ GTF_ATTRIBUTE_TOKEN_ALL+ GTF_ATTRIBUTE_TOKEN_SEP+ GTF_ATTRIBUTE_TOKEN_READS,
	predReadsSplitTag= GTF_ATTRIBUTE_TOKEN_PRED+ GTF_ATTRIBUTE_TOKEN_SEP+ GTF_ATTRIBUTE_TOKEN_TID+ GTF_ATTRIBUTE_TOKEN_SEP+ GTF_ATTRIBUTE_TOKEN_READS, 
	predReadsUniqTag= GTF_ATTRIBUTE_TOKEN_PRED+ GTF_ATTRIBUTE_TOKEN_SEP+ GTF_ATTRIBUTE_TOKEN_EXC+ GTF_ATTRIBUTE_TOKEN_SEP+ GTF_ATTRIBUTE_TOKEN_READS,
	lengthAllTag= GTF_ATTRIBUTE_LENGTH+ GTF_ATTRIBUTE_TOKEN_SEP+ GTF_ATTRIBUTE_TOKEN_ALL,
	lengthSplitTag= GTF_ATTRIBUTE_LENGTH+ GTF_ATTRIBUTE_TOKEN_SEP+ GTF_ATTRIBUTE_TOKEN_TID,
	lengthUniqTag= GTF_ATTRIBUTE_LENGTH+ GTF_ATTRIBUTE_TOKEN_SEP+ GTF_ATTRIBUTE_TOKEN_EXC,
	normReadsAllTag= GTF_ATTRIBUTE_TOKEN_BALANCED+ GTF_ATTRIBUTE_TOKEN_SEP+ GTF_ATTRIBUTE_TOKEN_ALL+ GTF_ATTRIBUTE_TOKEN_SEP+ GTF_ATTRIBUTE_TOKEN_READS,
	normReadsSplitTag= GTF_ATTRIBUTE_TOKEN_BALANCED+ GTF_ATTRIBUTE_TOKEN_SEP+ GTF_ATTRIBUTE_TOKEN_TID+ GTF_ATTRIBUTE_TOKEN_SEP+ GTF_ATTRIBUTE_TOKEN_READS, 
	normReadsUniqTag= GTF_ATTRIBUTE_TOKEN_BALANCED+ GTF_ATTRIBUTE_TOKEN_SEP+ GTF_ATTRIBUTE_TOKEN_EXC+ GTF_ATTRIBUTE_TOKEN_SEP+ GTF_ATTRIBUTE_TOKEN_READS;

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
			
			
			if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) {
				Log.message("\tappending relative measurements");
                Log.progressStart("progress");
			}
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
				if (ss[2].equals(GFF_FEATURE_FRAGMENT)|| ss[2].equals(GFF_FEATURE_PAIRED)) {
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
							if (ss[i].equals(obsReadsAllTag)) 
								appendReads[0][0]= Double.parseDouble(sss== null?target: sss[dim]);
							else if (ss[i].equals(obsReadsSplitTag)) 
								appendReads[0][1]= Double.parseDouble(sss== null?target: sss[dim]);
							else if (ss[i].equals(obsReadsUniqTag)) 
								appendReads[0][2]= Double.parseDouble(sss== null?target: sss[dim]);
							else if (ss[i].equals(predReadsAllTag)) 
								appendReads[1][0]= Double.parseDouble(sss== null?target: sss[dim]);
							else if (ss[i].equals(predReadsSplitTag)) 
								appendReads[1][1]= Double.parseDouble(sss== null?target: sss[dim]);
							else if (ss[i].equals(predReadsUniqTag)) 
								appendReads[1][2]= Double.parseDouble(sss== null?target: sss[dim]);
							else if (ss[i].equals(normReadsAllTag)) 
								appendReads[1][0]= Double.parseDouble(sss== null?target: sss[dim]);
							else if (ss[i].equals(normReadsSplitTag)) 
								appendReads[1][1]= Double.parseDouble(sss== null?target: sss[dim]);
							else if (ss[i].equals(normReadsUniqTag)) 
								appendReads[1][2]= Double.parseDouble(sss== null?target: sss[dim]);
							else if (ss[i].equals(lengthAllTag)) 
								appendLengths[0]= Double.parseDouble(sss== null?target: sss[dim]);
							else if (ss[i].equals(lengthSplitTag)) 
								appendLengths[1]= Double.parseDouble(sss== null?target: sss[dim]);
							else if (ss[i].equals(lengthUniqTag)) 
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
					if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
						System.err.println("[FAILED] Cannot move file.");
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
			String tag= GTF_ATTRIBUTES_BASE[i];
			tag+= GTF_ATTRIBUTE_TOKEN_SEP;
			for (int j = 0; j < 3; j++) {	// all, split, uniq
				String tag2= tag+ GTF_ATTRIBUTES_RESOLUTION[j];
				tag2+= GTF_ATTRIBUTE_TOKEN_SEP;
				for (int x = 0; x < 2; x++) {	// freq, cov
					int pos= (i* 6)+ (j* 2)+ x;
					if (lengths[j]< 0|| reads[i][j]< 0) {
						sb[pos]= null;
						continue;
					}
					String tag3= tag2+ GTF_ATTRIBUTES_MEASUREMENT[x+1];
					if (sb[pos].length()== 0) {
						sb[pos].append(" ");
						sb[pos].append(tag3);
						sb[pos].append(" "); // \"
					} else
						sb[pos].append(",");

					double rfreq= 0;
					//assert((nrReadsMapped== 0)== (reads[i][j]== 0));
					double base= nrMappingsPairsMappedWt+ nrMappingsPairsMappedWg+ nrMappingsPairsMappedCt+ nrMappingsPairsMappedCg;
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
		return MyFile.stripExtension(fileGTF.getName())+ "__"
					+ MyFile.stripExtension(fileBED.getName());
	}
	
	private String getNameProfile() {
		String fName= getCompositeFName();
		if (fName== null)
			return null;
		String sfx= FileHelper.getCompressionExtension(compressionProfiles);
		return fName+ SFX_PROFILES+ (sfx== null? "": Constants.DOT+ sfx);
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
		return getCompositeFName()+ Constants.DOT+ SFX_GTF;
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
		return s+ SFX_INSERTSIZE+ Constants.DOT+ "txt";
	}
	
	public File getFileISize() {
		
		if (fileISize == null) {
			fileISize = createTempFile(getNameISize(), null);
		}

		return fileISize;
		
	}
	
	private String getNameLP() {
		return getCompositeFName()+ SFX_LP;
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
		return s+ SFX_MAPPED+ Constants.DOT+ SFX_BED;
	}
	
	
	private String getNameNotMappedReads() {
		String s= getCompositeFName();
		if (s== null)
			return null;
		return s+ SFX_NOTMAPPED+ Constants.DOT+ SFX_BED;
	}
	public File getFileMappedReads() {
		if (fileMappedReads == null) {
			fileMappedReads= createTempFile(
					getNameMappedReads(),
					SFX_BED);
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

            Log.progressFinish(StringUtils.OK, true);

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
						String lenString= i< BIN_LEN.length? Integer.toString(BIN_LEN[i]): "big";
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
							if (strand== STRAND_ENABLED) {
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
	
	public static byte FORMAT_BED= 0, FORMAT_GTF= 1, FORMAT_SAM= 2;
	public static byte mapFileType= FORMAT_SAM;
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
	

	
	public static final byte MODE_LEARN= 0, MODE_RECONSTRUCT= 1;
	public static final byte MODE_COUNT= 0, MODE_CONSTRAINT= 1, MODE_COMPLETE= 2;
	TProfileFunction func= new TProfileFunction(this); 
	private Vector<Thread> threadPool= new Vector<Thread>();
	int maxThreads= 1;
	Vector<String> origLines= null;
	int checkGTFscanExons= 0;
	public boolean explore(byte mode) {

		nrSingleTranscriptLoci= 0;
		nrReadsLoci= 0;
		
		BEDobject2[] leftover= null;
		
		SyncIOHandler2 handler= new SyncIOHandler2(10* 1024* 1024);
		
		if (mode== MODE_LEARN) {
			nrReadsSingleLoci= 0;
			nrReadsSingleLociMapped= 0;
		} 
		
		//System.out.println(System.getProperty("java.library.path"));
		long t0= System.currentTimeMillis();
		try {
			
			//this.gtfReader= null;
			//GFFReader gtfReader= getGTFreader();
			gtfReader.reset();
			if (bedWrapper== null)
				bedWrapper= getBedReader();
			else
				bedWrapper.reset();
			
			if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) {
				if (mode== MODE_LEARN) 
					System.err.println("\n[LEARN] Scanning the input and getting the attributes.");
				else if (mode== MODE_RECONSTRUCT)
					System.err.println("\n[SOLVE] Deconvolving reads of overlapping transcripts.");
			}
			final String profiling= "profiling ", decomposing= "decomposing "; 

				if (mode== MODE_LEARN)
					Log.progressStart(profiling);
				else if (mode== MODE_RECONSTRUCT)
					Log.progressStart(decomposing);


			if (mode== MODE_LEARN) 
				gtfReader.setKeepOriginalLines(false);
			else if (mode== MODE_RECONSTRUCT)
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
				if (mode== MODE_RECONSTRUCT)
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
					else if (mode== MODE_LEARN)
						continue;	// performance for not reading beds
					
					BEDobject2[] beds= null;

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
					
					if (beds!= null&& beds.length> 0) {
						
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
*/						 
							readObjects+= beds.length;
//						if (beds.length> 0&& mode== MODE_RECONSTRUCT)
//							System.err.println(gene[i].toUCSCString()+ " "+ beds.length);
					}
					
//					if (i>= 1500)
//						System.err.println("read "+beds.length+" objects");
					
					if (mode== MODE_LEARN&& beds!= null&& beds.length> 0) {
//						if (Constants.progress!= null) 
//							Constants.progress.setString(profiling+ gene[i].getGeneID());
						solve(gene[i], beds, false);
					}
					else if (mode== MODE_RECONSTRUCT) {

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
								" cluster "+ gene[i].getGeneID()+
								", "+beds.length+" reads.");
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
			if (mode== MODE_LEARN) {

				if (pairedEnd&& outputISize)
					writeISizes();
				if (pairedEnd&& func.getTProfiles()!= null) {
					insertMinMax= getInsertMinMax();
				}

				outputMappingStats(System.err, System.currentTimeMillis()- t0);

				if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) {
					
					//System.err.println(" OK.");
/*					System.err.println("\tfirst round finished .. took "+((System.currentTimeMillis()- t0)/ 1000)+ " sec.\n\n\t"
							+ nrSingleTranscriptLoci+" single transcript loci\n\t"							
							+ getBedReader().getNrLines()+ " mappings\n\t"
							+ nrReadsSingleLoci+" mappings in and around these loci(+/-"+tolerance+"nt)\n\t"
							+ nrReadsSingleLociMapped+" mappings map to annotation\n\t"
							+ ((strand== STRAND_SPECIFIC)?(nrMappingsWrongStrandWS+ nrMappingsWrongStrandWA+ nrMappingsWrongStrandCS+ nrMappingsWrongStrandCA)+" mappings map to annotation in antisense direction,\n\t":"")
							//+ (pairedEnd?(nrReadsSingleLociPotentialPairs+ " mappings form potential pairs,\n\t"):"")
							+ (pairedEnd?(nrReadsSingleLociPairsMapped* 2)+" mappings in annotation-mapped pairs\n\t":"")
							//+ nrReadsSingleLociNoAnnotation+ " mappings do NOT match annotation,\n\t"
							//+ (uniform?"":func.profiles.size()+" profiles collected\n\t")
							+ readLenMin+ ","+ readLenMax+ " min/max read length\n\t"							
							+ (pairedEnd&& insertMinMax!= null?insertMinMax[0]+","+insertMinMax[1]+" min/max insert size\n\t":""));
*/							
					//nrUniqueReads= getBedReader().getNrUniqueLinesRead();
					//System.err.println("\ttotal lines in file "+nrUniqueReads);
					System.err.println();
				}
				
			} else if (mode== MODE_RECONSTRUCT) {
				while (threadPool.size()> 0&& threadPool.elementAt(0).isAlive())
					try {
						threadPool.elementAt(0).join();
					} catch (Exception e) {
						; //:)
					}
				getWriter().flush();
				getWriter().close();

				outputMappingStats(System.err, System.currentTimeMillis()- t0);
//					if (fileMappings!= null)
//						getSammy().close();
				
				//assert(nrUniqueReads==getBedReader().getNrUniqueLinesRead());	// take out for cheat
			}
			
		} catch (Exception e1) {
			e1.printStackTrace();
			return false;
		}

		return true;
	}
	
	private void outputMappingStats(PrintStream p, long time) {
		
		if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) {
			p.println();
			int allMapped= nrMappingsSingleWSt+ nrMappingsSingleWAt+ nrMappingsSingleCSt+ nrMappingsSingleCAt
							+ nrMappingsSingleWSg+ nrMappingsSingleWAg+ nrMappingsSingleCSg+ nrMappingsSingleCAg,
				allWatson= nrMappingsSingleWSt+ nrMappingsSingleWAt+ nrMappingsSingleWSg+ nrMappingsSingleWAg,
				allCrick= nrMappingsSingleCSt+ nrMappingsSingleCAt+ nrMappingsSingleCSg+ nrMappingsSingleCAg,
				allSense= nrMappingsSingleWSt+ nrMappingsSingleCSt+ nrMappingsSingleWSg+ nrMappingsSingleCSg,
				allAsense= nrMappingsSingleWAt+ nrMappingsSingleCAt+ nrMappingsSingleWAg+ nrMappingsSingleCAg,
				allWS= nrMappingsSingleWSt+ nrMappingsSingleWSg,
				allWA= nrMappingsSingleWAt+ nrMappingsSingleWAg,
				allCS= nrMappingsSingleCSt+ nrMappingsSingleCSg,
				allCA= nrMappingsSingleCAt+ nrMappingsSingleCAg,
				allWrongStrand= nrMappingsWrongStrandWSt+ nrMappingsWrongStrandWAt+ nrMappingsWrongStrandCSt+ nrMappingsWrongStrandCAt
								+ nrMappingsWrongStrandWSg+ nrMappingsWrongStrandWAg+ nrMappingsWrongStrandCSg+ nrMappingsWrongStrandCAg;
			if (pairedEnd) {
				// pairs mapped
				allMapped+= nrMappingsPairsMappedWt+ nrMappingsPairsMappedCt+ nrMappingsPairsMappedWg+ nrMappingsPairsMappedCg;
				allWatson+= nrMappingsPairsMappedWt+ nrMappingsPairsMappedWg;
				allWS+= nrMappingsPairsMappedWt+ nrMappingsPairsMappedWg;
				allWA+= nrMappingsPairsMappedWt+ nrMappingsPairsMappedWg;
				allCrick+= nrMappingsPairsMappedCt+ nrMappingsPairsMappedCg;
				allCS+= nrMappingsPairsMappedCt+ nrMappingsPairsMappedCg;
				allCA+= nrMappingsPairsMappedCt+ nrMappingsPairsMappedCg;
				allSense+= nrMappingsPairsMappedWt+ nrMappingsPairsMappedCt+ nrMappingsPairsMappedWg+ nrMappingsPairsMappedCg;
				allAsense+= nrMappingsPairsMappedWt+ nrMappingsPairsMappedCt+ nrMappingsPairsMappedWg+ nrMappingsPairsMappedCg;
				// pairs w/o tx evidence
				allMapped+= nrMappingsPairsWoTxEvidenceWt+ nrMappingsPairsWoTxEvidenceCt+ nrMappingsPairsWoTxEvidenceWg+ nrMappingsPairsWoTxEvidenceCg;
				allWatson+= nrMappingsPairsWoTxEvidenceWt+ nrMappingsPairsWoTxEvidenceWg;
				allWS+= nrMappingsPairsWoTxEvidenceWt+ nrMappingsPairsWoTxEvidenceWg;
				allWA+= nrMappingsPairsWoTxEvidenceWt+ nrMappingsPairsWoTxEvidenceWg;
				allCrick+= nrMappingsPairsWoTxEvidenceCt+ nrMappingsPairsWoTxEvidenceCg;
				allCS+= nrMappingsPairsWoTxEvidenceCt+ nrMappingsPairsWoTxEvidenceCg;
				allCA+= nrMappingsPairsWoTxEvidenceCt+ nrMappingsPairsWoTxEvidenceCg;
				allSense+= nrMappingsPairsWoTxEvidenceWt+ nrMappingsPairsWoTxEvidenceCt+ nrMappingsPairsWoTxEvidenceWg+ nrMappingsPairsWoTxEvidenceCg;
				allAsense+= nrMappingsPairsWoTxEvidenceWt+ nrMappingsPairsWoTxEvidenceCt+ nrMappingsPairsWoTxEvidenceWg+ nrMappingsPairsWoTxEvidenceCg;
				// pairs in wrong orientation
				allMapped+= nrMappingsPairsWrongOrientationWt+ nrMappingsPairsWrongOrientationCt+ nrMappingsPairsWrongOrientationWg+ nrMappingsPairsWrongOrientationCg;
				allWatson+= nrMappingsPairsWrongOrientationWt+ nrMappingsPairsWrongOrientationWg;
				allWS+= nrMappingsPairsWrongOrientationWt+ nrMappingsPairsWrongOrientationWg;
				allWA+= nrMappingsPairsWrongOrientationWt+ nrMappingsPairsWrongOrientationWg;
				allCrick+= nrMappingsPairsWrongOrientationCt+ nrMappingsPairsWrongOrientationCg;
				allCS+= nrMappingsPairsWrongOrientationCt+ nrMappingsPairsWrongOrientationCg;
				allCA+= nrMappingsPairsWrongOrientationCt+ nrMappingsPairsWrongOrientationCg;
				allSense+= nrMappingsPairsWrongOrientationWt+ nrMappingsPairsWrongOrientationCt+ nrMappingsPairsWrongOrientationWg+ nrMappingsPairsWrongOrientationCg;
				allAsense+= nrMappingsPairsWrongOrientationWt+ nrMappingsPairsWrongOrientationCt+ nrMappingsPairsWrongOrientationWg+ nrMappingsPairsWrongOrientationCg;
			}
			if (stranded) {
				allMapped+= nrMappingsWrongStrandWSt+ nrMappingsWrongStrandWAt+ nrMappingsWrongStrandCSt+ nrMappingsWrongStrandCAt
							+ nrMappingsWrongStrandWSg+ nrMappingsWrongStrandWAg+ nrMappingsWrongStrandCSg+ nrMappingsWrongStrandCAg;
				allWatson+= nrMappingsWrongStrandWSt+ nrMappingsWrongStrandWAt+ nrMappingsWrongStrandWSg+ nrMappingsWrongStrandWAg;
				allWS+= nrMappingsWrongStrandWSt+ nrMappingsWrongStrandWSg;
				allWA+= nrMappingsWrongStrandWAt+ nrMappingsWrongStrandWAg;
				allCrick+= nrMappingsWrongStrandCSt+ nrMappingsWrongStrandCAt+ nrMappingsWrongStrandCSg+ nrMappingsWrongStrandCAg;
				allCS+= nrMappingsWrongStrandCSt+ nrMappingsWrongStrandCSg;
				allCA+= nrMappingsWrongStrandCAt+ nrMappingsWrongStrandCAg;
				allSense+= nrMappingsWrongStrandWSt+ nrMappingsWrongStrandCSt+ nrMappingsWrongStrandWSg+ nrMappingsWrongStrandCSg;
				allAsense+= nrMappingsWrongStrandWAt+ nrMappingsWrongStrandCAt+ nrMappingsWrongStrandWAg+ nrMappingsWrongStrandCAg;
			}
			int all= getBedReader().getNrLines();
			DecimalFormat fmt = new DecimalFormat("###,###,###");
			String allStr = fmt.format(all);
			p.println("\treconstruction finished .. took "+(time/ 1000)+ " sec.\n\n");
			p.println("\t"+ allStr+ " ("+ StringUtils.append(' ', StringUtils.fprint(all * 100f / all, 2), 6, true)
					+ "%) mappings in file");
			p.println("\t"+ StringUtils.append(' ', fmt.format(allMapped), allStr.length(), true)
					+ " ("+ StringUtils.append(' ', StringUtils.fprint(allMapped * 100f / all, 2), 6, true)
					+ "%) mappings to annotation");
			if (Constants.verboseLevel> Constants.VERBOSE_NORMAL) {
				p.println("\t\tW "+ fmt.format(allWatson)+ "; C "+ fmt.format(allCrick)
						+ "; S "+ fmt.format(allSense)+ "; A "+ fmt.format(allAsense));
				p.println("\t\tWS "+ fmt.format(allWS)+ "; WA "+ fmt.format(allWA)+ "; CS "+ fmt.format(allCS)+ "; CA "+ fmt.format(allCA));
			}
			if (stranded&& (allWrongStrand)> 0) {
				p.println("\t"+ StringUtils.append(' ', fmt.format(allWrongStrand), allStr.length(), true)
						+" ("+ StringUtils.append(' ', StringUtils.fprint((allWrongStrand) * 100f / all, 2), 6, true)+ "%) mappings with wrong strand");
				if (Constants.verboseLevel> Constants.VERBOSE_NORMAL) {
					p.println("\t\tW "+ fmt.format(nrMappingsWrongStrandWSt+ nrMappingsWrongStrandWAt+ nrMappingsWrongStrandWSg+ nrMappingsWrongStrandWAg)+ "; C "+ fmt.format(nrMappingsWrongStrandCSt+ nrMappingsWrongStrandCAt+ nrMappingsWrongStrandCSg+ nrMappingsWrongStrandCAg)
							+ "; S "+ fmt.format(nrMappingsWrongStrandWSt+ nrMappingsWrongStrandCSt+ nrMappingsWrongStrandWSg+ nrMappingsWrongStrandCSg)+ "; A "+ fmt.format(nrMappingsWrongStrandWAt+ nrMappingsWrongStrandCAt+ nrMappingsWrongStrandWAg+ nrMappingsWrongStrandCAg));
					p.println("\t\tWS "+ fmt.format(nrMappingsWrongStrandWSt+ nrMappingsWrongStrandWSg)+ "; WA "+ fmt.format(nrMappingsWrongStrandWAt+ nrMappingsWrongStrandWAg)+ "; CS "
							+ fmt.format(nrMappingsWrongStrandCSt+ nrMappingsWrongStrandCSg)+ "; CA "+ fmt.format(nrMappingsWrongStrandCAt+ nrMappingsWrongStrandCAg));
					p.println("\t\tT "+ fmt.format(nrMappingsWrongStrandWSt+ nrMappingsWrongStrandWAt+ nrMappingsWrongStrandCSt+ nrMappingsWrongStrandCAt)
							+ "; G "+ fmt.format(nrMappingsWrongStrandWSg+ nrMappingsWrongStrandWAg+ nrMappingsWrongStrandCSg+ nrMappingsWrongStrandCAg));
					p.println("\t\tTW "+ fmt.format(nrMappingsWrongStrandWSt+ nrMappingsWrongStrandWAt)+ "; TC "+ fmt.format(nrMappingsWrongStrandCSt+ nrMappingsWrongStrandCAt)
							+ "; TS "+ fmt.format(nrMappingsWrongStrandWSt+ nrMappingsWrongStrandCSt)+ "; TA "+ fmt.format(nrMappingsWrongStrandWAt+ nrMappingsWrongStrandCAt));
					p.println("\t\tGW "+ fmt.format(nrMappingsWrongStrandWSg+ nrMappingsWrongStrandWAg)+ "; GC "+ fmt.format(nrMappingsWrongStrandCSg+ nrMappingsWrongStrandCAg)
							+ "; GS "+ fmt.format(nrMappingsWrongStrandWSg+ nrMappingsWrongStrandCSg)+ "; GA "+ fmt.format(nrMappingsWrongStrandWAg+ nrMappingsWrongStrandCAg));
				}				
			}
			if (pairedEnd) {
				int allWrongOrientation= nrMappingsPairsWrongOrientationWt+ nrMappingsPairsWrongOrientationCt+ nrMappingsPairsWrongOrientationWg+ nrMappingsPairsWrongOrientationCg;
				if (allWrongOrientation> 0) {
					p.println("\t"+ StringUtils.append(' ', fmt.format(allWrongOrientation), allStr.length(), true)
							+ " ("+ StringUtils.append(' ', StringUtils.fprint(allWrongOrientation * 100f / all, 2), 6, true)+"%) mappings show wrong pair orientation");
					if (Constants.verboseLevel> Constants.VERBOSE_NORMAL)  {
						p.println("\t\tW "+ fmt.format(nrMappingsPairsWrongOrientationWt+ nrMappingsPairsWrongOrientationWg)+ "; C "+ fmt.format(nrMappingsPairsWrongOrientationCt+ nrMappingsPairsWrongOrientationCg));
						p.println("\t\tTW "+ fmt.format(nrMappingsPairsWrongOrientationWt)+ "; TC "+ fmt.format(nrMappingsPairsWrongOrientationCt));
						p.println("\t\tGW "+ fmt.format(nrMappingsPairsWrongOrientationWg)+ "; GC "+ fmt.format(nrMappingsPairsWrongOrientationCg));
					}
				}
				int allPairsWoTx= nrMappingsPairsWoTxEvidenceWt+ nrMappingsPairsWoTxEvidenceCt+ nrMappingsPairsWoTxEvidenceWg+ nrMappingsPairsWoTxEvidenceCg; 
				if (allPairsWoTx> 0) {
					p.println("\t"+ StringUtils.append(' ', fmt.format(allPairsWoTx), allStr.length(), true)
							+ " ("+ StringUtils.append(' ', StringUtils.fprint(allPairsWoTx * 100f / all, 2), 6, true)+ "%) mappings without support for pair");
					if (Constants.verboseLevel> Constants.VERBOSE_NORMAL) {
						p.println("\t\tW "+ fmt.format(nrMappingsPairsWoTxEvidenceWt+ nrMappingsPairsWoTxEvidenceWg)+ "; C "+ fmt.format(nrMappingsPairsWoTxEvidenceCt+ nrMappingsPairsWoTxEvidenceCg));
						p.println("\t\tTW "+ fmt.format(nrMappingsPairsWoTxEvidenceWt)+ "; TC "+ fmt.format(nrMappingsPairsWoTxEvidenceCt));
						p.println("\t\tGW "+ fmt.format(nrMappingsPairsWoTxEvidenceWg)+ "; GC "+ fmt.format(nrMappingsPairsWoTxEvidenceCg));
					}
				}
				int allSingle= nrMappingsSingleWSt+ nrMappingsSingleWAt+ nrMappingsSingleCSt+ nrMappingsSingleCAt+
								nrMappingsSingleWSg+ nrMappingsSingleWAg+ nrMappingsSingleCSg+ nrMappingsSingleCAg;
				if (allSingle> 0) {
					p.println("\t"+ StringUtils.append(' ', fmt.format(allSingle), allStr.length(), true)
							+ " ("+ StringUtils.append(' ', StringUtils.fprint(allSingle * 100f / all, 2), 6, true)+ "%) mappings without mate");
					if (Constants.verboseLevel> Constants.VERBOSE_NORMAL) {
						p.println("\t\tW "+ fmt.format(nrMappingsSingleWSt+ nrMappingsSingleWAt+ nrMappingsSingleWSg+ nrMappingsSingleWAg)+ "; C "+ fmt.format(nrMappingsSingleCSt+ nrMappingsSingleCAt+ nrMappingsSingleCSg+ nrMappingsSingleCAg)
								+ "; S "+ fmt.format(nrMappingsSingleWSt+ nrMappingsSingleCSt+ nrMappingsSingleWSg+ nrMappingsSingleCSg)+ "; A "+ fmt.format(nrMappingsSingleWAt+ nrMappingsSingleCAt+ nrMappingsSingleWAg+ nrMappingsSingleCAg));
						p.println("\t\tWS "+ fmt.format(nrMappingsSingleWSt+ nrMappingsSingleWSg)+ "; WA "+ fmt.format(nrMappingsSingleWAt+ nrMappingsSingleWAg)+ "; CS "+ fmt.format(nrMappingsSingleCSt+ nrMappingsSingleCSg)+ "; CA "+ fmt.format(nrMappingsSingleCAt+ nrMappingsSingleCAg));
						p.println("\t\tTW "+ fmt.format(nrMappingsSingleWSt+ nrMappingsSingleWAt)+ "; TC "+ fmt.format(nrMappingsSingleCSt+ nrMappingsSingleCAt));
						p.println("\t\tGW "+ fmt.format(nrMappingsSingleWSg+ nrMappingsSingleWAg)+ "; GC "+ fmt.format(nrMappingsSingleCSg+ nrMappingsSingleCAg));
					}
				}
				int allPairs= nrMappingsPairsMappedWt+ nrMappingsPairsMappedCt+ nrMappingsPairsMappedWg+ nrMappingsPairsMappedCg;
				p.println("\t"+ StringUtils.append(' ', fmt.format(allPairs), allStr.length(), true)
						+" ("+ StringUtils.append(' ', StringUtils.fprint(allPairs * 100f / all, 2), 6, true)+ "%) mappings correctly mapped");
				if (Constants.verboseLevel> Constants.VERBOSE_NORMAL) {
					p.println("\t\tW "+ fmt.format(nrMappingsPairsMappedWt+ nrMappingsPairsMappedWg)+ "; C "+ fmt.format(nrMappingsPairsMappedCt+ nrMappingsPairsMappedCg));
					p.println("\t\tTW "+ fmt.format(nrMappingsPairsMappedWt)+ "; TC "+ fmt.format(nrMappingsPairsMappedCt));
					p.println("\t\tGW "+ fmt.format(nrMappingsPairsMappedWg)+ "; GC "+ fmt.format(nrMappingsPairsMappedCg));
				}
			}
			//+ nrMultiMaps+" mapped multiply.\n\n\t"
			p.println((outputGene?"\n\t"+ nrLoci+ " loci, "+ nrLociExp+ " detected":"")
					+ (outputTranscript?"\n\t"+nrTx+" transcripts, "+nrTxExp+" detected":"")
					+ (outputEvent?"\n\t"+ nrEvents+" ASevents of dimension "+eventDim+", "+nrEventsExp+" detected":"")
					+ "\n"
					//+ nrUnsolved+" unsolved systems."
					);
		}					

		
	}

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
		outputEvent= true, 
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
	private void solve(Gene gene, BEDobject2[] beds, boolean decompose) {
		
		// create LP and solve
		LocusSolver3 lsolver= new LocusSolver3();
		lsolver.init(gene, beds, decompose);
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
			lsolver.start(); //;
			
		} else
			lsolver.run();
	}
	
	private Thread getLastThread() {
		synchronized(FluxCapacitorNew.this.threadPool) {
			if (this.threadPool.size()> 0)
				return this.threadPool.get(this.threadPool.size()- 1);
			else 
				return null;
		}
	}

	
	public static final String GFF_FEATURE_JUNCTION= "junction", GFF_FEATURE_PAIRED= "paired",
		GFF_FEATURE_FRAGMENT= "fragment";
	
	public static final String GTF_ATTRIBUTE_PVAL= "falsification";
	
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
		nrLociUnderPredicted= 0,
		nrLociOverPredicted= 0,
		nrTx= 0, 
		nrEvents= 0, 
		nrEventsExp= 0, 
		nrTxExp= 0, 
		nrReadsLoci= 0, 
		nrReadsSingleLociMapped= 0, 
		nrReadsSingleLociPairsMapped= 0, 
		nrUnsolved= 0, 
		nrReadsSingleLoci= 0,
		nrSingleTranscriptLoci= 0,	// counted in learn AND decompose redundantly
		nrSingleTranscriptLearn= 0,	
		nrReadsSingleLociNoAnnotation= 0,
		
		nrMultiMaps= 0,
		// new mapping counters
		nrMappingsSingleWSt= 0, nrMappingsSingleWAt= 0, nrMappingsSingleCSt= 0, nrMappingsSingleCAt= 0,
		nrMappingsSingleWSg= 0, nrMappingsSingleWAg= 0, nrMappingsSingleCSg= 0, nrMappingsSingleCAg= 0,
		nrMappingsWrongStrandWSt= 0, nrMappingsWrongStrandWAt= 0, nrMappingsWrongStrandCSt= 0, nrMappingsWrongStrandCAt= 0,
		nrMappingsWrongStrandWSg= 0, nrMappingsWrongStrandWAg= 0, nrMappingsWrongStrandCSg= 0, nrMappingsWrongStrandCAg= 0,
		nrMappingsPairsMappedWt= 0, nrMappingsPairsMappedCt= 0,
		nrMappingsPairsMappedWg= 0, nrMappingsPairsMappedCg= 0,
		nrMappingsPairsWoTxEvidenceWt= 0, nrMappingsPairsWoTxEvidenceCt= 0,
		nrMappingsPairsWoTxEvidenceWg= 0, nrMappingsPairsWoTxEvidenceCg= 0,
		nrMappingsPairsWrongOrientationWt= 0, nrMappingsPairsWrongOrientationCt= 0,
		nrMappingsPairsWrongOrientationWg= 0, nrMappingsPairsWrongOrientationCg= 0;
	
	boolean 
		map= false, 
		decompose= false, 
		uniform= false;
	long nrReadsAll= 0;
	double costModelPar= Double.NaN;
	float[] costBounds= new float[] {0.95f, Float.NaN};	// how much of the original observation can be subs/add
	private static final int BIG= 999999999;
	public static final String VALUE_NA= "0";	// NA
	
	int[] profileBoundaries;
	private int getBinIdx(int len) {
		int p= Arrays.binarySearch(profileBoundaries, len);
		p= (p<0)?-p-1:p;
		return p;
	}
	
	BinVector isizeV; 
	synchronized void addInsertSize(int isize) {
		isizeV.incrTuple(isize);
	}
	
	private String getAttributeOF(double val, GraphLPsolver solver, int readCount) {

		StringBuilder sb= new StringBuilder(GTF_ATTRIBUTE_PVAL);
		sb.append(" \"");
		if (val > BIG|| solver== null) {
			if (val> BIG)
				System.currentTimeMillis();
			sb.append(" \""+VALUE_NA+ "\";" );
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

	private static final String FLOAT_STRING_0= "0.0";
	private Vector<Vector<Edge>> eeV= new Vector<Vector<Edge>>();
	
	static AtomicLong along= new AtomicLong((0L ^ 0x5DEECE66DL) & ((1L << 48) - 1));
	
	/**
	 * @deprecated
	 * @return
	 */
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
		double c= 0; 	// TODO was global
        double dbl= (l / (double)(1L << 53))*((c++% 2== 0)?100:1); 
        
		return Math.max(0.001, dbl);
	}
		

	private static final double GAUSS_FACTOR = 1d/ Math.sqrt(2d* Math.PI);
	
	public static final byte BYTE_0= (byte) 0, BYTE_1= (byte) 1, BYTE_MINUS_1= (byte) (-1);
	private static final String NULL_STR= "0";
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
			Transcript.setEdgeConfidenceLevel((byte) -1);	// trust no one
			
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

	public static int[] BIN_EXP= new int[] {10, 100};
	public static int[] BIN_LEN= new int[] {1000, 2000};
	int[][] profileStub, profileStubRev;
	/**
	 * profile managing the matrices
	 */
	Profile profile;

	public final static String SFX_LPOUT= ".lp";
	private BEDobject2[] readBedFile(Gene gene, int from, int to, byte mode) {
		
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
//		if (gene.getGeneID().equals("chr19:1609293-1652326C"))
//			System.currentTimeMillis();
		if (beds== null)
			return null;
		
		return beds;
		
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

	private static String FNAME_PROPERTIES= "capacitor.prop";
	static void readProperties() {
		String wrapper= System.getProperty(Constants.PROPERTY_KEY_WRAPPER_BASE);
		if (wrapper== null)
			wrapper= FNAME_PROPERTIES;
		else
			wrapper+= File.separator+ FNAME_PROPERTIES;
		File pFile= new File(wrapper);		
		Properties props= new Properties();
		try {
			props.load(new FileInputStream(pFile));
		} catch (Exception e) {
			; // :)
		}
	
		if (props.containsKey(PROPERTY_BUILD))
			version= (String) props.get(PROPERTY_BUILD);
		if (props.containsKey(PROPERTY_JDK)) {
			try {
				float v= Float.parseFloat((String) props.get(PROPERTY_JDK));
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

	public static int loadLibraries() {
		
		if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) 
			System.err.println("[PRE-CHECK] I am checking availability of the required lpsolve JNI libs.");
//		if (SystemInspector.checkRuntime())
//			miss= true;
//		if (SystemInspector.checkRuntimeCirco()) {
//			System.exit(-1);
//		}
	
		// check to load from java.library.path 
		VersionInfo lpVer= null;
		try {
			System.loadLibrary(LPSOLVE_LIB_NAME);
			System.loadLibrary(LPSOLVE_JNI_NAME);
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
		String nativePath= System.getProperty(CLI_LONG_LIB);
		if (nativePath== null) {
			nativePath= basePath+ File.separator+ SUBDIR_NATIVELIBS+ File.separator+ SUBDIR_LPSOLVE+ File.separator+ 
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
							shell= SHELL_BASH;
						else if (s.endsWith("csh"))
							shell= SHELL_CSH;
						else if (s.endsWith("ksh"))
							shell= SHELL_KSH;
					}
					p.destroy();
				} catch (Exception e) {
					; // :)
				}
				
				if (shell== SHELL_CSH)
					writer.write("setenv ");
				else
					writer.write("export ");
				
				if (SystemInspector.getOSgroup()== SystemInspector.OS_GROUP_MACOSX)
					writer.write("DY");
				writer.write("LD_LIBRARY_PATH");
				if (shell== SHELL_CSH)
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

}
