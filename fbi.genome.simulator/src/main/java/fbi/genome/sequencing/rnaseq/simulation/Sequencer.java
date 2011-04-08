package fbi.genome.sequencing.rnaseq.simulation;

import fbi.commons.ByteArrayCharSequence;
import fbi.commons.Log;
import fbi.commons.file.FileHelper;
import fbi.commons.io.IOHandler;
import fbi.commons.io.IOHandlerFactory;
import fbi.commons.thread.StoppableRunnable;
import fbi.commons.thread.ThreadedQWriter;
import fbi.genome.io.BufferedBACSReader;
import fbi.genome.io.Fasta;
import fbi.genome.io.UnixStreamSort;
import fbi.genome.io.UnixStreamSort2;
import fbi.genome.io.UnixStreamSort2.DesignatedHierarchicalFieldComparator;
import fbi.genome.io.gff.GFFReader;
import fbi.genome.io.rna.FMRD;
import fbi.genome.model.*;
import fbi.genome.model.bed.BEDobject;
import fbi.genome.model.bed.BEDobject2;
import fbi.genome.model.constants.Constants;
import fbi.genome.sequencing.rnaseq.simulation.error.ModelPool;

import java.io.*;
import java.util.*;
import java.util.zip.ZipEntry;
import java.util.zip.ZipException;
import java.util.zip.ZipFile;
import java.util.zip.ZipOutputStream;

//import io.Sammy;
//import net.sf.samtools.SAMFileWriter;
//import net.sf.samtools.SAMFileWriterFactory;

public class Sequencer implements StoppableRunnable {
    public final static String NAME_ANN= "Annotation", NAME_EXP= "Expression", NAME_RTFRAG= "Library", NAME_SEQ= "Sequencing";
	static final String SFX_FASTA= "fasta", SFX_FASTQ= "fastq", SFX_SAM= "sam";
	static final byte[] DAVID_CORRECT= {40,40,40,40,40,40,40,40,40,40,40,40,39,38,37,36,36,35,34,33},
		DAVID_WRONG= {40,38,36,34,32,30,28,26,24,22,20,18,16,14,12,10,8,6,4,2};
	public static final String TAG_LID= "l", TAG_TID= "t", TAG_MOL= "m", TAG_READ= "r", TAG_FRAGLENGTH= "f", TAG_SEP=";", TAG_EQ= "=";
	FluxSimulatorSettings settings;
	boolean stop= false;
	long totalReads= -1, nrOfFrags= 0;
	int cntLoci, cntTrpts, cntExons;
	Hashtable<String,int[][]> mapFrags;
	ModelPool babes= null;
	boolean multiThread= false;
	boolean outputSAM= false;
	IOHandler rw;
	int cntPlus, cntMinus;
	
    public Sequencer(FluxSimulatorSettings settings) {
		this.settings= settings;
	}

    class ZipperThread extends Thread {
    	
    	BufferedBACSReader in;
    	ZipOutputStream out;
    	ByteArrayCharSequence cs;
    	long totBytes;
    	int linesRec;
    	
    	public ZipperThread(InputStream istream, OutputStream ostream, int lineLength, long bytesExpected) {
			super("ZipperThread");
			in= new BufferedBACSReader(istream, lineLength* 10);
			out= new ZipOutputStream(ostream);
			cs= new ByteArrayCharSequence(lineLength);
			totBytes= bytesExpected;
		}
    	
    	@Override
    	public void run() {
			
			ByteArrayCharSequence lastID= null;
			linesRec= 0;
			ZipEntry ze= null;
			int cnt= 0;
			long currBytes= 0;
			int lastPerc= 0;
			while(in.readLine(cs)> 0&& !stop) {
				++linesRec;
				currBytes+= cs.length()+ 1;
                Log.progress(currBytes, totBytes);
				try {
					lastID= writeOutZip(cs, out, lastID);
				} catch (Exception e) {
					e.printStackTrace();
					break;
				}
				int pEnd= cs.p1- 1;
				try {
					out.write(cs.a, cs.start, pEnd- cs.start);		// only write the coordinates
					out.write(BYTE_NL);
				} catch (IOException e) {
					e.printStackTrace();
					break;
				}
				++cnt;
			}
			try {
				out.flush();
				out.close();
			} catch (IOException e) {
				if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
					e.printStackTrace();
				return;
			}

			
    	}
    }
    
    private static int ctrProcessors= 0;
    class Processor extends Thread {

    	Gene gene;
    	Object lock= new Object();
		BEDobject2 obj;
		ByteArrayCharSequence cs;
		boolean stop= false;

    	public Processor() {
    		super("Sequencing Processor "+ ++ctrProcessors);
    		cs= new ByteArrayCharSequence(128);
    		obj= new BEDobject2(128);
		}
    	
    	public void setGene(Gene g) {
			this.gene= g;
    	}
    	
    	public void close() {
    		stop= true;
    		while (isAlive())
        		interrupt();
	    		try {
					join();
				} catch (InterruptedException e) {
					; // :)
				}
    	}
    	
    	@Override
    	public void run() {
    		
    		while((!isStop())) {	// 20101201 killed && !stop, do not do that here, called by close() (gene loss)
    			if (isAlive()) {
		    		synchronized (lock) {
		        		while(!stop&& gene== null)
							try {
								lock.wait();
							} catch (InterruptedException e) {
								; 
							}
		    		}
    			}
	    		
	    		// alternative to rw
	    		//BufferedByteArrayReader buffy= new BufferedByteArrayReader();
	    		String baseID= null;
	    		if (gene== null)
	    			break;
	    		else
	    			baseID= gene.getGeneID()+ FluxSimulatorSettings.SEP_LOC_TID;
				for (int j = 0; (!isStop())&& j < gene.getTranscripts().length; j++) {	// 20101202: killed (!stop)&& , do not do that here, called by close() (gene loss)
					Transcript t= gene.getTranscripts()[j];
					String compID= baseID+ t.getTranscriptID();
					
					ZipEntry ze= zipHash.remove(compID);	// t.getChromosome()+ (char) Fragmenter.BYTE_SEP_LC_TX+ 
					if (ze== null)
						continue;	// not in frg file
					
					try {
						InputStream is= zzFile.getInputStream(ze);
						//rw.addStream(is);	// use inputstream for rw
						
						// 20101215 BACS reader not suitable
						//BufferedBACSReader buffy= new BufferedBACSReader(is);
						BufferedReader buffy= new BufferedReader(new InputStreamReader(is));
						int k= 0;

						//while(rw.readLine(is, cs)> 0) {	// (buffy.readLine(is, cs)> 0)
						//while(buffy.readLine(cs)> 0) {
						String s= null;
						while((s= buffy.readLine())!= null) {
							incrementFragCtr(t);
							cs.init(s);
							int fstart= cs.getTokenInt(0);
							int fend= cs.getTokenInt(1);
							
							double r= rnd.nextDouble();
							if (r< p) { 
								double q= rndFiftyFifty.nextDouble();
								obj.reset();
								if (settings.isPairedEnd()|| q< 0.5) {
									process(true, t, fstart, fend, k);
									++cntPlus;
								}
								if (settings.isPairedEnd()|| q> 0.5) {									
									process(false, t, fstart, fend, k);
									++cntMinus;
								}
							} 
							++k;
						}
						//rw.removeStream(is);
						buffy.close(); 	
						is.close();
						
					} catch (IOException e) {
						e.printStackTrace();
					}
				}
				this.gene= null;
    		}
    	}



		private void process(boolean left, Transcript t, int fstart, int fend, int k) throws IOException {
			
			byte absDir= (byte) (t.getStrand()>= 0?1: -1),
				antiDir= (byte) (t.getStrand()>= 0?-1: 1);
			int rLen= settings.getReadLength();
			int flen= fend- fstart+ 1;
			
			 ++totalReads;
			 if (flen< rLen)
				 ++cntTruncReads;
			 String id= t.getGene().getGeneID()+ FluxSimulatorSettings.SEP_LOC_TID+ t.getTranscriptID();
			 if (map.containsKey(id))
				 map.put(id, map.get(id)+ 1);
			 else 
				 map.put(id, long0);
			
			 
			 // bed object
			 if (left) {
				 createRead(obj, 
						 fstart, Math.min(fstart+ settings.getReadLength()- 1, fend), 	// start, end
						 t, null, k, absDir, 
						 fstart, fend, left);
			 } else {
				 createRead(obj, 
						 Math.max(fend- settings.getReadLength()+ 1, fstart), fend, 	// start, end
						 t, null, k, antiDir, 
						 fstart, fend, left);
			 }
			 
			 rw.writeLine(obj, oBed);
//			 try {
//				cs.append('\n');
//				oBed.write(cs.a, cs.start, cs.length());
//			} catch (IOException e) {
//				System.err.println("\n[ERROR] write error during processing");
//				e.printStackTrace();
//			}
			 
			 
			 // fasta seq
			 if (settings.isFastQ()&& Graph.overrideSequenceDirPath!= null) {
				 if (left) 
					 createQname(obj, cs, t, null, k, absDir, fstart, fend,
							 fstart,Math.min(fstart+ rLen- 1, fend));
				 else 
					 createQname(obj, cs, t, null, k, antiDir, fstart, fend,
							 Math.max(fend- rLen+ 1, fstart), fend);
				 
				 createQSeq(cs, obj, absDir, t.getStrand(), rLen, flen/*fstart, fend, t*/);
				 rw.writeLine(cs, oFasta);
//				 try {
//					cs.append('\n');
//					oFasta.write(cs.a, cs.start, cs.length());
//				 } catch (IOException e) {
//					System.err.println("\n[ERROR] write error during processing fasta");
//					e.printStackTrace();
//				 }	

			 } 
		}


    }
    
    
	private boolean convertToSAM() {
		/*Sammy sam= new Sammy(settings.seqFile, FileHelper.replaceSfx(settings.seqFile, Sammy.SFX_SAM), false);
		if (settings.fastQ)
			sam.setFastaQFile(getFASTAfile());
		sam.run();*/
		return true;
	}
	
	boolean init() {
		
		if (settings.getProFile()== null) {
			if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) 
				System.err.println("\t[OOPS] no input for sequencing");
			return false;
		} 
		if (settings.getSeqFile()== null) {
			if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) 
				System.err.println("\t[AIAIAI] no output for sequencing");
			return false;
		}
		if (settings.getSeqFile().exists()) {
			if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) {
				if (!FileHelper.checkForOverwrite(System.err, settings.getSeqFile()))
					return false;
			} else {
				settings.getSeqFile().delete();
			}			
		}
		
		
		if (settings.getReadLength()<= 0|| settings.getReadNr()<= 0
				|| (!settings.getProFile().exists())
				|| (!settings.getProFile().canRead())
				//|| (!settings.getSeqFile().getParentFile().canWrite()) // NPException for .pro file without path
				) {
			if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) {
				System.err.println("\t[TOOLITTLE] there are parameters missing: ");
				if ((!settings.getProFile().exists())|| (!settings.getProFile().canRead()))
					System.err.println("\t"+settings.getProFile().getAbsolutePath());
				if (!settings.getSeqFile().canWrite())
					System.err.println("\t"+settings.getSeqFile().getAbsolutePath());
				if (settings.getReadLength()<= 0)
					System.err.println("\t"+settings.PAR_SEQ_READ_LENGTH+" ");
				if (settings.getReadNr()<= 0)
					System.err.println("\t"+settings.PAR_SEQ_READ_NUMBER+" ");
				System.err.println("\n");
			}
			return false;
		} else {
			if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) {
				System.err.println("\t"+settings.PAR_PRO_FNAME+"\t"+settings.getProFile().getAbsolutePath());
				System.err.println("\t"+settings.PAR_SEQ_FNAME+"\t"+settings.getSeqFile().getAbsolutePath());
				System.err.println("\t"+settings.PAR_SEQ_READ_LENGTH+"\t"+settings.getReadLength());
				System.err.println("\t"+settings.PAR_SEQ_READ_NUMBER+"\t"+settings.getReadNr());
			}
		}
		return true;
	}

	double p= -1;
	public void run() {
		
		if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) 
			System.err.println("[SEQUENCING] getting the reads");
		
		if ((!init())|| stop) 
			return;		
		if ((!writeInitialFile())|| stop)
			return;
		sequence();
		
		//System.err.println("sequencing "+ (System.currentTimeMillis()- t0)/1000+ " sec. + "+cntPlus+", - "+cntMinus);
		
	}
	
	public boolean loadErrors() {
		
		// load model
		if (settings.getErrFile()!= null) {
//			String s= settings.getErrFile().getAbsolutePath();
			babes= ModelPool.read(settings.getErrFile(), settings);
			if (babes== null)
				return false;

            // check qualities for issued #48
            if(!babes.hasQualities() && settings.isFastQ()){
                settings.setFastQ(false);
                Log.warn("FastQ output requested, but the model does not support qualities. Disabled FastQ output!");
            }
			return true; 
		} 
		// else
		return false;
	}	

	public boolean loadStats() {
		if (!isFinished())
			return false;
		
		try {
            String s= "initializing sequencer ";
            Log.progressStart(s);

			totalReads= FileHelper.countLines(settings.getSeqFile().getAbsolutePath());
			
			
			// TODO OutOfMemoryError, do something else
			// Exception in thread "Blocking Thread" java.lang.OutOfMemoryError: Java heap space
			// at java.lang.String.substring(String.java:1770)
			// at java.util.StringTokenizer.nextToken(StringTokenizer.java:335)

/*			BufferedReader buffy= new BufferedReader(new FileReader(settings.getSeqFile()));
			totalReads= 0;
			HashSet hashTrp= new HashSet<String>(), hashLoc= new HashSet<String>();
			StringTokenizer st;
			long totBytes= settings.getSeqFile().length(), bytesRead= 0l;
			int perc= 0;
			for (String s; (s= buffy.readLine())!= null; ++totalReads) {
				bytesRead+= s.length()+ 1;
				if ((bytesRead* 10d/ totBytes)> perc) {
					if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) {
						if (Constants.progress!= null)
							Constants.progress.progress();
						else {
							System.err.print("*");
							System.err.flush();
						}
					}
					++perc;
				}
				st= new StringTokenizer(s, "\t"); 
				if (st.countTokens()< 4)
					continue;
				st.nextToken();
				st.nextToken();
				st.nextToken();
				st= new StringTokenizer(st.nextToken(), ";");	// TODO general name delimiter
				try {
					hashLoc.add(st.nextToken());
					hashTrp.add(st.nextToken());
				} catch (Exception e) {
					; // :)
				}
			}
			cntLoci= hashLoc.size();
			cntTrpts= hashTrp.size();
			buffy.close();
*/			
		} catch (Exception e) {
			return false;
		}
		
		return true;
	}
	
	Random rnd= new Random(), rndFiftyFifty= new Random();
	Hashtable<CharSequence,Long> map;
	private synchronized void incrementFragCtr(Transcript t) {
		++totalFrags;
		//hashTrp.add(t.getTranscriptID());
		//hashLoc.add(t.getGene().getGeneID());
	}
	
	int totalFrags= 0;
	Long long0= new Long(1);

	HashSet<String> hashTrp, hashLoc;
	File zipFile= null; //new File("N:\\tmp\\master.zip");
	ZipFile zzFile;
	Hashtable<CharSequence, ZipEntry> zipHash;
	void writeInitialFile_old() {
		try {
			String sss= "sequencing init";
            Log.progressStart(sss);

			BufferedReader buffy= new BufferedReader(new FileReader(settings.getFrgFile()));
			long bytesRead= 0, totBytes= settings.getFrgFile().length();			
			int cnt= 0, perc= 0;
			mapFrags= new Hashtable<String,int[][]>(settings.getProfiler().getIds().length);
			Vector<int[]> v= new Vector<int[]>();
			String[] token;
			String currentID= null;
			for (String s; 
				!stop&& (s= buffy.readLine())!= null; 
				++cnt, bytesRead+= s.length()+ System.getProperty("line.separator").length()) {
                Log.progress(bytesRead, totBytes);
				token= s.split(Fragmenter.FRG_FILE_TAB);
				
				if (!token[2].equals(currentID)) {
					if (v.size()> 0) {
						int[][] m= new int[v.size()][];
						for (int i = 0; i < m.length; i++) 
							m[i]= v.elementAt(i);
						v.removeAllElements();
						mapFrags.put(currentID, m);
					}
					currentID= token[2];
				}
				
				int[] p= new int[2];
				p[0]= Integer.parseInt(token[0]);
				p[1]= Integer.parseInt(token[1]);
				v.add(p);
			}
			buffy.close();
			
			nrOfFrags= cnt;
			p= settings.getReadNr()/ (double) cnt;

			if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) 
				System.err.println();
			
			if (!stop) {
				sequence();
			}
			
			//convertToSAM();
			
		} catch (Exception e) {
			return;
		}

	}
	
	int writeInitialFileSort_old(ByteArrayCharSequence cs) {
		// Sort
		try {
			String sss= "init sorting";
            Log.progressStart(sss);

			File inFile= settings.getFrgFile();
			InputStream is= new FileInputStream(inFile);
			File tmpOut= File.createTempFile("sim", "_sorted");
			final OutputStream os= new FileOutputStream(tmpOut);
			rw= IOHandlerFactory.getDefaultHandler();//new SyncIOHandler2(2);
			rw.addStream(is, 10* 1024);
			rw.addStream(os, 10* 1024);
			PipedOutputStream out = new PipedOutputStream();
			PipedInputStream in = new PipedInputStream(out);
			final BufferedWriter writer= new BufferedWriter(new OutputStreamWriter(
					out));
			fbi.genome.io.UnixStreamSort.DesignatedHierarchicalFieldComparator compi=
				new fbi.genome.io.UnixStreamSort.DesignatedHierarchicalFieldComparator(3);
			final UnixStreamSort sorter= new UnixStreamSort(in, compi);
			final byte[] bb= new byte[100];			
			Thread t= new Thread("Sorted Writer") {
				public void run() {
					while(!stop) {
						try {
							int len= sorter.getOutInStream().read(bb);
							if (len== -1)
								break;
							else
								rw.write(bb, 0, len, os);
						} catch (Exception e) {
							e.printStackTrace();
						}
					}
				};
			};
			// rw.start();
			
			t.start();
			sorter.start();
			BufferedOutputStream bout= new BufferedOutputStream(out);
			long totBytes= inFile.length(), currBytes= 0;
			int lastPerc= 0, cnt= 0;
			while(rw.readLine(is,cs)> 0) {
				++cnt;
				currBytes+= cs.length()+ 1;
                Log.progress(currBytes, totBytes);
				bout.write(cs.a, cs.start, cs.length());
				bout.write(BYTE_NL);
			}
			bout.close();			
			sorter.join();
			t.join();
			rw.close();
			//rw.join();
            Log.progressFinish();
			return cnt;
		} catch (Exception e) {
			e.printStackTrace();
			return -1;
		}
	}
	
	int writeInitialFileZip(ByteArrayCharSequence cs, File in, File out) {
		try {
			String sss= "init zipping";
            Log.progressStart(sss);

			FileInputStream inStram= new FileInputStream(in);
			rw= IOHandlerFactory.getDefaultHandler();//new SyncIOHandler2(2);
			rw.addStream(inStram);
//			if (FluxSimulatorSettings.optDisk)
//				rw.start();
			File zipF= new File(System.getProperty("java.io.tmpdir")+ File.separator
					+ "master.zip");
			ZipOutputStream zipOut= new ZipOutputStream(new BufferedOutputStream(new FileOutputStream(zipF)));
			//ByteArrayCharSequence cs= new ByteArrayCharSequence(50);
			ByteArrayCharSequence lastID= null;
			ZipEntry ze= null;
			int cnt= 0;
			long totBytes= in.length(), currBytes= 0;
			int lastPerc= 0;
			while(rw.readLine(inStram, cs)!= 0) {
				currBytes+= cs.length()+ 1;
                Log.progress(currBytes, totBytes);
				cs.resetFind();
				ByteArrayCharSequence id= cs.getToken(2);
				int pEnd= cs.p1- 1;
				if (!id.equals(lastID)) {
					try {
						zipOut.putNextEntry(new ZipEntry(id.toString()));
					} catch (ZipException e) {
						e.printStackTrace();
					}
					lastID= id.cloneCurrentSeq();
				}
				cs.append(BYTE_NL);
				zipOut.write(cs.a, cs.start, pEnd- cs.start);	// only write the coordinates
				++cnt;
			}
			rw.close();
			zipOut.close();
            Log.progressFinish();
			
			return cnt;
			
		} catch (Exception e) {
			e.printStackTrace();
			return -1;
		}
	}
	
	boolean writeInitialFile() {
		
		try {
			
			if (zipFile== null)
				zipFile= File.createTempFile("sim", "master.gz");
			
			ByteArrayCharSequence cs= new ByteArrayCharSequence(100);
			File inFile= settings.getFrgFile();
			nrOfFrags= sortAndZip(cs, inFile, zipFile);
			
			// create ZIP
			// hash entries
			if (zipFile== null|| !zipFile.exists())
				return false;
			zipHash= new Hashtable<CharSequence, ZipEntry>(settings.getProfiler().getIds().length);
			ZipFile zFile= new ZipFile(zipFile);
			Enumeration e= zFile.entries();
			ZipEntry ze;
			while (e.hasMoreElements()) {
				ze= (ZipEntry) e.nextElement();
				zipHash.put(ze.getName(), ze);
	
				// 20101215 BACS reader not good
/*				InputStream is= zFile.getInputStream(ze);
//				BufferedReader buffy= new BufferedReader(new InputStreamReader(is));
				BufferedBACSReader buffy= new BufferedBACSReader(ze.getSize(), is);
				System.err.print("\ntest "+ ze.getName());
				int k= 0;
				while(buffy.readLine(cs)> 0) {
//				while(buffy.readLine()!= null) {
					++k;
				}
				
				System.err.println("\t"+k+" lines.");
				
				is.close();
*/				
			}
			zFile.close();
			
			// stats
			p= settings.getReadNr()/ (double) nrOfFrags;

			return true;
			
		} catch (Throwable e) {
			if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
				e.printStackTrace();			
			return false;
		}
		
		
	}
	
	private Processor[] getProcessorPool(int nr) {
		if (processorPool == null) {
			processorPool = new Processor[nr];
			for (int i = 0; i < processorPool.length; i++) { 
				processorPool[i]= new Processor();
				if (nr> 1)
					processorPool[i].start();
			}
		}

		return processorPool;
	}
	
	Processor[] processorPool;
	FileOutputStream oBed, oFasta;
	boolean sequence() { 
		try {
			if (zipFile== null)
				return false;
			
			String sss= "sequencing";
            Log.progressStart(sss);

			GFFReader reader= new GFFReader(settings.getRefFile().getCanonicalPath());
			reader.setReadAheadLimit(500);
			reader.setSilent(true);
			reader.setStars(true);
			
			// vars
			Gene[] g;
			totalReads= 0;
			long totalMols= 0;
			totalFrags= 0;
			cntPolyA= 0;
			cntTruncReads= 0;
			
/*			Iterator<int[][]> iter= mapFrags.values().iterator();
			while(iter.hasNext())
				totalMols+= iter.next().length;
			for (int i = 0; i < settings.getProfiler().getMolecules().length; i++) 
				totalMols+= settings.getProfiler().getMolecules()[i];
*/			
			hashTrp= new HashSet<String>(); 
			hashLoc= new HashSet<String>();
//			ThreadedQWriter qwriter= null, qfasta= null;
//			BufferedWriter bwriter= null, bfasta= null;
			File tmpFile= File.createTempFile("flux",NAME_SEQ,settings.getTmpDir()), tmpFasta= null;
//			if (!multiThread)
//				bwriter= new BufferedWriter(new FileWriter(tmpFile));
			if (settings.isFastQ()&& settings.getGenDir()!= null) {
				tmpFasta= File.createTempFile("flux",NAME_SEQ,settings.getTmpDir());
				Graph.overrideSequenceDirPath= settings.getGenDir().getAbsolutePath();
//				if (!multiThread)
//					bfasta= new BufferedWriter(new FileWriter(tmpFasta));
			}
//			char[] seq= new char[settings.getReadLength()];
//			Random rndQual= new Random();
//			final char[] pA= new char[settings.getReadLength()], pT= new char[settings.getReadLength()];
//			for (int i = 0; i < pT.length; i++) {
//				pA[i]= 'a';
//				pT[i]= 't';
//			}
			//final String polyA= new String(pA), polyT= new String(pT);
//			byte[] sequals= null;
//			if (babes!= null) 
//				sequals= new byte[settings.getReadLength()];		

			// probability for fragment selection
			p= 1;
			if (settings.getReadNr()< nrOfFrags)
				p= settings.getReadNr()/ (double) nrOfFrags;	
			//System.err.println(settings.readNr+" reads, "+nrOfFrags+ " frags, p= "+p);
			map= new Hashtable<CharSequence,Long>(settings.getProfiler().getMolecules().length);
			
			// init IO
			zzFile= new ZipFile(zipFile);
			oBed= new FileOutputStream(tmpFile);
			
			// MARCELS PROBLEM
			// A fatal error has been detected by the Java Runtime Environment:
			// Internal Error (safepoint.cpp:237), pid=19287, tid=1098824016
			// Error: guarantee(PageArmed == 0,"invariant")
			// JRE version: 6.0_18-b07
			// Java VM: Java HotSpot(TM) 64-Bit Server VM (16.0-b13 mixed mode linux-amd64 )
			// If you would like to submit a bug report, please visit:
			// VM_Operation (0x00000000410e20a0): RevokeBias, mode: safepoint, requested by thread 0x00000000401d8000
			// Java Threads: ( => current thread )
			// 0x00000000401d8000 JavaThread "Sequencing Processor 1" [_thread_in_vm, id=19810,
			// also happened in rw thread 
			// workaround with -XX:-UseBias
			rw= IOHandlerFactory.getDefaultHandler();//new SyncIOHandler2(2);
			rw.addStream(oBed);
			if (tmpFasta!= null) {
				oFasta= new FileOutputStream(tmpFasta);
				rw.addStream(oFasta);
			}
			// GFFReader and Zipfile unfortunately not
//			if (FluxSimulatorSettings.optDisk)
//				rw.start();
			Processor[] processors= getProcessorPool(Math.min(settings.getMaxThreads(), 1));
			cntPlus= cntMinus= 0;
			for (reader.read(); !stop&& (g= reader.getGenes())!= null; reader.read()) {
				
				for (int i = 0; (!isStop())&& i < g.length; i++) {
					
					// process
					if (processors.length== 1) {
						processors[0].gene= g[i];
						processors[0].run();
					} else {
						boolean searching= true;
						while (searching&& !isStop())
							for (int j = 0; j < processors.length; j++) {
								synchronized (processors[j].lock) {
									if (processors[j].gene== null) {
										processors[j].gene= g[i];
										processors[j].lock.notify();
										searching= false;
										break;
									}
	//								while (processors[j].gene!= null)
	//									processors[j].lock.wait();
	//								processors[j].gene= g[i];
	//								processors[j].lock.notifyAll();
	//								searching= false;
	//								break;
								}
							}
					}
				}				
			}
			for (int i = 0; i < processors.length; i++) {
				processors[i].close();
				processors[i]= null;
			}
			this.processorPool= null;
			ctrProcessors= 0;
			zzFile.close();
			
			// stats
			cntLoci= hashLoc.size();
			cntTrpts= hashTrp.size();
			hashLoc= null;
			hashTrp= null;
			System.gc();

            Log.progressFinish();
            Log.message("");
            Log.message("\t"+ totalFrags+ " fragments found");
            Log.message("\t"+ totalReads+ " reads sequenced");
            Log.message("\t"+ cntPolyA+ " reads fall in poly-A tail");
            Log.message("\t"+ cntTruncReads+ " truncated reads");

			
			// I/O end
			rw.close();
//			oBed.flush();
//			oBed.close();
//			if (oFasta!= null) {
//				oFasta.flush();
//				oFasta.close();
//			}
			
			if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) 				
				System.err.println("\n\tMoving temporary BED file");
			FileHelper.move(tmpFile, settings.getSeqFile());
            Log.progressFinish();

			if (tmpFasta!= null) {
				if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
					System.err.println("\n\tCopying qFasta file");
				File fileFASTA= getFASTAfile();				
				FileHelper.move(tmpFasta, fileFASTA);
                Log.progressFinish();
			}
			
			if (!stop) 
				FluxSimulatorSettings.appendProfile(settings, FluxSimulatorSettings.PRO_COL_NR_SEQ, map);
		
			zipFile.delete();
			zipFile= null;
			return true;
			
		} catch (Exception e) {
			//if (Constants.verboseLevel>= Constants.VERBOSE_ERRORS)
				e.printStackTrace();
			if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
				System.err.println(" FAILED");
			zipFile.delete();
			zipFile= null;
			return false;
		}
	}
	
	private boolean hasQualities() {
		if (babes== null|| !babes.hasQualities())
			return false;
		return true;
	}
	
	private File getFASTAfile() {
		return FileHelper.replaceSfx(settings.getSeqFile(), "."+ (hasQualities()? SFX_FASTQ: SFX_FASTA));
	}

	/**
	 * @deprecated
	 * @param rndQual
	 * @param seq
	 * @param sequals
	 * @return
	 */
	private String createQual(Random rndQual, char[] seq, byte[] sequals) {
		 for (int l = 0; sequals!= null&& l < seq.length; l++) {
			 byte dbg= (byte) (sequals[l]+ 64);
			 seq[l]= (char) (sequals[l]+ ((ModelPool.qualLevels[1]== ModelPool.QLEVEL_ILLUMINA[1])?64:33));
//			 if (sequals[l]< 0)
//				 System.currentTimeMillis();
		 }
			 //(33+ DAVID_CORRECT[rndQual.nextInt(DAVID_CORRECT.length)]);
		return new String(seq);
	}
	
	private void createQual(ByteArrayCharSequence sequals) {
		byte[] a= sequals.a;
		for (int l = sequals.start; sequals!= null&& l < sequals.end; l++) {
			a[l]+= (ModelPool.qualLevels[1]== ModelPool.QLEVEL_ILLUMINA[1])? 64: 33;
		}
	}

	/**
	 * @deprecated
	 * @param obj
	 * @param absDir
	 * @param tDir
	 * @param seq
	 * @param sequals
	 * @return
	 */
	private String createQSeq(BEDobject obj, byte absDir, byte tDir, char[] seq, byte[] sequals) {
		 
		String s= obj.getChrom().equals(getPolyAobj().getChrom())?"":obj.readSequence();
//		 if (absDir< 0) {
//			 s= Graph.reverseSequence(s);
//			 s= Graph.complementarySequence(s);
//		 } 
		 s= s.toUpperCase();
		 try {
			 System.arraycopy(s.toCharArray(),0,seq,0,s.length());
		 } catch (ArrayIndexOutOfBoundsException e) {
			 e.printStackTrace();
		 }
        // c was always false
//		 for (int i = 0; c&& i < s.length(); i++)
//			 seq[i]= (c&& rnd.nextBoolean())? seq[i/2]: seq[i];
		 if (s.length()< seq.length) {
			 if (absDir== tDir)
				 Arrays.fill(seq,s.length(),seq.length,'a');
			 else {
				 System.arraycopy(seq, 0, seq, (seq.length- s.length()), s.length());
				 Arrays.fill(seq,0,(seq.length- s.length()),'t');
			 }
		 }
		 
		 if (babes!= null) 
			 babes.apply(seq, sequals);
		 return new String(seq);
	}

	private static final char DOT= '.', PLUS= '+', MINUS= '-';
	/**
	 * @deprecated
	 * @param obj
	 * @param t
	 * @param t2
	 * @param molNr
	 * @param absDir
	 * @param fragStart
	 * @param fragEnd
	 * @param readStart
	 * @param readEnd
	 * @return
	 */
	private String createQname(BEDobject obj, Transcript t, Transcript t2, 
			long molNr, byte absDir, int fragStart, int fragEnd, int readStart, int readEnd) {
		StringBuilder sb= new StringBuilder();
		if (babes== null|| !babes.hasQualities())
			sb.append(">");
		else
			sb.append("@");
		if (obj== null) { 
			sb.append(createReadName(t, t2, molNr, absDir, fragStart, fragEnd, readStart, readEnd));
			sb.append(FMRD.DELIM_FMRD[0]);
			sb.append(obj.getChrom());
			sb.append(FluxSimulatorSettings.DELIM_FMOLI);
			sb.append(DOT);
			sb.append(FluxSimulatorSettings.DELIM_FMOLI);
			sb.append(DOT);
			sb.append(FluxSimulatorSettings.DELIM_FMOLI);
			sb.append(absDir> 0?PLUS: MINUS); 
		} else { 
			sb.append(obj.getName());
			sb.append(FMRD.DELIM_FMRD[0]);
			sb.append(obj.getChrom());
			sb.append(FluxSimulatorSettings.DELIM_FMOLI);
			sb.append(Integer.toString(obj.getAbsoluteStart()));
			sb.append(FluxSimulatorSettings.DELIM_FMOLI);
			sb.append(Integer.toString(obj.getAbsoluteEnd()));
			sb.append(FluxSimulatorSettings.DELIM_FMOLI);
			sb.append(absDir> 0?PLUS: MINUS); 
			if (obj.getBlockCount()> 1) {
				sb.append(FluxSimulatorSettings.DELIM_FMOLI);
				sb.append(obj.getBlockSizes().toString());
				sb.append(FluxSimulatorSettings.DELIM_FMOLI);
				sb.append(obj.getBlockStarts().toString());
			}
		}
		
		return sb.toString();
	}

	public static final byte BYTE_TAB= '\t', BYTE_DOT= '.', BYTE_COMMA= ',', BYTE_0= '0', BYTE_1= '1', BYTE_PLUS= 43, BYTE_MINUS= '-', BYTE_GT= 62, BYTE_AT= 64, BYTE_NL= '\n';
	private void createQname(BEDobject2 obj2, ByteArrayCharSequence cs, Transcript t, Transcript t2, 
			long molNr, byte absDir, int fragStart, int fragEnd, int readStart, int readEnd) {
		
		byte[] a= obj2.a;
		cs.reset();
        int p1= obj2.getNameP1(), p2= obj2.getNameP2();
        cs.ensureLength(0, 1+ (p2-p1));
		byte[] b= cs.a;
		b[0]= (babes== null|| !babes.hasQualities())? BYTE_GT: BYTE_AT;
		++cs.end;
		assert(p1> 0&& p2> 0);
		System.arraycopy(a, p1, b, 1, p2- p1);
		cs.end+= p2- p1;
		cs.append(BYTE_NL);		
		
		
/*		StringBuilder sb= new StringBuilder();
		if (babes== null|| !babes.hasQualities())
			sb.append(">");
		else
			sb.append("@");
		if (obj== null) { 
			sb.append(createReadName(t, t2, molNr, absDir, fragStart, fragEnd, readStart, readEnd));
			sb.append(FMRD.DELIM_FMRD[0]);
			sb.append(obj.getChrom());
			sb.append(FluxSimulatorSettings.DELIM_FMOLI);
			sb.append(DOT);
			sb.append(FluxSimulatorSettings.DELIM_FMOLI);
			sb.append(DOT);
			sb.append(FluxSimulatorSettings.DELIM_FMOLI);
			sb.append(absDir> 0?PLUS: MINUS); 
		} else { 
			sb.append(obj.getName());
			sb.append(FMRD.DELIM_FMRD[0]);
			sb.append(obj.getChrom());
			sb.append(FluxSimulatorSettings.DELIM_FMOLI);
			sb.append(Integer.toString(obj.getAbsoluteStart()));
			sb.append(FluxSimulatorSettings.DELIM_FMOLI);
			sb.append(Integer.toString(obj.getAbsoluteEnd()));
			sb.append(FluxSimulatorSettings.DELIM_FMOLI);
			sb.append(absDir> 0?PLUS: MINUS); 
			if (obj.getBlockCount()> 1) {
				sb.append(FluxSimulatorSettings.DELIM_FMOLI);
				sb.append(obj.getBlockSizes().toString());
				sb.append(FluxSimulatorSettings.DELIM_FMOLI);
				sb.append(obj.getBlockStarts().toString());
			}
		}
		
		return sb.toString();
*/		
	}

	private static ByteArrayCharSequence writeOutZip(ByteArrayCharSequence cs, ZipOutputStream zout, ByteArrayCharSequence lastID) throws IOException {
		cs.resetFind();
		ByteArrayCharSequence id= cs.getToken(2);
		if (!id.equals(lastID)) {
			zout.putNextEntry(new ZipEntry(id.toString()));
			return id.cloneCurrentSeq();
		}
		return lastID;
	}
	
	public static final String BED_NAME_SEPARATOR= "_", BED_NAME_FORWARD= "1", BED_NAME_REVERSE= "2";
	private static final int[] int1array= new int[]{0};
	/**
	 * @deprecated
	 * @param t
	 * @param t2
	 * @param molNr
	 * @param absDir
	 * @param fragStart
	 * @param fragEnd
	 * @param readStart
	 * @param readEnd
	 * @return
	 */
	private String createReadName(Transcript t, Transcript t2, 
			long molNr, byte absDir, int fragStart, int fragEnd, int readStart, int readEnd) {
		
		// FURI
		String s= t.getGene().getGeneID()+ FluxSimulatorSettings.DELIM_FMOLI
			+ (t2== null? t.getTranscriptID(): t2.getTranscriptID())+ FluxSimulatorSettings.DELIM_FMOLI
			+ Long.toString(molNr+1)+ FluxSimulatorSettings.DELIM_FMOLI+ 
			+ (t2== null? t.getExonicLength(): t2.getExonicLength())+ FluxSimulatorSettings.DELIM_FMOLI
			+ Integer.toString(fragStart)+ FluxSimulatorSettings.DELIM_FMOLI+ Integer.toString(fragEnd)+ FluxSimulatorSettings.DELIM_FMOLI
			+ Integer.toString(readStart)+ FluxSimulatorSettings.DELIM_FMOLI+ Integer.toString(readEnd);
		
		if (settings.isPairedEnd()) { 
			s+= Character.toString(FMRD.DELIM_FMRD[0])+ FMRD.ID_PE+ (absDir== t.getStrand()?FMRD.PE_OPT[0]:FMRD.PE_OPT[1]);
		}
		
		return s;
	} 
	
	private BEDobject polyAobj;
	private BEDobject getPolyAobj() {
		if (polyAobj == null) {
			polyAobj = new BEDobject("polyA", (byte) 1);
		}

		return polyAobj;
	}
	
	private static final ByteArrayCharSequence CHR_POLYA= new ByteArrayCharSequence("polyA");
	
	private BEDobject2 createReadPolyA(BEDobject2 obj, int start, int end, Transcript t, Transcript t2, 
			long molNr, byte absDir, int fragStart, int fragEnd) {
		
		obj.reset();
		
		// intransparent and slower
/*		obj.setChromosome(CHR_POLYA);
		obj.setStart(0);
		obj.setEnd(end- start);
		FMRD.addReadName(obj, t, t2, molNr, absDir, fragStart, fragEnd, start, end, settings.isPairedEnd());
		obj.setStrand(absDir);
*/
		obj.append(CHR_POLYA);
		obj.append(BYTE_TAB);
		obj.append(0);
		obj.append(BYTE_TAB);
		obj.append(end- start);
		obj.append(BYTE_TAB);
		FMRD.appendReadName(obj, t, t2, molNr, absDir, fragStart, fragEnd, start, end, settings.isPairedEnd());
		obj.append(BYTE_TAB);
		obj.append(BYTE_0);
		obj.append(BYTE_TAB);
		obj.append((absDir>= 0? BYTE_PLUS: BYTE_MINUS));
		obj.append(BYTE_TAB);
		obj.append(BYTE_DOT);
		obj.append(BYTE_TAB);
		obj.append(BYTE_DOT);
		obj.append(BYTE_TAB);
		obj.append(BYTE_0);
		obj.append(BYTE_COMMA);
		obj.append(BYTE_0);
		obj.append(BYTE_COMMA);
		obj.append(BYTE_0);
		obj.append(BYTE_TAB);
		obj.append(BYTE_1);
		obj.append(BYTE_TAB);
		obj.append(end- start+ 1);
		obj.append(BYTE_TAB);
		obj.append(BYTE_0);
		
		return obj;
	}
	
	public final static byte[] BYTE_ARRAY_FROM_STRAND_TO_BLOCKS= new byte[] {'\t', '.', '\t', '.', '\t', '0', ',', '0', ',', '0'};
	                         
	private BEDobject2 createRead(BEDobject2 obj, int start, int end, Transcript t, Transcript t2, 
			long molNr, byte absDir, int fragStart, int fragEnd, boolean left) {
		
		int tlen= t.getExonicLength();
		if (start> tlen) {
			++cntPolyA;
			return createReadPolyA(obj, start, end, t, t2, molNr, absDir, fragStart, fragEnd);	// read in polyA tail
		}
		int offsStart= 0, offsEnd= 0;
		if (start< 1) {
			offsStart= start- 1; // t-coords 1-based, offs negative 
			start= 1;
		}
		if (end> tlen) {
			offsEnd= end- tlen;	// positive, pos's after annotated end
			end= tlen;
			if (!left)
				++cntPolyA;
		} else if (end< 1) {
			offsEnd= end- 1;	// negative, pos's before annotated start
			end= 1;
		}
		
		
		// bed boundaries
		byte strand= t.getStrand();
		int bedStart= Math.abs(t.getGenomicPosition(start-1)),
			bedEnd= Math.abs(t.getGenomicPosition(end-1));	// uncorrected positions, 1-based
		int idxExA= t.getExonIdx(strand* bedStart), 
			idxExB= t.getExonIdx(strand* bedEnd);	// exon indices
		if (idxExA== -1|| idxExB== -1) {
			System.err.println("[INCONSISTENT] strand "+ strand+ " idx "+ idxExA+", "+idxExB);
		}
		bedEnd= offsEnd>= 0? bedEnd+ (offsEnd* strand)
				: bedStart+ (offsEnd* strand);	// use original bedstart, before!
		bedStart+= offsStart* strand;			// correct out of range
		if (bedStart> bedEnd) {
			if (t.getStrand()>= 0)
				System.err.println("[INCONSISTENT] start "+ bedStart+", end "+bedEnd+", strand "+t.getStrand());
			int h= bedStart; bedStart= bedEnd; bedEnd= h;	// swap for neg strand
		}
		--bedStart; // lower the lower pos, 0-based

		// build object
		obj.reset();
		obj.append(t.getChromosome());
		obj.append(BYTE_TAB);
		obj.append(bedStart);
		obj.append(BYTE_TAB);
		obj.append(bedEnd);
		obj.append(BYTE_TAB);
		FMRD.appendReadName(obj, t, t2, 
				molNr, absDir, fragStart, fragEnd, 
				start, end, settings.isPairedEnd());
		obj.append(BYTE_TAB);
		obj.append(BYTE_0);
		obj.append(BYTE_TAB);
		obj.append(absDir>= 0? BYTE_PLUS: BYTE_MINUS);
		obj.append(BYTE_ARRAY_FROM_STRAND_TO_BLOCKS, 0,
				BYTE_ARRAY_FROM_STRAND_TO_BLOCKS.length);	// spare if there are blocks?
		if (idxExA== idxExB) {	// no blocks, spare?
			// spare?
			obj.append(BYTE_TAB);
			obj.append(BYTE_1);	// 1 block
			obj.append(BYTE_TAB);
			obj.append(bedEnd- bedStart);
			obj.append(BYTE_TAB);
			obj.append(BYTE_0);
			return obj;
		}
		
//		if (offsStart== -1)
//			System.currentTimeMillis();
		// blocks
		if (strand< 0) {
			int h= idxExA; idxExA= idxExB; idxExB= h;	// swap for iterating oder left->right
//			int h= bedStart; bedStart= bedEnd; bedEnd= h;
//			h= offsStart; offsStart= -offsEnd; offsEnd= -offsStart;
		}
		
		int nr= Math.abs(idxExB- idxExA)+ 1;
		obj.append(BYTE_TAB);
		obj.append(nr);		
		Exon[] ee= t.getExons();	// idx in trpt dir
		for (int i = idxExA ;; i+= (strand>= 0? 1: -1)) {			
			int gp= ((strand>= 0&& i== 0&& offsStart< 0)|| (strand< 0&& i== (ee.length- 1)&& offsEnd> 0))? bedStart+ 1
					: Math.max(Math.abs(ee[i].getStart()), bedStart+ 1);					
			int x= gp- (bedStart+ 1);
			obj.setNextBlockStart(x);
			x=  ((strand>= 0&& i== (ee.length- 1)&& offsEnd> 0)|| (strand<0&& i== 0&& offsStart< 0)? bedEnd
					: Math.min(Math.abs(ee[i].getEnd()), bedEnd))- gp+ 1;
			obj.setNextBlockSize(x);
			if (i== idxExB)
				break;
		}
		
		// DEBUG check
/*		 String[] ss= obj.getToken(BEDobject2.FN_BLOCK_SIZES).toString().split(",");
		 int sum= 0;
		 for (int i = 0; i < ss.length; i++) 
			sum+= Integer.parseInt(ss[i]);
		 if (sum!= settings.getReadLength()) {
			 System.currentTimeMillis();
		 }
*/		
		return obj;
	}
	
	public boolean setStop(boolean newStop) {
		if (newStop) 
			return setStop();
		this.stop= newStop;
		return true;
	}
	
	public boolean setStop() {
		this.stop= true;
		return true;
	}
	
	public boolean isStop() {
		return this.stop;
	}

    /**
     * Returns an error message if something is broken or missing and null if everything is fine
     *
     * @return message error message or null
     */
	public String isReady() {
        if(settings == null) return "No Setting specified!";
        File file = settings.getFrgFile();
        if(file == null || !file.exists()){
            if(file == null)return "No Fragmentation file specified! Check your parameters file!";
            if(!file.exists() || settings.getFrgFile().length() == 0) return "Fragmentation file " + file.getAbsolutePath() + " not found or empty. Make sure fragmentation was done !";
        }
		return null;
	}
	
	public boolean isFinished() {
		if (settings!= null&& settings.getSeqFile()!= null&& settings.getSeqFile().exists()&& settings.getSeqFile().length()> 0)
			return true;
		return false;
	}

	public int getCntLoci() {
		return cntLoci;
	}

	public void setCntLoci(int cntLoci) {
		this.cntLoci = cntLoci;
	}

	public int getCntTrpts() {
		return cntTrpts;
	}

	public void setCntTrpts(int cntTrpts) {
		this.cntTrpts = cntTrpts;
	}

	public int getCntExons() {
		return cntExons;
	}

	public void setCntExons(int cntExons) {
		this.cntExons = cntExons;
	}

	public long getTotalReads() {
		return totalReads;
	}

	public ModelPool getBabes() {
		return babes;
	}

	public boolean isBabeCompatible() {
		if (babes== null) 
			return true;
		
		if (babes.getReadLength()!= settings.getReadLength())
			return false;
		return true;
	}
	
	public void setBabes(ModelPool babes) {
		this.babes = babes;
	}

	private static final byte BYTE_a= 97, BYTE_t= 116;
	private int cntPolyA= 0, cntTruncReads= 0;
	private ByteArrayCharSequence createQSeq(ByteArrayCharSequence cs, BEDobject2 obj, byte absDir, byte tDir, int len, int flen/*int fstart, int fend, Transcript t*/) {

			//int flen= fend- fstart+ 1;
			cs.ensureLength(cs.end, len);
			int x;
			for (x= 0; x< CHR_POLYA.length(); ++x) 
				if (obj.charAt(x)!= CHR_POLYA.charAt(x))
					break;
			int fStart= cs.start, seqStart= cs.end;	// start of fasta obj
			if (x!= CHR_POLYA.length()) {		// not in pA, at least partially
				obj.readSequence(cs);	// not completely in poly-A
			 	cs.toUpperCase(seqStart, cs.end);			 	
			}
			
			// Issue 36 (20100516): sequencing poly-A in fragments < readLen
			// change from len (readLength) to Math.min(flen, len)
			int diff= Math.min(flen, len)- (cs.end- seqStart);	// fill trailing As
			if (diff> 0) {
				
				
				// prevent polyA in the middle of a transcript
				if (absDir== tDir)
					Arrays.fill(cs.a, cs.end, cs.end+ diff, BYTE_a);
				else {
					System.arraycopy(cs.a, fStart, cs.a, fStart+ diff, diff);
					Arrays.fill(cs.a, fStart, fStart+ diff, BYTE_t);
				}
				cs.end+= diff;
				
				//++cntPolyA; // count only reads that fully fall into pA, see createRead()
			}
			
			// create qual seq
			if (babes!= null) {
				//int mark= cs.end;
				babes.apply(cs, seqStart);
			}
			
			return cs;
	}

	Random rndQual= new Random();
	boolean sequence_old() { 
			try {
				String sss= "sequencing";
                Log.progressStart(sss);

				GFFReader reader= new GFFReader(settings.getRefFile().getCanonicalPath());
				reader.setReadAheadLimit(500);
				reader.setSilent(true);
				reader.setStars(true);
				Gene[] g;
				ThreadedFragFileReader fragReader= new ThreadedFragFileReader(new FileReader(settings.getFrgFile().getAbsolutePath()));
				fragReader.read();	// thread killed, currently intransparent
				BEDobject obj;
				totalReads= 0;
				long totalMols= 0, totalFrags= 0;
				Iterator<int[][]> iter= mapFrags.values().iterator();
				while(iter.hasNext())
					totalMols+= iter.next().length;
				for (int i = 0; i < settings.getProfiler().getMolecules().length; i++) 
					totalMols+= settings.getProfiler().getMolecules()[i];
				HashSet<String> hashTrp= new HashSet<String>(), hashLoc= new HashSet<String>();
				ThreadedQWriter qwriter= null, qfasta= null;
				BufferedWriter bwriter= null, bfasta= null;
				File tmpFile= File.createTempFile("flux",NAME_SEQ,settings.getTmpDir()), tmpFasta= null;
				if (!multiThread)
					bwriter= new BufferedWriter(new FileWriter(tmpFile));
				if (settings.isFastQ()&& settings.getGenDir()!= null) {
					tmpFasta= File.createTempFile("flux",NAME_SEQ,settings.getTmpDir());
					Graph.overrideSequenceDirPath= settings.getGenDir().getAbsolutePath();
					if (!multiThread)
						bfasta= new BufferedWriter(new FileWriter(tmpFasta));
				}
				char[] seq= new char[settings.getReadLength()];
				final char[] pA= new char[settings.getReadLength()], pT= new char[settings.getReadLength()];
				for (int i = 0; i < pT.length; i++) {
					pA[i]= 'a';
					pT[i]= 't';
				}
				//final String polyA= new String(pA), polyT= new String(pT);
				byte[] sequals= null;
				if (babes!= null) 
					sequals= new byte[settings.getReadLength()];		
	
				// probability for fragment selection
				p= 1;
				if (settings.getReadNr()< nrOfFrags)
					p= settings.getReadNr()/ (double) nrOfFrags;	
				//System.err.println(settings.readNr+" reads, "+nrOfFrags+ " frags, p= "+p);
				Hashtable<CharSequence,Long> map= new Hashtable<CharSequence,Long>(settings.getProfiler().getMolecules().length);
				Long long0= new Long(1);

				
				for (reader.read(); !stop&& (g= reader.getGenes())!= null; reader.read()) {
					
					if (multiThread) {
						if (qwriter!= null) {
							qwriter.flush();
							qwriter.close();
							try {
								qwriter.join();
							} catch (InterruptedException shit) {
								qwriter.setStop();
							}
						}
						if (qfasta!= null) {
							qfasta.flush();
							qfasta.close();
							try {
								qfasta.join();
							} catch (InterruptedException shit) {
								qfasta.setStop();
							}
						}
						qwriter= new ThreadedQWriter(tmpFile);
						qwriter.setAppend(true);
						qwriter.init();
						qwriter.start();
						if (tmpFasta!= null) {
							qfasta= new ThreadedQWriter(tmpFasta);
							qfasta.setAppend(true);
							qfasta.init();
							qfasta.start();
						}
					}
					
					for (int i = 0; (!isStop())&& i < g.length; i++) {
						for (int j = 0; (!isStop())&& j < g[i].getTranscripts().length; j++) {
							Transcript t= g[i].getTranscripts()[j];
							Transcript t2= null;
							int[][] m= mapFrags.remove(t.getTranscriptID());
							if (m!= null) {
								 for (int k = 0; k < m.length; k++, ++totalFrags) {
									 double r= rnd.nextDouble();
									 if (r< p) { 
										 hashTrp.add(t.getTranscriptID());
										 hashLoc.add(t.getGene().getGeneID());
	//									 if (totalReads% 1000== 0&& (totalReads* 10d/ settings.readNr> perc)) {
	//										 ++perc;
	//										 if (Constants.progress!= null)
	//											 Constants.progress.setValue(perc);
	//									 }
										 
										 double q= rndFiftyFifty.nextDouble();
										 int fragLen= (m[k][1]-m[k][0]+1);
										 if (settings.isPairedEnd()|| q< 0.5) {
											 byte absDir= (byte) (t.getStrand()> 0?1: -1);
											 ++totalReads;
											 if (map.containsKey(t.getTranscriptID()))
												 map.put(t.getTranscriptID(), map.get(t.getTranscriptID())+ 1);
											 else 
												 map.put(t.getTranscriptID(), long0);
											 obj= createRead(m[k][0],Math.min(m[k][0]+ settings.getReadLength()- 1, m[k][1]), t, t2, k, absDir, 
													 m[k][0], m[k][1]);
											 if (multiThread)
												 qwriter.add(obj);
											 else
												 bwriter.write(obj.toString()+"\n");
											 if (settings.isFastQ()) {
												 String qname= createQname(obj,t,t2,k,absDir,m[k][0],m[k][1],
														 m[k][0],Math.min(m[k][0]+ settings.getReadLength()- 1, m[k][1])), 
												 	qseq= createQSeq(obj, absDir, t.getStrand(), seq, sequals);
												 if (multiThread) {
													 qfasta.add(qname);
													 qfasta.add(qseq);
												 } else {
													 bfasta.write(qname+ "\n");
													 bfasta.write(qseq+ "\n");
												 }
												 if (sequals!= null&& babes.hasQualities()) {
													 if (multiThread) {
														 qfasta.add(Fasta.QFASTA_PLUS);
														 qfasta.add(createQual(rndQual, seq, sequals));
													 } else {
														 bfasta.write(Fasta.QFASTA_PLUS+ "\n");
														 bfasta.write(createQual(rndQual, seq, sequals)+ "\n");
													 }
												 }
											 } 
										 }
									 	 if (settings.isPairedEnd()|| q> 0.5) {
											 byte absDir= (byte) (t.getStrand()> 0?-1:1);
											 ++totalReads;
											 if (map.containsKey(t.getTranscriptID()))
												 map.put(t.getTranscriptID(), map.get(t.getTranscriptID())+ 1);
											 else 
												 map.put(t.getTranscriptID(), long0);
											 
											 obj= createRead(Math.max(m[k][1]- settings.getReadLength()+ 1, m[k][0]),m[k][1], t, t2, k, absDir,
													 m[k][0],m[k][1]);
											 if (multiThread)
												 qwriter.add(obj);
											 else
												 bwriter.write(obj+ "\n");
											 if (settings.isFastQ()) {
												 String qname= createQname(obj,t,t2,k,absDir,m[k][0],m[k][1], Math.max(m[k][1]- settings.getReadLength()+ 1, m[k][0]),m[k][1]), 
												 	qseq= createQSeq(obj,absDir, t.getStrand(), seq, sequals);
												 if (multiThread) {
													 qfasta.add(qname);
													 qfasta.add(qseq);
												 } else {
													 bfasta.write(qname+ "\n");
													 bfasta.write(qseq+ "\n");
												 }
												 if (sequals!= null&& babes.hasQualities()) {
													 if (multiThread) {
														 qfasta.add(Fasta.QFASTA_PLUS);
														 qfasta.add(createQual(rndQual, seq, sequals));
													 } else {
														 bfasta.write(Fasta.QFASTA_PLUS+ "\n");
														 bfasta.write(createQual(rndQual, seq, sequals)+ "\n");
													 }
												 }
											 }
											 
									 	 }
									 } 
								 }
								 
								 if (mapFrags.size()== 0) {
									 fragReader.read();
									 mapFrags= fragReader.getMapFrags();
									 i= -1;
									 break;
								 }
							}
						}
						if (mapFrags== null)
							break;
					}
					
				}
				fragReader.close();
				cntLoci= hashLoc.size();
				cntTrpts= hashTrp.size();
				hashLoc= null;
				hashTrp= null;
				System.gc();
				if (multiThread) {
					if (qwriter!= null) {
						qwriter.flush();
						qwriter.close();
						try {
							qwriter.join();
						} catch (Exception e) {
							qwriter.setStop();
						}
					}
					if (qfasta!= null) {
						qfasta.flush();
						qfasta.close();
						try {
							qfasta.join();
						} catch (InterruptedException shit) {
							qfasta.setStop();
						}
					}
				} else {
					bwriter.flush();
					bwriter.close();
					if (bfasta!= null) {
						bfasta.flush();
						bfasta.close();
					}
				}
				

                Log.progressStart("Moving temporary BED file");
				FileHelper.move(tmpFile, settings.getSeqFile());
                Log.progressFinish();

				if (tmpFasta!= null) {
                    Log.progressStart("Copying qFasta file");
					File fileFASTA= getFASTAfile();
					FileHelper.move(tmpFasta, fileFASTA);
                    Log.progressFinish();
				}
				
				if (!stop) {
					FluxSimulatorSettings.appendProfile(settings, FluxSimulatorSettings.PRO_COL_NR_SEQ, map);
				}
				
				return true;
				
			} catch (Exception e) {
				//if (Constants.verboseLevel>= Constants.VERBOSE_ERRORS)
					e.printStackTrace();
				if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
					System.err.println(" FAILED");
				return false;
			}
		}

	/**
	 * @deprecated
	 * @param start
	 * @param end
	 * @param t
	 * @param t2
	 * @param molNr
	 * @param absDir
	 * @param fragStart
	 * @param fragEnd
	 * @return
	 */
	private BEDobject createRead(int start, int end, Transcript t, Transcript t2, 
				long molNr, byte absDir, int fragStart, int fragEnd) {
			
			int tlen= t.getExonicLength();
	//		if (!dir)
	//			System.currentTimeMillis();
			
			String name= createReadName(t, t2, molNr, absDir, fragStart, fragEnd, start, end);
	
			if (start> tlen)
				return createReadPolyA(start, end, t, t2, molNr, absDir, fragStart, fragEnd);	// read in polyA tail		
			int offsStart= 0, offsEnd= 0;
			if (start< 1) {
				offsStart= start- 1; // neg offset
				start= 1;
			}
			if (end< 1) {
				offsEnd= end -1;	// also neg, can maybe occur
				end= 1;
			}
			if (end> tlen) {
				offsEnd= end- tlen;	// + 1, wrong length for single block units, killed
				end= tlen;
			}
			assert(offsStart<= 0);
			
			if (offsEnd< 0) {
				int bedStart= t.get5PrimeEdge()- 1+ offsStart, 
					bedEnd= t.get5PrimeEdge()+ offsEnd;
				BEDobject obj= new BEDobject(
						t.getChromosome(), 
						absDir,
						bedStart,
						bedEnd
				);
	//			obj.setBlockCount(1);
	//			obj.setBlockStarts(int1array);
	//			obj.setBlockSizes(new int[] {bedEnd- bedStart});
				obj.setName(name);
				return obj;
			}
			start= t.getGenomicPosition(start-1);
			end= t.getGenomicPosition(end-1);
			if (start> end) {
				int h= start; start= end; end= h;
			}
			// TODO efficiency
			DirectedRegion reg= new DirectedRegion(start, end, t.getStrand());
			reg.setChromosome(t.getChromosome());
			DirectedRegion[] regs= DirectedRegion.intersect(new DirectedRegion[] {reg}, t.getExons());
			BEDobject obj= new BEDobject(regs);
			obj.setStrand(absDir);
			obj.setName(name);
			
			if (offsStart< 0)
				obj.extend(true, Math.abs(offsStart));
	//		if (offsEnd< 0)
	//			obj.extend(false, Math.abs(offsEnd));
			
	//		if (offsStart< 0) 
	//			obj.trim(true, offsStart);
			
			// Do not include the poly-A region
	/*		if (offsEnd> 0) {
				if (obj.getBlockSizes().length== 1) {	// TODO debug, is wrong for that
					if (t.isForward()) {
						obj.setEnd(obj.getEnd()+ offsEnd);
						obj.getBlockSizes()[obj.getBlockSizes().length- 1]+= offsEnd;
					} else {
						obj.setStart(obj.getStart()- offsEnd);
						obj.getBlockSizes()[0]+= offsEnd;
					}
				}
			}
	*/		
			
	/*		int sum= 0;
			for (int i = 0; i < obj.getBlockSizes().length; i++) 
				sum+= obj.getBlockSizes()[i];
			if (sum!= settings.readLength)
				System.currentTimeMillis();
	*/			
	
			assert(obj.check());
			 
			return obj;
		}

	/**
	 * @deprecated
	 * @param start
	 * @param end
	 * @param t
	 * @param t2
	 * @param molNr
	 * @param absDir
	 * @param fragStart
	 * @param fragEnd
	 * @return
	 */
	private BEDobject createReadPolyA(int start, int end, Transcript t, Transcript t2, 
				long molNr, byte absDir, int fragStart, int fragEnd) {
			BEDobject obj= getPolyAobj();
			obj.setStrand(absDir);
			obj.setStart(0);
			obj.setEnd(end- start);
	//		obj.setBlockCount(1);
	//		obj.setBlockStarts(new int[] {0});
	//		obj.setBlockSizes(new int[] {end- start});
			obj.setName(createReadName(t, t2, molNr, absDir, fragStart, fragEnd, start, end));
			return obj;
		}

	int sortAndZip(ByteArrayCharSequence cs, File inFile, File outFile) {
		// Sort
		try {
			String sss= "initing";
            Log.progressStart(sss);

			// IO handler disk
			InputStream istream= new FileInputStream(inFile);
			OutputStream ostream= new FileOutputStream(outFile);
			rw= IOHandlerFactory.getDefaultHandler();//new SyncIOHandler2(2);
			rw.addStream(ostream, 10* 1024);
			
			// setup sorter / zipper pipes
			PipedOutputStream pipeOut= null;
			DesignatedHierarchicalFieldComparator compi= new DesignatedHierarchicalFieldComparator(3);
			ZipperThread zipper= null;
			UnixStreamSort2 sorter= null;
			if (settings.getMaxThreads()== 1) {
				sorter= new UnixStreamSort2(inFile, compi);
				sorter.run();
				
				istream= new FileInputStream(sorter.getOutFile());
				rw.addStream(istream, 10* 1024);
				
				
			} else {
				rw.addStream(istream, 10* 1024);
				pipeOut= new PipedOutputStream();
				PipedInputStream pipeIn = new PipedInputStream(pipeOut);
				sorter= new UnixStreamSort2(pipeIn, compi);
				sorter.setSilent(true);
				zipper= new ZipperThread(
						sorter.getOutInStream(), 
						ostream,
						100,
						inFile.length()
				);
				sorter.start();
				zipper.start();
			}
//			if (FluxSimulatorSettings.optDisk)
//				rw.start();

			// feed sorter
			OutputStream feedOut= pipeOut== null?new ZipOutputStream(ostream):new BufferedOutputStream(pipeOut);
			int cnt= 0, perc= 0;
			long totBytes= inFile.length(), currBytes= 0;
			ByteArrayCharSequence lastID= null;
			int linesRec= 0;
			while(rw.readLine(istream, cs)> 0&& !stop) {
				currBytes+= cs.length()+ 1;				
				++cnt;
				//System.err.println(cnt+ "\t"+ currBytes);
                Log.progress(currBytes, totBytes);
				//System.err.println("w: "+ cs.toString());
				boolean closeAll= false;
				try {  
					if (zipper== null) {
						lastID= writeOutZip(cs, (ZipOutputStream) feedOut, lastID);
						++linesRec;
					}
					feedOut.write(cs.a, cs.start, cs.length());
					feedOut.write(BYTE_NL);
				} catch (InterruptedIOException ex) {
					closeAll= true;
				} catch (IOException ex) {
					closeAll= true;
				} catch (Exception ex) {
					closeAll= true;
				}
				if (closeAll) {
					assert(stop);
					istream.close();
					ostream.close();
					break;
				}
			}

            Log.progressFinish();

			if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) {
				System.err.println("\t"+ cnt+ " lines submitted");
			}
			
			// wait for zipping
			sss= "zipping";
            Log.progressStart(sss);
			feedOut.close();
			if (sorter.isAlive())
				sorter.join();			
			if (zipper!= null) {
				zipper.join();
				linesRec= zipper.linesRec;
			}
			rw.close();
//			if (rw.isAlive())
//				rw.join();
			ostream.close();
			
			assert(cnt== zipper.linesRec);
			
			// end
			if (stop) 
				return -1;

            Log.progressFinish();
		    Log.message("\t"+ linesRec+ " lines zipped");

			return cnt;
			
		} catch (Exception e) {
			if (!(e instanceof InterruptedException)) {
				e.printStackTrace();
			}
			return -1;
		}
	}	
}
