package fbi.genome.sequencing.rnaseq.simulation;

import fbi.commons.ByteArrayCharSequence;
import fbi.commons.Log;
import fbi.commons.file.FileHelper;
import fbi.commons.file.ReverseFileReader;
import fbi.commons.thread.StoppableRunnable;
import fbi.genome.io.ThreadedBufferedByteArrayStream;
import fbi.genome.io.gff.GFFReader;
import fbi.genome.model.Gene;
import fbi.genome.model.commons.Distribution;
import fbi.genome.model.commons.IntVector;
import fbi.genome.model.constants.Constants;

import java.io.*;
import java.util.*;

public class Profiler implements StoppableRunnable {
	
	public static final byte STAT_NONE= 0, STAT_ANN= 4, STAT_RELFREQ= 5, STAT_MOL= 6;
	FluxSimulatorSettings settings;
	byte status= -1;
	ByteArrayCharSequence[] ids= null, locIDs;
	int[] len= null;
	long[] molecules= null;
	boolean[] cds= null;
	private GFFReader gffReader;
	boolean stop= false;
	int cntLoci= -1;
	float txLocAvg= -1, txLocMed= -1, txLoc1Q= -1, txLoc3Q= -1, txLocSTD= -1, lenMed= -1, lenAvg= -1, lenSTD= -1, len1Q= -1, len3Q= -1, lenMin= -1, lenMax= -1;
	HashSet<CharSequence> sfHi, sfMed, sfLo;
	Hashtable<ByteArrayCharSequence,int[]> mapLenExp;
	boolean annotationChecked= false;
	
	public Profiler(FluxSimulatorSettings settings) {
		this.settings= settings;
	}
	
	
	
	public void run() {
        Log.info("PROFILING", "I am assigning the expression profile");
		status= getStatus();
		if (status== STAT_NONE) {
			if (!readAnnotation())
				throw new RuntimeException("Error while reading annotations!");
			status= STAT_ANN;
			writeProfile();
		}
		
		if (!isReady()) 
			return;
		if (!profile())
			throw new RuntimeException("Error creating profile !");

		if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
			Log.info(null,"\tmolecules\t"+sumMol);
		
	}
	
	public boolean isFinishedReadAnnotation() {
		return (locIDs!= null&& ids!= null&& ids.length> 0&& cds!= null&& len!= null&& 
				ids.length== len.length&& ids.length== locIDs.length&& ids.length== cds.length); 
	}
	
	public boolean isFinishedExpression() {
		return (isFinishedReadAnnotation()&& molecules!= null&& ids.length== molecules.length); 
	}
	
	public boolean isReady() {
		
		// TODO make cells -> dec_par, nb_molecules
		if (settings== null|| Double.isNaN(settings.getExpDistrP1())|| Double.isNaN(settings.getExpDistrP2())|| Double.isNaN(settings.getDecDistrP1())|| settings.getNbMolecules()<= 0|| 
				settings.getProFile()== null) {	// || (!settings.getProFile().canWrite()) // fails on winOS
			if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) {
				if (settings== null) 
					System.err.println("\t[NOPARAMS] I have no parameters!");
				else {
					System.err.print("\t[TOOLITTLE] there are parameters missing: ");
					if (Double.isNaN(settings.getExpDistrP1()))
						System.err.print(settings.PAR_EXPR_K+" ");
					if (Double.isNaN(settings.getExpDistrP2()))
						System.err.print(settings.PAR_EXPR_X0+" ");
					if (Double.isNaN(settings.getDecDistrP1()))
						System.err.print(settings.PAR_DEC_P1+" ");
					if (settings.getNbMolecules()<= 0)
						System.err.print(settings.PAR_NB_MOLECULES+" ");
					if ((settings.getProFile()== null)|| (!FileHelper.canWrite(settings.getProFile())))
						System.err.print("valid workfile ");
					System.err.println("\n");
				}
			}
			return false;
		} else {
			if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) {
				if (settings.getNbCells()> 0)
					System.err.println("\t"+settings.PAR_NB_CELLS+"\t"+settings.getNbCells());
				System.err.println("\t"+settings.PAR_NB_MOLECULES+"\t"+settings.getNbMolecules());
				System.err.println("\t"+settings.PAR_EXPR_K+"\t"+settings.getExpDistrP1());
				System.err.println("\t"+settings.PAR_EXPR_X0+"\t"+settings.getExpDistrP2());
				System.err.println("\t"+settings.PAR_DEC_P1+"\t"+settings.getDecDistrP1());
				try {
					System.err.println("\t"+FluxSimulatorSettings.PAR_PRO_FNAME+"\t"+settings.getProFile().getCanonicalPath());
				} catch (IOException e) {
					;	// :)
				}
			}
		}
		return true;
	}
	
	byte getStatus() {
		
		if (settings.getProFile()== null|| !settings.getProFile().exists())
			return STAT_NONE;
		
		try {
			ReverseFileReader rreader= new ReverseFileReader(settings.getProFile().getCanonicalPath());
			String s= rreader.readLine();
			if (s== null)
				return STAT_NONE;
			
			String[] tokens= s.split("\\s");
			return (byte) tokens.length;
			
		} catch (Exception e) {
			//e.printStackTrace();
			return STAT_NONE;
		}
		
	}

	public boolean readAnnotation() {
		try {
			gffReader= null;
			GFFReader reader= getGFFReader();
			if (reader== null|| isStop())
				return false;
			
			//if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) 
			Constants.progress.setString("Reading reference annotation");
			reader.read();
			java.util.Vector<ByteArrayCharSequence> v= new java.util.Vector<ByteArrayCharSequence>(30000);
			Vector<ByteArrayCharSequence>vLoc= new java.util.Vector<ByteArrayCharSequence>(30000);
			Vector<Boolean> vBoo= new java.util.Vector<Boolean>(30000);
			IntVector lenV= new IntVector(30000); /*, txLoc= new IntVector(20000);*/
			HashMap<ByteArrayCharSequence,ByteArrayCharSequence> locMap= 
					new HashMap<ByteArrayCharSequence,ByteArrayCharSequence>();
			mapLenExp= new Hashtable<ByteArrayCharSequence,int[]>();
			if (Constants.progress!= null)
				Constants.progress.setValue(0);
			for (Gene[] g; (!stop)&& (g= reader.getGenes())!= null; reader.read()) {				
				for (int i = 0; (!stop)&& i < g.length; i++) {
					++cntLoci;
					int cntTrpt= 0;
					for (int j = 0; (!stop)&& j < g[i].getTranscripts().length; j++) {
						
						if (g[i].getTranscripts()[j].isCoding()&& (!settings.isLoadCoding()))
							continue;
						if ((!g[i].getTranscripts()[j].isCoding())&& (!settings.isLoadNoncoding()))
							continue;
						if (g[i].getTranscripts()[j].isCoding())
							vBoo.add(true);
						else
							vBoo.add(false);
						v.add(new ByteArrayCharSequence(g[i].getTranscripts()[j].getTranscriptID()));
						ByteArrayCharSequence locName= new ByteArrayCharSequence(g[i].getGeneID());
						if (locMap.containsKey(locName))
							locName= locMap.get(locName);
						vLoc.add(locName);
						int[] a= new int[2];
						a[0]= g[i].getTranscripts()[j].getExonicLength();
						lenV.add(a[0]);
					}
					//txLoc.add(cntTrpt);
				}
			}
			if (stop)
				return false;
			
			cntLoci= locMap.size();			
			ids= new ByteArrayCharSequence[v.size()];
			for (int i = 0; i < v.size(); i++) 
				ids[i]= v.elementAt(i);
			v= null;
			locIDs= new ByteArrayCharSequence[vLoc.size()];
			for (int i = 0; i < vLoc.size(); i++) 
				locIDs[i]= vLoc.elementAt(i);
			vLoc= null;
			cds= new boolean[vBoo.size()];
			for (int i = 0; i < vBoo.size(); i++) 
				cds[i]= vBoo.elementAt(i);
			vBoo= null;
			
			len= lenV.toIntArray();
			//int[] as= txLoc.toIntArray();
			lenV= null;
			//txLoc= null;
			System.gc();
			
			calcStats();
			
			if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
				Constants.progress.message("\tfound "+ids.length+" transcripts\n");
		
			System.gc();

			
			
		} catch (Exception e) {
			e.printStackTrace();
			return false;
		}
		return true;
	}
	
	public boolean writeProfile() {
		try {
//			if (settings.getProFile().exists())
//				return false;
			BufferedWriter writer= new BufferedWriter(new FileWriter(settings.getProFile()));
			for (int i = 0; i < ids.length; i++) {
				int x= i;
//				if (FluxSimulator.c)
//					x= ids.length- 1- i;
				writer.write(locIDs[x]+ settings.PRO_FILE_SEP+ ids[x]+ settings.PRO_FILE_SEP+ (cds[x]?settings.PRO_FILE_CDS:settings.PRO_FILE_NC)
						+settings.PRO_FILE_SEP+ Integer.toString(len[x])+ "\n");
			}
			writer.flush();
			writer.close();
		} catch (Exception e) {
			return false;
		}
		return true;
	}
	
	public int getLength(ByteArrayCharSequence id) {
		return mapLenExp.get(id)[0];
	}
	
	
	
	private void calcStats() {
		
		CharSequence lastLocID= null;
		int asCtr= 1;
		IntVector asV= new IntVector();
		cntLoci= 0;
		for (int i = 0; i < locIDs.length; i++) {
			if ((!locIDs[i].equals(lastLocID))) {
				++cntLoci;
				if (lastLocID!= null) {
					asV.add(asCtr);
					asCtr= 1;
				}
				lastLocID= locIDs[i];
			} else
				++asCtr;
		}
		
		int[] as= asV.toIntArray();
		Arrays.sort(as);
		Distribution dist= new Distribution(as);
		txLocMed= (float) dist.getMedian();
		txLocAvg= (float) dist.getMean();
		txLocSTD= (float) dist.getStandardDeviation();
		txLoc1Q= (float) dist.get1stQuart();
		txLoc3Q= (float) dist.get3rdQuart();
		as= null;
		dist= new Distribution(len.clone());
		System.gc();

		lenMin= (float) dist.getMin();
		lenMax= (float) dist.getMax();
		lenAvg= (float) dist.getMean();
		lenMed= (float) dist.getMedian();
		len1Q= (float) dist.get1stQuart();
		len3Q= (float) dist.get3rdQuart();
		lenSTD= (float) dist.getStandardDeviation();
	}

	private GFFReader getGFFReader() {
//		if (gffReader == null) {
			gffReader = new GFFReader(settings.getRefFile().getAbsolutePath());	
			if (gffReader== null)
				return null;
			try {
				if ((!annotationChecked)&& (!gffReader.isApplicable())) {
					File refFile= gffReader.createSortedFile();
					if (refFile== null)
						return null;
					settings.setRefFile(new File(settings.getProFile().getParent()+ File.separator+ refFile.getName()));
					if (!refFile.equals(settings.getRefFile())) {
						if (!FileHelper.move(refFile, settings.getRefFile(), Constants.progress))
							settings.setRefFile(refFile);
					}
					gffReader= new GFFReader(settings.getRefFile().getAbsolutePath());
				}
				annotationChecked= true;
				gffReader.setSilent(true);
				gffReader.setStars(true);
				
			} catch (Exception e) {
				return null;
			}
//		}
	
		return gffReader;
		
	}
	
	double sumMol= 0;
	public boolean profile() {
		
		if (Constants.progress!= null)			
			Constants.progress.setString("profiling");
		try {
			if (ids== null) {
				ids= new ByteArrayCharSequence[FileHelper.countLines(settings.getProFile().getCanonicalPath())];
				len= new int[ids.length];
			}
			double sumRF= 0; sumMol= 0;
			sfHi= null; sfMed= null; sfLo= null;
			molecules= new long[ids.length];
			double[] relFreq= new double[ids.length];
			
			if (status< STAT_RELFREQ) {	// generate ranks
				
				// generate random permutation of ranks
				Random r= new Random();
				for (int i = 0; i < molecules.length; i++) 
					molecules[i]= 1+ r.nextInt(molecules.length-1); //i+1;	// ranks
//				for (int k = molecules.length - 1; k > 0; k--) {
//				    int w = (int) Math.floor(r.nextDouble() * (k+1));
//				    long temp = molecules[w];
//				    molecules[w] = molecules[k];
//				    molecules[k] = temp;
//				}
				
				relFreq= new double[ids.length];
				for (int i = 0; (!stop)&& i < relFreq.length; i++) {
					double par= pareto(molecules[i], settings.getExpDistrP1(), settings.getExpDistrP2());
					double exp= exponential(molecules[i], settings.getDecDistrP1());
					sumRF+= (relFreq[i]= par* exp);
				}
				for (int i = 0; (!stop)&& i < relFreq.length; ++i) {	// normalize
					relFreq[i]/= sumRF;
					sumMol+= (molecules[i]= Math.round(relFreq[i]* settings.getNbMolecules()));
					if (molecules[i]> 0)
						mapLenExp.put(getCombinedID(i), new int[] {len[i], (int) molecules[i]});
				}
			}

		
			if (Constants.progress!= null) {
				Constants.progress.setString("OK");
				Constants.progress.finish();
			} else if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
				System.err.println(" OK");
			
			if (!isStop()) {
				Hashtable<CharSequence,Long> map= new Hashtable<CharSequence,Long>(settings.getProfiler().getMolecules().length);
				for (int i = 0; i < settings.getProfiler().getMolecules().length; i++)
					if (settings.getProfiler().getMolecules()[i]!= 0) {
						ByteArrayCharSequence locNtid= settings.getProfiler().getLocIDs()[i].cloneCurrentSeq();
						locNtid.append(Character.toString(FluxSimulatorSettings.SEP_LOC_TID));
						locNtid.append(settings.getProfiler().getIds()[i]);
						map.put(locNtid,
								settings.getProfiler().getMolecules()[i]);
					}
				if (!FluxSimulatorSettings.appendProfile(settings, FluxSimulatorSettings.PRO_COL_NR_MOL, map))
					return false;
			}
			
			return true;
			
		} catch (Exception e) {
			if (Constants.progress!= null) {
				Constants.progress.setString("FAILED");
				Constants.progress.finish();
			} else if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
				System.err.println(" FAILED");
			if (Constants.verboseLevel>= Constants.VERBOSE_ERRORS)
				e.printStackTrace();
			return false;
		}
	}

	public static double exponential(double rank, double par1) {
	//		double val= Math.exp(- (Math.pow(rank, 2)/ Math.pow(par1, 2))
	//				- (rank/par1));
			double val= Math.exp(- (Math.pow(rank/par1, 2))	// / 122000000
					- (rank/par1));	// 7000
	
			return val;
		}

	public static double pareto(double rank, double par1, double par2) {
	//		double val= par2/ Math.pow(rank, par1)
	//			* Math.exp(- (Math.pow(rank, 2)/ 122000000)
	//					- (rank/7000));
			double val= Math.pow(rank/par2, par1)/* 2731598d*/;	// par1= 0,6  par2= (41627d/ 2731598d)
			return val;
		}



	public boolean isStop() {
		return stop;
	}

	public boolean setStop(boolean stop) {
		if (stop)
			return setStop();

		this.stop = stop;
		if (gffReader!= null)
			return gffReader.setStop(stop);
		return true;		
	}

	public ByteArrayCharSequence[] getIds() {
		return ids;
	}

	public int[] getLen() {
		return len;
	}

	public long[] getMolecules() {
		return molecules;
	}

	public double getSumMol() {
		return sumMol;
	}



	public float getTxLocAvg() {
		return txLocAvg;
	}



	public void setTxLocAvg(float txLocAvg) {
		this.txLocAvg = txLocAvg;
	}



	public float getTxLocMed() {
		return txLocMed;
	}



	public void setTxLocMed(float txLocMed) {
		this.txLocMed = txLocMed;
	}



	public float getTxLoc1Q() {
		return txLoc1Q;
	}



	public void setTxLoc1Q(float txLoc1Q) {
		this.txLoc1Q = txLoc1Q;
	}



	public float getTxLoc3Q() {
		return txLoc3Q;
	}



	public void setTxLoc3Q(float txLoc3Q) {
		this.txLoc3Q = txLoc3Q;
	}



	public float getTxLocSTD() {
		return txLocSTD;
	}



	public void setTxLocSTD(float txLocSTD) {
		this.txLocSTD = txLocSTD;
	}



	public float getLenMed() {
		return lenMed;
	}



	public float getLenAvg() {
		return lenAvg;
	}



	public float getLenSTD() {
		return lenSTD;
	}



	public float getLen1Q() {
		return len1Q;
	}



	public float getLen3Q() {
		return len3Q;
	}



	public float getLenMin() {
		return lenMin;
	}



	public float getLenMax() {
		return lenMax;
	}



	public int getCntLoci() {
		return cntLoci;
	}


	boolean createXPRmaps() {
		if (molecules== null|| ids== null)
			return false;
		sfHi= new HashSet<CharSequence>();
		sfMed= new HashSet<CharSequence>();
		sfLo= new HashSet<CharSequence>();
		for (int i = 0; i < molecules.length; i++) {
			if (Math.ceil(molecules[i]/ (double) settings.getNbCells())>= FluxSimulatorSettings.SF_CELL_HI)
				sfHi.add(ids[i]);
			else if (Math.ceil(molecules[i]/ (double) settings.getNbCells())>= FluxSimulatorSettings.SF_CELL_MED)
				sfMed.add(ids[i]);
			else if (molecules[i]> 0)
				sfLo.add(ids[i]);
		}
		return true;
	}

	public HashSet<CharSequence> getSfHi() {
		if (sfHi == null) {
			createXPRmaps();
		}

		return sfHi;
	}



	public HashSet<CharSequence> getSfMed() {
		if (sfMed == null) {
			createXPRmaps();
		}

		return sfMed;
	}



	public HashSet<CharSequence> getSfLo() {
		if (sfLo == null) {
			createXPRmaps();
		}

		return sfLo;
	}
	
	public boolean loadStats() {
		if (settings.getProFile()== null|| (!settings.getProFile().exists()))
			return false;

		int lim= Integer.MAX_VALUE;	// last working token
		try {
			if (Constants.progress!= null) {
				Constants.progress.setString("initializing profiler ");
				Constants.progress.setValue(0);
			}
			int lines= FileHelper.countLines(settings.getProFile().getAbsolutePath());
			//String lineSep= FileHelper.getLineSeparator() // TODO
			ids= new ByteArrayCharSequence[lines]; 
			locIDs= new ByteArrayCharSequence[lines]; 
			len= new int[lines];
			molecules= new long[lines];
			cds= new boolean[lines];
			HashMap<ByteArrayCharSequence,ByteArrayCharSequence> locIDset= 
				new HashMap<ByteArrayCharSequence,ByteArrayCharSequence>();	// TODO MyHashSet.get(Object o)
			int ptr= -1, perc= 0;
			long bytesRead= 0, bytesTot= settings.getProFile().length();
			mapLenExp= new Hashtable<ByteArrayCharSequence, int[]>();
			
			BufferedInputStream istream= new BufferedInputStream(new FileInputStream(settings.getProFile()));
			ThreadedBufferedByteArrayStream buffy= 
				new ThreadedBufferedByteArrayStream(10* 1024, istream, true, false);
			ByteArrayCharSequence cs= new ByteArrayCharSequence(1024);
			for (buffy.readLine(cs); cs.end> 0; buffy.readLine(cs)) {
				
				cs.resetFind();
				
				bytesRead+= cs.length()+ 1;	// TODO fs
				if (ptr% 1000== 0&& bytesRead* 10d/ bytesTot> perc) {
					++perc;
					if (Constants.progress!= null)
						Constants.progress.progress();	// setValue(perc)
					else if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) { 
						System.err.print("*");
						System.err.flush();
					}
				}
				if (lim< 3)
					break;	// give up
				
				++ptr;
				int tok= 0;
				ByteArrayCharSequence x= cs.getToken(tok++);
				if (x== null) { 
					lim= tok- 2;
					continue;
				} else {
					if (locIDset.containsKey(x)) 
						locIDs[ptr]= locIDset.get(x);						
					else {
						locIDs[ptr]= x.cloneCurrentSeq();
						locIDset.put(locIDs[ptr],locIDs[ptr]);	
					}
				}
				
				if (lim< tok)
					continue;
				x= cs.getToken(tok++);
				if (x== null) {
					lim= tok- 2;
					continue;
				} else
					ids[ptr]= x.cloneCurrentSeq();
				
				if (lim< tok)
					continue;
				x= cs.getToken(tok++);
				if (x== null) {
					lim= tok- 2;
					continue;
				} else {
					if (x.equals(FluxSimulatorSettings.PRO_FILE_CDS))
						cds[ptr]= true;
					else if (x.equals(FluxSimulatorSettings.PRO_FILE_NC))
						cds[ptr]= false;
					else {
						lim= 1;
						continue;
					}
				}

				if (lim< tok)
					continue;
				int y= cs.getTokenInt(tok++);
				if (y== Integer.MIN_VALUE) {
					lim= tok- 2;
					continue;
				} else 
					len[ptr]= y;
				
				tok++;	// perc
				
				if (lim< tok)
					continue;
				y= cs.getTokenInt(tok++);
				if (y== Integer.MIN_VALUE) 
					lim= tok- 3;
				else {
					molecules[ptr]= y;
					if (molecules[ptr]> 0) {
							mapLenExp.put(getCombinedID(ptr),	// Issue32: ids[ptr]  
								new int[] {len[ptr], (int) molecules[ptr]});
					}
					lim= tok- 1;
				}
			}
			istream.close();
			buffy.close();
		} catch (Exception e) {
			lim= -1; // :)
			if (Constants.progress!= null)
				Constants.progress.finish();
			e.printStackTrace();
			return false;
		}
		
		if (lim< 4)
			molecules= null;
		if (lim< 3)
			len= null;
		if (lim< 2) {
			ids= null;
			locIDs= null;
		} else if (!isStop())
			calcStats();
		if (Constants.progress!= null)
			Constants.progress.finish();
		return true;
	}



	public ByteArrayCharSequence[] getLocIDs() {
		return locIDs;
	}



	public void setLocIDs(ByteArrayCharSequence[] locIDs) {
		this.locIDs = locIDs;
	}



	public boolean[] getCds() {
		return cds;
	}



	public void setCds(boolean[] cds) {
		this.cds = cds;
	}



	public void setMolecules(long[] molecules) {
		this.molecules = molecules;
	}



	public boolean setStop() {
		this.stop= true;
		if (gffReader!= null)
			return gffReader.setStop();
		return true;
	}



	public ByteArrayCharSequence getCombinedID(int i) {
		if (i< 0|| i>= getIds().length)
			return null;
		ByteArrayCharSequence locID= getLocIDs()[i];
		ByteArrayCharSequence tID= getIds()[i];
		
		ByteArrayCharSequence cs= new ByteArrayCharSequence(locID.length()+ tID.length()+ 1);
		System.arraycopy(locID.a, locID.start, cs.a, cs.end, locID.length());
		cs.end+= locID.length();
		cs.a[cs.end++]= FluxSimulatorSettings.SEP_LOC_TID;
		System.arraycopy(tID.a, tID.start, cs.a, cs.end, tID.length());
		cs.end+= tID.length();
		
		return cs;
	}
}
