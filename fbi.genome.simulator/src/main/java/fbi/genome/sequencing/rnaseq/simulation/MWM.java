package fbi.genome.sequencing.rnaseq.simulation;


import commons.file.FileHelper;
import fbi.genome.io.bed.BEDwrapper;
import fbi.genome.io.gff.GFFReader;
import fbi.genome.model.Gene;
import fbi.genome.model.Graph;
import fbi.genome.model.Transcript;
import fbi.genome.model.bed.BEDobject2;
import fbi.genome.model.constants.Constants;

import java.io.*;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Random;

public class MWM implements WeightMatrix {

	HashMap<String, Double> map= null;
	long sum= 0;
	int motifLen= -1;
	boolean plus;
	
	public static MWM create(File f, boolean plus) {
		try {
			BufferedReader buffy= new BufferedReader(new FileReader(f));
			MWM myMWM= new MWM();
			myMWM.plus= plus;
			myMWM.map= new HashMap<String, Double>(10000);	// TODO
			myMWM.sum= 0l;
			int ctr= 0;
			for (String s; (s= buffy.readLine())!= null;++ctr) {
				String[] ss= s.split("\\s");
				if ((plus&& !ss[0].equals("+"))|| ((!plus)&& !ss[0].equals("-")))
					continue;
				int mol= Integer.parseInt(ss[2]);
				if(myMWM.motifLen== -1)
					myMWM.motifLen= ss[1].length();
				else {
					if (myMWM.motifLen!= ss[1].length()) {
						if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) {
							System.err.println("\n\t[UH-OH] motif length not consistent:");
							System.err.println("\tfile\t"+ f.getAbsolutePath());
							System.err.println("\tline\t"+ (ctr+ 1));
							System.err.println("\tlength\t"+ ss[1].length()+"\t(before "+myMWM.motifLen+")");
						}
						return null;
					}
				}
				myMWM.map.put(ss[1], (double) mol);
				myMWM.sum+= mol;
			}
			
			// rel. freq
			Iterator<String> iter= myMWM.map.keySet().iterator();
			while(iter.hasNext()) {
				String k= iter.next();
				double v= myMWM.map.remove(k);
				v/= myMWM.sum;
				myMWM.map.put(k, v);
			}
			
			return myMWM;
		} catch (Exception e) {
			e.printStackTrace();
		}
		return null;
	}
	
	protected MWM(){
	}
	
	
	public float[] apply(CharSequence s) {
		float[] a= new float[s.length()];
		Arrays.fill(a, 0);
		for (int i = (plus? 0: motifLen- 1); i < s.length()- (plus? motifLen- 1: 0); i++) {
			int x= plus? i: i- motifLen+ 1;
			String k= s.subSequence(x, x+ motifLen).toString();
			if (!map.containsKey(k))
				a[i]= 0f;
			else
				a[i]= (float) map.get(k).doubleValue();
		}
		return a;
	}

	/*
		 * 
	/Users/micha/projects/simulator/data/wold08/mm_brain_spiked/SRR001356_uniq.bed
	/Users/micha/projects/simulator/data/wold08/mm_brain_spiked/SRR001357_uniq.bed
	/Users/micha/projects/simulator/data/wold08/mm_liver_spiked/SRR001358_uniq.bed
	/Users/micha/projects/simulator/data/wold08/mm_liver_spiked/SRR001359_uniq.bed
	/Users/micha/projects/simulator/data/wold08/mm_liver_spiked/SRR001360_uniq.bed
	/Users/micha/projects/simulator/data/wold08/mm_muscle_spiked/SRR001361_uniq.bed
	/Users/micha/projects/simulator/data/wold08/mm_muscle_spiked/SRR001362_uniq.bed
		 */
		/*
	/Users/micha/projects/simulator/data/wold10/mm_muscle_-24h/SRR037945_uniq.bed
	/Users/micha/projects/simulator/data/wold10/mm_muscle_-24h/SRR037946_uniq.bed
	/Users/micha/projects/simulator/data/wold10/mm_muscle_120h/SRR037950_uniq.bed
	/Users/micha/projects/simulator/data/wold10/mm_muscle_120h/SRR037951_uniq.bed
	/Users/micha/projects/simulator/data/wold10/mm_muscle_168h/SRR037952_uniq.bed
	/Users/micha/projects/simulator/data/wold10/mm_muscle_168h/SRR037953_uniq.bed
	/Users/micha/projects/simulator/data/wold10/mm_muscle_168h/SRR037954_uniq.bed
	/Users/micha/projects/simulator/data/wold10/mm_muscle_60h/SRR037947_uniq.bed
	/Users/micha/projects/simulator/data/wold10/mm_muscle_60h/SRR037948_uniq.bed
	/Users/micha/projects/simulator/data/wold10/mm_muscle_60h/SRR037949_uniq.bed
		 */
		
		/*
		 * /Users/micha/annotation/mm9_UCSCgenes_fromUCSC101210_sorted.gtf 
	/Users/micha/genomes/mm9/
	/Users/micha/projects/simulator/data/wold10/mm_muscle_-24h/SRR037945_uniq.bed
		 */
		public static void main(String[] args) {
			
	/*		
			if (args.length< 3) {
				System.err.println("Usage: PWM <GTF> <genomeDir> <BED1> [BEDi]");
				System.exit(-1);
			}
	*/
			args= new String[] {
					"/Users/micha/annotation/mm9_UCSCgenes_fromUCSC101210_sorted.gtf", 
					"/Users/micha/genomes/mm9/",
					"/Users/micha/projects/simulator/data/wold10/readData/SRR037945.bed_sorted"
			};
			
			
			Graph.overrideSequenceDirPath= args[1];
			System.err.println("processing.. "+ args[2]);
			args[0]= extractMatrices(args[0], args[2]);
		}

		/**
		 * @deprecated
		 * @param fileGTF
		 * @param fileBed
		 * @return
		 */
		private static String extractMatrices(String fileGTF, String fileBed) {
			
			Random rnd= new Random();
			String ref= fileGTF;
			try {
				
	//			HashMap<String, IntVector[]> mapExprWeights= null;
	//			if (fileExprWeights!= null)
	//				mapExprWeights= getExprWeights(fileExprWeights);
				
				int flank5= 30, flank3= 50;
				
				GFFReader anoReader= new GFFReader(fileGTF);
				if(false&& !anoReader.isApplicable()){
					System.err.println("\tsorting GTF file");
					File f= anoReader.createSortedFile();
					System.err.println("\tsorted file in "+ f.getAbsolutePath());
					ref= f.getAbsolutePath();
					anoReader= new GFFReader(f.getAbsolutePath());
				}
				System.err.println();
				
				File ff= new File(fileBed+ "_sorted");
				BEDwrapper bedReader= null;
				if (ff.exists()) {
					System.err.println("\tusing sorted file "+ ff.getName());
					fileBed= ff.getAbsolutePath();
					bedReader= new BEDwrapper(fileBed);
				} else {
					bedReader= new BEDwrapper(fileBed);
					if(false&& !bedReader.isApplicable()) {
						System.err.println("\tsorting BED file");
						File f= bedReader.sortBED(new File(fileBed));
						if (FileHelper.move(f, ff, null))
							fileBed= ff.getAbsolutePath();
						else
							fileBed= f.getAbsolutePath();
						System.err.println("\tsorted file in "+ fileBed);
						bedReader= new BEDwrapper(fileBed);
					}
				}
	
				System.err.print("\tprocessing ");
				System.err.flush();
				int cntTrpt= 0, cntTrptP= 0, cntTrptM= 0, cntTxWreads= 0, cntReads= 0;
				HashMap<String, Integer> mapSense= new HashMap<String, Integer>(1000), mapASense= new HashMap<String, Integer>(1000);  

				Gene[] genes= null;
				int invalidPlus= 0, invalidMinus= 0, cntPlus= 0, cntMinus= 0, txPlus= 0, txMinus= 0;
				String lastChr= null;
				byte lastStrand= 0;
				for (anoReader.read();(genes= anoReader.getGenes())!= null; anoReader.read()) {
					for (int i = 0; i < genes.length; i++) {
						if (genes[i].getTranscriptCount()> 1)	// non-AS genes
							continue;
						++cntTrpt;
						if (genes[i].getStrand()> 0)
							++cntTrptP;
						else
							++cntTrptM;
						Transcript t= genes[i].getTranscripts()[0];
						if (lastStrand== 0)
							lastStrand= genes[i].getStrand();
						else {
							if (genes[i].getStrand()!= lastStrand&& genes[i].getChromosome().equals(lastChr)) 
								bedReader.reset(lastChr);
						}
						lastChr= genes[i].getChromosome();
						lastStrand= genes[i].getStrand();

						// random sequences
						if (1== 1) {
							String s= t.getSplicedSequence().toUpperCase();
							if (s.length()<= 2*(flank5+ flank3+ 1))
								continue;
							for (int j = 0; j < 100; j++) {
								int p= rnd.nextInt(s.length());								
								boolean sens= j% 2== 0;
								while (((sens&& (p< flank5|| p+ flank3+ 1> s.length()))
										|| (p< flank3|| p+ flank5+ 1> s.length()))) {
									p= rnd.nextInt(s.length());
								}

								String k= null;
								if (sens) {
									k= s.substring(p- flank5, p+ flank3+ 1);
								} else {	// read asense to s
									k= s.substring(p- flank3, p+ flank5+ 1);
									k= Graph.reverseSequence(k);
									k= Graph.complementarySequence(k);
								}
								++cntReads;
								HashMap<String, Integer> map= sens? mapSense: mapASense;
								if (map.containsKey(k))
									map.put(k, map.get(k)+ 1);
								else
									map.put(k, 1);
							}
							continue;
						}
						
						BEDobject2[] beds= bedReader.read(t.getChromosome(), Math.abs(t.getStart()), Math.abs(t.getEnd()));
						if (beds== null|| beds.length< 100)
							continue;
						++cntTxWreads;
						if (t.getStrand()>0)
							++txPlus;
						else
							++txMinus;
						String s= t.getSplicedSequence().toUpperCase();
						int readNr= beds.length, invalid= 0;
						while(readNr> 100&& invalid!= readNr) {	// from: invalid+ (beds.length- readNr)!= beds.length
							int p= rnd.nextInt(beds.length);
							while(beds[p]== null)
								p= rnd.nextInt(beds.length);
							boolean sens= beds[p].getStrand()== t.getStrand();
							int tstart= t.getExonicPosition(beds[p].getStart()+ 1),
								tend= t.getExonicPosition(beds[p].getEnd());	// t-coordinates, 0-based
							if (tstart< 0|| tstart>= s.length()|| 
									((sens&& (tstart< flank5|| tstart+ flank3+ 1> s.length()))
									|| (tend< flank3|| tend+ flank5+ 1> s.length()))) {
								++invalid;
								if (beds[p].getStrand()> 0)
									++invalidPlus;
								else
									++invalidMinus;
								continue;
							}
							beds[p]= null;
							--readNr;
						}
						for (int j = 0; j < beds.length; j++) {
							
							if (beds[j]== null)
								continue;
							if (beds[j].getStrand()> 0)
								++cntPlus;
							else
								++cntMinus;
							// get t-coordinates
							int tstart= t.getExonicPosition(beds[j].getStart()+ 1),
								tend= t.getExonicPosition(beds[j].getEnd());	// t-coordinates, 0-based
							if (tstart< 0|| tstart>= s.length())
								continue;
							
							// count on subsequence
							boolean sens= beds[j].getStrand()== t.getStrand();
							HashMap<String, Integer> map= sens? mapSense: mapASense;
							String k= null;
							if (sens) {
								if (tstart< flank5|| tstart+ flank3+ 1> s.length())
									continue;
								k= s.substring(tstart- flank5, tstart+ flank3+ 1);
							} else {	// read asense to s
								if (tend< flank3|| tend+ flank5+ 1> s.length())
									continue;
								k= s.substring(tend- flank3, tend+ flank5+ 1);
								k= Graph.reverseSequence(k);
								k= Graph.complementarySequence(k);
							}							
							++cntReads;
							if (map.containsKey(k))
								map.put(k, map.get(k)+ 1);
							else
								map.put(k, 1);
	
						}
					}
				}
				System.err.println(" OK");
				System.err.println("\tFound "+ cntTrpt+ " non-AS, " + cntTrptP+ "+, "+ cntTrptM+ "-"+
						"\n\t"+ cntTxWreads+" of them with >=100 reads, " + txPlus+ "+, "+ txMinus+ "-"+
						"\n\t"+ invalidPlus+" inv+, "+ invalidMinus+ "inv-"+
						"\n\t"+ cntReads + " reads in total, "+ cntPlus+ "+, "+ cntMinus+ "-");
				System.err.println();
				
				// output
				BufferedWriter writer= new BufferedWriter(new FileWriter(fileBed+"_sense_random.kmer"));
				Iterator<String> iter= mapSense.keySet().iterator();
				while(iter.hasNext()) {
					String k= iter.next();
					writer.write(k+ "\t"+ mapSense.get(k)+ "\n");
				}
				writer.flush();
				writer.close();
				writer= new BufferedWriter(new FileWriter(fileBed+"_asense_random.kmer"));
				iter= mapASense.keySet().iterator();
				while(iter.hasNext()) {
					String k= iter.next();
					writer.write(k+ "\t"+ mapASense.get(k)+ "\n");
				}
				writer.flush();
				writer.close();
				System.err.println("wrote "+ fileBed+"_sense.kmer,\n\tand "+ fileBed+"_asense.kmer.");
				
			} catch (Exception e) {
				e.printStackTrace();
			}
			
			return ref;
		}
}
