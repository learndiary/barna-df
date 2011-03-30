package fbi.genome.sequencing.rnaseq.simulation;

import fbi.commons.file.FileHelper;
import fbi.genome.io.bed.BEDwrapper;
import fbi.genome.io.gff.GFFReader;
import fbi.genome.model.Gene;
import fbi.genome.model.Graph;
import fbi.genome.model.Transcript;
import fbi.genome.model.bed.BEDobject2;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Random;

public class KmerExtractor {

	/**
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
				
				int flank5= 5, flank3= 5, kmer= 6;				
				int mapSize= (int) (Math.pow(4, kmer)* 2);
				HashMap<CharSequence, double[]> mapSense= new HashMap<CharSequence, double[]>(mapSize),
												mapASense= new HashMap<CharSequence, double[]>(mapSize);
				
				
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
	
						BEDobject2[] beds= bedReader.read(t.getChromosome(), Math.abs(t.getStart()), Math.abs(t.getEnd()));
						if (beds== null|| beds.length< 100)
							continue;
						++cntTxWreads;
						if (t.getStrand()>0)
							++txPlus;
						else
							++txMinus;
						String s= t.getSplicedSequence().toUpperCase();
						int readNr= 0, invalid= 0;
						while(readNr< 100&& readNr+ invalid< beds.length) {	// from: invalid+ (beds.length- readNr)!= beds.length
							int p= rnd.nextInt(beds.length);
							while(beds[p]== null)
								p= rnd.nextInt(beds.length);
							boolean sens= beds[p].getStrand()== t.getStrand();
							int tstart= t.getExonicPosition(beds[p].getStart()+ 1),
								tend= t.getExonicPosition(beds[p].getEnd());	// t-coordinates, 0-based
							if (tstart== Integer.MIN_VALUE|| tend== Integer.MIN_VALUE
									|| tstart< 0|| tstart>= s.length() 
									|| ((sens&& (tstart< flank5|| tstart+ flank3+ 1> s.length()))
									|| (tend< flank3|| tend+ flank5+ 1> s.length()))) {
								++invalid;
								if (beds[p].getStrand()> 0)
									++invalidPlus;
								else
									++invalidMinus;
								beds[p]= null;
								continue;
							}
							beds[p]= null;
							++readNr;
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
							
							// count on subsequence
							boolean sens= beds[j].getStrand()== t.getStrand();
							HashMap<CharSequence, double[]> map= sens? mapSense: mapASense;
							for (int n = -flank5; n <= flank3; ++n) {
								String k= null;
								if (sens) {
									if (tstart== Integer.MIN_VALUE|| tstart+ n< 0|| tstart+ n+ 6>= s.length())
										continue;
									k= s.substring(tstart+ n, tstart+ n+ 6);
								} else {	// read asense to s
									if (tend== Integer.MIN_VALUE|| tend- n- 6< 0|| tend- n>= s.length())
										continue;
									k= s.substring(tend- n- 6+ 1, tend- n+ 1);
									k= Graph.reverseSequence(k);
									k= Graph.complementarySequence(k);
								}							
								++cntReads;
								if (map.containsKey(k)) 
									++map.get(k)[n+ flank5];
								else {
									double[] d= new double[flank3+ flank5+ 1];
									Arrays.fill(d, 0d);
									++d[n+flank5];
									map.put(k, d);
								}
							}
	
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
				output(mapSense, flank5+ 1, fileBed+"_sense_kmat");
				output(mapASense, flank5+ 1, fileBed+"_asense_kmat");
				System.err.println("wrote "+ fileBed+"_sense.kmat,\n\tand "+ fileBed+"_asense.kmat");
				
			} catch (Exception e) {
				e.printStackTrace();
			}
			
			return ref;
		}

	private static void output(HashMap<CharSequence, double[]> mapSense,
			int firstPos, String fName) {
		try {
			BufferedWriter writer= new BufferedWriter(new FileWriter(fName));
			writer.write("#\t"+ firstPos+ "\n");
			Object[] kk= mapSense.keySet().toArray();
			Arrays.sort(kk);
			for (int i = 0; i < kk.length; i++) {
				writer.write(kk[i].toString());	// TODO
				double[] d= mapSense.get(kk[i]);
				for (int k = 0; k < d.length; k++) 
					writer.write("\t"+ d[k]);
				writer.write("\n");
			}
			writer.flush();
			writer.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
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

}
