package fbi.genome.sequencing.rnaseq.simulation;

import fbi.commons.file.FileHelper;
import fbi.genome.io.bed.BEDwrapper;
import fbi.genome.io.gff.GFFReader;
import fbi.genome.model.Gene;
import fbi.genome.model.Graph;
import fbi.genome.model.Transcript;
import fbi.genome.model.bed.BEDobject2;
import fbi.genome.model.commons.IntVector;

import java.io.*;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Vector;

public class PWM implements WeightMatrix {
	double[][] pwm= null;
	// either sorted positions, or continuous range [min;max]
	int[] pos= null;	// 0-based, 0= ref.position (+1)
	int min= Integer.MAX_VALUE, max= Integer.MIN_VALUE;
	
	static final int MIN_DEFAULT= Integer.MAX_VALUE, MAX_DEFAULT= Integer.MIN_VALUE;
	
	public static PWM create(File f) throws Exception{
		BufferedReader buffy= new BufferedReader(new FileReader(f));
		Vector<String> v= new Vector<String>();		
		for (String s= null; (s= buffy.readLine())!= null;) {
			s= s.trim();
			if (s.length()== 0|| s.startsWith("#"))
				continue;
			v.add(s);
		}
		if (v.size()== 0)
			return null;
		
		// positions
		String[] ss= v.elementAt(0).split("\\s");
		int[] pos;
		if (ss.length== 1&& ss[0].contains(";")) {	// min, max
			int p= ss[0].indexOf(";");
			int min= Integer.parseInt(ss[0].substring(0, p)), max= Integer.parseInt(ss[0].substring(p+1));
			pos= new int[max- min+ 1];
			for (int i = min; i <= max; i++) 
				pos[i- min]= i;
		} else {
			pos= new int[ss.length];
			for (int i = 0; i < pos.length; i++) 
				pos[i]= Integer.parseInt(ss[i]);
		}
		
		// matrix
		CharSequence[] kmers= new CharSequence[v.size()- 1];
		double[][] a= new double[v.size()- 1][];
		for (int i = 1; i < v.size(); i++) {
			ss= v.elementAt(i).split("\\s");
			assert(ss.length- 1== pos.length);
			a[i-1]= new double[ss.length- 1];
			kmers[i-1]= ss[0];
			for (int k = 1; k < ss.length; k++) 
				a[i-1][k- 1]= Double.parseDouble(ss[k]);
		}
		
		return new PWM(kmers, pos, a);
	}

	public PWM(int min, int max, double[][] pwm) {
		this.min= min;
		this.max= max;
		setPWM(pwm);
	}

	public PWM(int[] pos, double[][] pwm) {
		set(pos, pwm);
		this.from= 0;
		this.to= pwm.length- 1;
	}
	
	CharSequence[] kmers;
	int kmerLen= 0;
	public PWM(CharSequence[] kmers, int[] pos, double[][] pwm) {
		//set(pos, pwm);
		this.kmers= kmers;
		kmerLen= kmers[0].length();
		this.pos= pos;
		this.pwm= pwm;
		sortKmers();
		sortPositions();
	}
	
	public void set(int[] pos2, double[][] pwm2) {
		this.pos= new int[pos2.length];
		System.arraycopy(pos2, 0, pos, 0, pos2.length);
		Arrays.sort(pos);
		this.pwm= new double[pwm2.length][];
		
		for (int i = 0; i < pwm2.length; i++) {
			int p= Arrays.binarySearch(pos, pos2[i]);
			pwm[p]= pwm2[i];
		}
		
		System.currentTimeMillis();
		
	}

	void setPWM(double[][] pwm) {
		
		// always norm to 1
		//if (pwm[0][0]> 1&& pwm[0][1]> 1) {
			for (int i = 0; i < pwm.length; i++) {
				double s= 0;
				for (int j = 0; j < pwm[i].length; j++) 
					s+= pwm[i][j];
				for (int j = 0; j < pwm[i].length; j++) { 
					pwm[i][j]= (pwm[i][j]/ s);
					//pwm[i][j]= Math.log(pwm[i][j]);
					//System.currentTimeMillis();
				}
			}
		//} 
		this.pwm= pwm;
		
	}

	/**
	 * @deprecated
	 */
	public float[] apply(CharSequence s) {
		if (pwm== null|| ((min== MIN_DEFAULT|| max== MAX_DEFAULT)&& pos== null)) 
			throw new RuntimeException("[PWM] Position-weight matrix not initialized.");
		float[] a= new float[s.length()];
		int[] ci= new int[s.length()];
		for (int i = 0; i < ci.length; i++) {
			char c= Character.toUpperCase(s.charAt(i));
			ci[i]= ((c=='A'|| c=='a')?0:((c=='C'|| c=='c')?1:((c=='G'|| c== 'g')?2:((c=='T'|| c== 't')?3:-1))));
		}
		float min= Float.MAX_VALUE, max= -1f;
		for (int i = 0; i < a.length; i++) {
			a[i]= apply(ci, i);
			min=Math.min(min, a[i]);
			max=Math.max(max, a[i]);
		}
		// stretch to [0,1]
/*		double d= max- min;
		for (int i = 0; i < a.length; i++) {
			double f= (a[i]- min)/ d;
			a[i]= (float) f;
		}
*/		
		return a;
	}
	

	int from, to;
	
	
	/**
	 * 
	 * @param ci s encoded as colnr of letter
	 * @param p position in s
	 * @return
	 */
	private float apply(int[] ci, int p) {
	
		double w= 1;
		// boundaries on s
/*		int lo= p+ (pos== null? min: pos[0]), up= p+ (pos== null? max: pos[pos.length- 1]);
		// boundaries on s
		int mlo= 0, mup= pwm.length;
		if (lo< 0) {
			mlo= -lo;
		}
		for (int i= 1;up> ci.length;++i) {
			--mup;			
			up= p+ (pos== null? max- i: pos[pos.length- 1- i]);
		}
*/		
		int c= 0;
		for (int i = from; i <= to; i++) {
			int pp= (pos== null)? p+ min+ i: p+ pos[i];
			if (pp< 0) 
				continue;
			if (pp>= ci.length) 
				break;
			
			++c;
			w*= pwm[i][ci[pp]];
		}
		//if (c< pwm.length) {			
			// int d= pwm.length- c;
		
			w= Math.pow(w, 1d/ c);
		
		// }
		return (float) w;
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
				"/Users/micha/annotation/spike_sequences.gtf", 
				"/Users/micha/genomes/spikes/",
				null
		};
		
		File dir= new File("/Users/micha/projects/simulator/data/wold10/");
		String[] files= dir.list();
		args= new String[30];
		args[0]= "/Users/micha/annotation/spike_sequences.gtf"; 
		args[1]= "/Users/micha/genomes/spikes/";
		int c=2;
		for (int i = 0; i < files.length; i++) {
			if (files[i].startsWith("exp")&& files[i].endsWith("bed"))
				args[c++]= dir.getAbsolutePath()+ File.separator+ files[i];
		}
		
		
		Graph.overrideSequenceDirPath= args[1];
		for (int i = 2; i < args.length; i++) {
			System.err.println("processing.. "+ args[i]);
			args[0]= extractMatrices(args[0], args[i]);
		}
	}

	private static String extractMatrices(String fileGTF, String fileBed) {
		
		String ref= fileGTF;
		try {
			
//			HashMap<String, IntVector[]> mapExprWeights= null;
//			if (fileExprWeights!= null)
//				mapExprWeights= getExprWeights(fileExprWeights);
			
			int flank5= 10, flank3= 20;
			
			GFFReader anoReader= new GFFReader(fileGTF);
			if(!anoReader.isApplicable()){
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
				if(!bedReader.isApplicable()) {
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
			int cntTrpt= 0, cntReads= 0;
			int[][] sense= new int[flank5+ flank3][], asense= new int[flank5+ flank3][];
			for (int i = 0; i < asense.length; i++) {
				sense[i]= new int[4];
				asense[i]= new int[4];
				for (int j = 0; j < asense[i].length; j++) {
					sense[i][j]= 0;
					asense[i][j]= 0;
				}
			}
			Gene[] genes= null;
			for (anoReader.read();(genes= anoReader.getGenes())!= null; anoReader.read()) {
				for (int i = 0; i < genes.length; i++) {
					if (genes[i].getTranscriptCount()> 1)	// non-AS genes
						continue;
					++cntTrpt;
					Transcript t= genes[i].getTranscripts()[0];
					BEDobject2[] beds= bedReader.read(t.getChromosome(), t.getStart(), t.getEnd());
					if (beds== null)
						continue;
					String s= t.getSplicedSequence().toUpperCase();
					for (int j = 0; j < beds.length; j++) {
						// get t-coordinates
						int tstart= t.getExonicPosition(beds[j].getStart()+ 1),
							tend= t.getExonicPosition(beds[j].getEnd());	// t-coordinates, 0-based
						if (tstart< 0|| tstart>= s.length())
							continue;
						
						// count on subsequence
						++cntReads;
						boolean sens= beds[j].getStrand()== t.getStrand();
						int[][] a= sens? sense: asense;
						if (sens) {	
							for (int k = 0; k < a.length; ++k) {
								int p= tstart- flank5+ k;
								if (p< 0)
									continue;
								if (p>= s.length())
									break;
								if (s.charAt(p)== 'A')
									++a[k][0];
								else if (s.charAt(p)== 'C')
									++a[k][1];
								else if (s.charAt(p)== 'G')
									++a[k][2];
								else if (s.charAt(p)== 'T')
									++a[k][3];
							}
							
						} else {	// read asense to s
							for (int k = 0; k < a.length; ++k) {
								int p= tend+ flank5- k;
								if (p>= s.length())
									continue;
								if (p< 0)
									break;
								// reverse complement?
/*								if (s.charAt(p)== 'A')
									++a[k][3];	// A -> count T
								else if (s.charAt(p)== 'C')
									++a[k][2];
								else if (s.charAt(p)== 'G')
									++a[k][1];
								else if (s.charAt(p)== 'T')
									++a[k][0];
*/								if (s.charAt(p)== 'A')
									++a[k][0];
								else if (s.charAt(p)== 'C')
									++a[k][1];
								else if (s.charAt(p)== 'G')
									++a[k][2];
								else if (s.charAt(p)== 'T')
									++a[k][3];
							}
						}

					}
				}
			}
			System.err.println(" OK");
			System.err.println("\tFound "+ cntTrpt+ " non-AS tx with "+ cntReads+ " reads.");
			System.err.println();
			
			// output
			BufferedWriter writer= new BufferedWriter(new FileWriter(fileBed+"_sense.pwm"));
			for (int i = 0; i < sense.length; i++) {
				int pos= (i>= flank5? i- flank5+ 1: i- flank5);
				writer.write(pos+ "\t"+ sense[i][0]+ "\t"+ sense[i][1]+ "\t"+ sense[i][2]+ "\t"+ sense[i][3]+ "\n");
			}
			writer.flush();
			writer.close();
			writer= new BufferedWriter(new FileWriter(fileBed+"_asense.pwm"));
			for (int i = 0; i < asense.length; i++) {
				int pos= (i>= flank5? i- flank5+ 1: i- flank5);
				writer.write(pos+ "\t"+ asense[i][0]+ "\t"+ asense[i][1]+ "\t"+ asense[i][2]+ "\t"+ asense[i][3]+ "\n");
			}
			writer.flush();
			writer.close();
			System.err.println("wrote "+ fileBed+"_sense.pwm, and "+ fileBed+"_asense.pwm.");
			
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		return ref;
	}

	private static HashMap<String, IntVector[]> getExprWeights(
			String fileExprWeights) {
		
		try {
			HashMap<String, IntVector[]> map= new HashMap<String, IntVector[]>();
			BufferedReader buffy= new BufferedReader(new FileReader(fileExprWeights));
			IntVector[] vStartEnd;
			for (String s= null; (s= buffy.readLine())!= null; ) {
				String[] ss= s.split("\\s");
				if (map.containsKey(ss[0]))
					vStartEnd= map.get(ss[0]);
				else {
					vStartEnd= new IntVector[2];
					vStartEnd[0]= new IntVector(1000);
					vStartEnd[1]= new IntVector(1000);
					map.put(ss[0], vStartEnd);
				}
				int start= Integer.parseInt(ss[1]);
			}
		} catch (Exception e) {
			// TODO: handle exception
		}
		
		return null;
	}

	public void setFrom(int from) {
		if (pos!= null) {
			int i1= Arrays.binarySearch(pos, from);
			if (i1>= 0)
				this.from= i1;
		}
	}

	public void setTo(int to) {
		if (pos!= null) {
			int i1= Arrays.binarySearch(pos, to);
			if (i1>= 0)
				this.to= i1;
		}
	}
	
	private void sortKmers() {
		CharSequence[] kmers2= kmers.clone();
		Arrays.sort(kmers2);
		
		for (int i = 0; i < kmers.length; i++) {
			int p= Arrays.binarySearch(kmers2, kmers[i]);
			if (p== i)
				continue;
			assert(p>= 0);
			
			double[] h= pwm[p];
			pwm[p]= pwm[i];
			pwm[i]= h;
			CharSequence c= kmers[p];
			kmers[p]= kmers[i];
			kmers[i]= c;
			
			--i;
		}
	}
	
	private void sortPositions() {
		int[] pos2= pos.clone();
		Arrays.sort(pos2);
		
		for (int i = 0; i < pos.length; i++) {
			int p= Arrays.binarySearch(pos2, pos[i]);
			if (p== i)
				continue;
			assert(p>= 0);
			
			for (int j = 0; j < pwm.length; j++) {	// all kmers
				double h= pwm[j][p];
				pwm[j][p]= pwm[j][i];
				pwm[j][i]= h;
			}
			int d= pos[p];
			pos[p]= pos[i];
			pos[i]= d;
			
			--i;
		}
	}
	
	public void invert() {
		for (int i = 0; i < pos.length; i++) 
			pos[i]= -pos[i]- (kmerLen- 1);
		sortPositions();
		for (int i = 0; i < kmers.length; i++) 
			kmers[i]= Graph.invertSequence(kmers[i].toString());// TODO work on charsequence
		sortKmers();
	}
	
	public void multiply() {

		for (int j = 0; j < pwm[0].length; ++j) {	// positions
			for (int i = 0; i < pwm.length; ++i) {	// kmers
				if (pwm[i][j]< 1)
					pwm[i][j]= -Math.log(pwm[i][j])*2;
				else
					pwm[i][j]*= 2;
			}
		}
	}
	
	public void makePDF() {
		for (int j = 0; j < pwm[0].length; ++j) {	// positions
			double colSum= 0d;
			for (int i = 0; i < pwm.length; ++i) 	// kmers
				colSum+= pwm[i][j];
			//assert(colSum!= 0);
			if (colSum== 0)
				return; // throw some runtime exception
			for (int i = 0; i < pwm.length; ++i) 
				pwm[i][j]/= colSum;
		}
	}
	
	
	/**
	 * 
	 * @param ci s encoded as colnr of letter
	 * @param p position in s
	 * @return
	 */
	public double apply(CharSequence s, int p) {
		
		double w= 0; // sum logs, numerical stability
			// 1; for multiplying probabilities
		int c= 0;
		for (int i = 0; i < pos.length; i++) {
			int pp= p+ pos[i];
			if (pp< 0) 
				continue;
			if (pp+ kmerLen>= s.length()) 
				break;
			
			++c;
			CharSequence seq= s.subSequence(pp, pp+ kmerLen);
			int kmerPos= Arrays.binarySearch(kmers, seq);

            if(kmerPos >= 0){
			    w+= Math.log(pwm[kmerPos][i]);
            }
			//w*= pwm[i][kmerPos];
		}
		if (c== 0) 
			return 0;	// avoid NaN
		
		// weight by information
		w= Math.exp(w/ c);
		//w= Math.pow(w, 1d/ c);

		assert(!(Double.isNaN(w)|| Double.isInfinite(w)|| w< 0));			
		return w;
	}

	/**
	 * @deprecated TO DO
	 * @param fileGTF
	 * @param fileBed
	 * @return
	 */
	private static String extractDinucleotideMatrices(String fileGTF, String fileBed) {
			
			String ref= fileGTF;
			try {
				int flank5= 10, flank3= 20;
				
				GFFReader anoReader= new GFFReader(fileGTF);
				if(!anoReader.isApplicable()){
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
					if(!bedReader.isApplicable()) {
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
				int cntTrpt= 0, cntReads= 0;
				int[][] sense= new int[flank5+ flank3][], asense= new int[flank5+ flank3][];
				int[][] diSense= new int[flank5+ flank3][], diASense= new int[flank5+ flank3][];
				for (int i = 0; i < asense.length; i++) {
					sense[i]= new int[4];
					asense[i]= new int[4];
					diSense[i]= new int[16];
					diASense[i]= new int[16];
					for (int j = 0; j < asense[i].length; j++) {
						sense[i][j]= 0;
						asense[i][j]= 0;
					}
					for (int j = 0; j < diASense.length; j++) {
						diSense[i][j]= 0;
						diASense[i][j]= 0;
					}
				}
				Gene[] genes= null;
				for (anoReader.read();(genes= anoReader.getGenes())!= null; anoReader.read()) {
					for (int i = 0; i < genes.length; i++) {
						if (genes[i].getTranscriptCount()> 1)
							continue;
						++cntTrpt;
						Transcript t= genes[i].getTranscripts()[0];
						BEDobject2[] beds= bedReader.read(t.getChromosome(), t.getStart(), t.getEnd());
						if (beds== null)
							continue;
						String s= t.getSplicedSequence().toUpperCase();
						for (int j = 0; j < beds.length; j++) {
							// get t-coordinates
							int tstart= t.getExonicPosition(beds[j].getStart()+ 1),
								tend= t.getExonicPosition(beds[j].getEnd());	// t-coordinates, 0-based
							if (tstart< 0|| tstart>= s.length())
								continue;
							
							// count on subsequence
							++cntReads;
							boolean sens= beds[j].getStrand()== t.getStrand();
							int[][] a= sens? sense: asense;
							int[][] b= sens? diSense: diASense;
							if (sens) {	
								for (int k = 0; k < a.length; ++k) {
									int p= tstart- flank5+ k;
									if (p< 0)
										continue;
									if (p>= s.length())
										break;
									if (s.charAt(p)== 'A') {
										++a[k][0];
										if (p> 0) {
											
										}
									} else if (s.charAt(p)== 'C')
										++a[k][1];
									else if (s.charAt(p)== 'G')
										++a[k][2];
									else if (s.charAt(p)== 'T')
										++a[k][3];
								}
								
							} else {	// read asense to s
								for (int k = 0; k < a.length; ++k) {
									int p= tend+ flank5- k;
									if (p>= s.length())
										continue;
									if (p< 0)
										break;
									// reverse complement?
	/*								if (s.charAt(p)== 'A')
										++a[k][3];	// A -> count T
									else if (s.charAt(p)== 'C')
										++a[k][2];
									else if (s.charAt(p)== 'G')
										++a[k][1];
									else if (s.charAt(p)== 'T')
										++a[k][0];
	*/								if (s.charAt(p)== 'A')
										++a[k][0];
									else if (s.charAt(p)== 'C')
										++a[k][1];
									else if (s.charAt(p)== 'G')
										++a[k][2];
									else if (s.charAt(p)== 'T')
										++a[k][3];
								}
							}
	
						}
					}
				}
				System.err.println(" OK");
				System.err.println("\tFound "+ cntTrpt+ " non-AS tx with "+ cntReads+ " reads.");
				System.err.println();
				
				// output
				BufferedWriter writer= new BufferedWriter(new FileWriter(fileBed+"_sense.pwm"));
				for (int i = 0; i < sense.length; i++) {
					int pos= (i>= flank5? i- flank5+ 1: i- flank5);
					writer.write(pos+ "\t"+ sense[i][0]+ "\t"+ sense[i][1]+ "\t"+ sense[i][2]+ "\t"+ sense[i][3]+ "\n");
				}
				writer.flush();
				writer.close();
				writer= new BufferedWriter(new FileWriter(fileBed+"_asense.pwm"));
				for (int i = 0; i < asense.length; i++) {
					int pos= (i>= flank5? i- flank5+ 1: i- flank5);
					writer.write(pos+ "\t"+ asense[i][0]+ "\t"+ asense[i][1]+ "\t"+ asense[i][2]+ "\t"+ asense[i][3]+ "\n");
				}
				writer.flush();
				writer.close();
				System.err.println("wrote "+ fileBed+"_sense.pwm, and "+ fileBed+"_asense.pwm.");
				
			} catch (Exception e) {
				e.printStackTrace();
			}
			
			return ref;
		}
	
	/**
	 * Returns the maximum product of the matrix
	 * @return
	 */
	double maxP= -1;
	public double getMaximumP() {
		if (maxP < 0) {
			maxP= 1d;
			for (int i = 0; i < pwm[0].length; i++) {
				double maxT= -1;
				for (int j = 0; j < pwm.length; j++) {
					if (pwm[j][i]> maxT)
						maxT= pwm[j][i];
				}
				maxP*= maxT;
			}
		}

		return maxP;
	}
	
}
