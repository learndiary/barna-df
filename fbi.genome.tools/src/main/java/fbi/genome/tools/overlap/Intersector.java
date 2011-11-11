package fbi.genome.tools.overlap;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Future;
import java.util.zip.GZIPInputStream;

import org.cyclopsgroup.jcli.ArgumentProcessor;
import org.cyclopsgroup.jcli.annotation.Cli;
import org.cyclopsgroup.jcli.annotation.Option;
import org.jfree.util.Log;

import fbi.commons.ByteArrayCharSequence;
import fbi.commons.Execute;
import fbi.commons.flux.FluxTool;
import fbi.commons.options.HelpPrinter;
import fbi.genome.io.BufferedBACSReader;
import fbi.genome.io.FileHelper;
import fbi.genome.model.bed.BEDobject2;

@Cli(name = "isect", description = "Intersector")
public class Intersector implements FluxTool<Void>{

	public static void main(String[] args) {
		
		File parFile= null;
		
		try {
			parFile= new File("/Users/micha/projects/demassy/download_new/cisgenome/bed/overlap.par"); 
				
//				FileHelper.createTempFile(
//					Intersector.class.getSimpleName(), 
//					"tst_overlap.par");
			BufferedWriter writer= new BufferedWriter(new FileWriter(parFile));
//			writer.write("/Users/micha/projects/demassy/download_new/tst2.txt.peak_peak.bed_clean\n");
//			writer.write("/Users/micha/projects/demassy/download_new/tst3.txt.peak_peak.bed_clean\n");
			String dir= "/Users/micha/projects/demassy/download_new/cisgenome/bed/";
			writer.write(dir+ "B6_200511_input.txt.peak_sorted.bed\n");
			writer.write(dir+ "B6_230611_input.txt.peak_sorted.bed\n");
			writer.write(dir+ "RJ2_190511_input.txt.peak_sorted.bed\n");
			writer.write(dir+ "RJ2_220611_input.txt.peak_sorted.bed\n");
			writer.close();
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(-1);
		}
		
		File outFile= new File("/Users/micha/projects/demassy/download_new/cisgenome/bed/overlap.bed");
		
		Execute.initialize(2);

		Intersector myIsector= new Intersector();
		myIsector.parFile= parFile;
		myIsector.outFile= outFile;

		Future<Void> captain= Execute.getExecutor().submit(myIsector);
		try {
			captain.get();
		} catch (InterruptedException e) {
			e.printStackTrace();
		} catch (ExecutionException e) {
			e.printStackTrace();
		}
		Execute.shutdown();

	}

    @Option(name = "p", longName = "parameter", description = "specify parameter file", displayName = "file", required = true)
    public void setParameters(File file) {
        this.parFile = file;
    }
	
    @Option(name = "o", longName = "output", description = "specify output file", displayName = "file", required = true)
    public void setOutput(File file) {
        this.outFile = file;
    }
	
	/**
	 * Pairwisely overlaps blocks in a file
	 * @param one the first file
	 * @param oneCode code of the first file
	 * @param two the second file
	 * @param twoCode code of the second file
	 * @param out the output file
	 */
	public void intersect(File one, int oneCode, File two, int twoCode, File out) {
		
		BufferedBACSReader buffy1= null, buffy2= null;
		BufferedWriter writer= null;
		try {
			ByteArrayCharSequence buf1= new ByteArrayCharSequence(160);
			ByteArrayCharSequence buf2= new ByteArrayCharSequence(160);
			if (one.getName().endsWith(".gz")) {
				buffy1= new BufferedBACSReader(new GZIPInputStream(new FileInputStream(one)));
				buffy2= new BufferedBACSReader(new GZIPInputStream(new FileInputStream(two)));
			} else {
				buffy1= new BufferedBACSReader(new FileInputStream(one));
				buffy2= new BufferedBACSReader(new FileInputStream(two));
			}
			writer= new BufferedWriter(new FileWriter(out));
			
			// TODO assume sorted files
			BEDobject2 bed1= null, bed2= null, bedTmp= null;
			while(buf1== null|| buf1.length()== 0|| buf1.startsWith("track"))
				buf1= buffy1.readLine(buf1);
			while(buf2== null|| buf2.length()== 0|| buf2.startsWith("track"))
				buf2= buffy2.readLine(buf2);
			if (buf1!= null)
				bed1= new BEDobject2(buf1);
			if (buf2!= null)
				bed2= new BEDobject2(buf2);
			while (bed1!= null|| bed2!= null) {

				if (bed1== null) {
					write(bed2, twoCode, writer);
					bed2= null;
				} else if (bed2== null) {
					write(bed1, oneCode, writer);
					bed1= null;
				} else {
					ByteArrayCharSequence chr1= bed1.getChr();
					ByteArrayCharSequence chr2= bed2.getChr();
					int cc= chr1.compareTo(chr2);
					if (cc> 0) {
						write(bed2, twoCode, writer);
						bed2= null;
					} else if (cc< 0) {
						write(bed1, oneCode, writer);
						bed1= null;
					} else {
						
						int start1= bed1.getStart();
						int end1= bed1.getEnd();
						int start2= bed2.getStart();
						int end2= bed2.getEnd();

						// one before the other
						if (end1<= start2) {
							write(bed1, oneCode, writer);
							bed1= null;
						} else if (end2<= start1) {
							write(bed2, twoCode, writer);
							bed2= null;
						} else {	// intersect
							// overlap
							int minStart= (start1< start2? start1: start2);
							int maxEnd= (end1< end2? end2: end1);

							// obs: there are max. 3 blocks produced, 
							// less iff start1== start2, or end1== end2
							
							// block1 (minStart,maxStart)
							int start= minStart;
							int end= (start1== minStart? start2: start1);
							if (start!= end) 
								bedTmp= write(start, end++, 
										start== start1? bed1: bed2, 
										start== start1? oneCode: twoCode, 
										null, -1, writer);
							// block2 (maxStart,minEnd)
							start= end;
							end= (end1== maxEnd? end2: end1);
							if (start!= end)
								bedTmp= write(start, end++, bed1, oneCode, bed2, twoCode, writer);
							// block3 (minEnd,maxEnd)
							start= end;
							end= (end1== maxEnd? end1: end2);
							if (start!= end)
								bedTmp= write(start, end++, 
										maxEnd== end1? bed1: bed2, 
										maxEnd== end1? oneCode: twoCode, 
										null, twoCode, writer);

							// consume completely divided,
							// set other to rest
							if (end1< end2) {
								bed1= null;
								bed2= bedTmp;
							} if (end2< end1) {
								bed2= null;
								bed1= bedTmp;
							} else {
								bed1= null;
								bed2= null;
							}
							
						}
						
					}
				}

				// read
				if (bed1== null) {
					buf1= buffy1.readLine(buf1);
					if (buf1!= null)
						bed1= new BEDobject2(buf1);
				}
				if (bed2== null) {
					buf2= buffy2.readLine(buf2);
					if (buf2!= null)
						bed2= new BEDobject2(buf2);
				}

			}
			
			
		} catch (Exception e) {
			e.printStackTrace();
		} finally {
			if (buffy1!= null)
				try {
					buffy1.close();
				} catch (IOException e) {
					e.printStackTrace();
				}	
			if (buffy2!= null)
				try {
					buffy2.close();
				} catch (IOException e) {
					e.printStackTrace();
				}	
			if (writer!= null)
				try {
					writer.close();
				} catch (IOException e) {
					e.printStackTrace();
				}	
		}
	}

	private BEDobject2 write(int start, int end, BEDobject2 bed1, int oneCode,
			BEDobject2 bed2, int twoCode, BufferedWriter writer) {
		
		BEDobject2 bed= new BEDobject2();
		bed.setChromosome(bed1.getChr());
		bed.setStart(start);
		bed.setEnd(end);
		int code1= oneCode> 0? oneCode: Integer.parseInt(bed1.getName().toString());
		int code2= 0;
		if (bed2!= null)
			code2= twoCode> 0? twoCode: Integer.parseInt(bed2.getName().toString());
		int code= code1+ code2;
		bed.setName(Integer.toString(code));
		int len= end- start;
		int score= (int) (bed1.getScore()* len/ (float) bed1.getLength());
		if (bed2!= null) {
			score+= (bed2.getScore()* len/ (float) bed2.getLength());
			score/= 2;	// arit. average of rel. weight for both overlapping parts
		}

		bed.setScore(score> 1000? 1000: score);
		bed.setStrand((byte) 1);
		
		write(bed, -1, writer);
		
		return bed;
	}

	private char[] charBuf= new char[160];
	private ByteArrayCharSequence ctrBuf= new ByteArrayCharSequence(10);
	private void write(BEDobject2 bed, int ocode, BufferedWriter writer) {
		
		if (ocode> 0)
			bed.setName(Integer.toString(ocode));	// TODO make BACS to work with single chars
		charBuf= bed.toCharArray(charBuf);		
		try {
			writer.write(charBuf, 0, bed.length());
			writer.write("\n");
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	/**
	 * Parameter file with a list of files that are to be intersect.
	 */
	protected File parFile;
	
	/**
	 * Output file to which the intersected bed is written.
	 */
	protected File outFile;
	
	/**
	 * Iterates files in the list provided, performing pairwise 
	 * overlaps.
	 */
	public Void call() throws Exception {

		BufferedReader buffy= null;
		try {
			buffy= new BufferedReader(new FileReader(parFile));
			
			String s= buffy.readLine();
			File f0= new File(s);
			s= buffy.readLine();
			File f= new File(s);
			File tmpFile= FileHelper.createTempFile(this.getClass().getSimpleName(), ".bed");
			int ctr= 1;
			intersect(f0, ctr, f, (ctr*2), tmpFile);
			
			for (s= null; (s= buffy.readLine())!= null;) {
				
				f0= tmpFile;
				f= new File(s);
				tmpFile= FileHelper.createTempFile(this.getClass().getSimpleName(), ".bed");
				ctr*= 2;
				
				intersect(f0, -1, f, (ctr* 2), tmpFile);

				f0.delete();
			}
			
			FileHelper.move(tmpFile, outFile);
			
		} catch (Exception e) {
			Log.error(e);
		} finally {
			if (buffy!= null)
				buffy.close();
		}
		
		return null;
	}

	/**
	 * Checks that file list, and each file in the list exists.
	 */
	public boolean validateParameters(HelpPrinter printer,
			ArgumentProcessor toolArguments) {
		
		if (parFile== null|| (!parFile.exists())|| (!parFile.canRead()))
			return false;
		BufferedReader buffy= null;
		try {
			buffy= new BufferedReader(new FileReader(parFile));
			int ctr= 0;
			for (String s= null; (s= buffy.readLine())!= null; ++ctr) {
				File f= new File(s);
				if ((!f.exists())|| (!f.canRead())) {
					Log.error("File inaccessible "+ f.getAbsolutePath());
					return false;
				}
			}
			if (ctr< 2)
				Log.error("Cannot intersect less than two files, you provided only "+ ctr);
			return true;
		} catch (Exception e) {
			Log.error(e);
			return false;
		} finally {
			if (buffy!= null)
				try {
					buffy.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
		}
	}

	
}
