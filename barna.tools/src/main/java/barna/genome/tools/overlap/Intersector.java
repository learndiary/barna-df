/*
 * Copyright (c) 2010, Micha Sammeth
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * The names of its contributors may be not used to endorse or promote
 *       products derived from this software without specific prior written
 *       permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL MICHA SAMMETH BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

package barna.genome.tools.overlap;

import barna.commons.ByteArrayCharSequence;
import barna.commons.Execute;
import barna.commons.cli.jsap.JSAPParameters;
import barna.commons.launcher.FluxTool;
import barna.io.BufferedBACSReader;
import barna.io.FileHelper;
import barna.model.bed.BEDobject2;
import com.martiansoftware.jsap.JSAPResult;
import com.martiansoftware.jsap.Parameter;
import org.jfree.util.Log;

import java.io.*;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Future;
import java.util.zip.GZIPInputStream;


public class Intersector implements FluxTool<Void> {

	public static byte MODE_INTERSECT= 1;
	public static byte MODE_OVERLAP= 2;
	
	public static void main(String[] args) {
		
		File parFile= new File("/Users/micha/projects/demassy/download_new/cisgenome/inter_over/overlap_may_b50_w10_c3_june_b50_w10_c2_pas");

/*		
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
*/		
		File outFile= new File("/Users/micha/projects/demassy/download_new/cisgenome/inter_over/overlap_may_b50_w10_c3_june_b50_w10_c2_pas_isect.bed");
		
		Execute.initialize(2);

		Intersector myIsector= new Intersector();
		myIsector.parFile= parFile;
		myIsector.outFile= outFile;
		myIsector.mode= MODE_INTERSECT;
		myIsector.ctr= 1;

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

	byte mode= -1;
	
    public void setParameters(File file) {
        this.parFile = file;
    }
	
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
		
		BEDobject2[] beds= new BEDobject2[2];
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

//				if (two.getName().contains("RJ2_220611_input")) { 
//					System.currentTimeMillis();
//					if (bed2!= null&& 
//						((bed2.getChr().equals("chr1")
//						&& bed2.getStart()>= 186430000
//						&& bed2.getEnd()<=186450000))
//						|| (!bed2.getChr().equals("chr1")))
//					System.currentTimeMillis();
//				}
				
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
						} else {
							beds[0]= bed1;
							beds[1]= bed2;
							if (mode== MODE_INTERSECT)
								bedTmp= intersect(beds, start1, end1, oneCode, 
										start2, end2, twoCode, writer);
							else if (mode== MODE_OVERLAP)
								bedTmp= overlap(beds, start1, end1, oneCode, 
										start2, end2, twoCode, writer);
							bed1= beds[0];
							bed2= beds[1];
							
							// consume completely joint,
							// set other to rest
							if (end1< end2) {
								bed1= null;
								bed2= bedTmp;
							} if (end2< end1) {
								bed2= null;
								bed1= bedTmp;
							} else {	// Lemma, in this case there exists no 3rd block
								bed1= null;
								bed2= null;
							}

						}
						
					} // end: on same chromosome
				} // end: overlapping objects

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
	
	private String toBinaryString(int d) {
		StringBuilder b= new StringBuilder(Integer.toBinaryString(d));
		for (int i = b.length(); i < parFileNr; i++) 
			b.insert(0, '0');
		b.reverse();
		return b.toString();
	}
	
	/**
	 * Creates and writes bed object, generated from one or two parents. 
	 * @param start 1st coordinate of the bed line to be created
	 * @param end 2nd coordinate of the bed line to be created
	 * @param bed1 1st parent object
	 * @param oneCode 1st parent's code
	 * @param bed2 2nd parent object
	 * @param twoCode 2nd parent's code
	 * @param writer the writer to which new bed object is written, or
	 * <code>null</code> if no output is desired
	 * @return the created bed object
	 */
	private BEDobject2 write(int start, int end, BEDobject2 bed1, int oneCode,
			BEDobject2 bed2, int twoCode, BufferedWriter writer) {
		
		BEDobject2 bed= new BEDobject2();
		bed.setChromosome(bed1.getChr());
		bed.setStart(start);
		bed.setEnd(end);
//		int code1= oneCode> 0? oneCode: Integer.parseInt(bed1.getName().toString());
		int code1= oneCode;
		try {
			String s= bed1.getName().toString();
			int p= s.indexOf(':');
			if (p> 0&& s.indexOf(':', p+1)< 0)
				code1= Integer.parseInt(s.substring(0, p));
		} catch (Exception e) {
			; // :)
		}
		int code2= 0;
		if (bed2!= null) {
//			code2= twoCode> 0? twoCode: Integer.parseInt(bed2.getName().toString());
			code2= twoCode;
			try {
				String s= bed2.getName().toString();
				int p= s.indexOf(':');
				if (p> 0&& s.indexOf(':', p+1)< 0)
					code2= Integer.parseInt(s.substring(0, p));
			} catch (Exception e) {
				; //:)
			}
				

		}
		
		// bin.join
		int code= (code1| code2);
		
		if (code< 0|| code> 256)
			System.currentTimeMillis();
		bed.setName(Integer.toString(code)+ ":"+ toBinaryString(code));
		int len= end- start;
		int score= (int) (bed1.getScore()* len/ (float) bed1.getLength());
		if (bed2!= null) {
			score+= (bed2.getScore()* len/ (float) bed2.getLength());
			score/= 2;	// arit. average of rel. weight for both overlapping parts
		}

		bed.setScore(score> 1000? 1000: score);
		bed.setStrand((byte) 1);
		
		if (writer!= null)
			write(bed, -1, writer);
		
		return bed;
	}

	private char[] charBuf= new char[160];
	private ByteArrayCharSequence ctrBuf= new ByteArrayCharSequence(10);
	/**
	 * Write bed object to the given writer
	 * @param bed the bed object
	 * @param ocode code of the bed object, or <code>&leq;1</code> if not specified
	 * @param writer the writer
	 */
	private void write(BEDobject2 bed, int ocode, BufferedWriter writer) {

		// try to preserve suitable code in input
		int code= ocode;
		try {
			String s= bed.getName().toString();
			int p= s.indexOf(':');
			if (p> 0&& s.indexOf(':', p+1)< 0)
				code= Integer.parseInt(s.substring(0, p));
		} catch (Exception e) {
			; // :)
		}
		if (code> 0)
			bed.setName(Integer.toString(code)+ ":"+ toBinaryString(ocode));	// TODO make BACS to work with single chars
		else
			System.currentTimeMillis();
		
		String s= bed.getName().toString();
		int p= s.indexOf(':');
		int ccode= Integer.parseInt(s.substring(0, p));
		if (ccode< 0|| ccode> 256)
			System.currentTimeMillis();

		// write
		charBuf= bed.toCharArray(charBuf);		
		try {
			writer.write(charBuf, 0, bed.length());
			writer.write(barna.commons.system.OSChecker.NEW_LINE);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	/**
	 * Parameter file with a list of files that are to be intersect.
	 */
	protected File parFile;
	
	/**
	 * Number of file entries in the parameter file.
	 */
	protected int parFileNr= 0;
	
	/**
	 * Output file to which the intersected bed is written.
	 */
	protected File outFile;
	
	/**
	 * Counter to start with
	 */
	int ctr= 1;
	
	/**
	 * Iterates files in the list provided, performing pairwise 
	 * overlaps.
	 */
	public Void call() throws Exception {

		BufferedReader buffy= null;
        long ll = FileHelper.countLines(parFile);
        if(ll > Integer.MAX_VALUE) throw new RuntimeException(ll + " value > Integer.MAX_VALUE");
        parFileNr= (int) ll;
		try {
			buffy= new BufferedReader(new FileReader(parFile));
			
			String s= buffy.readLine();
			File f0= new File(s);
			s= buffy.readLine();
			File f= new File(s);
			File tmpFile= FileHelper.createTempFile(this.getClass().getSimpleName(), ".bed");
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


    @Override
    public String getName() {
        return "isect";
    }

    @Override
    public String getDescription() {
        return "Intersector";
    }

    @Override
    public List<Parameter> getParameter() {
        ArrayList<Parameter> parameters = new ArrayList<Parameter>();
        parameters.add(JSAPParameters.flaggedParameter("parameter", 'p').type(File.class).help("specify parameter file").valueName("file").required().get());
        parameters.add(JSAPParameters.flaggedParameter("output", 'o').type(File.class).help("specify output file").valueName("output").required().get());
        return parameters;

    }

    @Override
    public boolean validateParameter(JSAPResult args) {
        setParameters(args.getFile("parameter"));
        setOutput(args.getFile("output"));
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

	protected BEDobject2 intersect(BEDobject2[] beds, int start1, int end1, int oneCode, 
			int start2, int end2, int twoCode, BufferedWriter writer) {
		
		// overlap
		int minStart= (start1< start2? start1: start2);
		int maxEnd= (end1< end2? end2: end1);

		// obs: there are max. 3 blocks produced, 
		// less iff start1== start2, or end1== end2

		// all block boundaries
		int startI1= minStart, 
			endI1= (start1== minStart? start2: start1),
			startI2= endI1,
			endI2= (end1== maxEnd? end2: end1),
			startI3= endI2,
			endI3= (end1== maxEnd? end1: end2);
		
		// block1 (minStart,maxStart)
		BEDobject2 bedTmp= null;
		if (startI1!= endI1) 
			bedTmp= write(startI1, endI1, 
					startI1== start1? beds[0]: beds[1], 
					startI1== start1? oneCode: twoCode, 
					null, -1, writer);	// block before intersection always written
		// block2 (maxStart,minEnd)
		if (startI2!= endI2)
			bedTmp= write(startI2, endI2, beds[0], oneCode, beds[1], twoCode, writer);	// intersection block always written
		// block3 (minEnd,maxEnd)
		if (startI3!= endI3)
			bedTmp= write(startI3, endI3, 
					maxEnd== end1? beds[0]: beds[1], 
					maxEnd== end1? oneCode: twoCode, 
					null, twoCode, null);	// never write out last block

		return bedTmp;
	}

	protected BEDobject2 overlap(BEDobject2[] beds, int start1, int end1, int oneCode, 
			int start2, int end2, int twoCode, BufferedWriter writer) {
		
		// overlap
		int minStart= (start1< start2? start1: start2);
		int maxEnd= (end1< end2? end2: end1);

		BEDobject2 bedTmp= write(minStart, maxEnd, 
				beds[0], oneCode, beds[1], twoCode, null);
	
		// Lemma, if both end at same position, 
		// their join is not needed for further overlaps
		if (end1== end2) 
			write(bedTmp, -1, writer);
		
		
		return bedTmp;
	}

	
}
