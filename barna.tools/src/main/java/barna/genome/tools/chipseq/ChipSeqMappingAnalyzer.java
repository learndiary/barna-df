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

package barna.genome.tools.chipseq;

import barna.commons.ByteArrayCharSequence;
import barna.commons.Execute;
import barna.commons.cli.jsap.JSAPParameters;
import barna.commons.launcher.Tool;
import barna.commons.log.Log;
import barna.commons.log.Log.Level;
import barna.flux.capacitor.reconstruction.Kernel;
import barna.io.FileHelper;
import barna.io.bed.BEDMappingIteratorDisk;
import barna.io.bed.BEDReader;
import barna.model.bed.BEDobject2;
import barna.model.commons.IntVector;
import barna.model.rna.UniversalReadDescriptor;
import barna.model.rna.UniversalReadDescriptor.Attributes;
import com.martiansoftware.jsap.JSAPResult;
import com.martiansoftware.jsap.Parameter;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Future;

/**
 * Iterates a mapping file chromosome by chromosome and
 * gets the distance statistics for the mappings.
 * 
 * @author Micha Sammeth (gmicha@gmail.com)
 *
 */
public class ChipSeqMappingAnalyzer implements Tool<int[]> {

	public static void main(String[] args) {
		
		Execute.initialize(2);
		
		File f= new File(			
			"/Users/micha/projects/demassy/download/IP5300109chrall_sorted.bed"
			//"/Users/micha/projects/demassy/download_new/B6+K4me3+200511_sorted_chrY.bed"
		);
				
		ChipSeqMappingAnalyzer myRun= new ChipSeqMappingAnalyzer();
		myRun.fileInput= f;
//		myRun.descriptor= new UniversalReadDescriptor();
//		myRun.descriptor.init(UniversalReadDescriptor.getDescriptor(UniversalReadDescriptor.DESCRIPTORID_PAIRED));
				
		Future<int[]> captain= Execute.getExecutor().submit(myRun);
		int[] distr= null;
		try {
			distr= captain.get();
		} catch (InterruptedException e) {
			e.printStackTrace();
		} catch (ExecutionException e) {
			e.printStackTrace();
		}

		if (distr!= null)
			System.err.println("peak "+ distr[1]+ ", bounds=["+ distr[0]+ ","+ distr[2]+ "]");
		Execute.shutdown();
	}
	
	/**
	 * Rules to parse the read ID.
	 */
	UniversalReadDescriptor descriptor= null;
	
	/**
	 * Default input file.
	 */
	File fileInput;
	
	/**
	 * Default output file.
	 */
	File fileOutput;
	
	/**
	 * The program run parameters.
	 */
	ChipSeqSettings parameters;
	
	int maxInsertSize= 1000;
	
	/**
	 * Empty constructor to comply with Tool implementation.
	 */
	public ChipSeqMappingAnalyzer() {		
	}
	
	/**
	 * Creates an instance based on the given mapping file with a 
	 * specific read descriptor for parsing IDs of mapped reads.
	 * @param inputFile file with mappings
	 * @param descriptor parsing instructions for read IDs
	 */
	public ChipSeqMappingAnalyzer(File inputFile, File outputFile, UniversalReadDescriptor descriptor) {
		this();
		this.fileInput= inputFile;
		this.fileOutput= outputFile;
		this.descriptor= descriptor;
	}
	
	public int[] call() throws Exception {
		
		// copy run parameters
		if (parameters!= null) {
			if (fileInput== null)
				fileInput= parameters.get(ChipSeqSettings.FILE_INPUT);
			if (descriptor== null) {
				descriptor= new UniversalReadDescriptor(parameters.get(ChipSeqSettings.READ_DESCRIPTOR));
			}
			if (fileOutput== null){
				fileOutput= parameters.get(ChipSeqSettings.FILE_OUTPUT);
                if(fileOutput == null || fileOutput.isDirectory()){
                    throw new RuntimeException("You have to specify the FILE_OUTPUT parameter");
                }
            }
		}
		if (descriptor== null) {
			descriptor= new UniversalReadDescriptor(UniversalReadDescriptor.DESCRIPTORID_SIMPLE);
		}
			
		
		// sort
		File sortedInput= new BEDReader(fileInput).getSortedFile(null,    //TODO bad
                ((descriptor != null && descriptor.isPaired()) ? BEDReader.COMPARATOR_PAIRED_END : null));
		System.gc();
		Thread.yield();
		
		// get distance distribution
		int[] a= getDistances(sortedInput, null);
		output2(a);
		
		// find peaks
		int[] bounds= findMainPeak(a);
		output(bounds);
		
		return bounds;
	}

	protected static File getFileOutput(File inputFile) {
		File f= new File(FileHelper.append(inputFile.getAbsolutePath(), 
				"_peaks", true, ".txt"));
		return f;
	}

	private void output(int[] bounds) {

		BufferedWriter writer= null;
		try {
			File f= fileOutput;
			if (f== null)
				f= new File(FileHelper.append(fileInput.getAbsolutePath(), 
						"_peaks", true, ".txt"));
			System.err.println("writing peak descriptions to "+ f.getAbsolutePath());
			writer= new BufferedWriter(new FileWriter(f));
			int mid= (bounds.length- 1)/ 2;
			writer.write(ChipSeqSettings.DISTO_PEAK_MAIN_MAX.getName()+ "\t"+bounds[mid]+ barna.commons.system.OSChecker.NEW_LINE);
			if (mid> 0&& mid< bounds.length- 1) {	// flanks
				writer.write(ChipSeqSettings.DISTO_PEAK_MAIN_LEFT.getName()+ "\t"+ bounds[mid- 1]+ barna.commons.system.OSChecker.NEW_LINE);
				writer.write(ChipSeqSettings.DISTO_PEAK_MAIN_RIGHT.getName()+ "\t"+ bounds[mid+ 1]+ barna.commons.system.OSChecker.NEW_LINE);
			}
			if (mid> 1&& mid< bounds.length- 2) {	// secondary peaks
				writer.write(ChipSeqSettings.DISTO_PEAK_SEC_LEFT.getName()+ "\t"+ bounds[mid- 2]+ barna.commons.system.OSChecker.NEW_LINE);
				writer.write(ChipSeqSettings.DISTO_PEAK_SEC_RIGHT.getName()+ "\t"+ bounds[mid+ 2]+ barna.commons.system.OSChecker.NEW_LINE);
			}
		
		} catch (Exception e) {
			throw new RuntimeException(e);
		} finally {
			try {
				if (writer!= null)
					writer.close();
			} catch (IOException e) {
				throw new RuntimeException(e);
			}
		}
	}

	private void output2(int[] a) {
		if (parameters!= null&& parameters.get(ChipSeqSettings.OUTPUT2)) {
			File f= null;
			if (parameters.get(ChipSeqSettings.FILE_OUTPUT2)== null)
				f= new File(FileHelper.append(fileInput.getAbsolutePath(), 
						"_distogram", true, ".txt"));
			else
				f= parameters.get(ChipSeqSettings.FILE_OUTPUT2);
			
			BufferedWriter writer= null;
			try {
				System.err.println("writing distogram to "+ f.getAbsolutePath());
				writer= new BufferedWriter(new FileWriter(f));
				for (int i = 0; i < a.length; i++) 
					writer.write(Integer.toString(i)+ "\t"+ Integer.toString(a[i])+ barna.commons.system.OSChecker.NEW_LINE);
			} catch (Exception e) {
				throw new RuntimeException(e);
			} finally {
				if (writer!= null)
					try {
						writer.close();
					} catch (IOException e) {
						throw new RuntimeException(e);
					}
			}
		}
	}

	/**
	 * @deprecated doesn't seem promising
	 */
	int[] findMainPeakDerivative(int[] a) {
		
		// smoothen
		float[] b= smoothen(a, 30);
		
		float max= -1;
		int argmax= -1;
		// init
		for (int i = 0; i < b.length; i++) {
			if (b[i]> max) {
				max= b[i];
				argmax= i;
			}
		}

		for (int i = 0; i < b.length- 1; i++) 
			b[i]= b[i+1]- b[i];
		for (int i = 0; i < b.length- 2; i++) 
			b[i]= b[i+1]- b[i];
		
		return null;
	}

	int[] getDistances(File f, int[] a) {
		
		BEDMappingIteratorDisk iter= new BEDMappingIteratorDisk(f);      //
		BEDobject2 bed1= null, bed2= null;
		int ctr= 0;
		boolean pairedEnd= descriptor.isPaired();
		Attributes att1= null, att2= null;
		IntVector multiPairs= null;
		long max= 100000000; // HiSeq lane mappings
		Log.progressStart("Distogram");
		while(iter.hasNext()) {
			++ctr;
			Log.progress(ctr, max);
			
			bed1= new BEDobject2(iter.next());
			if (pairedEnd) {
				att1= descriptor.getAttributes(bed1.getName(), att1);
				if (att1.flag!= 1)
					continue;	// start at mate /1 for paired-end
			} else {
				if (bed1.getStrand()< 0)
					continue;	// only start at + for single mappings
			}
			ByteArrayCharSequence chr1= bed1.getChr();
			iter.mark();
			
			// skip
			while(iter.hasNext()) {	
				bed2= new BEDobject2(iter.next());
				if (pairedEnd) {
					if (!bed2.getName().equals(bed1.getName()))
						break;	// start when other read ID
				} else {	// single read
					if (Math.abs(bed2.getStart()- bed1.getStart())> maxInsertSize
							|| !bed2.getChr().equals(chr1))
						break;	// bad: stop
					if (bed2.getStrand()< 0)
						break;	// good: start when first - mapping
				}
			}
			
			// scan
			if (pairedEnd) {
				if (multiPairs== null)
					multiPairs= new IntVector(10);
				else
					multiPairs.removeAll();
			}
			for(; bed2!= null; bed2= (iter.hasNext()? new BEDobject2(iter.next()): null)) {	// connect
				if (!bed2.getChr().equals(chr1))
					break;
				if (pairedEnd) {	
					att2= descriptor.getAttributes(bed2.getName(), att2);
					if (att2.flag!= 2|| !att1.id.equals(att2.id))
						break;	// stop when not same mate family anymore
				} else {	
					if (bed2.getStrand()> 0)
						continue;	// skip + mappings
				}
				
				int dist= 0;
				if (pairedEnd) {
					if (bed1.getStart()< bed2.getStart()) {
						if (bed1.getStrand()< 0|| bed2.getStrand()> 0)
							continue;	// wrong orientation
						dist= bed2.getEnd()- bed1.getStart();
					} else { 
						if (bed1.getStrand()> 0|| bed2.getStrand()< 0)
							continue;	// wrong orientation
						dist= bed1.getEnd()- bed2.getStart();
					}
					
					if (dist< maxInsertSize)
						multiPairs.add(dist);
				} else {
					dist= bed2.getEnd()- bed1.getStart();
					if (dist> maxInsertSize)
						break;
					a= add(dist, a);
				}
				
			}
			// add best paired-end distance
			if (pairedEnd&& multiPairs.size()> 0) {
				int min= Integer.MAX_VALUE;
				for (int i = 0; i < multiPairs.length; i++) {
					if (multiPairs.get(i)< min)
						min= multiPairs.get(i);
				}
				a= add(min, a);	// take min distance pairing
			}
			
			iter.reset();
		}

		Log.progressFinish();
		return a;
	}

	private int[] add(int dist, int[] a) {
		if (a== null|| dist>= a.length) {
			int[] b= new int[dist+ 1];
			if (a!= null)
				System.arraycopy(a, 0, b, 0, a.length);
			else
				Arrays.fill(b, 0);
			a= b;
		}
		++a[dist];
		return a;
	}

	float[] smoothen(int[] a, int w) {
			
		float[] b= new float[a.length];
		double[] triweight= Kernel.getTriweightKernel(w, 1.09);
		// assume 0-values beyond chr edges to smoothen
		for (int i = 0; i < a.length; i++) {
			double sum= 0;
			for (int j = -w; j < 0; j++)
				sum+= (i+ j< 0? 0: triweight[j+ w]* a[i+ j]);
			for (int j = 0; j <= w; j++) 
				sum+= (i+ j>= a.length? 0: triweight[j+ w]* a[i+ j]);
			
			b[i]= (float) sum;
		}
			
		return b;
	}

	/**
	 * from up to down
	 * @param a
	 */
	int[] findMainPeak(int[] a) {
		
		// smoothen
		float[] b= smoothen(a, 30);
		
		float max= -1;
		int argmax= -1;
		// init
		for (int i = 0; i < b.length; i++) {
			if (b[i]> max) {
				max= b[i];
				argmax= i;
			}
		}
		
		int x1= 0, x2= argmax, x3= argmax, x4= a.length- 1;
		while(x2> x1|| x4> x3) {
			
			// lower bound
			if (x2> x1) {
				float y= b[--x2];
				for (int i = x2- 1; i >= x1; --i) {
					if (b[i]>= y) {
						x1= i;
						break;
					}
				}
			}
			
			// upper bound
			if (x3< x4) {
				float y= b[++x3];
				for (int i = x3+ 1; i <= x4; ++i) {
					if (b[i]>= y) {
						x4= i;
						break;
					}
				}
			}
		}
	
		return new int[] {x2, argmax, x3};
	}

	public void setParameters(File file) {
		try {
			this.parameters= ChipSeqSettings.createSettings(file);
		} catch (Exception e) {
			throw new RuntimeException(e);
		}
	}

    @Override
    public String getName() {
        return "mapdist";
    }

    @Override
    public String getDescription() {
        return "Mapping Distribution";
    }

    @Override
    public String getLongDescription() {
        return null;
    }


    @Override
    public List<Parameter> getParameter() {
        ArrayList<Parameter> parameters = new ArrayList<Parameter>();
        parameters.add(JSAPParameters.flaggedParameter("parameter", 'p').type(File.class).help("specify parameter file").valueName("file").required().get());
        return parameters;

    }

    @Override
    public boolean validateParameter(JSAPResult args) {
        setParameters(args.getFile("parameter"));
		if (parameters== null) {
			System.err.println("I have no parameter file and I want to scream!");
			return false;
		}
		
		if (parameters.get(ChipSeqSettings.FILE_INPUT)== null|| !parameters.get(ChipSeqSettings.FILE_INPUT).canWrite()) {
			System.err.println("Invalid input file "+ parameters.get(ChipSeqSettings.FILE_INPUT).getAbsolutePath());
			return false;
		} else
			fileInput= parameters.get(ChipSeqSettings.FILE_INPUT);
		
		if (parameters.get(ChipSeqSettings.READ_DESCRIPTOR)!= null) {
			descriptor= UniversalReadDescriptor.createTestDescriptor();
			try {
				descriptor.init(parameters.get(ChipSeqSettings.READ_DESCRIPTOR));
			} catch (RuntimeException e) {
				System.err.println("Invalid read descriptor "+ parameters.get(ChipSeqSettings.READ_DESCRIPTOR));
				if (Log.getLogLevel()== Level.DEBUG)
					e.printStackTrace();
				return false;
			}
		}
			
		return true;
	}
	
}
