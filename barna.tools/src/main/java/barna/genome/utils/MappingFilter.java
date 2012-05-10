package barna.genome.utils;/*
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

import barna.commons.ByteArrayCharSequence;
import barna.commons.Execute;
import barna.commons.cli.jsap.JSAPParameters;
import barna.commons.launcher.FluxTool;
import barna.commons.log.Log;
import barna.io.BEDMappingIteratorDisk;
import barna.io.BEDMappingIteratorDisk;
import barna.io.FileHelper;
import barna.io.bed.BEDFileReader;
import barna.io.bed.BEDFileReader;
import barna.io.rna.UniversalReadDescriptor;
import barna.io.rna.UniversalReadDescriptor.Attributes;
import barna.model.bed.BEDobject2;
import com.martiansoftware.jsap.JSAPResult;
import com.martiansoftware.jsap.Parameter;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Future;

/**
 * Class to filter mappings according to custom criteria.
 * @author Micha Sammeth (gmicha@gmail.com)
 *
 */

public class MappingFilter implements FluxTool<Void> {

	public static void main(String[] args) {
		
		Execute.initialize(2);

		MappingFilter myFilter= new MappingFilter();
		myFilter.setInputFile(new File(
				"/Users/micha/projects/demassy/download_new/B6+K4me3+200511_sorted_chrY_testchrx.bed"
				//"/Users/micha/projects/demassy/download_new/transfer/B6+K4me3+200511.sorted.bed"
		));
		myFilter.setDescriptor(UniversalReadDescriptor.DESCRIPTORID_PAIRED);
		myFilter.setInsertMin(60);
		myFilter.setInsertMax(250);
		
		Future<Void> captain= Execute.getExecutor().submit(myFilter);
		try {
			captain.get();
		} catch (InterruptedException e) {
			e.printStackTrace();
		} catch (ExecutionException e) {
			e.printStackTrace();
		}
		Execute.shutdown();
	}
	
	/**
	 * Default suffix for output files
	 */
	public static String SFX_FILTERED= "_filtered";
	
	/**
	 * Creates the default output file.
	 * @param in the input file
	 * @return output file
	 */
	public static File getDefaultOutputFile(File in) {
		String s= FileHelper.append(in.getAbsolutePath(), SFX_FILTERED, 
				false, null);
		return new File(s);
	}
	
	File fileIn;
	File fileOut;
	UniversalReadDescriptor descriptor;
	int insertMin= 170;
	int insertMax= 230;
	
	public Void call() throws Exception {

		// output
		if (fileOut== null)
			fileOut= getDefaultOutputFile(fileIn);
		Log.info("Output in "+ fileOut.getAbsolutePath());
		
		// sort
		File sortedInput= fileIn;
		sortedInput= new BEDFileReader(fileIn).getSortedFile(null,
                ((descriptor != null && descriptor.isPaired()) ? BEDFileReader.COMPARATOR_PAIRED_END : null));

		// doit
		filter(sortedInput, fileOut);
		
		return null;
	}
	
	

	protected void filter(File sortedInput, File fileOut) {
		
		try {
			BEDMappingIteratorDisk iter= new BEDMappingIteratorDisk(sortedInput);
			BufferedWriter writer= new BufferedWriter(new FileWriter(fileOut));

			BEDobject2 bed1, bed2, bedLeft, bedRight; 
			Attributes at1= null, at2= null;
			CharSequence bucketID= null;
			HashSet<ByteArrayCharSequence> map= new HashSet<ByteArrayCharSequence>(10);
			
			Log.progressStart("Filtering");
			long fsize= sortedInput.length();
			long ctr= 0l;
			
			while(iter.hasNext()) {
				
				bed1= new BEDobject2(iter.next());
				
				ctr+= bed1.length()+ 1;
				Log.progress(ctr, fsize);

				at1= descriptor.getAttributes(bed1.getName(), at1);
				if (at1.flag!= 1)
					continue;
				
				boolean written1= false;
				iter.mark();
				
				// output bucket
				if (bucketID== null)
					bucketID= at1.id;
				else if (!at1.id.equals(bucketID)) {
					Iterator<ByteArrayCharSequence> iter2= map.iterator();
					while(iter2.hasNext()) {
						writer.write(iter2.next().toCharArray());
						writer.write("\n");
					}
					writer.flush();
					map.clear();
					bucketID= at1.id;
				}
				
				while (iter.hasNext()) {
					bed2= new BEDobject2(iter.next());
					at2= descriptor.getAttributes(bed2.getName(), at2);
					if (!at1.id.equals(at2.id))
						break;
					
					if (bed1.getStart()< bed2.getStart()) {
						bedLeft= bed1;
						bedRight= bed2;
					} else {
						bedLeft= bed2;
						bedRight= bed1;
					}
					if (!(bedLeft.getStrand()> 0&& bedRight.getStrand()< 0))
						continue;	// orientation check
					
					int dist= bedRight.getEnd()- bedLeft.getStart();					
					if (dist> insertMin&& dist< insertMax) {
						if (!written1) {
							writer.write(bed1.toCharArray());
							writer.write("\n");
							written1= true;
						}
						// don't break for + --> - - constellation
						map.add(bed2);
					}
					
				} // end inner hasNext()

				iter.reset();
				
			} // end outer iter.next()

			// last bucket
			Iterator<ByteArrayCharSequence> iter2= map.iterator();
			while(iter2.hasNext()) {
				writer.write(iter2.next().toCharArray());
				writer.write("\n");
			}
			writer.flush();
			map.clear();

			writer.close();
			
			Log.progressFinish("Done.", true);
			
		} catch (Exception e) {
			throw new RuntimeException(e);
		}
		
	}



    @Override
    public String getName() {
        return "mfilter";
    }

    @Override
    public String getDescription() {
        return "Mapping filter";
    }

    @Override
    public List<Parameter> getParameter() {
        ArrayList<Parameter> parameters = new ArrayList<Parameter>();

        parameters.add(JSAPParameters.flaggedParameter("input", 'i').type(File.class).help("set input mappings file (BED)").valueName("bed").required().get());
        parameters.add(JSAPParameters.flaggedParameter("descriptor", 'd').help("set read descriptor").valueName("descriptor").required().get());
        parameters.add(JSAPParameters.flaggedParameter("insertmin").help("set insert size minimum").type(Integer.class).defaultValue("170").get());
        parameters.add(JSAPParameters.flaggedParameter("insertmax").help("set insert size maximum").type(Integer.class).defaultValue("230").get());

        return parameters;

    }

    @Override
    public boolean validateParameter(JSAPResult args) {
        setInputFile(args.getFile("input"));
        setDescriptor(args.getString("descriptor"));
        setInsertMax(args.getInt("insertmax"));
        setInsertMin(args.getInt("insertmin"));

		if (fileIn== null|| (!fileIn.exists())|| (!fileIn.canRead())
				|| descriptor== null)
			return false;
		return true;
	}

	
    public void setInputFile(File file) {
        this.fileIn= file;
    }


    public void setInsertMin(int imin) {
        this.insertMin= imin;
    }


    public void setInsertMax(int imax) {
        this.insertMax= imax;
    }

    public void setDescriptor(String descriptor) {
    	this.descriptor= new UniversalReadDescriptor();
    	try {
    		this.descriptor.init(UniversalReadDescriptor.getDescriptor(descriptor));
    	} catch (Exception e) {
    		e.printStackTrace();
    		this.descriptor= null;
    	}
    }
	

}
