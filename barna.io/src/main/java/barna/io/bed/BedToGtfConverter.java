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

package barna.io.bed;

//import gphase.LaVista;

import barna.commons.utils.ArrayUtils;
import barna.io.FileHelper;
import barna.io.gtf.GTFwrapper;
import barna.model.bed.BEDobject;
import barna.model.gff.GFFObject;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.HashMap;
import java.util.Vector;

public class BedToGtfConverter {

	static barna.model.commons.MyFile inFile, outFile;
	barna.model.commons.MyFile inputFile, outputFile;
	
	static final String errorMsg = "usage: BedToGtfConverter <inputFile> [outputFileName]\n\n" + "where\n" + "<inputFile>\ta BED file\n" + "[outpuFileName] an optional name for the output GTF (default is <inputFileName>.gtf" + "\n\nmicha, may 07";
	
	static void parseArguments(String[] args) {
		if (args== null|| args.length< 1) {
			System.err.println(errorMsg);
			System.exit(-1);
		}
		
		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-out")) {
				outFile= new barna.model.commons.MyFile(args[++i]);
				continue;
			}
			inFile= new barna.model.commons.MyFile(args[i]);
		}
		
		if (!inFile.exists()) {
			System.err.println(errorMsg);
			System.exit(-1);
		}
		if (outFile== null) {
			outFile= new barna.model.commons.MyFile(
					FileHelper.getPathOnly(inFile)+ File.separator+
					inFile.getFileNameOnly()+ ".gtf");
		}
	}
	
	public static void main(String[] args) {
		long t0= System.currentTimeMillis();
		parseArguments(args);
		BedToGtfConverter conv= new BedToGtfConverter(inFile, outFile);
		conv.convert();
		System.out.println("took "+(System.currentTimeMillis()- t0)/1000+" sec.");
	}
	
	public BedToGtfConverter(barna.model.commons.MyFile inFile, barna.model.commons.MyFile outFile) {
		this.inputFile= inFile;
		this.outputFile= outFile;
	}
	
	public void convert() {
		try {
			System.out.println("Reading..");
			BEDFileReader reader= new BEDFileReader(inputFile.getAbsolutePath());
			reader.read(0);
			BEDobject[] beds= reader.getBeds();
			Vector gtfV= new Vector(beds.length);
			System.out.println("Converting..");
			String src= inputFile.getFileNameOnly();
			HashMap<String,Object> map= new HashMap<String,Object>(beds.length,1f);
			new File(FileHelper.getAbsolutePathWithoutExtension(outputFile)+"_error.txt").delete();
			for (int i = 0; i < beds.length; i++) {
				
				GFFObject[] gtfs= GFFObject.fromBed(beds[i]);
				for (int j = 0; gtfs!= null&&  j < gtfs.length; j++) {
					gtfs[j].setSeqname(beds[i].getChrom().toString());
					gtfs[j].setSource(src);					
					if (initIDs(gtfs[j],beds[i],map))
						gtfV.add(gtfs[j]);
					
				}
				map.put(beds[i].getName().toString(),beds[i].getName());
			}
			
			System.out.println("Writing..");
			GFFObject[] gtfs= (GFFObject[]) ArrayUtils.toField(gtfV);
			GTFwrapper writer= new GTFwrapper(outputFile.getAbsolutePath());
			writer.setChromosomeWise(false);
			// TODO rewrite
//			writer.setGtfObj(gtfs);
//			writer.write();
			
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		
	}

	private boolean initIDs(GFFObject gtf, BEDobject bed, HashMap<String,Object> map) {
		String evnt_tid= bed.getName().toString();
		if (evnt_tid.contains("-")) {	// TIDs are awked in..
			String[] tokens= evnt_tid.split("-");
			gtf.addAttribute(GFFObject.TRANSCRIPT_ID_TAG, tokens[0]);
			evnt_tid= tokens[1];
		}
		if (map.get(evnt_tid)!= null) {
			try {
				BufferedWriter errorWriter= new BufferedWriter(new FileWriter(FileHelper.getAbsolutePathWithoutExtension(outputFile)+"_error.txt", true));
				errorWriter.write(evnt_tid+"\tdouble entry\n");
				errorWriter.flush();
				errorWriter.close();
			} catch (Exception e) {
				e.printStackTrace();
			}
			return false;
		} 

		String[] ev_exon= evnt_tid.split("_");
		if (ev_exon.length!= 2) {
			try {
				BufferedWriter errorWriter= new BufferedWriter(new FileWriter(FileHelper.getAbsolutePathWithoutExtension(outputFile)+"_error.txt", true));
				errorWriter.write(evnt_tid+"\twrong format <> 2 '_' delimited fields\n");
				errorWriter.flush();
				errorWriter.close();
			} catch (Exception e) {
				e.printStackTrace();
			}
			return false;
		}
		String gid= ev_exon[0];
		for (int k = 1; k < ev_exon.length- 1; k++) 
			gid+="_"+ ev_exon[k];			
		
		// TODO check whether still needed
//		gtf.addAttribute(LaVista.GTF_ATTRIBUTE_GROUP_ID, gid);	// event					
		
		if (ev_exon[1].equals("all"))
			gtf.setFeature("event");
		else
			gtf.setFeature("probe_"+ev_exon[ev_exon.length-1]);
		
		return true;

	}
}
