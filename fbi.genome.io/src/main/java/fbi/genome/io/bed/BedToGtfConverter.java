package fbi.genome.io.bed;

//import gphase.LaVista;

import fbi.commons.tools.ArrayUtils;
import fbi.genome.io.gff.GFFReader;
import fbi.genome.model.bed.BEDobject;
import fbi.genome.model.gff.GFFObject;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.HashMap;
import java.util.Vector;

public class BedToGtfConverter {

	static fbi.genome.model.commons.MyFile inFile, outFile;
	fbi.genome.model.commons.MyFile inputFile, outputFile;
	
	static final String errorMsg = "usage: BedToGtfConverter <inputFile> [outputFileName]\n\n" + "where\n" + "<inputFile>\ta BED file\n" + "[outpuFileName] an optional name for the output GTF (default is <inputFileName>.gtf" + "\n\nmicha, may 07";
	
	static void parseArguments(String[] args) {
		if (args== null|| args.length< 1) {
			System.err.println(errorMsg);
			System.exit(-1);
		}
		
		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-out")) {
				outFile= new fbi.genome.model.commons.MyFile(args[++i]);
				continue;
			}
			inFile= new fbi.genome.model.commons.MyFile(args[i]);
		}
		
		if (!inFile.exists()) {
			System.err.println(errorMsg);
			System.exit(-1);
		}
		if (outFile== null) {
			outFile= new fbi.genome.model.commons.MyFile(
					inFile.getPathOnly()+ File.separator+
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
	
	public BedToGtfConverter(fbi.genome.model.commons.MyFile inFile, fbi.genome.model.commons.MyFile outFile) {
		this.inputFile= inFile;
		this.outputFile= outFile;
	}
	
	public void convert() {
		try {
			System.out.println("Reading..");
			BEDwrapper reader= new BEDwrapper(inputFile.getAbsolutePath());
			reader.read(0);
			BEDobject[] beds= reader.getBeds();
			Vector gtfV= new Vector(beds.length);
			System.out.println("Converting..");
			String src= inputFile.getFileNameOnly();
			HashMap<String,Object> map= new HashMap<String,Object>(beds.length,1f);
			new File(outputFile.getAbsolutePathWithoutExtension()+"_error.txt").delete();
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
			GFFReader writer= new GFFReader(outputFile.getAbsolutePath());
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
				BufferedWriter errorWriter= new BufferedWriter(new FileWriter(outputFile.getAbsolutePathWithoutExtension()+"_error.txt", true));
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
				BufferedWriter errorWriter= new BufferedWriter(new FileWriter(outputFile.getAbsolutePathWithoutExtension()+"_error.txt", true));
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
