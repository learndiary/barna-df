/*
 * This file is part of the Flux Library.
 *
 * The code of the Flux Library may be freely distributed and modified under the terms of the
 * European Union Public Licence (EUPL) published on the web site <http://www.osor.eu/eupl/european-union-public-licence-eupl-v.1.1>.
 * Copyright for the code is held jointly by the individual authors, who should be listed
 * in @author doc comments. According to Article 5 and Article 11 of the EUPL, publications that
 * include results produced by the Flux Library are liable to reference the Work,
 * see the Flux Library homepage <http://flux.sammeth.net> for more information.
 */

package fbi.genome.sequencing.rnaseq.reconstruction;

import fbi.commons.Log;
import fbi.genome.io.FileHelper;
import fbi.genome.io.rna.UniversalReadDescriptor;
import fbi.genome.model.constants.Constants;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;

/**
 * Flux Capacitor settings
 */
public class FluxCapacitorParameters {

	public static final String PAR_YES= "YES"; 
	public static final String PAR_NO= "NO"; 
	public static final String PAR_ANNOTATION_FILE= "ANNOTATION_FILE"; 
	public static final String PAR_SORTED_ANNOTATION_FILE= "SORTED_ANNOTATION_FILE";
	public static final String PAR_MAPPING_FILE= "MAPPING_FILE";
	public static final String PAR_COVERAGE_FILE= "COVERAGE_FILE";
	public static final String PAR_SORTED_MAPPINGS_FILE= "SORTED_MAPPINGS_FILE";
	public static final String PAR_STDOUT_FILE= "STDOUT_FILE";
	public static final String PAR_FORCE= "FORCE";	
	public static final String PAR_STDERR_FILE= "STDERR_FILE"; 
	public static final String PAR_TMP_DIR= "TMP_DIR";	
	public static final String PAR_COPY_INPUT= "COPY_INPUT";	
	public static final String PAR_READ_NUMBER= "READ_NUMBER"; 
	public static final String PAR_READ_DESCRIPTOR= "READ_DESCRIPTOR";	 
	public static final String PAR_ANNOTATION_MAPPING= "ANNOTATION_MAPPING"; 
	public static final String ANNOTATION_MAPPING_PAIRED= "PAIRED", 
								ANNOTATION_MAPPING_STRANDED= "STRANDED", 
								ANNOTATION_MAPPING_SINGLE= "SINGLE",
								ANNOTATION_MAPPING_COMBINED= "COMBINED";
	public static final String PAR_PROFILE_FILE= "PROFILE_FILE";
	public static final String PAR_MAPPED_READ_FILE= "MAPPED_READ_FILE";
	public static final String PAR_NOT_MAPPED_READ_FILE= "NOT_MAPPED_READ_FILE";
	public static final String PAR_INSERT_SIZE_FILE= "INSERT_SIZE_FILE";
	public static final String PAR_LP_ZIP_FILE= "LP_ZIP_FILE";
	public static final String PAR_SORT_IN_RAM= "SORT_IN_RAM";
	
	public boolean check() {
		
		if (pairedEnd&& !descriptor.isPaired()) { 
			if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
				System.err.println("[DESCRIPTOR] For paired-end quantification, a read descriptor describing the pairs is required.");
			return false;
		}
		if (stranded&& !descriptor.isStranded()) { 
			if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
				System.err.println("[DESCRIPTOR] For stranded quantification, a read descriptor describing sense/antisense annotation is required.");
			return false;
		} 


		return true;
	}
	
	public static FluxCapacitorParameters create(File parFile) {
		try {
			FluxCapacitorParameters pars= new FluxCapacitorParameters();
			BufferedReader buffy= new BufferedReader(new FileReader(parFile));
			for (String line= null; (line= buffy.readLine())!= null;) {
				if (line.trim().length()== 0)
					continue;
				String[] ll= line.split("\\s");
				if (ll[0].equalsIgnoreCase(PAR_ANNOTATION_FILE)) {
					pars.fileAnnotation= new File(ll[1]);
					if (!pars.fileAnnotation.exists()) {
						notexist(ll[1]);
						buffy.close();
						return null;
					}
				} else if (ll[0].equalsIgnoreCase(PAR_MAPPING_FILE)) {
					pars.fileMappings= new File(ll[1]);
					if (!pars.fileMappings.exists()) {
						notexist(ll[1]);
						buffy.close();
						return null;
					}
				} else if (ll[0].equalsIgnoreCase(PAR_STDERR_FILE)) {
					pars.fileStdErr= new File(ll[1]);
					if (!pars.fileStdErr.exists()) {
						notexist(ll[1]);
						buffy.close();
						return null;
					}
				} else if (ll[0].equalsIgnoreCase(PAR_STDOUT_FILE)) {
					pars.fileStdOut= new File(ll[1]);
					File parent= pars.fileStdOut.getParentFile();
					if (!parent.exists()) {
						notexist(parent.getAbsolutePath());
						buffy.close();
						return null;
					}
				} else if (ll[0].equalsIgnoreCase(PAR_SORTED_MAPPINGS_FILE)) {
					pars.fileMappingsSorted= new File(ll[1]);
					File parent= pars.fileMappingsSorted.getParentFile();
					if (!parent.exists()) {
						notexist(parent.getAbsolutePath());
						buffy.close();
						return null;
					}
				} else if (ll[0].equalsIgnoreCase(PAR_COVERAGE_FILE)) {
					pars.fileCoverage= new File(ll[1]);
					File parent= pars.fileCoverage.getParentFile();
					if (!parent.exists()) {
						notexist(parent.getAbsolutePath());
						buffy.close();
						return null;
					}
				} else if (ll[0].equalsIgnoreCase(PAR_SORTED_ANNOTATION_FILE)) {
					pars.fileAnnotationSorted= new File(ll[1]);
					File parent= pars.fileAnnotationSorted.getParentFile();
					if (!parent.exists()) {
						notexist(parent.getAbsolutePath());
						buffy.close();
						return null;
					}
				} else if (ll[0].equalsIgnoreCase(PAR_PROFILE_FILE)) {
					pars.fileProfile= new File(ll[1]);	// is created if non-existent
				} else if (ll[0].equalsIgnoreCase(PAR_MAPPED_READ_FILE)) {
					pars.fileMapped= new File(ll[1]);
					if (!pars.fileMapped.exists()) {
						notexist(ll[1]);
						buffy.close();
						return null;
					}
				} else if (ll[0].equalsIgnoreCase(PAR_NOT_MAPPED_READ_FILE)) {
					pars.fileNotmapped= new File(ll[1]);
					if (!pars.fileNotmapped.exists()) {
						notexist(ll[1]);
						buffy.close();
						return null;
					}
				} else if (ll[0].equalsIgnoreCase(PAR_INSERT_SIZE_FILE)) {
					pars.fileInsert= new File(ll[1]);
					if (!pars.fileInsert.exists()) {
						notexist(ll[1]);
						buffy.close();
						return null;
					}
				} else if (ll[0].equalsIgnoreCase(PAR_LP_ZIP_FILE)) {
					if (PAR_LP_ZIP_FILE.endsWith(FileHelper.SFX_ZIP))
						pars.fileLPzip= new File(ll[1]);
					else
						pars.fileLPzip= new File(ll[1]+ ".zip");
					if (!pars.fileLPzip.getParentFile().exists()) {
						notexist(ll[1]);
						buffy.close();
						return null;
					}					
				} else if (ll[0].equalsIgnoreCase(PAR_TMP_DIR)) {
					File f= new File(ll[1]);
					if (!f.exists()) {
						notexist(ll[1]);
						buffy.close();
						return null;
					}
					System.setProperty(Constants.PROPERTY_TMPDIR, f.getAbsolutePath());
				} else if (ll[0].equalsIgnoreCase(PAR_FORCE)) {
					pars.force= true;
				} else if (ll[0].equalsIgnoreCase(PAR_COPY_INPUT)) {
					pars.ioInTemp= true;
				} else if (ll[0].equalsIgnoreCase(PAR_ANNOTATION_MAPPING)) {
					String[] lll= ll[1].split(",");
					for (int i = 0; i < lll.length; i++) {
						if (lll[i].equalsIgnoreCase(ANNOTATION_MAPPING_PAIRED))
							pars.pairedEnd= true;
						else if (lll[i].equalsIgnoreCase(ANNOTATION_MAPPING_STRANDED))
							pars.stranded= true;
						else if (lll[i].equalsIgnoreCase(ANNOTATION_MAPPING_COMBINED)) {
							pars.stranded= true;
							pars.pairedEnd= true;
						} else if (lll[i].equalsIgnoreCase(ANNOTATION_MAPPING_SINGLE)) {
							pars.pairedEnd= false;
						} else {
							if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) {
								System.err.println("[BAD] No valid annotation mapping "+ ll[1]);
							}
							buffy.close();
							return null;
						}
					}
				} else if (ll[0].equalsIgnoreCase(PAR_READ_NUMBER)) {
					pars.readNr= Long.parseLong(ll[1]);
					if (pars.readNr< 0) {
						if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) {
							System.err.println("[BAD] No valid read number "+ ll[1]);
						}
						buffy.close();
						return null;
					}
				} else if (ll[0].equalsIgnoreCase(PAR_READ_DESCRIPTOR)) {
					pars.descriptor= new UniversalReadDescriptor();
					try {
						pars.descriptor.init(ll[1]);
					} catch (Exception e) {
						if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) {
							System.err.println("[BAD] No valid read descriptor "+ ll[1]);
						}
						buffy.close();
						return null;
					}
				} else if (ll[0].equalsIgnoreCase(PAR_SORT_IN_RAM)) {
					if (ll[1].equalsIgnoreCase(PAR_YES))
						pars.sortInRam= true;
					else if (ll[1].equalsIgnoreCase(PAR_NO))
						pars.sortInRam= false;
					else
						Log.error("Invalid value "+ ll[1]+ ", please use "+ PAR_YES+ " or "+
								PAR_NO+ " to specify a value for "+ PAR_SORT_IN_RAM);
				} else {
					if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) {
						System.err.println("[UNKNOWN] I couldnt understand parameter "+ ll[0]);
					}
					buffy.close();
					return null;
				}
			}
			buffy.close();
			return pars;
		} catch (Exception e) {
			if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
				e.printStackTrace();
			return null;
		}
	}
	
	
	private static void notexist(String string) {
		if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP) {
			System.err.println("[INVALID] no path to "+ string);
		}
	}


	File fileAnnotation, fileAnnotationSorted, fileMappings, fileMappingsSorted, 
		fileProfile, fileMapped, fileNotmapped, fileInsert, fileLPzip, fileStdErr, fileStdOut,
		fileCoverage;
	boolean force= false, ioInTemp= false, inputOverwrite= false, pairedEnd= false, stranded= false;
	boolean sortInRam= false;
	long readNr= -1;
	UniversalReadDescriptor descriptor;
	
	public FluxCapacitorParameters() {
		
	}
	
}
