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

package fbi.genome.io;

import fbi.genome.model.gff.GFFObject;
import fbi.genome.io.gtf.GTFwrapper;
import fbi.genome.model.Gene;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

public class DNA2RNA {

	public static void main(String[] args) {
		DNA2RNA myConverter= new DNA2RNA(
				"c:\\annotations\\hg18_EnsemblGenes_fromUCSC090615_sorted.gtf", 
				"m:\\projects\\solexa\\publication\\interspecies\\interspecies2\\vingron\\1160342s_tableS9b_sorted.gtf", 
				"m:\\projects\\solexa\\publication\\interspecies\\interspecies2\\vingron\\1160342s_tableS9b_exonic.gtf"
		);
		myConverter.dna2rna();
	}
	
	File annFile, inFile, outFile;
	public DNA2RNA(String annFileName, String inFileName, String outFileName) {
		this.annFile= new File(annFileName);
		this.inFile= new File(inFileName);
		this.outFile= new File(outFileName);
	}
	
	
	GTFwrapper annReader, inReader;
	private GTFwrapper getAnnReader() {
		if (annReader == null) {
			annReader = new GTFwrapper(annFile.getAbsolutePath());
			if (!annReader.isApplicable()) {
				annFile= annReader.sort();
				annReader = new GTFwrapper(annFile.getAbsolutePath());
			}
				
			annReader.setReadGene(true);
		}

		return annReader;
	}
	private GTFwrapper getInReader() {
		if (inReader == null) {
			inReader = new GTFwrapper(inFile.getAbsolutePath());
			if (!inReader.isApplicable()) {
				inFile= inReader.sort();
				inReader = new GTFwrapper(inFile.getAbsolutePath());
			}
			inReader.setReadGene(false);
			inReader.setReadGTF(true);
			inReader.setReadAllLines();
			inReader.setSilent(true);
		}

		return inReader;
	}
	
	
	Gene[] getGenes() {
		try {
			getAnnReader().read();
		} catch (Exception e) {			
			e.printStackTrace();
		}
		return getAnnReader().getGenes();
	}
	GFFObject[] getObjects() {
		try {
			getInReader().read();
		} catch (Exception e) {			
			e.printStackTrace();
		}
		return getInReader().getGtfObj();
	}
	
	BufferedWriter buffy;
	private BufferedWriter getBuffy() {
		if (buffy == null) {
			try {
				buffy = new BufferedWriter(new FileWriter(outFile));
			} catch (IOException e) {
				e.printStackTrace();
			}
		}

		return buffy;
	}
	
	void write(String s) {
		try {
			getBuffy().write(s+ "\n");
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	void dna2rna() {
		try {
			while (true) {
				GFFObject[] obs= getObjects();
				Gene[] genes= getGenes();
				int i= 0, j = 0;
				while (true) {
					while (obs!= null&& i< obs.length&& genes!= null&& j < genes.length) {
						int c= genes[j].getChromosome().compareTo(obs[i].getChromosome());
						if (c== 0) 
							c= genes[j].getStrand()- obs[i].getStrand();
						
						if (c> 0) { // gene has to move
							++j;
							continue;
						} else if (c< 0) {	// obj has to move
							++i;
							continue;
						}
						int end= Math.max(obs[i].getStrand()* obs[i].getStart(),
								obs[i].getStrand()* obs[i].getEnd());
						if (genes[j].get5PrimeEdge()> end) {
							System.err.println("skipped obj "+ obs[i].getTranscriptID());
							++i;
							continue;
						}
						for (int m = 0; m < genes[j].getTranscripts().length; m++) {
							int k= 0;
							for (k = i; k < obs.length; k++) {
								end= Math.max(obs[k].getStrand()* obs[k].getStart(), obs[k].getStrand()* obs[k].getEnd());
								if (end> genes[j].get3PrimeEdge())
									break;
								if (genes[j].getTranscripts()[m].
										getTranscriptID().equals(obs[k].getTranscriptID())) {
									int genStart= obs[i].getStart()* obs[i].getStrand(),
										genEnd= obs[i].getEnd()* obs[i].getStrand();
									int txStart= genes[j].getTranscripts()[m].getExonicPosition(genStart),
										txEnd= genes[j].getTranscripts()[m].getExonicPosition(genEnd);
									obs[i].setStart(txStart);
									obs[i].setEnd(txEnd);
									write(obs[i].toString());
								}
							}
							i= k;
						}
						if (i!= obs.length)
							++j;
					}
					if (i== obs.length) {
						obs= getObjects();
						i= 0;
						--j;
					}
					if (j== genes.length) {
						genes= getGenes();
						j= 0;
					}
				}
			}			
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
