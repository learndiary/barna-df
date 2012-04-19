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

package barna.io;

import barna.io.gtf.GTFwrapper;
import barna.model.Gene;
import barna.model.gff.GFFObject;

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
