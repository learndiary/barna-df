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

package barna.io.gtf;

import barna.io.FileHelper;
import barna.model.*;
import barna.model.gff.GFFObject;

import java.io.*;
import java.util.HashMap;

/**
 * Reads locus by locus and
 * - determines orientation of mRNAs and ESTs
 * - annotates CDSs on mRNAs
 * @author micha
 *
 */
public class GFFreAnnotator {

	public static final String COMPLEMENT= "complement(";
	
	class CDS {
		int start, stop, qstart, qstop, qsize;
		boolean startCodon, stopCodon, complement, sense;
		@Override
		public String toString() {
			StringBuilder sb= new StringBuilder();
			if (complement)
				sb.append("complement(");
			if (!startCodon)
				sb.append("<");			
			sb.append(start);
			sb.append("..");
			if (!stopCodon)
				sb.append(">");
			sb.append(stop);
			if (complement)
				sb.append(")");
			return sb.toString();
		}
	}
	
	File fileAnnotation, fileCDS, fileGenomeDir, fileOutput;
	
	public GFFreAnnotator(File annotation, File cds, File genomeDir, File outFile) {
		this.fileAnnotation= annotation;
		this.fileCDS= cds;
		this.fileGenomeDir= genomeDir;
		this.fileOutput= outFile;
	}
	
	public static void main(String[] args) {
		GFFreAnnotator myAnnotator= new GFFreAnnotator(new File(args[0]), new File(args[1]), new File(args[2]), new File(args[3]));
		myAnnotator.run();
	}
	
	//#hg17.all_mrna.strand   hg17.all_mrna.qName     hg17.all_mrna.tName     hg17.all_mrna.blockSizes        hg17.all_mrna.qStarts   hg17.all_mrna.tStarts hg17.cds.name
	private int fieldStrand= 0, fieldQName= 1, fieldTName= 2, fieldBSizes= 3, fieldQStarts= 4, fieldtStarts= 5, fieldCDS= 6;
	private void initFieldNr(String s) {
		String[] ss= s.split("\t");
		for (int i = 0; i < ss.length; i++) {
			if (ss[i].equalsIgnoreCase("hg17.all_mrna.strand"))
				fieldStrand= i;
			else if (ss[i].equalsIgnoreCase("hg17.all_mrna.qName"))
				fieldQName= i;
			else if (ss[i].equalsIgnoreCase("hg17.all_mrna.tName"))
				fieldTName= i;
			else if (ss[i].equalsIgnoreCase("hg17.all_mrna.blockSizes"))
				fieldBSizes= i;
			else if (ss[i].equalsIgnoreCase("hg17.all_mrna.qStarts"))
				fieldQStarts= i;
			else if (ss[i].equalsIgnoreCase("hg17.all_mrna.tStarts"))
				fieldtStarts= i;
			else if (ss[i].equalsIgnoreCase("hg17.cds.name"))
				fieldCDS= i;
		}
	}
	
	
	/*
	 * #hg17.all_mrna.strand   hg17.all_mrna.qName     hg17.all_mrna.qSize     hg17.all_mrna.qStart    hg17.all_mrna.qEnd      hg17.all_mrna.tName   hg17.all_mrna.tStart    hg17.cds.name
+       AL831987        4633    1       4611    chr1    67102660        join(248..316,316..1848)
+       AY124190        2554    0       2554    chr1    67102662        234..1769
+       AY124189        2629    0       2629    chr1    67102662        234..1844
	 */
	HashMap<String, CDS> getCDSmap() {
		HashMap<String, CDS> map= new HashMap<String, GFFreAnnotator.CDS>(200000, 1f);
		int cntTotal= 0, cntComplete= 0, cnt3trunc= 0, cnt5trunc= 0, cnt53trunc= 0, 
			cntComplement= 0, cntFrag= 0, cntNA= 0;
		System.err.print("reading CDS file..");
		System.err.flush();
		try {
			BufferedReader buffy= new BufferedReader(new FileReader(fileCDS.getAbsolutePath()));
			String tid, cds;
			for(String s= null; (s= buffy.readLine())!= null;) {
				if (s.startsWith("#")) {
					if (map.size()== 0) {
						initFieldNr(s);
					}
					continue;
				}
				++cntTotal;
				String[] pp= s.split("\t");
				// 0=strand, 1=qName, 2=qSize, 3=qStart, 4=qEnd, 5=tName, 6=tStart, 7=name
/*				int p1= s.indexOf("\t"), p2= s.indexOf("\t", p1+ 1),
					p3= s.indexOf("\t", p2+ 1), p4= s.indexOf("\t", p3+ 1), p5= s.indexOf("\t", p4+ 1),
					p6= s.indexOf(" ", p5+ 1);
				if (p6> 0) 
					s= s.substring(0, p6);
*/				
				//tid= new String(s.substring(p1+ 1, p5));	// clone for map
				tid= pp[1]+ "\t"+ pp[5]+ "\t"+ pp[6]; // tid+ chr+ chr_start
				cds= pp[pp.length- 1];
				// "9517 isolated from HL-60 cills including"
				int p= cds.indexOf(" ");
				if (p> 0) {
					for (int i = p; i < cds.length(); ++i) {
						if (Character.isLetter(cds.charAt(i))) {
							cds= cds.substring(0, p);
							break;
						}
					}
				}
				if (cds.startsWith("n/a")) {
					++cntNA;
					continue;
				}
				if (cds.contains("join")) {	// fragmented, skip
					++cntFrag;
					continue;
				}
				
				CDS cdsObj= new CDS();
				cdsObj.qsize= Integer.parseInt(pp[2]);
				cdsObj.qstart= Integer.parseInt(pp[3]);
				cdsObj.qstop= Integer.parseInt(pp[4]);
				if (s.startsWith("+"))
					cdsObj.sense= true;
				else
					cdsObj.sense= false;
				if (cds.startsWith(COMPLEMENT)) {
					cds= cds.substring(COMPLEMENT.length(), cds.length()- 1);
					cdsObj.complement= true;
					++cntComplement;
				} else
					cdsObj.complement= false;
				if (cds.startsWith("<")) {
					cds= cds.substring(1);
					cdsObj.startCodon= false;
				} else {
					cdsObj.startCodon= true;
				}
				int sep= cds.indexOf("..");
				cdsObj.start= Integer.parseInt(cds.substring(0, sep));
				if (cds.charAt(sep+ 2)== '>') {
					cdsObj.stopCodon= false;
					cdsObj.stop= Integer.parseInt(cds.substring(sep+ 3));
					if (cdsObj.startCodon)
						++cnt3trunc;
					else 
						++cnt53trunc;
				} else {
					cdsObj.stopCodon= true;
					try {
						cdsObj.stop= Integer.parseInt(cds.substring(sep+ 2));
					} catch (NumberFormatException ex) {
						System.err.println(barna.commons.system.OSChecker.NEW_LINE+ s);
						ex.printStackTrace();
					}
					if (cdsObj.startCodon) 
						++cntComplete;
					else
						++cnt5trunc;
				}
				
				map.put(tid, cdsObj);
			}
			
			buffy.close();
			System.err.println("done.");
			System.err.println("total\t"+cntTotal);
			System.err.println("complete\t"+cntComplete);
			System.err.println("5trunc\t"+cnt5trunc);
			System.err.println("3trunc\t"+cnt3trunc);
			System.err.println("53trunc\t"+cnt53trunc);
			System.err.println();
			System.err.println("complement\t"+cntComplement+"(contained in above)");
			System.err.println("join\t"+cntComplement+"(discarded)");
			System.err.println("n/a\t"+cntNA+"(nothint to do)");
			System.err.println();
			return map;
			
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(-1);
		}
		return null;
	}
	
	public static boolean checkAcceptableIntronSeq(String donSeq, String accSeq) {
//		   if (($seq1 =~ /^GT/i && $seq2 =~ /AG$/i)
//		   ||($seq1 =~ /^GC/i && $seq2 =~ /AG$/i)
//		   ||($seq1 =~ /^ATATC/i && $seq2 =~ /AG$/i)
//		   ||($seq1 =~ /^ATATC/i && $seq2 =~ /AC$/i)
//		   ||($seq1 =~ /^ATATC/i && $seq2 =~ /AT$/i)
//		   ||($seq1 =~ /^GTATC/i && $seq2 =~ /AT$/i)
//		   ||($seq1 =~ /^ATATC/i && $seq2 =~ /AA$/i) 

		if (donSeq.length()!= 5|| accSeq.length()!= 2) {
			System.err.println("checkAcceptableIntronSeq mismatch");
			return false;
		}
		
		String donSeq2= donSeq.substring(0,2);
		if (!((donSeq.startsWith("G")|| donSeq.startsWith("A"))&& accSeq.startsWith("A")))
			return false;
		if (accSeq.equals("AG")) {
			if (donSeq2.equals("GT")|| donSeq2.equals("GC")|| donSeq.equals("GTATC"))
				return true;
			return false;
		} else {
			if (donSeq.equals("ATATC"))
				return true;
			if (donSeq.equals("GTATC")&& accSeq.equals("AT"))
				return true;
			return false;
		}
	}
	
	public static int checkAcceptableIntron(SpliceSite donor, SpliceSite acceptor) {
		String donSeq= Graph.readSequence(donor, 0, 3).toUpperCase();
		String accSeq= Graph.readSequence(acceptor, 3, 0).toUpperCase();
		
		if (checkAcceptableIntronSeq(donSeq, accSeq.substring(3)))
			return 1;
		String h= donSeq;
		donSeq= Graph.reverseSequence(Graph.complementarySequence(accSeq));
		accSeq= Graph.reverseSequence(Graph.complementarySequence(h));
		if (checkAcceptableIntronSeq(donSeq, accSeq.substring(3)))
			return -1;
		return 0;
	}
	
	public void run() {
		
		if ((!fileGenomeDir.exists())|| (!fileGenomeDir.isDirectory())) {
			System.err.println("invalid genome dir "+ fileGenomeDir.getAbsolutePath());
			System.exit(-1);
		}
		Graph.overrideSequenceDirPath= fileGenomeDir.getAbsolutePath();
		
		HashMap<String, CDS> cdsMap= null;
		if (fileCDS!= null&& fileCDS.exists())
			cdsMap= getCDSmap();
		
		GTFwrapper reader= new GTFwrapper(fileAnnotation.getAbsolutePath());
		reader.setNoIDs(null);
		if (!reader.isApplicable()) {
			File tmpGTF= reader.sort();
			fileAnnotation= new File(
				fileAnnotation.getAbsolutePath()+"_sorted"
			);
			if (FileHelper.move(tmpGTF, this.fileAnnotation, null)) {
				reader= new GTFwrapper(fileAnnotation.getAbsolutePath());
				System.err.println("\n\tmoved sorted file to: "+ fileAnnotation.getAbsolutePath());
			} else
				reader= new GTFwrapper(tmpGTF.getAbsolutePath());
		}
		Gene[] genes= null;
		int cntSwap= 0, cntCDSano= 0, cntCDSdis= 0, cntGapAli= 0;
		System.err.print("re-annotating..");
		System.err.flush();
		try {
			BufferedWriter writer= new BufferedWriter(new FileWriter(fileOutput.getAbsolutePath()));
			for(reader.read(); (genes= reader.getGenes())!= null; reader.read()) {
				for (int i = 0; i < genes.length; i++) {
					int ref= 0, nRef= 0;
					for (int j = 0; j < genes[i].getTranscripts().length; j++) {
						Transcript t= genes[i].getTranscripts()[j];
						if(t.getTranscriptID().equals("BC156865_dup1"))
							System.currentTimeMillis();
						byte source= t.getSourceType();
						if (source== Transcript.ID_SRC_REFSEQ) {
							++ref;
							continue;
						}
						
						// determine directionality
						int ssSense= 0, ssAnti= 0;
						boolean turnAround= false;
						for (int k = 1; k < t.getExons().length; k++) {
							SpliceSite donor= new SpliceSite(t.getExons()[k-1].get3PrimeEdge(), 
									SpliceSite.TYPE_DONOR, genes[i]),
									acceptor= new SpliceSite(t.getExons()[k].get5PrimeEdge(),
									SpliceSite.TYPE_ACCEPTOR, genes[i]);
							int dir= checkAcceptableIntron(donor, acceptor);
							if (dir== 1) 
								++ssSense;
							else if (dir== -1)
								++ssAnti;
						}
						if (ssAnti> ssSense) {
							++cntSwap;
							turnAround= true;
							Gene g2= new Gene(null);
							g2.setChromosome(t.getChromosome());
							g2.setStrand((byte) ((-1)* t.getStrand()));
							Transcript t2= new Transcript(t.getTranscriptID());							
							t2.setStrand((byte) ((-1)* t.getStrand()));
							t2.setSource(t.getSource());
							t2.setGene(g2);
							for (int k = 0; k < t.getExons().length; k++) {
								Exon e1= t.getExons()[t.getExons().length- 1- k];
								SpliceSite acc1= e1.getAcceptor(),
											don1= e1.getDonor();
//								SpliceSite acc2= don1.isDonor()?new SpliceSite((-1)* don1.getPos(), SpliceSite.TYPE_ACCEPTOR, g2)
//									:(don1.isSoftEdge()?new SpliceSite((-1)* don1.getPos(), SpliceSite.TYPE_SOFT_START, g2)
//										:new SpliceSite((-1)* don1.getPos(), SpliceSite.TYPE_HARD_START, g2));
//								SpliceSite don2= acc1.isAcceptor()?new SpliceSite((-1)* acc1.getPos(), SpliceSite.TYPE_DONOR, g2)
//								:(acc1.isSoftEdge()?new SpliceSite((-1)* acc1.getPos(), SpliceSite.TYPE_SOFT_END, g2)
//									:new SpliceSite((-1)* acc1.getPos(), SpliceSite.TYPE_HARD_END, g2));

								//Exon ee= new Exon(t2,"",acc2.getPos(),don2.getPos());
								Exon e2= new Exon(t2,"",(-1)* don1.getPos(), (-1)* acc1.getPos());
								t2.addExon(e2);
								//exons2[k]= t.getExons()[exons2.length- 1- k];
							}
							t= t2;
						} else
							++nRef;
						
						// annotate CDS
						if (cdsMap!= null) {
							String key= t.getTranscriptID()+ "\t"+ t.getChromosome()+ "\t"+ Integer.toString(Math.abs(t.getStart())- 1);// tid+ chr+ chr_start
							CDS cds= cdsMap.get(key);
							if (source== Transcript.ID_SRC_MRNA&& cds!= null) {
								byte strandTx= t.getStrand(), 
									strandAli= cds.sense?(byte)1:(byte)-1,
									strandCDS= cds.complement?(byte)-1:(byte)1;
								if (turnAround)
									strandTx*= -1;
								// annotate?
								if (strandTx* strandAli== strandCDS) {
									int tlen= t.getExonicLength();
									Translation tra= new Translation(t);
									int estart= cds.start, 
										estop= cds.stop;
									if (turnAround&& cds.complement) {
										int h= estart;									
										estart= cds.qsize- estop;
										estop= cds.qsize- h;
										h= cds.qstart;
										cds.qstart= cds.qstop;
										cds.qstop= h;
									}
									if (estart- 1< cds.qstart) {
										cds.startCodon= false;
										estart= cds.qstart+ 1;
									}
									if (estop- 1> cds.qstop) {
										cds.stopCodon= false;
										estop= cds.qstop;
									}
									// gaps in alignment
									//if (estop- cds.qstart- 1>= tlen)
									//	estop= tlen+ cds.qstart;
	//								if (cds.qstop- cds.qstart+ 1!= tlen) {
	//									int x= t.getGenomicPosition(estop);
	//									estop= tlen- (cds.qstop- estop);
	//									int y= t.getGenomicPosition(estop);
	//									System.currentTimeMillis();
	//									cds= null;
	//								}
									int start= t.getGenomicPosition(estart- cds.qstart- 1),
									stop= t.getGenomicPosition(estop- cds.qstart- 1);
	
									if (stop- start== estop- estart) {
										tra.set5PrimeEdge(start);
										tra.set3PrimeEdge(stop);
										if (cds.startCodon)
											tra.setCodonStart(new SpliceSite(start, SpliceSite.TYPE_CODON_START, t.getGene()));
										if (cds.stopCodon)
											tra.setCodonStop(new SpliceSite(stop, SpliceSite.TYPE_CODON_STOP, t.getGene()));
										t.setTranslations(new Translation[]{tra});
										++cntCDSano;
									} else 
										++cntGapAli;
								} else
									++cntCDSdis;
							}
							
							// output
	/*						if (!turnAround) {
								if (cds== null)
									continue;
								else {
									//System.out.println("annotated: "+ cds);
									continue;
								}
							} else
								System.out.println("swapped:");
	*/						
							if (turnAround&& cds!= null)
								System.currentTimeMillis();
						}
						
						GFFObject[] obj= GFF.fromTranscript(t, true, true, true, true, false);
						for (int k = 0; k < obj.length; k++) {
							writer.write(obj[k].toString()+ barna.commons.system.OSChecker.NEW_LINE);
						}
					}
				}
			}
			writer.flush();
			writer.close();
			
			System.err.println("finished.");
			System.err.println("swapped\t"+ cntSwap);
			System.err.println("annotated CDS\t"+cntCDSano);
			System.err.println("discarded CDS\t"+(cntCDSdis+ cntGapAli));
			System.err.println("\t\tdisregarded "+ cntCDSdis);
			System.err.println("\t\tgapped ali "+ cntGapAli);
		} catch (Exception ex) {
			ex.printStackTrace();
			System.exit(-1);
		}
	}
	
}
