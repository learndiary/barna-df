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

package barna.model.tools;

import barna.model.Transcript;
import barna.model.Translation;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.Vector;

public class NMDSimulator {
	public static int MIN_DIST_NC_EJC= 50;
	public static int MAX_DIST_NC_EJC= 55;
	public static int MIN_ORF_LENGTH_AA= 35;

	public static void main(String[] args) {
		//testNMD();
		//testNMD_noncoding();
		//outputNMDEvents();
	}
	
	Transcript trpt= null;
	
	public NMDSimulator(Transcript t) {
		this.trpt= t;
	}
	
	/**
	 * The stop codon must be in the last exon or no further than 50bp`
	 * from the end of the penultimate exon. [HAVANA]
	 *
	 * @param tln a translation
	 * @param minDistNt a minimum distance in nucleotides
	 * @return <code>true</code> if the translation frame is terminating more than the provided distance upstream of
	 * the exon junction complex, <code>false</code> otherwise
	 */
	public boolean isTerminatingUpstreamOfEJC(Translation tln, int minDistNt) {
		
		int critPoint= tln.get3PrimeEdge(); 
		int[] prematStops= tln.getPrematureStops();	// genomic positions
		if (prematStops.length> 0&& prematStops[0]< critPoint)
			critPoint= prematStops[0];
		critPoint= trpt.getGenomicPosition(trpt.getExonicPosition(critPoint)+ minDistNt);
		
		int i;
		for (i = 0; i < trpt.getExons().length; i++) 
			if (trpt.getExons()[i].get5PrimeEdge()> critPoint)	// compare 5' !!
				break;
				
		if (i< trpt.getExons().length)
			return true;
		return false;
	}

	/**
	 * Never annotate an ATG starting internal of another CDS &gt; 35 aa upstream
	 * of the ATG as is subject to NMD. [HAVANA]
	 * 
	 * @param trans a translation
	 * @param maxDistAA a minimum distance in amino acids
	 * @deprecated HAVANA seems to mean something different by this phrase,
	 * see <code>hasUsORF</code>
	 * @return <code>true</code> if the CDS shows an internal ATG, <code>false</code> otherwise
	 */
	public boolean isInternalATG(Translation trans, int maxDistAA) {
		int maxDistNt= maxDistAA* 3;
		int startPos= trpt.getExonicPosition(trans.get5PrimeEdge());
		if (startPos< maxDistNt)
			return false;
		int frame= startPos% 3;
		String seq= trpt.getSplicedSequence();
		seq= seq.substring(frame);	// to search in the correct frame
		int limit= startPos- frame- maxDistNt;
		int i;	// find last inframe start before ATG (outside the safe distance)
		int usATGPos= 0;
		for (i = 0; (i+3) < limit; i+=3) 
			if (seq.substring(i, i+3).equalsIgnoreCase(Translation.START_CODON))
				usATGPos= i;
		if (usATGPos== 0)
			return false;	// no upstream atg >35 aa found
		
		seq= seq.substring(usATGPos, startPos- frame);
		for (i = 0; (i+3) < seq.length(); i+=3) 
			if (seq.substring(i, i+3).equalsIgnoreCase(Translation.STOP_CODONS[0])||
				seq.substring(i, i+3).equalsIgnoreCase(Translation.STOP_CODONS[1])||
				seq.substring(i, i+3).equalsIgnoreCase(Translation.STOP_CODONS[2]))
					return false;	// inframe stop found, ATG is not internal to other ORF
		return true;
	}

	/**
	 * Never annotate an ATG starting internal of another CDS &gt; 35 aa upstream
	 * of the ATG as is subject to NMD. [HAVANA]
	 * 
	 * @param trans a translation
	 * @param minSizeAA minimal size in amino acids
	 * @return <code>true</code> if there is an upstream open reading frame with the given attributes,
	 * <code>false</code> otherwise
	 */
	public boolean hasUsORF(Translation trans, int minSizeAA) {
		int minSizeNt= minSizeAA* 3;
		int startPos= trans.getTranscript().getExonicPosition(trans.get5PrimeEdge());
		if (startPos< minSizeNt)
			return false;	// cannot have a us ORF of minSize
	
		Translation[] uORFs= trans.getUsORF();
		int cntORF= 0;
		Vector usORFV= new Vector();
		for (int i = 0; uORFs!= null&& i < uORFs.length; i++) {
			if (uORFs[i].getSplicedLength()< minSizeNt)	// bit redundant since already not predicted, but important when predMinSize<> minSize here
				continue;
			++cntORF;
			usORFV.add(uORFs[i]);
		}
		
		return (cntORF> 0);
	}
	
	void checkATG() {
		Translation trans= trpt.getTranslations()[0];
//		int[] ncPos= trans.getNCPosAA();
//		
//		if (ncPos== null|| ncPos.length< 1)
//			return false;
		
		PrintStream p= null;
		try {
			FileOutputStream f= new FileOutputStream("_nonATGseqs.txt", true);
			p= new PrintStream(f);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		String seq= trpt.getCDSSequenceNt().toUpperCase();
		if (!seq.startsWith("ATG")) {
			p.println(">"+trans.getTranscript().getTranscriptID()+
					" "+ trans.getTranscript().getChromosome()+
					" "+ trans.getStrand());
			p.println(seq);
		} else {
			if (trans.getStrand()< 0)
				System.out.println(trans.getTranscript().getTranscriptID());
		}
		p.flush();
		p.close();

	}
	
	public boolean isNMD(Translation trln) {
		boolean ejcDS= isTerminatingUpstreamOfEJC(trln, MIN_DIST_NC_EJC);
		boolean atgUS= false;	//hasUsORF(trln, MIN_ORF_LENGTH_AA);
		return (ejcDS|| atgUS);
	}
	public boolean isNMD() {
		Translation trln= trpt.findLongestORF();
		return isNMD(trln);
	}
}
