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

package barna.flux.simulator;

import barna.io.gtf.GTFwrapper;
import barna.model.Gene;
import barna.model.Transcript;
import barna.model.bed.BEDobject2;
import org.junit.Test;

import java.io.File;
import java.util.Random;

import static junit.framework.Assert.assertEquals;

public class SequencerTest {

    @Test
	public void testReadGeneration() throws Exception {
		
		int nrTests= 100; // nr. of tests per transcript
		int maxBoundary= 100; // nt added at beginning/end of tx
		Random rand= new Random();
		
		GTFwrapper reader= new GTFwrapper(new File(
				getClass().getResource("/sacCer2_sgdGene_fromUCSC110329.gtf").getFile()).getAbsolutePath());
		Gene[] genes;
		for (reader.read(); (genes=reader.getGenes())!= null; reader.read()) {
			for (int i = 0; i < genes.length; i++) {
				Gene g= genes[i];
				int n= g.getTranscriptCount();
				for (int j = 0; j < n; j++) {
					Transcript t= g.getTranscripts()[j];
					int tlen= t.getExonicLength();						
					for (int k = 0; k < nrTests; k++) {
						int maxLeft= maxBoundary;
						if (t.getStrand()>= 0)
							maxLeft= Math.min(maxLeft, t.get5PrimeEdge());
						int maxRight= maxBoundary;
						if (t.getStrand()< 0)
							maxRight= Math.min(maxRight, Math.abs(t.get3PrimeEdge()));
						int boundLeft= (maxLeft== 0? 0: rand.nextInt(maxLeft));
						int boundRight= (maxRight== 0? 0: rand.nextInt(maxRight));
						
						byte absDir= t.getStrand();
						int fragStart= rand.nextInt(tlen+ boundLeft+ boundRight)- boundLeft;
						int fragEnd= fragStart+ rand.nextInt(tlen+ boundRight- fragStart);
						int fragLen= fragEnd- fragStart+ 1;
						if (fragLen<= 1)
							continue;
						int readLen= 1+ rand.nextInt(fragLen- 1);
						boolean left= rand.nextBoolean();
						int start= (left? fragStart: fragEnd- readLen+ 1);
						int end= (left? fragStart+ readLen- 1: fragEnd);
						assert(start<= end);
						BEDobject2 obj= new BEDobject2();
						int polyA= Sequencer.createRead(obj, start, end, t, k,
                                absDir, fragStart, fragEnd, left, 1);
						
						assertEquals(readLen, obj.getLength());
						int bedLeft= t.getGenomicPosition(start);
						int bedRight= t.getGenomicPosition(end);
						if (absDir< 0) {
							int h= bedLeft;
							bedLeft= -bedRight;
							bedRight= -h;
						}
						--bedLeft;
						
						if (polyA> 0) {
							if (polyA== readLen) {
								assertEquals(0, obj.getStart());
								assertEquals(readLen, obj.getEnd());
							}
						} else {
							assertEquals(bedLeft, obj.getStart());
							assertEquals(bedRight, obj.getEnd());
						}						
						obj.getLength();
					}
				}
			}
		}
	}
	
}
