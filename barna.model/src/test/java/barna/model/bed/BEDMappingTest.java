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

package barna.model.bed;

import barna.commons.ByteArrayCharSequence;
import org.junit.Test;

import static org.junit.Assert.assertEquals;


public class BEDMappingTest {
	
	@Test
	public void testInit() {
		String chr= "chrX";
		int start= 123;
		int end= 456;
		String name= "test1";
		int score= 1;
		String strand= "-";
		int thickStart= 0;
		int thickEnd= 0;
		String col= "0,0,0";
		int blockNr= 2;
		int blockSizes1= 111, blockSizes2= 222;
		int blockStart1= 0, blockStart2= 111;
		
		String bedLine= chr+ "\t"+ start+ "\t"+ end+ "\t"+name+ "\t"+ score+ "\t"+ strand+
			"\t"+ thickStart+ "\t"+ thickEnd+ "\t"+ col+ "\t"+ blockNr+ "\t"+ 
			blockSizes1+ ","+ blockSizes2+ "\t"+ blockStart1+ ","+ blockStart2;
		ByteArrayCharSequence bacs= new ByteArrayCharSequence(bedLine);
		BEDMapping bed2= new BEDMapping(bacs);

		assertEquals(start, bed2.getStart());
		assertEquals(end, bed2.getEnd());
		assertEquals(chr, bed2.getChromosome().toString());
		assertEquals(score, bed2.getScore());
		assertEquals((byte) -1, bed2.getStrand());
		assertEquals(blockNr, bed2.getBlockCount());
		assertEquals(blockSizes1, bed2.getNextBlockSize());
		assertEquals(blockStart1, bed2.getNextBlockStart());
		assertEquals(blockSizes2, bed2.getNextBlockSize());
		assertEquals(blockStart2, bed2.getNextBlockStart());
	}
}
