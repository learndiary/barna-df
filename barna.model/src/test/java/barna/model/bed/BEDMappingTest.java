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
import barna.model.rna.UniversalReadDescriptor;
import org.junit.BeforeClass;
import org.junit.Test;

import static org.junit.Assert.assertEquals;


public class BEDMappingTest {

    private static BEDMapping bed;

	@BeforeClass
    public static void init() {
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
        bed= new BEDMapping(bacs);
    }

    @Test
	public void testInit() {
		assertEquals(123, bed.getStart());
		assertEquals(456, bed.getEnd());
		assertEquals("chrX", bed.getChromosome().toString());
		assertEquals(1, bed.getScore());
		assertEquals((byte) -1, bed.getStrand());
		assertEquals(2, bed.getBlockCount());
		assertEquals(111, bed.getNextBlockSize());
		assertEquals(0, bed.getNextBlockStart());
		assertEquals(222, bed.getNextBlockSize());
		assertEquals(111, bed.getNextBlockStart());
	}

    @Test
    public void testName() throws Exception {
        UniversalReadDescriptor d = new UniversalReadDescriptor();
        d.init(UniversalReadDescriptor.DESCRIPTORID_PAIRED);
        bed.descriptor = d;
        assertEquals("test1", bed.getName(true).toString());
        assertEquals("test1", bed.getName(false).toString());
        bed.setName("test1/1");
        assertEquals("test1/1", bed.getName(true).toString());
        assertEquals("test1", bed.getName(false).toString());
    }

    @Test (expected=UnsupportedOperationException.class)
    public void testSequence() throws Exception {
        bed.getSequence();
    }

    @Test (expected=UnsupportedOperationException.class)
    public void testCigar() throws Exception {
        bed.getCigar();
    }

    @Test
    public void testMateFlag() throws Exception {
        bed.setDescriptor(UniversalReadDescriptor.getDefaultDescriptor());
        assertEquals(bed.getMateFlag(),0);
        UniversalReadDescriptor d = new UniversalReadDescriptor();
        d.init(UniversalReadDescriptor.DESCRIPTORID_PAIRED);
        ByteArrayCharSequence oldName = bed.getName(true);
        // test mate 1
        bed.setDescriptor(d);
        bed.setName(oldName+"/1");
        assertEquals(bed.getMateFlag(),1);
        // test mate 2
        bed.setName(oldName+"/2");
        assertEquals(bed.getMateFlag(),2);
        bed.setName(oldName);
    }
}
