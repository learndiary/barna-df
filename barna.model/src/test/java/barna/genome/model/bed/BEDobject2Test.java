package barna.genome.model.bed;

import barna.commons.ByteArrayCharSequence;
import org.junit.Test;

import static junit.framework.Assert.assertEquals;


public class BEDobject2Test {
	
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
		BEDobject2 bed2= new BEDobject2(bacs);

		assertEquals(start, bed2.getStart());
		assertEquals(end, bed2.getEnd());
		assertEquals(chr, bed2.getChr().toString());
		assertEquals(score, bed2.getScore());
		assertEquals((byte) -1, bed2.getStrand());
		assertEquals(blockNr, bed2.getBlockCount());
		assertEquals(blockSizes1, bed2.getNextBlockSize());
		assertEquals(blockStart1, bed2.getNextBlockStart());
		assertEquals(blockSizes2, bed2.getNextBlockSize());
		assertEquals(blockStart2, bed2.getNextBlockStart());
	}
}
