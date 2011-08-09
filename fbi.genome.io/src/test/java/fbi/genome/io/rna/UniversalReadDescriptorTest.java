package fbi.genome.io.rna;

import static org.junit.Assert.*;

import org.junit.Test;

import fbi.genome.io.rna.UniversalReadDescriptor.Attributes;

public class UniversalReadDescriptorTest {

	
	String simRead1= "chr1:4797974-4836816W:NM_008866:2:2433:548:757:548:S/1",
	 	simRead2= "chr1:4797974-4836816W:NM_008866:2:2433:548:757:548:S/2",
	 	simRead3= "chr1:4797974-4836816W:NM_008866:2:2433:548:757:548:A/1",
	 	simRead4= "chr1:4797974-4836816W:NM_008866:2:2433:548:757:548:A/2";
	
	String barnaID1= "BILLIEHOLIDAY:5:100:1000:1190/1s";
	String barnaID2= "BILLIEHOLIDAY:5:100:1000:1190/1a";
	String barnaID3= "BILLIEHOLIDAY:5:100:1000:1190/2";
	String barnaID4= "BILLIEHOLIDAY:5:100:1000:1190/2a";
	String barnaID5= "BILLIEHOLIDAY:5:100:1000:1190/1";
	
	String oldCshlID1= "BILLIEHOLIDAY:5:100:1000:1190/1_strand2";
	String oldCshlID2= "BILLIEHOLIDAY:5:100:1000:1190/1_strand0";
	
	String newCshlCombined= "MARILYN_0005:7:1:2804:1011#0/1";

	@Test
	public void testSimpleDescriptor() {
		UniversalReadDescriptor descriptor= new UniversalReadDescriptor();
		assertTrue(descriptor.init(UniversalReadDescriptor.getDescriptor(UniversalReadDescriptor.DESCRIPTORID_SIMPLE)));
	}
	
	@Test
	public void testSimpleID() {
		UniversalReadDescriptor descriptor= new UniversalReadDescriptor();
		descriptor.init(UniversalReadDescriptor.getDescriptor(UniversalReadDescriptor.DESCRIPTORID_SIMPLE));
		Attributes a= 
			descriptor.getAttributes(simRead1, null);
		assertEquals(a.id, simRead1);
		
		a= descriptor.getAttributes(simRead2, a);
		assertEquals(a.id, simRead2);
		
		a= descriptor.getAttributes(simRead3, a);
		assertEquals(a.id, simRead3);
		
		a= descriptor.getAttributes(simRead4, a);
		assertEquals(a.id, simRead4);
	}

	@Test
	public void testPairedDescriptor() {
		UniversalReadDescriptor descriptor= new UniversalReadDescriptor();
		assertTrue(descriptor.init(UniversalReadDescriptor.getDescriptor(UniversalReadDescriptor.DESCRIPTORID_PAIRED)));
	}
	
	@Test
	public void testPairedID() {
		UniversalReadDescriptor descriptor= new UniversalReadDescriptor();
		descriptor.init(UniversalReadDescriptor.getDescriptor(UniversalReadDescriptor.DESCRIPTORID_PAIRED));
		Attributes a= 
			descriptor.getAttributes(simRead1, null);
		assertEquals(a.id, simRead1.substring(0, simRead1.lastIndexOf('/')));
		assertEquals(a.flag, 1);
		
		a= descriptor.getAttributes(simRead2, a);
		assertEquals(a.id, simRead2.substring(0, simRead2.lastIndexOf('/')));
		assertEquals(a.flag, 2);
		
		a= descriptor.getAttributes(simRead3, a);
		assertEquals(a.id, simRead3.substring(0, simRead3.lastIndexOf('/')));
		assertEquals(a.flag, 1);
		
		a= descriptor.getAttributes(simRead4, a);
		assertEquals(a.id, simRead4.substring(0, simRead4.lastIndexOf('/')));
		assertEquals(a.flag, 2);
		
	}
	
	
	@Test
	public void testSimulatorDescriptor() {
		UniversalReadDescriptor descriptor= new UniversalReadDescriptor();
		assertTrue(descriptor.init(UniversalReadDescriptor.getDescriptor(UniversalReadDescriptor.DESCRIPTORID_SIMULATOR)));
	}
	
	@Test
	public void testSimulatorID() {
		UniversalReadDescriptor descriptor= new UniversalReadDescriptor();
		descriptor.init(UniversalReadDescriptor.getDescriptor(UniversalReadDescriptor.DESCRIPTORID_SIMULATOR));
		Attributes a= 
			descriptor.getAttributes(simRead1, null);
		assertEquals(a.id, simRead1.substring(0, simRead1.lastIndexOf(":548:")));
		assertEquals(a.strand, 1);
		assertEquals(a.flag, 1);
		
		a= descriptor.getAttributes(simRead2, a);
		assertEquals(a.id, simRead2.substring(0, simRead2.lastIndexOf(":548:")));
		assertEquals(a.strand, 1);
		assertEquals(a.flag, 2);
		
		a= descriptor.getAttributes(simRead3, a);
		assertEquals(a.id, simRead3.substring(0, simRead3.lastIndexOf(":548:")));
		assertEquals(a.strand, 2);
		assertEquals(a.flag, 1);
		
		a= descriptor.getAttributes("chr1:4797974-4836816W:NM_008866:2:2433:548:757:548:A/2", a);
		assertEquals(a.id, simRead4.substring(0, simRead4.lastIndexOf(":548:")));
		assertEquals(a.strand, 2);
		assertEquals(a.flag, 2);
		
	}
	
	@Test
	public void testBarnaDescriptor() {
		UniversalReadDescriptor descriptor= new UniversalReadDescriptor();
		assertTrue(descriptor.init(UniversalReadDescriptor.getDescriptor(UniversalReadDescriptor.DESCRIPTORID_BARNA)));
	}
	
	@Test
	public void testBarnaID() {
		UniversalReadDescriptor descriptor= new UniversalReadDescriptor();
		descriptor.init(UniversalReadDescriptor.getDescriptor(UniversalReadDescriptor.DESCRIPTORID_BARNA));
		Attributes a= 
			descriptor.getAttributes(barnaID1, null);
		assertEquals(a.id, barnaID1.substring(0, barnaID1.lastIndexOf("/")));
		assertEquals(a.flag, 1);
		assertEquals(a.strand, 1);
		
		a= descriptor.getAttributes(barnaID2, a);
		assertEquals(a.id, barnaID2.substring(0, barnaID2.lastIndexOf("/")));
		assertEquals(a.flag, 1);
		assertEquals(a.strand, 2);
		
		a= descriptor.getAttributes(barnaID3, a);
		assertEquals(a.id, barnaID3.substring(0, barnaID3.lastIndexOf("/")));
		assertEquals(a.flag, 2);
		assertEquals(a.strand, 0);
		
		a= descriptor.getAttributes(barnaID4, a);
		assertEquals(a.id, barnaID4.substring(0, barnaID4.lastIndexOf("/")));
		assertEquals(a.flag, 2);
		assertEquals(a.strand, 2);
		
		a= descriptor.getAttributes(barnaID5, a);
		assertEquals(a.id, barnaID5.substring(0, barnaID5.lastIndexOf("/")));
		assertEquals(a.flag, 1);
		assertEquals(a.strand, 0);
		
	}
	

	@Test
	public void testCSHLoldDescriptor() {
		UniversalReadDescriptor descriptor= new UniversalReadDescriptor();
		assertTrue(descriptor.init(UniversalReadDescriptor.getDescriptor(UniversalReadDescriptor.DESCRIPTORID_MATE_STRAND_CSHL)));
	}
	
	@Test
	public void testCSHLoldID() {
		UniversalReadDescriptor descriptor= new UniversalReadDescriptor();
		descriptor.init(UniversalReadDescriptor.getDescriptor(UniversalReadDescriptor.DESCRIPTORID_MATE_STRAND_CSHL));
		Attributes a= 
			descriptor.getAttributes(oldCshlID1, null);
		assertEquals(a.id, oldCshlID1.substring(0, oldCshlID1.lastIndexOf("/")));
		assertEquals(a.flag, 1);
		assertEquals(a.strand, 2);
		
		a= descriptor.getAttributes(oldCshlID2, null);
		assertEquals(a.id, oldCshlID2.substring(0, oldCshlID2.lastIndexOf("/")));
		assertEquals(a.flag, 1);
		assertEquals(a.strand, 0);
	}

	@Test
	public void testCSHLcombiDescriptor() {
		UniversalReadDescriptor descriptor= new UniversalReadDescriptor();
		assertTrue(descriptor.init(UniversalReadDescriptor.getDescriptor(UniversalReadDescriptor.DESCRIPTORID_MATE1_SENSE)));
	}
	
	@Test
	public void testCSHLcombiID() {
		UniversalReadDescriptor descriptor= new UniversalReadDescriptor();
		descriptor.init(UniversalReadDescriptor.getDescriptor(UniversalReadDescriptor.DESCRIPTORID_MATE1_SENSE));
		Attributes a= 
			descriptor.getAttributes(newCshlCombined, null);
		assertEquals(a.id, newCshlCombined.substring(0, newCshlCombined.lastIndexOf("/")));
		assertEquals(a.flag, 1);
		assertEquals(a.strand, 1);
		
	}
}
