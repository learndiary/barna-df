package fbi.genome.io.rna;

import static org.junit.Assert.*;

import org.junit.Test;

import fbi.genome.io.rna.UniversalReadDescriptor.Attributes;

public class UniversalReadDescriptorTest {

	String descriptorIDmissing= "";
	String descriptorIDoptional= 
		UniversalReadDescriptor.SYMBOL_OPT_LEFT
		+UniversalReadDescriptor.TAG_ID
		+UniversalReadDescriptor.SYMBOL_OPT_RIGHT;
	String descriptorIDnoLeft= 
		UniversalReadDescriptor.SYMBOL_OPT_LEFT
		+UniversalReadDescriptor.TAG_STRAND
		+UniversalReadDescriptor.SYMBOL_OPT_RIGHT
		+UniversalReadDescriptor.SYMBOL_TAG_LEFT
		+UniversalReadDescriptor.TAG_ID
		+UniversalReadDescriptor.SYMBOL_TAG_RIGHT;
	String descriptorIDnoRight= 
		UniversalReadDescriptor.SYMBOL_TAG_LEFT
		+UniversalReadDescriptor.TAG_ID
		+UniversalReadDescriptor.SYMBOL_TAG_RIGHT
		+UniversalReadDescriptor.SYMBOL_OPT_LEFT
		+UniversalReadDescriptor.TAG_STRAND
		+UniversalReadDescriptor.SYMBOL_OPT_RIGHT;
	String descriptorPairOptional= 
		UniversalReadDescriptor.SYMBOL_TAG_LEFT
		+UniversalReadDescriptor.TAG_ID
		+UniversalReadDescriptor.SYMBOL_TAG_RIGHT
		+UniversalReadDescriptor.SYMBOL_OPT_LEFT
		+UniversalReadDescriptor.TAG_PAIR
		+UniversalReadDescriptor.SYMBOL_OPT_RIGHT;
	String descriptorBracketMismatch= 
		UniversalReadDescriptor.SYMBOL_TAG_LEFT
		+UniversalReadDescriptor.TAG_ID
		+UniversalReadDescriptor.SYMBOL_OPT_RIGHT;
	
	String simRead1= "chr1:4797974-4836816W:NM_008866:2:2433:548:757:548:S/1";
	String simRead2= "chr1:4797974-4836816W:NM_008866:2:2433:548:757:548:S/2";
	String simRead3= "chr1:4797974-4836816W:NM_008866:2:2433:548:757:548:A/1";
	String simRead4= "chr1:4797974-4836816W:NM_008866:2:2433:548:757:548:A/2";
	
	String barnaID1= "BILLIEHOLIDAY:5:100:1000:1190/1s";
	String barnaID2= "BILLIEHOLIDAY:5:100:1000:1190/1a";
	String barnaID3= "BILLIEHOLIDAY:5:100:1000:1190/2";
	String barnaID4= "BILLIEHOLIDAY:5:100:1000:1190/2a";
	String barnaID5= "BILLIEHOLIDAY:5:100:1000:1190/1";
	
	String oldCshlID1= "BILLIEHOLIDAY:5:100:1000:1190/1_strand2";
	String oldCshlID2= "BILLIEHOLIDAY:5:100:1000:1190/1_strand0";
	
	String newCshlCombined= "MARILYN_0005:7:1:2804:1011#0/1";


	@Test
	public void testInvalidDescriptor() {
		UniversalReadDescriptor descriptor= new UniversalReadDescriptor();
		try {
			descriptor.init(descriptorIDmissing);
			fail(descriptorIDmissing);
		} catch (Exception e) {
			; //:)
		}
		try {
			descriptor.init(descriptorIDoptional);
			fail(descriptorIDoptional);
		} catch (Exception e) {
			; //:)
		}
		try {
			descriptor.init(descriptorIDnoLeft);
			fail(descriptorIDnoLeft);
		} catch (Exception e) {
			; //:)
		}
		try {
			descriptor.init(descriptorIDnoRight);
			fail(descriptorIDnoRight);
		} catch (Exception e) {
			; //:)
		}
		try {
			descriptor.init(descriptorPairOptional);
			fail(descriptorPairOptional);
		} catch (Exception e) {
			; //:)
		}
	}
	
	@Test
	public void testSimpleDescriptor() {
		UniversalReadDescriptor descriptor= new UniversalReadDescriptor();
		try {
			descriptor.init(UniversalReadDescriptor.getDescriptor(UniversalReadDescriptor.DESCRIPTORID_SIMPLE));	
		} catch (Exception e) {
			fail(e.getMessage());
		}
		
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
		try {
			descriptor.init(UniversalReadDescriptor.getDescriptor(UniversalReadDescriptor.DESCRIPTORID_PAIRED));
		} catch (Exception e) {
			fail(e.getMessage());
		}
		
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
		try {
			descriptor.init(UniversalReadDescriptor.getDescriptor(UniversalReadDescriptor.DESCRIPTORID_SIMULATOR));
		} catch (Exception e) {
			fail(e.getMessage());
		}
		
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
		
		a= descriptor.getAttributes(simRead4, a);
		assertEquals(a.id, simRead4.substring(0, simRead4.lastIndexOf(":548:")));
		assertEquals(a.strand, 2);
		assertEquals(a.flag, 2);

		// alternative control
		descriptor= new UniversalReadDescriptor();
		descriptor.init(
				UniversalReadDescriptor.SYMBOL_TAG_LEFT+
				UniversalReadDescriptor.TAG_ID+
				UniversalReadDescriptor.SYMBOL_TAG_RIGHT+
				":???:"+
				UniversalReadDescriptor.SYMBOL_TAG_LEFT+
				UniversalReadDescriptor.TAG_STRAND+
				UniversalReadDescriptor.SYMBOL_TAG_RIGHT+
				UniversalReadDescriptor.SYMBOL_SET_LEFT+
				"S,A"+
				UniversalReadDescriptor.SYMBOL_SET_RIGHT+
				"/"+
				UniversalReadDescriptor.SYMBOL_TAG_LEFT+
				UniversalReadDescriptor.TAG_PAIR+
				UniversalReadDescriptor.SYMBOL_TAG_RIGHT
		);
		a= descriptor.getAttributes(simRead4, a);
		assertEquals(a.id, simRead4.substring(0, simRead4.lastIndexOf(":548:")));
		assertEquals(a.strand, 2);
		assertEquals(a.flag, 2);

		// negative control
		descriptor= new UniversalReadDescriptor();
		descriptor.init(
				UniversalReadDescriptor.SYMBOL_TAG_LEFT+
				UniversalReadDescriptor.TAG_ID+
				UniversalReadDescriptor.SYMBOL_TAG_RIGHT+
				":??:"+
				UniversalReadDescriptor.SYMBOL_TAG_LEFT+
				UniversalReadDescriptor.TAG_STRAND+
				UniversalReadDescriptor.SYMBOL_TAG_RIGHT+
				UniversalReadDescriptor.SYMBOL_SET_LEFT+
				"S,A"+
				UniversalReadDescriptor.SYMBOL_SET_RIGHT+
				"/"+
				UniversalReadDescriptor.SYMBOL_TAG_LEFT+
				UniversalReadDescriptor.TAG_PAIR+
				UniversalReadDescriptor.SYMBOL_TAG_RIGHT
		);
		a= descriptor.getAttributes(simRead4, a);
		assertNull(a);
	}
	
	@Test
	public void testBarnaDescriptor() {
		UniversalReadDescriptor descriptor= new UniversalReadDescriptor();
		try {
			descriptor.init(UniversalReadDescriptor.getDescriptor(UniversalReadDescriptor.DESCRIPTORID_BARNA));
		} catch (Exception e) {
			fail(e.getMessage());
		}
		
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
		try {
			descriptor.init(UniversalReadDescriptor.getDescriptor(UniversalReadDescriptor.DESCRIPTORID_MATE1_SENSE));			
		} catch (Exception e) {
			fail(e.getMessage());
		}
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
		try {
			descriptor.init(UniversalReadDescriptor.getDescriptor(UniversalReadDescriptor.DESCRIPTORID_MATE1_SENSE));
		} catch (Exception e) {
			fail(e.getMessage());
		}
		
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
