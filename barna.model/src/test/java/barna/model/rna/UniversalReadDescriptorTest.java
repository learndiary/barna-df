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

package barna.model.rna;

import barna.model.rna.UniversalReadDescriptor.Attributes;
import org.junit.Test;

import static org.junit.Assert.*;

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
	
	static String simRead1= "chr1:4797974-4836816W:NM_008866:2:2433:548:757:S/1";
	static String simRead2= "chr1:4797974-4836816W:NM_008866:2:2433:548:757:S/2";
	static String simRead3= "chr1:4797974-4836816W:NM_008866:2:2433:548:757:A/1";
	static String simRead4= "chr1:4797974-4836816W:NM_008866:2:2433:548:757:A/2";
	
	static String barnaID1= "BILLIEHOLIDAY:5:100:1000:1190/1s";
	static String barnaID2= "BILLIEHOLIDAY:5:100:1000:1190/1a";
	static String barnaID3= "BILLIEHOLIDAY:5:100:1000:1190/2";
	static String barnaID4= "BILLIEHOLIDAY:5:100:1000:1190/2a";
	static String barnaID5= "BILLIEHOLIDAY:5:100:1000:1190/1";
	
	static String oldCshlID1= "BILLIEHOLIDAY:5:100:1000:1190/1_strand2";
	static String oldCshlID2= "BILLIEHOLIDAY:5:100:1000:1190/1_strand0";
	
	static String newCshlCombined= "MARILYN_0005:7:1:2804:1011#0/1";

	static String genevaRead1= "@HWI-ST661:130:C037KACXX:6:1101:1483:2196 1:N:0:CAGATC";
	static String genevaRead2= "@HWI-ST661:130:C037KACXX:6:1101:1483:2196 2:N:0:CAGATC";
	

	@Test
	public void testInvalidDescriptor() {
		UniversalReadDescriptor descriptor= UniversalReadDescriptor.createTestDescriptor();
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
		try {
            UniversalReadDescriptor descriptor= new UniversalReadDescriptor(
                    UniversalReadDescriptor.getDescriptor(UniversalReadDescriptor.DESCRIPTORID_SIMPLE));
		} catch (Exception e) {
			fail(e.getMessage());
		}
		
	}
	
	@Test
	public void testSimpleID() {
        UniversalReadDescriptor descriptor= new UniversalReadDescriptor(
                UniversalReadDescriptor.getDescriptor(UniversalReadDescriptor.DESCRIPTORID_SIMPLE));
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
		try {
            UniversalReadDescriptor descriptor= new UniversalReadDescriptor(
                    UniversalReadDescriptor.getDescriptor(UniversalReadDescriptor.DESCRIPTORID_PAIRED));
		} catch (Exception e) {
			fail(e.getMessage());
		}
		
	}
	
	@Test
	public void testPairedID() {

        UniversalReadDescriptor descriptor= new UniversalReadDescriptor(
                UniversalReadDescriptor.getDescriptor(UniversalReadDescriptor.DESCRIPTORID_PAIRED));
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
		try {
            UniversalReadDescriptor descriptor= new UniversalReadDescriptor(
                    UniversalReadDescriptor.getDescriptor(UniversalReadDescriptor.DESCRIPTORID_SIMULATOR));
		} catch (Exception e) {
			fail(e.getMessage());
		}
		
	}

    @Test
    public void testSenseDescriptor() {
        try {
            UniversalReadDescriptor descriptor= new UniversalReadDescriptor(
                    UniversalReadDescriptor.getDescriptor(UniversalReadDescriptor.DESCRIPTORID_SENSE));
        } catch (Exception e) {
            fail(e.getMessage());
        }
    }

    @Test
    public void testSenseID() {
        UniversalReadDescriptor descriptor= new UniversalReadDescriptor(
                UniversalReadDescriptor.getDescriptor(UniversalReadDescriptor.DESCRIPTORID_SENSE));
        Attributes a=
                descriptor.getAttributes(simRead1, null);
        assertEquals(a.id, simRead1);
        assertEquals(1, a.strand);

        a= descriptor.getAttributes(simRead2, a);
        assertEquals(a.id, simRead2);
        assertEquals(1, a.strand);

        a= descriptor.getAttributes(simRead3, a);
        assertEquals(a.id, simRead3);
        assertEquals(1, a.strand);

        a= descriptor.getAttributes(simRead4, a);
        assertEquals(a.id, simRead4);
        assertEquals(1, a.strand);
    }

    @Test
    public void testAntisenseID() {
        UniversalReadDescriptor descriptor= new UniversalReadDescriptor(UniversalReadDescriptor.getDescriptor(UniversalReadDescriptor.DESCRIPTORID_ANTISENSE));
        Attributes a=
                descriptor.getAttributes(simRead1, null);
        assertEquals(a.id, simRead1);
        assertEquals(2, a.strand);

        a= descriptor.getAttributes(simRead2, a);
        assertEquals(a.id, simRead2);
        assertEquals(2, a.strand);

        a= descriptor.getAttributes(simRead3, a);
        assertEquals(a.id, simRead3);
        assertEquals(2, a.strand);

        a= descriptor.getAttributes(simRead4, a);
        assertEquals(a.id, simRead4);
        assertEquals(2, a.strand);
    }

    @Test
    public void testAntisenseDescriptor() {
        UniversalReadDescriptor descriptor= UniversalReadDescriptor.createTestDescriptor();
        try {
            descriptor.init(UniversalReadDescriptor.getDescriptor(UniversalReadDescriptor.DESCRIPTORID_ANTISENSE));
        } catch (Exception e) {
            fail(e.getMessage());
        }
    }

    @Test
	public void testSimulatorID() {
        UniversalReadDescriptor descriptor= new UniversalReadDescriptor(
                UniversalReadDescriptor.getDescriptor(UniversalReadDescriptor.DESCRIPTORID_SIMULATOR));
		Attributes a= 
			descriptor.getAttributes(simRead1, null);
		assertEquals(a.id, simRead1.substring(0, simRead1.lastIndexOf(":")));
		assertEquals(a.strand, 1);
		assertEquals(a.flag, 1);
		
		a= descriptor.getAttributes(simRead2, a);
		assertEquals(a.id, simRead2.substring(0, simRead2.lastIndexOf(":")));
		assertEquals(a.strand, 1);
		assertEquals(a.flag, 2);
		
		a= descriptor.getAttributes(simRead3, a);
		assertEquals(a.id, simRead3.substring(0, simRead3.lastIndexOf(":")));
		assertEquals(a.strand, 2);
		assertEquals(a.flag, 1);
		
		a= descriptor.getAttributes(simRead4, a);
		assertEquals(a.id, simRead4.substring(0, simRead4.lastIndexOf(":")));
		assertEquals(a.strand, 2);
		assertEquals(a.flag, 2);

		// alternative control, used for question mark operator before
		descriptor= new UniversalReadDescriptor(
				UniversalReadDescriptor.SYMBOL_TAG_LEFT+
				UniversalReadDescriptor.TAG_ID+
				UniversalReadDescriptor.SYMBOL_TAG_RIGHT+
				":"+
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
		assertEquals(a.id, simRead4.substring(0, simRead4.lastIndexOf(":")));
		assertEquals(a.strand, 2);
		assertEquals(a.flag, 2);

		// negative control
		descriptor= new UniversalReadDescriptor(
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
		UniversalReadDescriptor descriptor= UniversalReadDescriptor.createTestDescriptor();
		try {
			descriptor.init(UniversalReadDescriptor.getDescriptor(UniversalReadDescriptor.DESCRIPTORID_BARNA));
		} catch (Exception e) {
			fail(e.getMessage());
		}
		
	}
	
	@Test
	public void testBarnaID() {
		UniversalReadDescriptor descriptor= new UniversalReadDescriptor(UniversalReadDescriptor.DESCRIPTORID_BARNA);
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
		try {
            UniversalReadDescriptor descriptor= new UniversalReadDescriptor(UniversalReadDescriptor.getDescriptor(UniversalReadDescriptor.DESCRIPTORID_MATE1_SENSE));
		} catch (Exception e) {
			fail(e.getMessage());
		}
	}
	
	@Test
	public void testCSHLoldID() {

        UniversalReadDescriptor descriptor= new UniversalReadDescriptor(UniversalReadDescriptor.getDescriptor(UniversalReadDescriptor.DESCRIPTORID_MATE_STRAND_CSHL));
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
		UniversalReadDescriptor descriptor= UniversalReadDescriptor.createTestDescriptor();
		try {
			descriptor.init(UniversalReadDescriptor.getDescriptor(UniversalReadDescriptor.DESCRIPTORID_MATE1_SENSE));
		} catch (Exception e) {
			fail(e.getMessage());
		}
		
	}
	
	@Test
	public void testCSHLcombiID() {
		UniversalReadDescriptor descriptor= new UniversalReadDescriptor(UniversalReadDescriptor.DESCRIPTORID_MATE1_SENSE);
		Attributes a=
			descriptor.getAttributes(newCshlCombined, null);
		assertEquals(a.id, newCshlCombined.substring(0, newCshlCombined.lastIndexOf("/")));
		assertEquals(a.flag, 1);
		assertEquals(a.strand, 1);
		
	}
	
	@Test
	public void testToString() {
		try {
			UniversalReadDescriptor descriptor= UniversalReadDescriptor.createTestDescriptor();
			
			String expr= UniversalReadDescriptor.getDescriptor(UniversalReadDescriptor.DESCRIPTORID_SIMPLE);
			descriptor.init(expr);
			assertEquals(expr, descriptor.toString());

			descriptor= UniversalReadDescriptor.createTestDescriptor();
			expr= UniversalReadDescriptor.getDescriptor(UniversalReadDescriptor.DESCRIPTORID_PAIRED);
			descriptor.init(expr);
			assertEquals(expr, descriptor.toString());

			descriptor= UniversalReadDescriptor.createTestDescriptor();
			expr= UniversalReadDescriptor.getDescriptor(UniversalReadDescriptor.DESCRIPTORID_STRAND_MATE);
			descriptor.init(expr);
			assertEquals(expr, descriptor.toString());
			
			descriptor= UniversalReadDescriptor.createTestDescriptor();
			expr= UniversalReadDescriptor.getDescriptor(UniversalReadDescriptor.DESCRIPTORID_MATE_STRAND_CSHL);
			descriptor.init(expr);
			assertEquals(expr, descriptor.toString());

			descriptor= UniversalReadDescriptor.createTestDescriptor();
			expr= UniversalReadDescriptor.getDescriptor(UniversalReadDescriptor.DESCRIPTORID_MATE1_SENSE);
			descriptor.init(expr);
			assertEquals(expr, descriptor.toString());

			descriptor= UniversalReadDescriptor.createTestDescriptor();
			expr= UniversalReadDescriptor.getDescriptor(UniversalReadDescriptor.DESCRIPTORID_MATE2_SENSE);
			descriptor.init(expr);
			assertEquals(expr, descriptor.toString());

			descriptor= UniversalReadDescriptor.createTestDescriptor();
			expr= UniversalReadDescriptor.getDescriptor(UniversalReadDescriptor.DESCRIPTORID_BARNA);
			descriptor.init(expr);
			assertEquals(expr, descriptor.toString());
			
			descriptor= UniversalReadDescriptor.createTestDescriptor();
			expr= UniversalReadDescriptor.getDescriptor(UniversalReadDescriptor.DESCRIPTORID_SIMULATOR);
			descriptor.init(expr);
			assertEquals(expr, descriptor.toString());


		} catch (Exception e) {
			fail(e.getMessage());
		}
	}
	
	@Test
	public void testGenevaID() {
		
		// @HWI-ST661:130:C037KACXX:6:1101:1483:2196 2:N:0:CAGATC
		String descriptorGeneva= 
				UniversalReadDescriptor.SYMBOL_TAG_LEFT
				+UniversalReadDescriptor.TAG_ID
				+UniversalReadDescriptor.SYMBOL_TAG_RIGHT
				+" "
				+UniversalReadDescriptor.SYMBOL_TAG_LEFT
				+UniversalReadDescriptor.TAG_PAIR
				+UniversalReadDescriptor.SYMBOL_TAG_RIGHT;

		
		UniversalReadDescriptor descriptor= new UniversalReadDescriptor(descriptorGeneva);

		Attributes a= null;
		a= descriptor.getAttributes(genevaRead1, a);
		String id= genevaRead1.split(" ")[0];
		assertEquals(id, a.id);
		assertEquals(1, a.flag);
		
		a= descriptor.getAttributes(genevaRead2, a);
		id= genevaRead2.split(" ")[0];
		assertEquals(id, a.id);
		assertEquals(2, a.flag);
		
	}
}
