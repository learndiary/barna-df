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

package barna.io.bed;

import barna.commons.ByteArrayCharSequence;
import barna.commons.Execute;
import barna.io.BufferedIteratorDisk;
import barna.io.rna.UniversalReadDescriptor;
import barna.model.bed.BEDobject2;
import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.*;
import java.util.Comparator;

import static junit.framework.Assert.*;

public class BEDiteratorDiskTest {

	@BeforeClass
	public static void init() {
		Execute.initialize(5);
	}
	
	@AfterClass
	public static void shutdown() {
		Execute.shutdown();
	}
	
	String[] beds= new String[] {
			"chr1\t6349539\t6378672\tchr1:6349412-6384812W:NM_001195732:58:858:129:214:S/2\t0\t+\t.\t.\t0,0,0\t56,19\t0,29114",
			"chr1\t6349550\t6378683\tchr1:6349412-6384812W:NM_001195732:18:858:129:214:A/2\t0\t-\t.\t.\t0,0,0\t45,30\t0,29103",
			"chr1\t6349539\t6378672\tchr1:6349412-6384812W:NM_001195732:18:858:129:214:S/1\t0\t+\t.\t.\t0,0,0\t56,19\t0,29114",
			"chr1\t6349550\t6378683\tchr1:6349412-6384812W:NM_001195732:58:858:129:214:A/1\t0\t-\t.\t.\t0,0,0\t45,30\t0,29103"
	};
	
	String[] bedSorted= new String[] {
			"chr1\t6349539\t6378672\tchr1:6349412-6384812W:NM_001195732:18:858:129:214:S/1\t0\t+\t.\t.\t0,0,0\t56,19\t0,29114",
			"chr1\t6349550\t6378683\tchr1:6349412-6384812W:NM_001195732:18:858:129:214:A/2\t0\t-\t.\t.\t0,0,0\t45,30\t0,29103",
			"chr1\t6349550\t6378683\tchr1:6349412-6384812W:NM_001195732:58:858:129:214:A/1\t0\t-\t.\t.\t0,0,0\t45,30\t0,29103",
			"chr1\t6349539\t6378672\tchr1:6349412-6384812W:NM_001195732:58:858:129:214:S/2\t0\t+\t.\t.\t0,0,0\t56,19\t0,29114"
	};
	@Test
	public void streamSortTest() {
		
		try {
			
			PipedInputStream pin= new PipedInputStream();
			PipedOutputStream pout= new PipedOutputStream(pin);
			
			UniversalReadDescriptor descriptor= new UniversalReadDescriptor();
			descriptor.init(UniversalReadDescriptor.getDescriptor(UniversalReadDescriptor.DESCRIPTORID_SIMULATOR));
			Comparator<CharSequence> comp= new BEDDescriptorComparator(descriptor);
			//Arrays.sort(beds, comp);
			//for (int i = 0; i < beds.length; i++) {
			//	System.out.println(beds[i]);
			//}
			File tmpFile= File.createTempFile(getClass().getSimpleName(), "bed");
			tmpFile.deleteOnExit();
			BufferedIteratorDisk biter= new BufferedIteratorDisk(pin, tmpFile, comp);
			biter.init();
			
			OutputStreamWriter writer= new OutputStreamWriter(pout);
			for (int i = 0; i < beds.length; i++) {
				writer.write(beds[i]);
				writer.write('\n');
			}
			writer.flush();
			writer.close();
			
			int i= 0;
			while(biter.hasNext()) {
				ByteArrayCharSequence cs= biter.next();
				BEDobject2 obj= new BEDobject2(cs);
				// do not exchange arguments, String.equalsTo() checks OID
				assertEquals(obj, bedSorted[i++]);
			}
			assertTrue(i== 4);
			tmpFile.delete();
			
		} catch (Exception e) {
			throw new RuntimeException(e); 
		}
	}
	
	private void test(String markedLine, int linesAhead, int charsAhead, boolean canJump) {
		File x= new File(getClass().getResource("/pulp_fiction_zed.txt").getFile());
		try {
			BufferedReader buffy= new BufferedReader(new FileReader(x));
			int refL= 0, refC= 0;
			boolean count= false;
			for(String s= null; (s= buffy.readLine())!= null;) {
				if (count) {
					++refL;
					refC+= s.length();
				}
				if (s.startsWith(markedLine))
					count= true;
			}
			buffy.close();
			
			UniversalReadDescriptor descriptor= new UniversalReadDescriptor();
			descriptor.init(UniversalReadDescriptor.getDescriptor(UniversalReadDescriptor.DESCRIPTORID_SIMULATOR));
			BufferedIteratorDisk biter= new BufferedIteratorDisk(x);
			biter.init();
			
			int tstL= 0, tstC= 0;
			count= false;
			boolean first= true;
			while(biter.hasNext()) {
				ByteArrayCharSequence cseq= biter.next();
				BEDobject2 cs= new BEDobject2(cseq);
				if (count) {
					++tstL;
					tstC+= cs.length();
				}
				if (first&& tstL> linesAhead&& tstC> charsAhead) { 
					biter.reset();
					tstL= 0;
					tstC= 0;
					first= false;
				}
				if (cs.startsWith(markedLine)&& !count) {
					count= true;
					biter.mark();
				}
			}
			if (canJump) {
				assertEquals(refL, tstL);
				assertEquals(refC, tstC);
			} else {
				assertFalse(tstL== refL);
				assertFalse(tstC== refC);
			}
		} catch (Exception e) {
			e.printStackTrace();
			fail();
		}
    }
	
	@Test
    public void testMarkResetSweep(){
    	String markedLine= "It's a chopper";
    	test(markedLine, 0, 1025, true);
    }
	
	@Test
    public void testMarkResetMem(){
    	String markedLine= "It's a chopper";
    	test(markedLine, 5, 0, true);
    }

}
