package fbi.genome.io.bed;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.OutputStreamWriter;
import java.io.PipedInputStream;
import java.io.PipedOutputStream;
import java.util.Arrays;
import java.util.Comparator;

import org.junit.Test;
import org.junit.BeforeClass;
import org.junit.AfterClass;

import fbi.commons.ByteArrayCharSequence;
import fbi.commons.Execute;
import fbi.genome.io.BufferedBACSReader;
import fbi.genome.io.rna.UniversalReadDescriptor;
import fbi.genome.model.bed.BEDobject2;
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
			BEDiteratorDisk biter= new BEDiteratorDisk(pin, false, comp);
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
				BEDobject2 obj= biter.next();
				// do not exchange arguments, String.equalsTo() checks OID
				assertEquals(obj, bedSorted[i++]);
			}
			assertTrue(i== 4);
			
		} catch (Exception e) {
			; // :)
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
			Comparator<CharSequence> comp= new BEDDescriptorComparator(descriptor);
			BEDiteratorDisk biter= new BEDiteratorDisk(x, true);
			biter.init();
			
			int tstL= 0, tstC= 0;
			count= false;
			boolean first= true;
			while(biter.hasNext()) {
				BEDobject2 cs= biter.next();
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
