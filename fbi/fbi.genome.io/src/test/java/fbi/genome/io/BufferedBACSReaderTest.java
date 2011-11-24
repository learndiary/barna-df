package fbi.genome.io;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;

import org.junit.Test;
import static junit.framework.Assert.*;

import fbi.commons.ByteArrayCharSequence;

public class BufferedBACSReaderTest {

	@Test
    public void testReadNormal(){
		File x= new File(getClass().getResource("/pulp_fiction_zed.txt").getFile());
		try {
			BufferedReader buffy= new BufferedReader(new FileReader(x));
			int refL= 0, refC= 0;
			for(String s= null; (s= buffy.readLine())!= null;) {
				++refL;
				refC+= s.length();
			}
			buffy.close();
			
			BufferedBACSReader bacs= new BufferedBACSReader(new FileInputStream(x));
			int tstL= 0, tstC= 0;
			for(ByteArrayCharSequence cs= null; (cs= bacs.readLine(cs))!= null;) {
				++tstL;
				tstC+= cs.length();
			}
			
			assertEquals(tstL, refL);
			assertEquals(tstC, refC);
		} catch (Exception e) {
			e.printStackTrace();
			fail();
		}
		
	}
	
	@Test
    public void testMarkResetBuffer(){
    	String markedLine= "Zed's dead";
    	test(markedLine, 1, -1, true);
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
			
			BufferedBACSReader bacs= new BufferedBACSReader(new FileInputStream(x));
			int tstL= 0, tstC= 0;
			count= false;
			boolean first= true;
			for(ByteArrayCharSequence cs= null; (cs= bacs.readLine(cs))!= null;) {
				if (count) {
					++tstL;
					tstC+= cs.length();
				}
				if (first&& tstL> linesAhead&& tstC> charsAhead) { 
					long b= bacs.reset();
					if (canJump)
						assertTrue(b< 0);
					else
						assertTrue(b>= 0);
					tstL= 0;
					tstC= 0;
					first= false;
				}
				if (cs.startsWith(markedLine)&& !count) {
					count= true;
					bacs.mark();
				}
				if(cs.startsWith("Baby-love"))
					System.currentTimeMillis();
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
    	test(markedLine, 0, 1025, false);
    }

}
