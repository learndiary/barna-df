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

package barna.io;

import barna.commons.ByteArrayCharSequence;
import org.junit.Test;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;

import static junit.framework.Assert.*;

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
