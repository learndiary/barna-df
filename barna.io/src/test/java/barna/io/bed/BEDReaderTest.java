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

import barna.commons.Execute;
import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.File;
import java.io.IOException;

import static junit.framework.Assert.*;

/**
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class BEDReaderTest {

    private static File testfile;

    @BeforeClass
    public static void setUp(){
        testfile = new File(BEDReaderTest.class.getResource("/test1.bed").getFile());
        Execute.initialize(4);
    }

    @AfterClass
    public static void tearDown() throws Exception {
        Execute.shutdown();
    }


    @Test
    public void testScanFile(){
        BEDReader wrapper = new BEDReader(new File(getClass().getResource("/test.bed").getFile()));
        wrapper.scanFile();

        //scanFileReadLines= 0;
        //countAll= 0; countEntire= 0; countSplit= 0; countReads= 0;
        assertEquals(17, wrapper.nrUniqueLinesRead );
        assertEquals(17, wrapper.countAll );
        assertEquals(5, wrapper.countSplit );
        assertEquals(12, wrapper.countEntire );
        assertEquals(17, wrapper.countReads );
    }

    @Test
    public void testReadDescriptorWithSpace(){
        BEDReader wrapper = new BEDReader(testfile.getAbsolutePath());
        wrapper.scanFile();

        //scanFileReadLines= 0;
        //countAll= 0; countEntire= 0; countSplit= 0; countReads= 0;
        assertEquals(1000, wrapper.nrUniqueLinesRead );
        assertEquals(1000, wrapper.countAll );
        assertEquals(75, wrapper.countSplit );
        assertEquals(925, wrapper.countEntire );
        assertEquals(972, wrapper.countReads );
    }

    @Test
    public void testIsApplicable() throws IOException {
        BEDReader wrapper = new BEDReader(testfile.getAbsolutePath());
        assertFalse(wrapper.isApplicable());
        File bedtest = File.createTempFile("bedtest", ".bed");
        wrapper.sort(bedtest);
        wrapper = new BEDReader(bedtest);
        assertTrue(wrapper.isApplicable());
        bedtest.delete();
    }
}
