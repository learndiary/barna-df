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

package barna.flux.simulator.error;

import barna.model.Qualities;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.File;
import java.io.IOException;

import static junit.framework.Assert.*;

/**
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class MapFileReaderTest {

    private static File testfile;

    @BeforeClass
    public static void setUp() {
        testfile = new File(MapFileReaderTest.class.getResource("/test.map").getFile());
    }

    @Test
    public void testReadOne() {
        MapFileReader reader = new MapFileReader(testfile, Qualities.Technology.Phred);
        try {
            Read read1 = reader.parseNext(false);
            assertNotNull(read1);

            assertEquals("HWUSI-EAS627_1:2:1:4:1299/1", read1.getName());
            assertEquals("ATTTTNNTCAAAAACTTTGTCTTTTTTTCTTTCCTCCCCTAAATTTTCCCCAATTTAAATTTTTCCCCCAGGGGTC", read1.getSequence());
            assertEquals(76, read1.getLength());
        } catch (IOException e) {
            e.printStackTrace();
            fail();
        } finally {
            reader.close();
        }
    }

    @Test
    public void testReadMappings() {
        MapFileReader reader = new MapFileReader(testfile, Qualities.Technology.Phred);
        try {
            Read read1 = reader.parseNext(false);
            assertNotNull(read1);

            read1 = reader.parseNext(false);
            // second read
            // should have one perfect mapping without missmatches
            assertNotNull(read1.getMappings());
            assertEquals(1, read1.getMappings().size());
            assertEquals("chr8:R11218885", read1.getMappings().get(0).getName());

            assertNull(read1.getMappings().get(0).getMissmatches());

            read1 = reader.parseNext(false);
            // third read
            // should have one mapping with 8 missmatches
            assertNotNull(read1.getMappings());
            assertEquals(3, read1.getMappings().size());
            assertEquals(8, read1.getMappings().get(0).getMissmatches().size());
            assertEquals("chr12:F10906850", read1.getMappings().get(0).getName());


        } catch (IOException e) {
            e.printStackTrace();
            fail();
        } finally {
            reader.close();
        }
    }
}
