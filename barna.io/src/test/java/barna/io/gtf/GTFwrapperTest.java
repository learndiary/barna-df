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

package barna.io.gtf;

import org.junit.Test;

import java.io.File;

import static junit.framework.Assert.*;

public class GTFwrapperTest {

    @Test
    public void testThatTheWrapperWorksWithGzippedFiles() throws Exception {
        File gzippedGtf = new File(getClass().getResource("/gzipped-gtf.gtf.gz").getFile());
        assertTrue(gzippedGtf.exists());

        GTFwrapper wrapperFile = new GTFwrapper(gzippedGtf);
        assertFalse(wrapperFile.isApplicable());

        File sorted = wrapperFile.sort();
        assertNotNull(sorted);
        assertTrue(new GTFwrapper(sorted).isApplicable());
        sorted.delete();
    }
    @Test
    public void testThatNormaltabsForFirstFieldsAreApplied() throws Exception {
        File gzippedGtf = new File(getClass().getResource("/valid-tab-space.gtf").getFile());
        assertTrue(gzippedGtf.exists());
        GTFwrapper wrapperFile = new GTFwrapper(gzippedGtf);
        File sorted = wrapperFile.sort();
        assertNotNull(sorted);
        assertTrue(new GTFwrapper(sorted).isApplicable());
        sorted.delete();
    }
    @Test
    public void testThatOnlySpaceSeparatedFields() throws Exception {
        File gzippedGtf = new File(getClass().getResource("/only-space.gtf").getFile());
        assertTrue(gzippedGtf.exists());
        GTFwrapper wrapperFile = new GTFwrapper(gzippedGtf);
        File sorted = wrapperFile.sort();
        assertNotNull(sorted);
        assertTrue(new GTFwrapper(sorted).isApplicable());
        sorted.delete();
    }
    @Test
    public void testMixedFileWithTabsAndSpaces() throws Exception {
        File gzippedGtf = new File(getClass().getResource("/mixed-space.gtf").getFile());
        assertTrue(gzippedGtf.exists());
        GTFwrapper wrapperFile = new GTFwrapper(gzippedGtf);
        File sorted = wrapperFile.sort();
        assertNotNull(sorted);
        assertTrue(new GTFwrapper(sorted).isApplicable());
        sorted.delete();
    }
    @Test
    public void testWithTabsInTheGTFFieldsAtTheEnd() throws Exception {
        File gzippedGtf = new File(getClass().getResource("/tab-fields.gtf").getFile());
        assertTrue(gzippedGtf.exists());
        GTFwrapper wrapperFile = new GTFwrapper(gzippedGtf);
        File sorted = wrapperFile.sort();
        assertNotNull(sorted);
        assertTrue(new GTFwrapper(sorted).isApplicable());
        sorted.delete();
    }
    @Test
    public void testThatInvalidGTFFailsWithReport() throws Exception {
        File gzippedGtf = new File(getClass().getResource("/invalid.gtf").getFile());
        assertTrue(gzippedGtf.exists());
        GTFwrapper wrapperFile = new GTFwrapper(gzippedGtf);
        try{
            File sorted = wrapperFile.sort();
            fail();
        }catch (Exception error){
        }
    }
}
