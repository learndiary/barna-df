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

package barna.flux.simulator;

import barna.commons.parameters.ParameterException;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.File;

import static junit.framework.Assert.*;

/**
 * Test flus simulator paramter file loading
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class FluxSimulatorSettingsTest {

    private static File par1;

    @BeforeClass
    public static void setUp(){
        par1 = new File(FluxSimulatorTest.class.getResource("/simulator_v10.par").getFile());
    }

    @Test
    public void testMinimal(){
        try {
            FluxSimulatorSettings s1 = FluxSimulatorSettings.createSettings(new File(FluxSimulatorTest.class.getResource("/minimal.par").getFile()));
            assertTrue(s1.get(FluxSimulatorSettings.LOAD_CODING));
            assertTrue(s1.get(FluxSimulatorSettings.LOAD_NONCODING));
            assertEquals(2.0d, s1.get(FluxSimulatorSettings.POLYA_SHAPE), 0.00001);

        } catch (Exception e) {
            e.printStackTrace();
            fail();
        }
    }
    @Test
    public void testTmpDirWritableChecks(){
        try {
            FluxSimulatorSettings s1 = FluxSimulatorSettings.createSettings(new File(FluxSimulatorTest.class.getResource("/test_tmp_dir.par").getFile()));
            fail();
        } catch (Exception e) {
            if(! (e instanceof ParameterException)) fail();
            StringBuilder tmpDir = new StringBuilder();
            tmpDir.append(new File(".").listRoots()[0]);
            tmpDir.append("some");
            tmpDir.append(File.separator);
            tmpDir.append("unknown");
            tmpDir.append(File.separator);
            tmpDir.append("directory");

            assertEquals("The temp-directory " + tmpDir + " does not exist or is not writable!", e.getMessage());
        }
    }

    @Test
    public void testPar1(){
        try {
            FluxSimulatorSettings s1 = FluxSimulatorSettings.createSettings(par1);

            assertNotNull(s1.get(FluxSimulatorSettings.REF_FILE));
            assertEquals(new File(FluxSimulatorTest.class.getResource("/spike_sequences.gtf").getFile()).getAbsolutePath(), s1.get(FluxSimulatorSettings.REF_FILE).getAbsolutePath());

            assertNotNull(s1.get(FluxSimulatorSettings.PRO_FILE));
            assertEquals("simulator_v10.pro", s1.get(FluxSimulatorSettings.PRO_FILE).getName());

            assertNotNull(s1.get(FluxSimulatorSettings.SEQ_FILE));
            assertEquals("simulator_v10.bed", s1.get(FluxSimulatorSettings.SEQ_FILE).getName());

            assertNotNull(s1.get(FluxSimulatorSettings.LIB_FILE));
            assertEquals("simulator_v10.lib", s1.get(FluxSimulatorSettings.LIB_FILE).getName());

            assertEquals(25d, s1.get(FluxSimulatorSettings.TSS_MEAN), 0.000001);

        } catch (Exception e) {
            e.printStackTrace();
            fail();
        }
    }
}
