/*
 * Copyright (c) 2012, Micha Sammeth, Thasso Griebel, Emilio Palumbo
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *      * Redistributions of source code must retain the above copyright
 *        notice, this list of conditions and the following disclaimer.
 *      * Redistributions in binary form must reproduce the above copyright
 *        notice, this list of conditions and the following disclaimer in the
 *        documentation and/or other materials provided with the distribution.
 *      * The names of its contributors may be not used to endorse or promote
 *        products derived from this software without specific prior written
 *        permission.
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

import barna.commons.log.Log;
import org.junit.*;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.security.DigestInputStream;
import java.security.MessageDigest;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.fail;

/**
 * Test adding a SEED and ensures that we create the same results for two independent runs
 */
public class FluxSimulatorBarna306Test {

    private static File pro1;
    private static File pro2;
    private static File lib1;
    private static File lib2;
    private static File seq1;
    private static File seq2;
    private static File settings;

    @BeforeClass
    public static void setUp() throws Exception {
        pro1 = File.createTempFile("BARNA-test-306", ".pro");
        pro2 = File.createTempFile("BARNA-test-306", ".pro");
        lib1 = File.createTempFile("BARNA-test-306", ".lib");
        lib2 = File.createTempFile("BARNA-test-306", ".lib");
        seq1 = File.createTempFile("BARNA-test-306", ".seq");
        seq2 = File.createTempFile("BARNA-test-306", ".seq");

        // disable any questions
        Log.setInteractive(false);
        // the setting file
        settings = new File(FluxSimulatorBarna306Test.class.getResource("/simulator-BARNA-306.par").getFile());
    }

    @AfterClass
    public static void tearDown() throws Exception {
        pro1.deleteOnExit();
        pro2.deleteOnExit();
        lib1.deleteOnExit();
        lib2.deleteOnExit();
        seq1.deleteOnExit();
        seq2.deleteOnExit();
    }

    private static byte[] md5(File file){
        InputStream is = null;
        MessageDigest md = null;
        try {
            md = MessageDigest.getInstance("MD5");
            is = new FileInputStream(file);
            is = new DigestInputStream(is, md);
            while(is.read() != -1){}
        }catch (Exception e){
            e.printStackTrace();
            fail();
        }finally {
            try {
                is.close();
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
        return md.digest();
    }


    @Test
    public void testProfilerEqual(){

        try {
            // setup output files
            File pro1 = File.createTempFile("BARNA-test-306", ".pro");
            File pro2 = File.createTempFile("BARNA-test-306", ".pro");

            pro1.deleteOnExit();
            pro2.deleteOnExit();

            // first run
            SimulationPipeline pipeline = new SimulationPipeline();
            pipeline.setFile(settings);
            pipeline.getSettings().set(FluxSimulatorSettings.SEED, 1L);
            pipeline.getSettings().set(FluxSimulatorSettings.PRO_FILE, pro1);
            pipeline.setExpression(true);
            pipeline.call();

            // second run
            pipeline = new SimulationPipeline();
            pipeline.setFile(settings);
            pipeline.getSettings().set(FluxSimulatorSettings.SEED, 1L);
            pipeline.getSettings().set(FluxSimulatorSettings.PRO_FILE, pro2);
            pipeline.setExpression(true);
            pipeline.call();
            assertArrayEquals(md5(pro1), md5(pro2));
        } catch (Exception e) {
            e.printStackTrace();
            fail();
        }

    }

    @Test
    public void testLibraryEquals(){
        try {
            // setup output files

            // first run
            SimulationPipeline pipeline = new SimulationPipeline();
            pipeline.setFile(settings);
            pipeline.getSettings().set(FluxSimulatorSettings.SEED, 1L);
            pipeline.getSettings().set(FluxSimulatorSettings.PRO_FILE, pro1);
            pipeline.getSettings().set(FluxSimulatorSettings.LIB_FILE, lib1);
            pipeline.getSettings().set(FluxSimulatorSettings.SEQ_FILE, seq1);
            pipeline.setExpression(true);
            pipeline.setLibrary(true);
//            pipeline.setSequence(true);
            pipeline.call();

            // second run
            pipeline = new SimulationPipeline();
            pipeline.setFile(settings);
            pipeline.getSettings().set(FluxSimulatorSettings.SEED, 1L);
            pipeline.getSettings().set(FluxSimulatorSettings.PRO_FILE, pro2);
            pipeline.getSettings().set(FluxSimulatorSettings.LIB_FILE, lib2);
            pipeline.getSettings().set(FluxSimulatorSettings.SEQ_FILE, seq2);
            pipeline.setExpression(true);
            pipeline.setLibrary(true);
//            pipeline.setSequence(true);
            pipeline.call();

            assertArrayEquals(md5(pro1), md5(pro2));
            assertArrayEquals(md5(lib1), md5(lib2));
        } catch (Exception e) {
            e.printStackTrace();
            fail();
        }

    }

    @Test
    public void testFullEqualRun(){
        try {
            // setup output files

            // first run
            SimulationPipeline pipeline = new SimulationPipeline();
            pipeline.setFile(settings);
            pipeline.getSettings().set(FluxSimulatorSettings.SEED, 1L);
            pipeline.getSettings().set(FluxSimulatorSettings.PRO_FILE, pro1);
            pipeline.getSettings().set(FluxSimulatorSettings.LIB_FILE, lib1);
            pipeline.getSettings().set(FluxSimulatorSettings.SEQ_FILE, seq1);
            pipeline.setExpression(true);
            pipeline.setLibrary(true);
            pipeline.setSequence(true);
            pipeline.call();

            // second run
            pipeline = new SimulationPipeline();
            pipeline.setFile(settings);
            pipeline.getSettings().set(FluxSimulatorSettings.SEED, 1L);
            pipeline.getSettings().set(FluxSimulatorSettings.PRO_FILE, pro2);
            pipeline.getSettings().set(FluxSimulatorSettings.LIB_FILE, lib2);
            pipeline.getSettings().set(FluxSimulatorSettings.SEQ_FILE, seq2);
            pipeline.setExpression(true);
            pipeline.setLibrary(true);
            pipeline.setSequence(true);
            pipeline.call();

            assertArrayEquals(md5(pro1), md5(pro2));
            assertArrayEquals(md5(lib1), md5(lib2));
            assertArrayEquals(md5(seq1), md5(seq2));
        } catch (Exception e) {
            e.printStackTrace();
            fail();
        }

    }
}
