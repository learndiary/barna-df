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

import barna.commons.Execute;
import barna.commons.log.Log;
import barna.io.FileHelper;
import barna.model.Graph;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

import java.io.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import static org.junit.Assert.*;

/**
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class SmallSimulationsTest {


    @Before
    public void setUp() throws Exception {
        Execute.initialize(4);
    }

    @After
    public void tearDown() throws Exception {
        Execute.shutdown();
    }

    @Test
    public void simpleSim1() {

        // disable any questions
        Log.setInteractive(false);
        // the setting file
        File settings = new File(getClass().getResource("/simulator_v10.par").getFile());

        SimulationPipeline pipeline = new SimulationPipeline();
        pipeline.setFile(settings);
        pipeline.setExpression(true);
        pipeline.setLibrary(true);
        pipeline.setSequence(true);

        try {
            pipeline.call();
        } catch (Exception e) {
            e.printStackTrace();
            fail();
        }

        // check read count
        assertEquals(pipeline.getSequencer().getTotalReads(), FileHelper.countLines(pipeline.getSettings().get(FluxSimulatorSettings.SEQ_FILE)));
    }
    @Test
    public void testBarna215() {

        // disable any questions
        Log.setInteractive(false);
        // the setting file
        File settings = new File(getClass().getResource("/simulator-BARNA-215.par").getFile());

        SimulationPipeline pipeline = new SimulationPipeline();
        pipeline.setFile(settings);
        pipeline.setExpression(true);
        pipeline.setLibrary(true);
        pipeline.setSequence(true);

        try {
            pipeline.call();
        } catch (Exception e) {
            //e.printStackTrace();
            fail();
        }

        // check read count
        assertEquals(pipeline.getSequencer().getTotalReads(), FileHelper.countLines(pipeline.getSettings().get(FluxSimulatorSettings.SEQ_FILE)));
    }

    @Test
    public void testRunForCustomErrorModelBarna106() {

        // disable any questions
        Log.setInteractive(false);
        // the setting file
        File settings = new File(getClass().getResource("/simulator-BARNA-106.par").getFile());

        SimulationPipeline pipeline = new SimulationPipeline();
        pipeline.setFile(settings);
        pipeline.setExpression(true);
        pipeline.setLibrary(true);
        pipeline.setSequence(true);

        try {
            pipeline.call();
        } catch (Exception e) {
            e.printStackTrace();
            fail();
        }

        // check read count
        assertEquals(pipeline.getSequencer().getTotalReads(), FileHelper.countLines(pipeline.getSettings().get(FluxSimulatorSettings.SEQ_FILE)));
    }

    @Test
    public void testRunWithCustomErrorModelAndParameters() {

        // disable any questions
        Log.setInteractive(false);
        // the setting file
        File settings = new File(getClass().getResource("/human100.par").getFile());

        SimulationPipeline pipeline = new SimulationPipeline();
        pipeline.setFile(settings);
        pipeline.setExpression(true);
        pipeline.setLibrary(true);
        pipeline.setSequence(true);

        try {
            pipeline.call();
        } catch (Exception e) {
            e.printStackTrace();
            fail();
        }

        // check read count
        assertEquals(pipeline.getSequencer().getTotalReads(), FileHelper.countLines(pipeline.getSettings().get(FluxSimulatorSettings.SEQ_FILE)));
        // we should have a fastq file
        assertTrue(new File(settings.getParentFile(), "human100.fastq").exists());
    }

    @Test
    public void testRunWithoutCustomErrorModelAndFastAOutput() {
        Graph.fileSep = null;
        // disable any questions
        Log.setInteractive(false);
        // the setting file
        File settings = new File(getClass().getResource("/human100-2.par").getFile());

        SimulationPipeline pipeline = new SimulationPipeline();
        pipeline.setFile(settings);
        pipeline.setExpression(true);
        pipeline.setLibrary(true);
        pipeline.setSequence(true);

        try {
            pipeline.call();
        } catch (Exception e) {
            e.printStackTrace();
            fail();
        }

        // check read count
        assertEquals(pipeline.getSequencer().getTotalReads(), FileHelper.countLines(pipeline.getSettings().get(FluxSimulatorSettings.SEQ_FILE)));
        // we should have a fastq file
        assertTrue(new File(settings.getParentFile(), "human100-2.fasta").exists());
    }

    @Test
    public void testRunWithoutDefaultErrorModelAndFastQOutput() {

        // disable any questions
        Log.setInteractive(false);
        // the setting file
        File settings = new File(getClass().getResource("/default-error-model.par").getFile());

        SimulationPipeline pipeline = new SimulationPipeline();
        pipeline.setFile(settings);
        pipeline.setExpression(true);
        pipeline.setLibrary(true);
        pipeline.setSequence(true);

        try {
            pipeline.call();
        } catch (Exception e) {
            e.printStackTrace();
            fail();
        }

        // check read count
        assertEquals(pipeline.getSequencer().getTotalReads(), FileHelper.countLines(pipeline.getSettings().get(FluxSimulatorSettings.SEQ_FILE)));
        // we should have a fastq file
        assertTrue(new File(settings.getParentFile(), "default-error-model.fastq").exists());
    }

    @Test
    public void testMinimalProfile() {
        // disable any questions
        Log.setInteractive(false);
        // the setting file
        File settings = new File(getClass().getResource("/minimal.par").getFile());

        SimulationPipeline pipeline = new SimulationPipeline();
        pipeline.setFile(settings);
        pipeline.setExpression(true);
        pipeline.setLibrary(false);
        pipeline.setSequence(false);

        try {
            pipeline.call();
        } catch (Exception e) {
            e.printStackTrace();
            fail();
        }

        // sum can vara
        //assertEquals(10000, pipeline.getProfiler().getSumMol());
    }


    @Test
    public void testMinimal() {
        // disable any questions
        Log.setInteractive(false);
        // the setting file
        File settings = new File(getClass().getResource("/minimal.par").getFile());

        SimulationPipeline pipeline = new SimulationPipeline();
        pipeline.setFile(settings);
        pipeline.setExpression(true);
        pipeline.setLibrary(true);
        pipeline.setSequence(true);

        try {
            pipeline.call();
        } catch (Exception e) {
            e.printStackTrace();
            fail();
        }

        // check read count
        assertEquals(pipeline.getSequencer().getTotalReads(), FileHelper.countLines(pipeline.getSettings().get(FluxSimulatorSettings.SEQ_FILE)));
    }

    @Test
    public void testStatisticsOutput() {
        // disable any questions
        Log.setInteractive(false);
        final StringWriter stringWriter = new StringWriter();
        PrintStream printStream = new PrintStream(new OutputStream() {
            @Override
            public void write(int b) throws IOException {
                stringWriter.write(b);
            }
        });
        Log.outputStream = printStream;
        Log.logStream = printStream;
        // the setting file
        File settings = new File(getClass().getResource("/test_output_stats.par").getFile());

        SimulationPipeline pipeline = new SimulationPipeline();
        pipeline.setFile(settings);
        pipeline.setExpression(true);
        pipeline.setLibrary(true);
        pipeline.setSequence(true);

        try {
            pipeline.call();
        } catch (Exception e) {
            e.printStackTrace();
            fail();
        }
        try {
            stringWriter.close();
            String allLines = stringWriter.toString();
            BufferedReader reader = new BufferedReader(new StringReader(allLines));
            long fragments = -1;
            long reads = -1;

            Pattern fragP = Pattern.compile("\\s+(\\d+).*fragments found.*");
            Pattern readP = Pattern.compile("\\s+(\\d+).*reads sequenced.*");
            for (String l = reader.readLine(); l != null; l = reader.readLine()) {
                Matcher m = fragP.matcher(l);
                System.out.println(l);
                if (m.matches()) {
                    fragments = Long.parseLong(m.group(1));
                }
                m = readP.matcher(l);
                if (m.matches()) {
                    reads = Long.parseLong(m.group(1));
                }
            }
            assertTrue(reads != -1);
            assertTrue(fragments != -1);
            //assertTrue(fragments != reads);
        } catch (IOException e) {
            e.printStackTrace();
            fail();
        }


        // check read count
        //assertEquals(pipeline.getSequencer().getTotalReads(), FileHelper.countLines(pipeline.getSettings().get(FluxSimulatorSettings.SEQ_FILE)));
    }

    @Test
    public void testpattern() throws Exception {
        Pattern fragP = Pattern.compile("\\s+(\\d+) fragments found");
        Matcher m = fragP.matcher("  123 fragments found");
        assertTrue(m.matches());
    }
}
