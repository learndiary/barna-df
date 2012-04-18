package barna.flux.simulator;

import barna.commons.Execute;
import barna.commons.log.Log;
import barna.io.FileHelper;
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
    public void simpleSim1(){

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
    public void testRunForCustomErrorModelBarna106(){

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
    public void testRunWithCustomErrorModelAndParameters(){

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
    public void testRunWithoutCustomErrorModelAndFastAOutput(){

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
    public void testMinimalProfile(){
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
    public void testMinimal(){
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
    public void testStatisticsOutput(){
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
            long fragments =-1;
            long reads =-1;

            Pattern fragP = Pattern.compile("\\s+(\\d+).*fragments found.*");
            Pattern readP = Pattern.compile("\\s+(\\d+).*reads sequenced.*");
            for(String l = reader.readLine(); l != null; l = reader.readLine()){
                Matcher m = fragP.matcher(l);
                System.out.println(l);
                if(m.matches()){
                    fragments = Long.parseLong(m.group(1));
                }
                m = readP.matcher(l);
                if(m.matches()){
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
        System.out.println(m.group(1));

    }
}
