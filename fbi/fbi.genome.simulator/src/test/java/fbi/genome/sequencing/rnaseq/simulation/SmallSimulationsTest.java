package fbi.genome.sequencing.rnaseq.simulation;

import fbi.commons.Execute;
import fbi.commons.Log;
import fbi.genome.io.FileHelper;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

import java.io.File;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.fail;

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

}
