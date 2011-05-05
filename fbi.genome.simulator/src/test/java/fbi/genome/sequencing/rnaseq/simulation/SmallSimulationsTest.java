package fbi.genome.sequencing.rnaseq.simulation;

import fbi.commons.Log;
import fbi.commons.file.FileHelper;
import org.junit.Test;

import java.io.File;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.fail;

/**
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class SmallSimulationsTest {

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
