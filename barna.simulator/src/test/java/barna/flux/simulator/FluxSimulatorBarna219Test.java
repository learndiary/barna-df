package barna.flux.simulator;

import barna.commons.log.Log;
import org.junit.Test;

import java.io.File;

import static org.junit.Assert.fail;

public class FluxSimulatorBarna219Test {

    @Test
    public void testBarna219(){
        // disable any questions
        Log.setInteractive(false);
        // the setting file
        File settings = new File(getClass().getResource("/BARNA-219-test-data/parameters.par").getFile());

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
    }
}
