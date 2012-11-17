package barna.flux.simulator;

import barna.commons.log.Log;
import barna.flux.simulator.distributions.AbstractDistribution;
import barna.flux.simulator.distributions.Distributions;
import barna.flux.simulator.distributions.EmpiricalDistribution;
import barna.flux.simulator.distributions.NormalDistribution;
import barna.flux.simulator.fragmentation.Fragmenter;
import org.junit.Test;

import java.io.File;

import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.fail;

public class FluxSimulatorBarna235Test {

    @Test
    public void testBarna235(){
        // disable any questions
        Log.setInteractive(false);
        // the setting file
        File settings = new File(getClass().getResource("/simulator-BARNA-235.par").getFile());

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
    }
}
