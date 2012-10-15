package barna.flux.simulator.fragmentation;

import barna.commons.log.Log;
import barna.flux.simulator.SimulationPipeline;
import barna.flux.simulator.distributions.AbstractDistribution;
import barna.flux.simulator.distributions.Distributions;
import barna.flux.simulator.distributions.EmpiricalDistribution;
import barna.flux.simulator.distributions.NormalDistribution;
import org.junit.Test;

import java.io.File;

import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.fail;

public class FluxSimulatorBarna235Test {

    @Test
    public void testDistributionInitialization(){
        long currentNumberOfFragments = 1000;
        String tmpFilePath =new File(getClass().getResource("/BARNA-235-test-data/frags.txt").getFile()).getAbsolutePath();
        AbstractDistribution dist = Distributions.parseDistribution("N(300,30)");
        assertNotNull(dist);
        int GEL_NB_BINS_LENGTH = 100;
        EmpiricalDistribution filterDist = new EmpiricalDistribution((NormalDistribution) dist, GEL_NB_BINS_LENGTH, 4d );
        EmpiricalDistribution originalDist = Fragmenter.readFragmentSizeDistribution(tmpFilePath, currentNumberOfFragments, filterDist.getMin(), filterDist.getMax());
        assertNotNull(filterDist);
        filterDist.normalizeToPrior(originalDist);
    }
}
