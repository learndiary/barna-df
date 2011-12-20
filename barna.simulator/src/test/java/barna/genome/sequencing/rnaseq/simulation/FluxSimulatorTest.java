package barna.genome.sequencing.rnaseq.simulation;

import barna.commons.launcher.Flux;
import barna.commons.launcher.FluxTool;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.File;
import java.util.List;

import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

/**
 * Test basic functions
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class FluxSimulatorTest {

    private static File sacCer2_sgdGene_fromUCSC110329;

    @BeforeClass
    public static void setup(){
        sacCer2_sgdGene_fromUCSC110329 = new File(FluxSimulatorTest.class.getResource("/sacCer2_sgdGene_fromUCSC110329.gtf").getFile());
    }

    @Test
    public void testFindTools(){
        List<FluxTool> tools = Flux.findTools();
        assertNotNull(tools);
        assertTrue(tools.size() > 2);
    }
}
