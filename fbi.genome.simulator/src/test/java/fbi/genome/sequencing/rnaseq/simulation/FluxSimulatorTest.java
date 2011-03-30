package fbi.genome.sequencing.rnaseq.simulation;

import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.File;

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
    public void testSimpleRun1(){
        FluxSimulator.main(new String[]{});
    }


}
