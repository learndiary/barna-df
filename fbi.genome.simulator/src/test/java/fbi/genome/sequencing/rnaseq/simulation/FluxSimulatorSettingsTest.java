package fbi.genome.sequencing.rnaseq.simulation;

import org.junit.BeforeClass;
import org.junit.Test;

import java.io.File;

import static junit.framework.Assert.*;

/**
 * Test flus simulator paramter file loading
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class FluxSimulatorSettingsTest {

    private static File par1;

    @BeforeClass
    public static void setUp(){
        par1 = new File(FluxSimulatorTest.class.getResource("/simulator_v10.par").getFile());
    }

    @Test
    public void testPar1(){
        try {
            FluxSimulatorSettings s1 = FluxSimulatorSettings.createSettings(par1);

            assertNotNull(s1.get(FluxSimulatorSettings.REF_FILE));
            assertEquals(FluxSimulatorTest.class.getResource("/spike_sequences.gtf").getFile(), s1.get(FluxSimulatorSettings.REF_FILE).getAbsolutePath());

            assertNotNull(s1.get(FluxSimulatorSettings.PRO_FILE));
            assertEquals("simulator_v10.pro", s1.get(FluxSimulatorSettings.PRO_FILE).getName());

            assertNotNull(s1.get(FluxSimulatorSettings.SEQ_FILE));
            assertEquals("simulator_v10.bed", s1.get(FluxSimulatorSettings.SEQ_FILE).getName());

            assertNotNull(s1.get(FluxSimulatorSettings.LIB_FILE));
            assertEquals("simulator_v10.lib", s1.get(FluxSimulatorSettings.LIB_FILE).getName());

            assertEquals(0d, s1.get(FluxSimulatorSettings.TSS_MEAN), 0.000001);

        } catch (Exception e) {
            e.printStackTrace();
            fail();
        }
    }
}
