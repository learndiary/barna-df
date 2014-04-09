package barna.flux.capacitor;

import barna.flux.capacitor.reconstruction.FluxCapacitorSettings;
import barna.flux.capacitor.reconstruction.PreprocessorTest;
import barna.flux.capacitor.utils.FluxCapacitorRunner;
import org.junit.Test;

import java.io.File;
import java.util.HashMap;
import java.util.Map;

/**
 * Tests based on RNA-Seq simulations.
 * User: micha
 */
public class SimulationTest extends GenericTest {

    // TODO apache.commons.math3 lib
    // http://commons.apache.org/proper/commons-math/javadocs/api-3.2/index.html

    @Test
    public void testSimulation01() throws Exception {

        if (1== 1)
            return; // PreprocessorTest.MAPPINGS_PSORT to be uploaded to artifactory

        Map pars = new HashMap();
        pars.put(FluxCapacitorSettings.ANNOTATION_FILE.getName(), PreprocessorTest.GENCODE_12);
        pars.put(FluxCapacitorSettings.MAPPING_FILE.getName(), PreprocessorTest.MAPPINGS_PSORT);
        pars.put(FluxCapacitorSettings.NO_FILE_CHECK.getName(), "TRUE");

        //pars.put(FluxCapacitorSettings.PROFILE_FILE.getName(), );

        File parFile = FluxCapacitorRunner.createTestDir(currentTestDirectory, pars);
        FluxCapacitorRunner.runCapacitor(parFile, null);
    }
}
