package barna.flux.capacitor;

import barna.commons.Execute;
import barna.flux.capacitor.profile.MappingStats;
import barna.flux.capacitor.reconstruction.FluxCapacitorSettings.AnnotationMapping;
import barna.flux.capacitor.utils.FluxCapacitorRunner;
import barna.io.FileHelper;
import org.junit.*;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.HashMap;
import java.util.Map;

import static junit.framework.Assert.*;

public class FluxCapacitorReadsOutputTest {

	final File GTF_SORTED= new File(getClass().getResource("/hg_havana_chr7_small_sorted.gtf").getFile());
	final File BED_SORTED= new File(getClass().getResource("/hg_chr7_small_sorted.bed").getFile());

	@BeforeClass
	public static void initExecuter() {
		Execute.initialize(2);
	}
	
	@AfterClass
	public static void shutdownExecuter() {
		Execute.shutdown();
	}
    File currentTestDirectory = null;
    @Before
    public void setUpTest() throws Exception {
        currentTestDirectory = FileHelper.createTempDir("FluxCapacitorIntegration", "", null);
    }

    @After
    public void cleanup(){
        if(currentTestDirectory != null){
            FileHelper.rmDir(currentTestDirectory);
        }
    }
	@Test
	public void testStasAreWrittenAndContainValidData() throws Exception {
        File statsFile = new File(currentTestDirectory, FileHelper.append(FluxCapacitorRunner.DEFAULT_OUTPUT_FILE.toString(), ".stats", true, ""));
        File proFile = new File(currentTestDirectory, FileHelper.append(FluxCapacitorRunner.DEFAULT_OUTPUT_FILE.toString(), ".profiles", true, ""));

        Map pars = new HashMap();
        pars.put("ANNOTATION_FILE", GTF_SORTED);
        pars.put("MAPPING_FILE", BED_SORTED);
        pars.put("ANNOTATION_MAPPING", AnnotationMapping.PAIRED);
        pars.put("READ_DESCRIPTOR", "CASAVA18");
        pars.put("STATS_FILE", statsFile);
        pars.put("PROFILE_FILE", proFile);

        File parFile = FluxCapacitorRunner.createTestDir(currentTestDirectory, pars);
        String[] params = {"--profile", "-p", parFile.getAbsolutePath()};

        MappingStats stats = FluxCapacitorRunner.runCapacitor(parFile, params);

        assertNotNull(stats);
        assertTrue(statsFile.exists());

        File outFile = new File(currentTestDirectory, FluxCapacitorRunner.DEFAULT_OUTPUT_FILE.toString());
        assertFalse(outFile.exists());
        assertTrue(proFile.exists());
        BufferedReader reader = new BufferedReader(new FileReader(proFile));
        String l = null;
        while((l = reader.readLine()) != null ){
            System.err.println(l);
        }
	}
}
