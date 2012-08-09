package barna.flux.capacitor;

import barna.commons.Execute;
import barna.flux.capacitor.integrationtest.FluxCapacitorRunner;
import barna.flux.capacitor.reconstruction.FluxCapacitor;
import barna.flux.capacitor.reconstruction.FluxCapacitorSettings.AnnotationMapping;
import barna.flux.capacitor.reconstruction.FluxCapacitorStats;
import barna.io.FileHelper;
import com.google.gson.GsonBuilder;
import org.junit.*;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.Future;

import static junit.framework.Assert.assertNotNull;
import static junit.framework.Assert.assertTrue;

public class FluxCapacitorReadsOutputTest {

	final File GTF_SORTED= new File(getClass().getResource("/hg_havana_chr7_small_sorted.gtf").getFile());
	final File BED_SORTED= new File(getClass().getResource("/hg_chr7_small_sorted.bed").getFile());

    protected FluxCapacitorStats runCapacitor(File parFile) throws Exception{
		FluxCapacitor capacitor= new FluxCapacitor();
		capacitor.setFile(parFile);
		Future<FluxCapacitorStats> captain= Execute.getExecutor().submit(capacitor);
        FluxCapacitorStats stats = captain.get();
        return stats;
	}

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
        File statsFile = new File(currentTestDirectory, "stats.txt");

        Map pars = new HashMap();
        pars.put("ANNOTATION_FILE", GTF_SORTED);
        pars.put("MAPPING_FILE", BED_SORTED);
        pars.put("ANNOTATION_MAPPING", AnnotationMapping.PAIRED);
        pars.put("READ_DESCRIPTOR", "CASAVA18");
        pars.put("STATS_FILE", statsFile);

        File parFile = FluxCapacitorRunner.createTestDir(currentTestDirectory, pars);

        FluxCapacitorStats stats = runCapacitor(parFile);

        assertNotNull(stats);
        assertTrue(statsFile.exists());

        FluxCapacitorStats loaded = new GsonBuilder().create().fromJson(new FileReader(statsFile), FluxCapacitorStats.class);
        File outFile = new File(currentTestDirectory, FluxCapacitorRunner.DEFAULT_OUTPUT_FILE.toString());
        assertTrue(outFile.exists());
        BufferedReader reader = new BufferedReader(new FileReader(outFile));
        String l = null;
        while((l = reader.readLine()) != null ){
            System.out.println(l);
            assertTrue(l.contains("reads"));
        }
	}
}
