package barna.flux.capacitor;

import barna.commons.Execute;
import barna.flux.capacitor.reconstruction.FluxCapacitorSettings;
import barna.flux.capacitor.reconstruction.FluxCapacitorSettings.AnnotationMapping;
import barna.flux.capacitor.reconstruction.FluxCapacitorStats;
import barna.flux.capacitor.utils.FluxCapacitorRunner;
import barna.io.FileHelper;
import barna.io.rna.UniversalReadDescriptor;
import com.google.gson.GsonBuilder;
import junit.framework.Assert;
import org.junit.*;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.EnumSet;
import java.util.HashMap;
import java.util.Map;

import static junit.framework.Assert.*;

public class FluxCapacitorTest {

    static final int SORTED = -1, UNSORT_GTF = 8, UNSORT_BED = 10;
    final File GTF_MM9_SORTED = new File(getClass().getResource("/mm9_chr1_chrX_sorted.gtf").getFile());
    final File GTF_MM9_UNSORTED = new File(getClass().getResource("/mm9_chr1_chrX.gtf").getFile());
    final File BED_MM9_SORTED = new File(getClass().getResource("/mm9_chr1_chrX_sorted.bed").getFile());
    final File BED_MM9_UNSORTED = new File(getClass().getResource("/mm9_chr1_chrX.bed").getFile());
    final File BED_MM9_SORTED_NO_CHR1 = new File(getClass().getResource("/mm9_chr1_chrX_sorted_no_chr1.bed").getFile());
    final File GTF_HG_SORTED = new File(getClass().getResource("/gencode_v12_hg_chr22_24030323-24041363.gtf").getFile());
    final File BED_HG_SORTED = new File(getClass().getResource("/test_hg_chr22_24030323-24041363.bed").getFile());

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
        currentTestDirectory = FileHelper.createTempDir("FluxCapacitorUnitTest", "", null);
    }

    @After
    public void cleanup(){
        if(currentTestDirectory != null){
            FileHelper.rmDir(currentTestDirectory);
        }
    }

    @Test
    public void testIOflatSortedWritableGTFflatSortedWritableBEDnoKeep() throws Exception {
        Map pars = new HashMap();
        pars.put("ANNOTATION_FILE", GTF_MM9_SORTED);
        pars.put("MAPPING_FILE", BED_MM9_SORTED);
        pars.put("ANNOTATION_MAPPING", AnnotationMapping.PAIRED);
        pars.put("READ_DESCRIPTOR", UniversalReadDescriptor.DESCRIPTORID_SIMULATOR);

        File parFile = FluxCapacitorRunner.createTestDir(currentTestDirectory, pars);

        FluxCapacitorRunner.runCapacitor(parFile);

        // check
        assertTrue(GTF_MM9_SORTED.exists());
        assertTrue(BED_MM9_SORTED.exists());
        assertTrue(new File(currentTestDirectory, FluxCapacitorRunner.DEFAULT_OUTPUT_FILE.toString()).exists());
        assertTrue(new File(currentTestDirectory, FluxCapacitorRunner.DEFAULT_PARAMETER_FILE.toString()).exists());
    }

    @Test
    public void testIOflatSortedWritableGTFflatUnsortedWritableBEDKeep() throws Exception {
        Map pars = new HashMap();
        pars.put("ANNOTATION_FILE", GTF_MM9_SORTED);
        pars.put("MAPPING_FILE", BED_MM9_UNSORTED);
        pars.put("ANNOTATION_MAPPING", AnnotationMapping.PAIRED);
        pars.put("READ_DESCRIPTOR", "SIMULATOR");
        pars.put("KEEP_SORTED", "tmp_sorted");

        File parFile = FluxCapacitorRunner.createTestDir(currentTestDirectory,pars);

        FluxCapacitorRunner.runCapacitor(parFile);

        // check
        assertTrue(GTF_MM9_SORTED.exists());
        assertTrue(BED_MM9_UNSORTED.exists());
        assertTrue(new File(currentTestDirectory, "tmp_sorted" + File.separator + BED_MM9_SORTED.getName()).exists());
        assertTrue(new File(currentTestDirectory, FluxCapacitorRunner.DEFAULT_OUTPUT_FILE.toString()).exists());
        assertTrue(new File(currentTestDirectory, FluxCapacitorRunner.DEFAULT_PARAMETER_FILE.toString()).exists());
    }

    @Test
    public void testIOflatUnsortedWritableGTFflatSortedWritableBEDKeep() throws Exception {
        Map pars = new HashMap();
        pars.put("ANNOTATION_FILE", GTF_MM9_UNSORTED);
        pars.put("MAPPING_FILE", BED_MM9_SORTED);
        pars.put("ANNOTATION_MAPPING", AnnotationMapping.PAIRED);
        pars.put("READ_DESCRIPTOR", "SIMULATOR");
        pars.put("KEEP_SORTED", "tmp_sorted");

        File parFile = FluxCapacitorRunner.createTestDir(currentTestDirectory,pars);

        FluxCapacitorRunner.runCapacitor(parFile);

        // check
        assertTrue(GTF_MM9_SORTED.exists());
        assertTrue(BED_MM9_UNSORTED.exists());
        assertTrue(new File(currentTestDirectory, "tmp_sorted" + File.separator + GTF_MM9_SORTED.getName()).exists());
        assertTrue(new File(currentTestDirectory, FluxCapacitorRunner.DEFAULT_OUTPUT_FILE.toString()).exists());
        assertTrue(new File(currentTestDirectory, FluxCapacitorRunner.DEFAULT_PARAMETER_FILE.toString()).exists());
    }

    @Test
    public void testIOgzippedSortedWritableGTFflatSortedWritableBEDnoKeep() throws Exception {
        Map pars = new HashMap();
        pars.put("ANNOTATION_FILE", GTF_MM9_SORTED);
        pars.put("MAPPING_FILE", BED_MM9_SORTED);
        pars.put("ANNOTATION_MAPPING", AnnotationMapping.PAIRED);
        pars.put("READ_DESCRIPTOR", "SIMULATOR");

        File parFile = FluxCapacitorRunner.createTestDir(currentTestDirectory,pars);

        FluxCapacitorStats stats = FluxCapacitorRunner.runCapacitor(parFile);

        // check
        assertTrue(GTF_MM9_SORTED.exists());
        assertTrue(BED_MM9_SORTED.exists());
        assertTrue(new File(currentTestDirectory, FluxCapacitorRunner.DEFAULT_OUTPUT_FILE.toString()).exists());
        assertTrue(new File(currentTestDirectory, FluxCapacitorRunner.DEFAULT_PARAMETER_FILE.toString()).exists());
    }

    @Test
    public void testNoDecompose() throws Exception {
        Map pars = new HashMap();
        pars.put("ANNOTATION_FILE", GTF_MM9_SORTED);
        pars.put("MAPPING_FILE", BED_MM9_SORTED);
        pars.put("ANNOTATION_MAPPING", AnnotationMapping.PAIRED);
        pars.put("READ_DESCRIPTOR", "SIMULATOR");
        pars.put("NO_DECOMPOSE", true);

        File parFile = FluxCapacitorRunner.createTestDir(currentTestDirectory,pars);

        FluxCapacitorStats stats = FluxCapacitorRunner.runCapacitor(parFile);
    }

    @Test
    public void testSJCount() throws Exception {

        Map pars = new HashMap();
        pars.put("ANNOTATION_FILE", GTF_MM9_SORTED);
        pars.put("MAPPING_FILE", BED_MM9_SORTED);
        pars.put("ANNOTATION_MAPPING", AnnotationMapping.PAIRED);
        pars.put("READ_DESCRIPTOR", "PAIRED");
        pars.put("NO_DECOMPOSE", true);
        pars.put("COUNT_ELEMENTS", EnumSet.of(FluxCapacitorSettings.CountElements.SPLICE_JUNCTIONS));

        File parFile = FluxCapacitorRunner.createTestDir(currentTestDirectory,pars);

        FluxCapacitorRunner.runCapacitor(parFile);
    }

    @Test
    public void testIntronsCount() throws Exception {

        Map pars = new HashMap();
        pars.put("ANNOTATION_FILE", GTF_MM9_SORTED);
        pars.put("MAPPING_FILE", BED_MM9_SORTED);
        pars.put("ANNOTATION_MAPPING", AnnotationMapping.PAIRED);
        pars.put("READ_DESCRIPTOR", "PAIRED");
        pars.put("NO_DECOMPOSE", true);
        pars.put("COUNT_ELEMENTS", EnumSet.of(FluxCapacitorSettings.CountElements.INTRONS));

        File parFile = FluxCapacitorRunner.createTestDir(currentTestDirectory,pars);

        FluxCapacitorRunner.runCapacitor(parFile);
    }

    @Test
    public void testAllCounters() throws Exception {

        Map pars = new HashMap();
        pars.put("ANNOTATION_FILE", GTF_MM9_SORTED);
        pars.put("MAPPING_FILE", BED_MM9_SORTED);
        pars.put("ANNOTATION_MAPPING", AnnotationMapping.PAIRED);
        pars.put("READ_DESCRIPTOR", "PAIRED");
        pars.put("NO_DECOMPOSE", true);
        pars.put("COUNT_ELEMENTS", EnumSet.allOf(FluxCapacitorSettings.CountElements.class));

        File parFile = FluxCapacitorRunner.createTestDir(currentTestDirectory,pars);

        FluxCapacitorRunner.runCapacitor(parFile);
    }

    @Test
    public void testDeconvolveAndCount() throws Exception {
        Map pars = new HashMap();
        pars.put("ANNOTATION_FILE", GTF_MM9_SORTED);
        pars.put("MAPPING_FILE", BED_MM9_SORTED);
        pars.put("ANNOTATION_MAPPING", AnnotationMapping.PAIRED);
        pars.put("READ_DESCRIPTOR", "PAIRED");
        pars.put("COUNT_ELEMENTS", EnumSet.allOf(FluxCapacitorSettings.CountElements.class));

        File parFile = FluxCapacitorRunner.createTestDir(currentTestDirectory,pars);

        FluxCapacitorRunner.runCapacitor(parFile);
    }

    @Test
    public void testReadsStranded() throws Exception {
        Map pars = new HashMap();
        pars.put("ANNOTATION_FILE", GTF_MM9_SORTED);
        pars.put("MAPPING_FILE", BED_MM9_SORTED);
        pars.put("ANNOTATION_MAPPING", AnnotationMapping.STRANDED);
        pars.put("READ_DESCRIPTOR", UniversalReadDescriptor.DESCRIPTORID_SENSE);

        File parFile = FluxCapacitorRunner.createTestDir(currentTestDirectory,pars);

        FluxCapacitorRunner.runCapacitor(parFile);

        // check
    }

    @Test
    public void testWithSortInRAM() throws Exception {
        Map pars = new HashMap();
        pars.put("ANNOTATION_FILE", GTF_MM9_SORTED);
        pars.put("MAPPING_FILE", BED_MM9_SORTED);
        pars.put("ANNOTATION_MAPPING", AnnotationMapping.PAIRED);
        pars.put("READ_DESCRIPTOR", "SIMULATOR");
        pars.put("SORT_IN_RAM", true);

        File parFile = FluxCapacitorRunner.createTestDir(currentTestDirectory,pars);

        FluxCapacitorRunner.runCapacitor(parFile);

        // check
    }

    @Test
    public void testStasAreWrittenAndContainValidData() throws Exception {
        File statsFile = new File(currentTestDirectory, "stats.txt");

        Map pars = new HashMap();
        pars.put("ANNOTATION_FILE", GTF_MM9_SORTED);
        pars.put("MAPPING_FILE", BED_MM9_SORTED);
        pars.put("ANNOTATION_MAPPING", AnnotationMapping.PAIRED);
        pars.put("READ_DESCRIPTOR", "SIMULATOR");
        pars.put("STATS_FILE", statsFile);

        File parFile = FluxCapacitorRunner.createTestDir(currentTestDirectory,pars);

        FluxCapacitorStats stats = FluxCapacitorRunner.runCapacitor(parFile);

        assertNotNull(stats);
        assertTrue(statsFile.exists());

        FluxCapacitorStats loaded = new GsonBuilder().create().fromJson(new FileReader(statsFile), FluxCapacitorStats.class);
        assertNotNull(loaded);

        assertEquals(loaded.getLociSingle(), stats.getLociSingle());
        assertEquals(loaded.getLociExp(), stats.getLociExp());
        assertEquals(loaded.getTxExp(), stats.getTxExp());
        assertEquals(loaded.getEventsExp(), stats.getEventsExp());
        assertEquals(loaded.getMappingsSingle(), stats.getMappingsSingle());
        assertEquals(loaded.getMappingsSinglePairs(), stats.getMappingsSinglePairs());
        assertEquals(loaded.getMappingsSinglePairsMapped(), stats.getMappingsSinglePairsMapped());
        assertEquals(loaded.getMappingsTotal(), stats.getMappingsTotal());
        assertEquals(loaded.getMappingsMapped(), stats.getMappingsMapped());
        assertEquals(loaded.getMappingsPairsNa(), stats.getMappingsPairsNa());
        assertEquals(loaded.getMappingsPairsWo(), stats.getMappingsPairsWo());
        assertEquals(loaded.getMappingsNotSens(), stats.getMappingsNotSens());
    }

    @Test
    public void testRPKMwithNrReadsMapped() throws Exception {
        File firstResults = new File(currentTestDirectory, "output/firstResults.gtf");

        Map pars = new HashMap();
        pars.put("ANNOTATION_FILE", GTF_MM9_SORTED);
        pars.put("MAPPING_FILE", BED_MM9_SORTED);
        pars.put("ANNOTATION_MAPPING", AnnotationMapping.PAIRED);
        pars.put("READ_DESCRIPTOR", "SIMULATOR");
        pars.put("STDOUT_FILE", firstResults);


        File parFile = FluxCapacitorRunner.createTestDir(currentTestDirectory,pars);

        FluxCapacitorRunner.runCapacitor(parFile);

        pars = new HashMap();
        pars.put("ANNOTATION_FILE", GTF_MM9_SORTED);
        pars.put("MAPPING_FILE", BED_MM9_SORTED);
        pars.put("ANNOTATION_MAPPING", AnnotationMapping.PAIRED);
        pars.put("READ_DESCRIPTOR", "SIMULATOR");
        pars.put("NR_READS_MAPPED", 7893);

        parFile = FluxCapacitorRunner.createTestDir(currentTestDirectory,pars);

        FluxCapacitorRunner.runCapacitor(parFile);

        BufferedReader b1 =null, b2 = null;
        try {
            b1 = new BufferedReader(new FileReader(firstResults));
            b2 = new BufferedReader(new FileReader(new File(currentTestDirectory,FluxCapacitorRunner.DEFAULT_OUTPUT_FILE.toString())));
            String s1, s2;
            while ((s1 = b1.readLine()) != null && (s2 = b2.readLine()) != null) {
                System.err.println(s1);
                String[] ss = s1.split("\\s");
                if (ss[1].equals("transcript")) {
                    if (ss[9].contains("NM_001159750"))
                        assertEquals(ss[ss.length - 1], "244929.484375");
                    else if (ss[9].contains("NM_001159751"))
                        assertEquals(ss[ss.length - 1], "32835.675781");
                    else if (ss[9].contains("NM_011541"))
                        assertEquals(ss[ss.length - 1], "77404.234375");
                    else if (ss[9].contains("NM_019397"))
                        assertEquals(ss[ss.length - 1], "27483.478516");
                    else
                        Assert.fail("Unknown Transcript ID: " + ss[9]);
                }

                assertEquals(s1, s2);
            }
            assertFalse(b1.ready());
            assertFalse(b2.ready());
        } catch (Exception e) {
            throw e;
        } finally {
            if (b1!=null)
                b1.close();
            if (b2!=null)
                b2.close();
        }
    }

    @Test
    public void testWrongReadDescriptor() throws Exception {
        Map pars = new HashMap();
        pars.put("ANNOTATION_FILE", GTF_MM9_SORTED);
        pars.put("MAPPING_FILE", BED_MM9_SORTED);
        pars.put("ANNOTATION_MAPPING", AnnotationMapping.PAIRED);
        pars.put("READ_DESCRIPTOR", UniversalReadDescriptor.DESCRIPTORID_MATE_STRAND_CSHL);

        File parFile = FluxCapacitorRunner.createTestDir(currentTestDirectory,pars);

        String msg = "";
        try {
            FluxCapacitorRunner.runCapacitor(parFile);
        } catch (Exception e) {
            msg = e.getMessage();
        }
        assertTrue(msg.contains("incompatible with read IDs"));
    }

    /**
     * A test to guarantee correct handling of loci without reads.
     */
    @Test
    public void testLocusWithoutReads() throws Exception {
        Map pars = new HashMap();
        pars.put("ANNOTATION_FILE", GTF_MM9_SORTED);
        pars.put("MAPPING_FILE", BED_MM9_SORTED_NO_CHR1);
        pars.put("ANNOTATION_MAPPING", AnnotationMapping.PAIRED);
        pars.put("READ_DESCRIPTOR", "SIMULATOR");

        File parFile = FluxCapacitorRunner.createTestDir(currentTestDirectory,pars);

        FluxCapacitorRunner.runCapacitor(parFile);
    }


    @Test
    public void testIOinsertSizes() throws Exception {
        File insFile = FileHelper.replaceSfx(new File(currentTestDirectory, FluxCapacitorRunner.DEFAULT_OUTPUT_FILE.toString()), "_ins.txt");

        Map pars = new HashMap();
        pars.put("ANNOTATION_FILE", GTF_MM9_SORTED);
        pars.put("MAPPING_FILE", BED_MM9_SORTED);
        pars.put("ANNOTATION_MAPPING", AnnotationMapping.PAIRED);
        pars.put("READ_DESCRIPTOR", "SIMULATOR");
        pars.put("INSERT_FILE", insFile);

        File parFile = FluxCapacitorRunner.createTestDir(currentTestDirectory,pars);

        FluxCapacitorStats stats = FluxCapacitorRunner.runCapacitor(parFile);

        assertNotNull(stats);
        assertEquals(1, stats.getLociSingle());
        assertEquals(0, stats.getLociExp());
        assertEquals(4, stats.getTxExp());
        assertEquals(0, stats.getEventsExp());
        assertEquals(566, stats.getMappingsSingle());
        assertEquals(586, stats.getMappingsSinglePairs());
        assertEquals(283, stats.getMappingsSinglePairsMapped());
        assertEquals(8005, stats.getMappingsTotal());
        assertEquals(8184, stats.getMappingsMapped());
        assertEquals(0, stats.getMappingsPairsNa());
        assertEquals(208, stats.getMappingsPairsWo());
        assertEquals(0, stats.getMappingsNotSens());

        // check
        BufferedReader buffy2 = null;
        try {
            buffy2 = new BufferedReader(new FileReader(insFile));
            String s = null;
            while ((s = buffy2.readLine()) != null) {
                String[] ss = s.split("\\s");
                assertEquals(9, ss.length);
            }
        } catch (Exception e) {
            throw e;
        } finally {
            if (buffy2!=null)
                buffy2.close();
        }
    }

    @Test
    public void testOutputProfiles() throws Exception {
        File proFile = new File(currentTestDirectory, FileHelper.append(FluxCapacitorRunner.DEFAULT_OUTPUT_FILE.toString(), "_profiles", true, "txt"));

        Map pars = new HashMap();
        pars.put("ANNOTATION_FILE", GTF_MM9_SORTED);
        pars.put("MAPPING_FILE", BED_MM9_SORTED);
        pars.put("ANNOTATION_MAPPING", AnnotationMapping.PAIRED);
        pars.put("READ_DESCRIPTOR", "SIMULATOR");
        pars.put("PROFILE_FILE", proFile);

        File parFile = FluxCapacitorRunner.createTestDir(currentTestDirectory,pars);

        FluxCapacitorRunner.runCapacitor(parFile);

        BufferedReader b1 = null;
        try {
            b1 = new BufferedReader(new FileReader(proFile));
            String s1;
            while ((s1 = b1.readLine()) != null) {
                System.err.println(s1);
            }
        } catch (Exception e) {
            throw e;
        } finally {
            if (b1!=null)
                b1.close();
        }
    }

    @Test
    public void testReadsPerTranscript() throws Exception {
        Map pars = new HashMap();
        pars.put("ANNOTATION_FILE", GTF_HG_SORTED);
        pars.put("MAPPING_FILE", BED_HG_SORTED);
        pars.put("ANNOTATION_MAPPING", AnnotationMapping.PAIRED);
        pars.put("READ_DESCRIPTOR", "CASAVA18");

        File parFile = FluxCapacitorRunner.createTestDir(currentTestDirectory,pars);

        FluxCapacitorRunner.runCapacitor(parFile);
    }
}
