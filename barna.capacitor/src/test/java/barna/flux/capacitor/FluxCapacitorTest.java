package barna.flux.capacitor;

import barna.commons.Execute;
import barna.flux.capacitor.profile.MappingStats;
import barna.flux.capacitor.reconstruction.FluxCapacitorSettings;
import barna.flux.capacitor.reconstruction.FluxCapacitorSettings.AnnotationMapping;
import barna.flux.capacitor.utils.FluxCapacitorRunner;
import barna.io.FileHelper;
import barna.io.rna.UniversalReadDescriptor;
import com.google.gson.GsonBuilder;
import junit.framework.Assert;
import org.junit.*;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.*;

import static junit.framework.Assert.*;

public class FluxCapacitorTest {

    static final int SORTED = -1, UNSORT_GTF = 8, UNSORT_BED = 10;
    final File GTF_MM9_SORTED = new File(getClass().getResource("/mm9_chr1_chrX_sorted.gtf").getFile());
    final File ENSEMBLE_SORTED = new File(getClass().getResource("/ensemble.gtf").getFile());
    final File BED_ENSEMBLE = new File(getClass().getResource("/15000first_lines.bed").getFile());
    final File GTF_MM9_UNSORTED = new File(getClass().getResource("/mm9_chr1_chrX.gtf").getFile());
    final File BED_MM9_SORTED = new File(getClass().getResource("/mm9_chr1_chrX_sorted.bed").getFile());
    final File BED_MM9_UNSORTED = new File(getClass().getResource("/mm9_chr1_chrX.bed").getFile());
    final File BED_MM9_SORTED_NO_CHR1 = new File(getClass().getResource("/mm9_chr1_chrX_sorted_no_chr1.bed").getFile());
    final File BED_MM9_PROFILE = new File(getClass().getResource("/mm9_chr1_chrX.profile").getFile());
    final File GTF_HG_SORTED = new File(getClass().getResource("/gencode_v12_hg_chr22_24030323-24041363.gtf").getFile());
    final File BED_HG_SORTED = new File(getClass().getResource("/test_hg_chr22_24030323-24041363.bed").getFile());
    final File BED_HG_PROFILE = new File(getClass().getResource("/test_hg_chr22_24030323-24041363.profile").getFile());

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
        pars.put("PROFILE_FILE", BED_MM9_PROFILE);
        pars.put("READ_DESCRIPTOR", UniversalReadDescriptor.DESCRIPTORID_SIMULATOR);

        File parFile = FluxCapacitorRunner.createTestDir(currentTestDirectory, pars);

        FluxCapacitorRunner.runCapacitor(parFile, null);

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
        pars.put("PROFILE_FILE", BED_MM9_PROFILE);
        pars.put("READ_DESCRIPTOR", "SIMULATOR");
        pars.put("KEEP_SORTED", "tmp_sorted");

        File parFile = FluxCapacitorRunner.createTestDir(currentTestDirectory, pars);

        FluxCapacitorRunner.runCapacitor(parFile, null);

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
        pars.put("PROFILE_FILE", BED_MM9_PROFILE);
        pars.put("READ_DESCRIPTOR", "SIMULATOR");
        pars.put("KEEP_SORTED", "tmp_sorted");

        File parFile = FluxCapacitorRunner.createTestDir(currentTestDirectory,pars);

        FluxCapacitorRunner.runCapacitor(parFile, null);

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
        pars.put("PROFILE_FILE", BED_MM9_PROFILE);
        pars.put("READ_DESCRIPTOR", "SIMULATOR");

        File parFile = FluxCapacitorRunner.createTestDir(currentTestDirectory,pars);

        MappingStats stats = FluxCapacitorRunner.runCapacitor(parFile, null);

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
        pars.put("PROFILE_FILE", BED_MM9_PROFILE);
        pars.put("READ_DESCRIPTOR", "SIMULATOR");
        pars.put("NO_DECOMPOSE", true);

        File parFile = FluxCapacitorRunner.createTestDir(currentTestDirectory,pars);

        MappingStats stats = FluxCapacitorRunner.runCapacitor(parFile, null);
    }

    @Test
    public void testSJCount() throws Exception {

        Map pars = new HashMap();
        pars.put("ANNOTATION_FILE", GTF_MM9_SORTED);
        pars.put("MAPPING_FILE", BED_MM9_SORTED);
        pars.put("PROFILE_FILE", BED_MM9_PROFILE);
        pars.put("READ_DESCRIPTOR", "PAIRED");
        pars.put("NO_DECOMPOSE", true);
        pars.put("COUNT_ELEMENTS", EnumSet.of(FluxCapacitorSettings.CountElements.SPLICE_JUNCTIONS));

        File parFile = FluxCapacitorRunner.createTestDir(currentTestDirectory,pars);

        FluxCapacitorRunner.runCapacitor(parFile, null);
    }

    @Test
    public void testIntronsCount() throws Exception {

        Map pars = new HashMap();
        pars.put("ANNOTATION_FILE", GTF_MM9_SORTED);
        pars.put("MAPPING_FILE", BED_MM9_SORTED);
        pars.put("PROFILE_FILE", BED_MM9_PROFILE);
        pars.put("READ_DESCRIPTOR", "PAIRED");
        pars.put("NO_DECOMPOSE", true);
        pars.put("COUNT_ELEMENTS", EnumSet.of(FluxCapacitorSettings.CountElements.INTRONS));

        File parFile = FluxCapacitorRunner.createTestDir(currentTestDirectory,pars);

        FluxCapacitorRunner.runCapacitor(parFile, null);
    }

    @Test
    public void testAllCounters() throws Exception {

        Map pars = new HashMap();
        pars.put("ANNOTATION_FILE", GTF_MM9_SORTED);
        pars.put("MAPPING_FILE", BED_MM9_SORTED);
        pars.put("PROFILE_FILE", BED_MM9_PROFILE);
        pars.put("READ_DESCRIPTOR", "PAIRED");
        pars.put("NO_DECOMPOSE", true);
        pars.put("COUNT_ELEMENTS", EnumSet.allOf(FluxCapacitorSettings.CountElements.class));

        File parFile = FluxCapacitorRunner.createTestDir(currentTestDirectory,pars);

        FluxCapacitorRunner.runCapacitor(parFile, null);
    }

    @Test
    public void testDeconvolveAndCount() throws Exception {
        Map pars = new HashMap();
        pars.put("ANNOTATION_FILE", GTF_MM9_SORTED);
        pars.put("MAPPING_FILE", BED_MM9_SORTED);
        pars.put("PROFILE_FILE", BED_MM9_PROFILE);
        pars.put("READ_DESCRIPTOR", "PAIRED");
        pars.put("COUNT_ELEMENTS", EnumSet.allOf(FluxCapacitorSettings.CountElements.class));

        File parFile = FluxCapacitorRunner.createTestDir(currentTestDirectory,pars);

        FluxCapacitorRunner.runCapacitor(parFile, null);
    }

    @Test
    public void testReadsStranded() throws Exception {
        Map pars = new HashMap();
        pars.put("ANNOTATION_FILE", GTF_MM9_SORTED);
        pars.put("MAPPING_FILE", BED_MM9_SORTED);
        pars.put("PROFILE_FILE", BED_MM9_PROFILE);
        pars.put("READ_DESCRIPTOR", UniversalReadDescriptor.DESCRIPTORID_SENSE);

        File parFile = FluxCapacitorRunner.createTestDir(currentTestDirectory,pars);

        FluxCapacitorRunner.runCapacitor(parFile, null);

        // check
    }

    @Test
    public void testWithSortInRAM() throws Exception {
        Map pars = new HashMap();
        pars.put("ANNOTATION_FILE", GTF_MM9_SORTED);
        pars.put("MAPPING_FILE", BED_MM9_SORTED);
        pars.put("PROFILE_FILE", BED_MM9_PROFILE);
        pars.put("READ_DESCRIPTOR", "SIMULATOR");
        pars.put("SORT_IN_RAM", true);

        File parFile = FluxCapacitorRunner.createTestDir(currentTestDirectory,pars);

        FluxCapacitorRunner.runCapacitor(parFile, null);

        // check
    }

    @Test
    public void testStasAreWrittenAndContainValidData() throws Exception {
        File statsFile = new File(currentTestDirectory, FileHelper.append(FluxCapacitorRunner.DEFAULT_OUTPUT_FILE.toString(), ".stats", true, ""));
        File proFile = new File(currentTestDirectory, FileHelper.append(FluxCapacitorRunner.DEFAULT_OUTPUT_FILE.toString(), ".profiles", true, ""));

        Map pars = new HashMap();
        pars.put("ANNOTATION_FILE", GTF_MM9_SORTED);
        pars.put("MAPPING_FILE", BED_MM9_SORTED);
        pars.put("ANNOTATION_MAPPING", AnnotationMapping.PAIRED);
        pars.put("READ_DESCRIPTOR", "SIMULATOR");
        pars.put("PROFILE_FILE", proFile);
        pars.put("STATS_FILE", statsFile);

        File parFile = FluxCapacitorRunner.createTestDir(currentTestDirectory,pars);
        String[] params = {"--profile", "-p", parFile.getAbsolutePath()};

        MappingStats stats = FluxCapacitorRunner.runCapacitor(parFile, params);
        stats = FluxCapacitorRunner.runCapacitor(parFile, null);

        assertNotNull(stats);
        assertTrue(statsFile.exists());

        MappingStats loaded = new GsonBuilder().create().fromJson(new FileReader(statsFile), MappingStats.class);
        assertNotNull(loaded);

        assertEquals(loaded.getSingleTxLoci(), stats.getSingleTxLoci());
        assertEquals(loaded.getLociExp(), stats.getLociExp());
        assertEquals(loaded.getTxsExp(), stats.getTxsExp());
        assertEquals(loaded.getEventsExp(), stats.getEventsExp());
        assertEquals(loaded.getReadsSingleTxLoci(), stats.getReadsSingleTxLoci());
        assertEquals(loaded.getMappingsSingleTxLoci(), stats.getMappingsSingleTxLoci());
        assertEquals(loaded.getMappingPairs(), stats.getMappingPairs());
        assertEquals(loaded.getMappingsTotal(), stats.getMappingsTotal());
        assertEquals(loaded.getMappingsMapped(), stats.getMappingsMapped());
        assertEquals(loaded.getMappingPairsNoTx(), stats.getMappingPairsNoTx());
        assertEquals(loaded.getPairsWrongOrientation(), stats.getPairsWrongOrientation());
        assertEquals(loaded.getMappingsWrongStrand(), stats.getMappingsWrongStrand());
    }

    @Test
    public void testRPKMwithNrReadsMapped() throws Exception {
        File firstResults = new File(currentTestDirectory, "output/firstResults.gtf");

        Map pars = new HashMap();
        pars.put("ANNOTATION_FILE", GTF_MM9_SORTED);
        pars.put("MAPPING_FILE", BED_MM9_SORTED);
        pars.put("PROFILE_FILE", BED_MM9_PROFILE);
        pars.put("ANNOTATION_MAPPING", AnnotationMapping.PAIRED);
        pars.put("READ_DESCRIPTOR", "SIMULATOR");
        pars.put("STDOUT_FILE", firstResults);


        File parFile = FluxCapacitorRunner.createTestDir(currentTestDirectory,pars);

        FluxCapacitorRunner.runCapacitor(parFile, null);

        pars = new HashMap();
        pars.put("ANNOTATION_FILE", GTF_MM9_SORTED);
        pars.put("MAPPING_FILE", BED_MM9_SORTED);
        pars.put("PROFILE_FILE", BED_MM9_PROFILE);
        pars.put("READ_DESCRIPTOR", "SIMULATOR");
        pars.put("NR_READS_MAPPED", 7893);

        parFile = FluxCapacitorRunner.createTestDir(currentTestDirectory,pars);

        FluxCapacitorRunner.runCapacitor(parFile, null);

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
        pars.put("PROFILE_FILE", BED_MM9_PROFILE);
        pars.put("ANNOTATION_MAPPING", AnnotationMapping.PAIRED);
        pars.put("READ_DESCRIPTOR", UniversalReadDescriptor.DESCRIPTORID_MATE_STRAND_CSHL);

        File parFile = FluxCapacitorRunner.createTestDir(currentTestDirectory,pars);

        String msg = "";
        try {
            FluxCapacitorRunner.runCapacitor(parFile, null);
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
        pars.put("PROFILE_FILE", BED_MM9_PROFILE);
        pars.put("ANNOTATION_MAPPING", AnnotationMapping.PAIRED);
        pars.put("READ_DESCRIPTOR", "SIMULATOR");

        File parFile = FluxCapacitorRunner.createTestDir(currentTestDirectory,pars);

        FluxCapacitorRunner.runCapacitor(parFile, null);
    }


    @Test
    public void testIOinsertSizes() throws Exception {
        File insFile = FileHelper.replaceSfx(new File(currentTestDirectory, FluxCapacitorRunner.DEFAULT_OUTPUT_FILE.toString()), "_ins.txt");

        Map pars = new HashMap();
        pars.put("ANNOTATION_FILE", GTF_MM9_SORTED);
        pars.put("MAPPING_FILE", BED_MM9_SORTED);
        pars.put("PROFILE_FILE", BED_MM9_PROFILE);
        pars.put("ANNOTATION_MAPPING", AnnotationMapping.PAIRED);
        pars.put("READ_DESCRIPTOR", "SIMULATOR");
        pars.put("INSERT_FILE", insFile);

        File parFile = FluxCapacitorRunner.createTestDir(currentTestDirectory,pars);

        MappingStats stats = FluxCapacitorRunner.runCapacitor(parFile, null);

        assertNotNull(stats);
        assertEquals(1, stats.getSingleTxLoci());
        assertEquals(2, stats.getLociExp());
        assertEquals(4, stats.getTxsExp());
        assertEquals(0, stats.getEventsExp());
        assertEquals(566, stats.getReadsSingleTxLoci());
        assertEquals(283, stats.getMappingsSingleTxLoci());
        assertEquals(586, stats.getMappingPairsSingleTxLoci());
        assertEquals(8005, stats.getMappingsTotal());
        assertEquals(8184, stats.getMappingsMapped());
        assertEquals(0, stats.getMappingPairsNoTx());
        assertEquals(192, stats.getPairsWrongOrientation());
        assertEquals(0, stats.getMappingsWrongStrand());

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
    public void testProfileOldWay() throws Exception {

        Map pars = new HashMap();
        pars.put("ANNOTATION_FILE", GTF_MM9_SORTED);
        pars.put("MAPPING_FILE", BED_MM9_SORTED);
        pars.put("ANNOTATION_MAPPING", AnnotationMapping.PAIRED);
        pars.put("READ_DESCRIPTOR", "SIMULATOR");

        File parFile = FluxCapacitorRunner.createTestDir(currentTestDirectory,pars);

        MappingStats stats = FluxCapacitorRunner.runCapacitor(parFile, null);

        assertNotNull(stats);
        assertEquals(1, stats.getSingleTxLoci());
        assertEquals(2, stats.getLociExp());
        assertEquals(4, stats.getTxsExp());
        assertEquals(0, stats.getEventsExp());
        assertEquals(566, stats.getReadsSingleTxLoci());
        assertEquals(283, stats.getMappingsSingleTxLoci());
        assertEquals(586, stats.getMappingPairsSingleTxLoci());
        assertEquals(8005, stats.getMappingsTotal());
        assertEquals(8184, stats.getMappingsMapped());
        assertEquals(0, stats.getMappingPairsNoTx());
        assertEquals(192, stats.getPairsWrongOrientation());
        assertEquals(0, stats.getMappingsWrongStrand());

    }

    @Test
    public void testOutputProfiles() throws Exception {
        File proFile = new File(currentTestDirectory, FileHelper.append(FluxCapacitorRunner.DEFAULT_OUTPUT_FILE.toString(), ".profiles", true, ""));

        Map pars = new HashMap();
        pars.put("ANNOTATION_FILE", GTF_MM9_SORTED);
        pars.put("MAPPING_FILE", BED_MM9_SORTED);
        pars.put("ANNOTATION_MAPPING", AnnotationMapping.PAIRED);
        pars.put("READ_DESCRIPTOR", "SIMULATOR");
        pars.put("PROFILE_FILE", proFile);

        File parFile = FluxCapacitorRunner.createTestDir(currentTestDirectory,pars);
        String[] params = {"--profile", "-p", parFile.getAbsolutePath()};

        FluxCapacitorRunner.runCapacitor(parFile, params);

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
    public void testCreatingCoverageFileWhereNoCoverageIsFound() throws Exception {
        // TEST for BARNA-269
        File proFile = new File(currentTestDirectory, FileHelper.append(FluxCapacitorRunner.DEFAULT_OUTPUT_FILE.toString(), "_profiles", true, "txt"));
        File coveragefile = File.createTempFile("coverage", ".test");
        coveragefile.deleteOnExit();
        Map pars = new HashMap();
        pars.put("ANNOTATION_FILE", ENSEMBLE_SORTED);
        pars.put("MAPPING_FILE", BED_ENSEMBLE);
        pars.put("ANNOTATION_MAPPING", AnnotationMapping.PAIRED);
        pars.put("READ_DESCRIPTOR", "PAIRED");
        pars.put("PROFILE_FILE", proFile);
        pars.put("COVERAGE_STATS", "true");
        pars.put("COVERAGE_FILE", coveragefile.getAbsolutePath());

        File parFile = FluxCapacitorRunner.createTestDir(currentTestDirectory,pars);
        FluxCapacitorRunner.runCapacitor(parFile);
        assertEquals(0, coveragefile.length());

    }
    @Test
    public void testCreatingCoverageFile() throws Exception {
        // TEST for BARNA-269
        File proFile = new File(currentTestDirectory, FileHelper.append(FluxCapacitorRunner.DEFAULT_OUTPUT_FILE.toString(), "_profiles", true, "txt"));
        File coveragefile = File.createTempFile("coverage", ".test");
        coveragefile.deleteOnExit();
        Map pars = new HashMap();
        pars.put("ANNOTATION_FILE", GTF_MM9_SORTED);
        pars.put("MAPPING_FILE", BED_MM9_SORTED);
        pars.put("ANNOTATION_MAPPING", AnnotationMapping.PAIRED);
        pars.put("READ_DESCRIPTOR", "PAIRED");
        pars.put("PROFILE_FILE", proFile);
        pars.put("COVERAGE_STATS", "true");
        pars.put("COVERAGE_FILE", coveragefile.getAbsolutePath());

        File parFile = FluxCapacitorRunner.createTestDir(currentTestDirectory,pars);
        FluxCapacitorRunner.runCapacitor(parFile);
        BufferedReader b1 = null;
        List<String> lines = new ArrayList<String>();
        try {
            b1 = new BufferedReader(new FileReader(coveragefile));
            String s1;
            while ((s1 = b1.readLine()) != null) {
                lines.add(s1);
            }
        } catch (Exception e) {
            throw e;
        } finally {
            if (b1!=null)
                b1.close();
        }
        assertEquals(1, lines.size());
        assertEquals("chrX:162960939-163023648C\tNM_019397\tCDS\t2700\t0\t0.0\t0\tNaN", lines.get(0));

    }

    @Test
    public void testReadsPerTranscript() throws Exception {
        Map pars = new HashMap();
        pars.put("ANNOTATION_FILE", GTF_HG_SORTED);
        pars.put("MAPPING_FILE", BED_HG_SORTED);
        pars.put("PROFILE_FILE", BED_HG_PROFILE);
        pars.put("ANNOTATION_MAPPING", AnnotationMapping.PAIRED);
        pars.put("READ_DESCRIPTOR", "PAIRED");

        File parFile = FluxCapacitorRunner.createTestDir(currentTestDirectory,pars);

        FluxCapacitorRunner.runCapacitor(parFile, null);
    }

}
