package barna.flux.capacitor;

import barna.commons.Execute;
import barna.commons.system.OSChecker;
import barna.flux.capacitor.profile.BiasProfiler;
import barna.flux.capacitor.profile.MappingStats;
import barna.flux.capacitor.profile.Profile;
import barna.flux.capacitor.reconstruction.FluxCapacitor;
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

    static {FluxCapacitor.DEBUG= false;}

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
    final File BED_HG_JUNCTION = new File(getClass().getResource("/chr1_329653_320881_junction.bed").getFile());
    final File GTF_HG_JUNCTION = new File(getClass().getResource("/chr1_329653_320881_junction.gtf").getFile());
    final File BAM_HG_MULTI = new File(getClass().getResource("/single_multimap.bam").getFile());
    final File GTF_HG_MULTI = new File(getClass().getResource("/single_multimap.gtf").getFile());

    @BeforeClass
    public static void initExecuter() {
        //Force en-US locale to use "." as the decimal separator in Windows OS
        if (OSChecker.isWindows()) {
            Locale.setDefault(new Locale("en", "US"));
        }
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
        pars.put(FluxCapacitorSettings.ANNOTATION_FILE.getName(), GTF_MM9_SORTED);
        pars.put(FluxCapacitorSettings.MAPPING_FILE.getName(), BED_MM9_SORTED);
        pars.put(FluxCapacitorSettings.PROFILE_FILE.getName(), BED_MM9_PROFILE);
        pars.put(FluxCapacitorSettings.READ_DESCRIPTOR.getName(), UniversalReadDescriptor.DESCRIPTORID_SIMULATOR);

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
        pars.put(FluxCapacitorSettings.ANNOTATION_FILE.getName(), GTF_MM9_SORTED);
        pars.put(FluxCapacitorSettings.MAPPING_FILE.getName(), BED_MM9_UNSORTED);
        pars.put(FluxCapacitorSettings.PROFILE_FILE.getName(), BED_MM9_PROFILE);
        pars.put(FluxCapacitorSettings.READ_DESCRIPTOR.getName(), UniversalReadDescriptor.DESCRIPTORID_SIMULATOR);
        pars.put(FluxCapacitorSettings.KEEP_SORTED.getName(), "tmp_sorted");

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
        pars.put(FluxCapacitorSettings.ANNOTATION_FILE.getName(), GTF_MM9_UNSORTED);
        pars.put(FluxCapacitorSettings.MAPPING_FILE.getName(), BED_MM9_SORTED);
        pars.put(FluxCapacitorSettings.PROFILE_FILE.getName(), BED_MM9_PROFILE);
        pars.put(FluxCapacitorSettings.READ_DESCRIPTOR.getName(), UniversalReadDescriptor.DESCRIPTORID_SIMULATOR);
        pars.put(FluxCapacitorSettings.KEEP_SORTED.getName(), "tmp_sorted");

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
        pars.put(FluxCapacitorSettings.ANNOTATION_FILE.getName(), GTF_MM9_SORTED);
        pars.put(FluxCapacitorSettings.MAPPING_FILE.getName(), BED_MM9_SORTED);
        pars.put(FluxCapacitorSettings.PROFILE_FILE.getName(), BED_MM9_PROFILE);
        pars.put(FluxCapacitorSettings.READ_DESCRIPTOR.getName(), UniversalReadDescriptor.DESCRIPTORID_SIMULATOR);

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
        pars.put(FluxCapacitorSettings.ANNOTATION_FILE.getName(), GTF_MM9_SORTED);
        pars.put(FluxCapacitorSettings.MAPPING_FILE.getName(), BED_MM9_SORTED);
        pars.put(FluxCapacitorSettings.PROFILE_FILE.getName(), BED_MM9_PROFILE);
        pars.put(FluxCapacitorSettings.READ_DESCRIPTOR.getName(), UniversalReadDescriptor.DESCRIPTORID_SIMULATOR);
        pars.put(FluxCapacitorSettings.DECONVOLUTE.getName(), false);

        File parFile = FluxCapacitorRunner.createTestDir(currentTestDirectory,pars);

        MappingStats stats = FluxCapacitorRunner.runCapacitor(parFile, null);
    }

    @Test
    public void testSJCount() throws Exception {

        Map pars = new HashMap();
        pars.put(FluxCapacitorSettings.ANNOTATION_FILE.getName(), GTF_MM9_SORTED);
        pars.put(FluxCapacitorSettings.MAPPING_FILE.getName(), BED_MM9_SORTED);
        pars.put(FluxCapacitorSettings.PROFILE_FILE.getName(), BED_MM9_PROFILE);
        pars.put(FluxCapacitorSettings.READ_DESCRIPTOR.getName(), UniversalReadDescriptor.DESCRIPTORID_PAIRED);
        pars.put(FluxCapacitorSettings.DECONVOLUTE.getName(), false);
        pars.put(FluxCapacitorSettings.COUNT_ELEMENTS.getName(), EnumSet.of(FluxCapacitorSettings.CountElements.SPLICE_JUNCTIONS));

        File parFile = FluxCapacitorRunner.createTestDir(currentTestDirectory,pars);

        FluxCapacitorRunner.runCapacitor(parFile, null);
    }

    @Test
    public void testIntronsCount() throws Exception {

        Map pars = new HashMap();
        pars.put(FluxCapacitorSettings.ANNOTATION_FILE.getName(), GTF_MM9_SORTED);
        pars.put(FluxCapacitorSettings.MAPPING_FILE.getName(), BED_MM9_SORTED);
        pars.put(FluxCapacitorSettings.PROFILE_FILE.getName(), BED_MM9_PROFILE);
        pars.put(FluxCapacitorSettings.READ_DESCRIPTOR.getName(), UniversalReadDescriptor.DESCRIPTORID_PAIRED);
        pars.put(FluxCapacitorSettings.DECONVOLUTE.getName(), false);
        pars.put(FluxCapacitorSettings.COUNT_ELEMENTS.getName(), EnumSet.of(FluxCapacitorSettings.CountElements.INTRONS));

        File parFile = FluxCapacitorRunner.createTestDir(currentTestDirectory,pars);

        FluxCapacitorRunner.runCapacitor(parFile, null);
    }

    @Test
    public void testAllCounters() throws Exception {

        Map pars = new HashMap();
        pars.put(FluxCapacitorSettings.ANNOTATION_FILE.getName(), GTF_MM9_SORTED);
        pars.put(FluxCapacitorSettings.MAPPING_FILE.getName(), BED_MM9_SORTED);
        pars.put(FluxCapacitorSettings.PROFILE_FILE.getName(), BED_MM9_PROFILE);
        pars.put(FluxCapacitorSettings.READ_DESCRIPTOR.getName(), UniversalReadDescriptor.DESCRIPTORID_PAIRED);
        pars.put(FluxCapacitorSettings.DECONVOLUTE.getName(), false);
        pars.put(FluxCapacitorSettings.COUNT_ELEMENTS.getName(), EnumSet.allOf(FluxCapacitorSettings.CountElements.class));

        File parFile = FluxCapacitorRunner.createTestDir(currentTestDirectory,pars);

        FluxCapacitorRunner.runCapacitor(parFile, null);
    }

    @Test
    public void testDeconvolveAndCount() throws Exception {
        Map pars = new HashMap();
        pars.put(FluxCapacitorSettings.ANNOTATION_FILE.getName(), GTF_MM9_SORTED);
        pars.put(FluxCapacitorSettings.MAPPING_FILE.getName(), BED_MM9_SORTED);
        pars.put(FluxCapacitorSettings.PROFILE_FILE.getName(), BED_MM9_PROFILE);
        pars.put(FluxCapacitorSettings.READ_DESCRIPTOR.getName(), UniversalReadDescriptor.DESCRIPTORID_PAIRED);
        pars.put(FluxCapacitorSettings.COUNT_ELEMENTS.getName(), EnumSet.allOf(FluxCapacitorSettings.CountElements.class));

        File parFile = FluxCapacitorRunner.createTestDir(currentTestDirectory,pars);

        FluxCapacitorRunner.runCapacitor(parFile, null);
    }

    @Test
    public void testReadsStranded() throws Exception {
        Map pars = new HashMap();
        pars.put(FluxCapacitorSettings.ANNOTATION_FILE.getName(), GTF_MM9_SORTED);
        pars.put(FluxCapacitorSettings.MAPPING_FILE.getName(), BED_MM9_SORTED);
        pars.put(FluxCapacitorSettings.PROFILE_FILE.getName(), BED_MM9_PROFILE);
        pars.put(FluxCapacitorSettings.READ_DESCRIPTOR.getName(), UniversalReadDescriptor.DESCRIPTORID_SENSE);

        File parFile = FluxCapacitorRunner.createTestDir(currentTestDirectory,pars);

        FluxCapacitorRunner.runCapacitor(parFile, null);

        // check
    }

    @Test
    public void testWithSortInRAM() throws Exception {
        Map pars = new HashMap();
        pars.put(FluxCapacitorSettings.ANNOTATION_FILE.getName(), GTF_MM9_SORTED);
        pars.put(FluxCapacitorSettings.MAPPING_FILE.getName(), BED_MM9_SORTED);
        pars.put(FluxCapacitorSettings.PROFILE_FILE.getName(), BED_MM9_PROFILE);
        pars.put(FluxCapacitorSettings.READ_DESCRIPTOR.getName(), UniversalReadDescriptor.DESCRIPTORID_SIMULATOR);
        pars.put(FluxCapacitorSettings.SORT_IN_RAM.getName(), true);

        File parFile = FluxCapacitorRunner.createTestDir(currentTestDirectory,pars);

        FluxCapacitorRunner.runCapacitor(parFile, null);

        // check
    }

    @Test
    public void testStasAreWrittenAndContainValidData() throws Exception {
        File statsFile = new File(currentTestDirectory, FileHelper.append(FluxCapacitorRunner.DEFAULT_OUTPUT_FILE.toString(), ".stats", true, null));
        File proFile = new File(currentTestDirectory, FileHelper.append(FluxCapacitorRunner.DEFAULT_OUTPUT_FILE.toString(), ".profile", true, null));

        Map pars = new HashMap();
        pars.put(FluxCapacitorSettings.ANNOTATION_FILE.getName(), GTF_MM9_SORTED);
        pars.put(FluxCapacitorSettings.MAPPING_FILE.getName(), BED_MM9_SORTED);
        pars.put(FluxCapacitorSettings.ANNOTATION_MAPPING.getName(), AnnotationMapping.PAIRED);
        pars.put(FluxCapacitorSettings.READ_DESCRIPTOR.getName(), UniversalReadDescriptor.DESCRIPTORID_SIMULATOR);
        pars.put(FluxCapacitorSettings.PROFILE_FILE.getName(), proFile);
        pars.put(FluxCapacitorSettings.STATS_FILE.getName(), statsFile);

        File parFile = FluxCapacitorRunner.createTestDir(currentTestDirectory,pars);
        String[] params = {"--profile", "-p", parFile.getAbsolutePath()};

        MappingStats stats = FluxCapacitorRunner.runCapacitor(parFile, params);
        File f= (File) pars.get(FluxCapacitorSettings.STDOUT_FILE.getName());
        f.delete(); // to avoid confirmation check
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
        assertEquals(loaded.getMappingPairsSingleTxLoci(), stats.getMappingPairsSingleTxLoci());
        assertEquals(loaded.getMappingsTotal(), stats.getMappingsTotal());
        assertEquals(loaded.getMappingsMapped(), stats.getMappingsMapped());
        assertEquals(loaded.getMappingPairsNoTx(), stats.getMappingPairsNoTx());
        assertEquals(loaded.getPairsWrongOrientation(), stats.getPairsWrongOrientation());
        assertEquals(loaded.getMappingsWrongStrand(), stats.getMappingsWrongStrand());
    }

    @Test
    public void testWrongReadDescriptor() throws Exception {
        Map pars = new HashMap();
        pars.put(FluxCapacitorSettings.ANNOTATION_FILE.getName(), GTF_MM9_SORTED);
        pars.put(FluxCapacitorSettings.MAPPING_FILE.getName(), BED_MM9_SORTED);
        pars.put(FluxCapacitorSettings.PROFILE_FILE.getName(), BED_MM9_PROFILE);
        pars.put(FluxCapacitorSettings.ANNOTATION_MAPPING.getName(), AnnotationMapping.PAIRED);
        pars.put(FluxCapacitorSettings.READ_DESCRIPTOR.getName(), UniversalReadDescriptor.DESCRIPTORID_MATE_STRAND_CSHL);

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
        pars.put(FluxCapacitorSettings.ANNOTATION_FILE.getName(), GTF_MM9_SORTED);
        pars.put(FluxCapacitorSettings.MAPPING_FILE.getName(), BED_MM9_SORTED_NO_CHR1);
        pars.put(FluxCapacitorSettings.PROFILE_FILE.getName(), BED_MM9_PROFILE);
        pars.put(FluxCapacitorSettings.ANNOTATION_MAPPING.getName(), AnnotationMapping.PAIRED);
        pars.put(FluxCapacitorSettings.READ_DESCRIPTOR.getName(), UniversalReadDescriptor.DESCRIPTORID_SIMULATOR);

        File parFile = FluxCapacitorRunner.createTestDir(currentTestDirectory,pars);

        FluxCapacitorRunner.runCapacitor(parFile, null);
    }


    @Test
    public void testIOinsertSizes() throws Exception {
        File insFile = FileHelper.replaceSfx(new File(currentTestDirectory, FluxCapacitorRunner.DEFAULT_OUTPUT_FILE.toString()), "_ins.txt");

        Map pars = new HashMap();
        pars.put(FluxCapacitorSettings.ANNOTATION_FILE.getName(), GTF_MM9_SORTED);
        pars.put(FluxCapacitorSettings.MAPPING_FILE.getName(), BED_MM9_SORTED);
        pars.put(FluxCapacitorSettings.PROFILE_FILE.getName(), BED_MM9_PROFILE);
        pars.put(FluxCapacitorSettings.ANNOTATION_MAPPING.getName(), AnnotationMapping.PAIRED);
        pars.put(FluxCapacitorSettings.READ_DESCRIPTOR.getName(), UniversalReadDescriptor.DESCRIPTORID_SIMULATOR);
        pars.put(FluxCapacitorSettings.INSERT_FILE.getName(), insFile);

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
        assertEquals(208, stats.getPairsWrongOrientation());
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
    public void testProfile() throws Exception {

        Map pars = new HashMap();
        pars.put(FluxCapacitorSettings.ANNOTATION_FILE.getName(), GTF_MM9_SORTED);
        pars.put(FluxCapacitorSettings.MAPPING_FILE.getName(), BED_MM9_SORTED);
        pars.put(FluxCapacitorSettings.PROFILE_FILE.getName(), BED_MM9_PROFILE);
        pars.put(FluxCapacitorSettings.ANNOTATION_MAPPING.getName(), AnnotationMapping.PAIRED);
        pars.put(FluxCapacitorSettings.READ_DESCRIPTOR.getName(), UniversalReadDescriptor.DESCRIPTORID_SIMULATOR);

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
        assertEquals(208, stats.getPairsWrongOrientation());
        assertEquals(0, stats.getMappingsWrongStrand());

    }

    @Test
    public void testProfileOldWay() throws Exception {

        Map pars = new HashMap();
        pars.put(FluxCapacitorSettings.ANNOTATION_FILE.getName(), GTF_MM9_SORTED);
        pars.put(FluxCapacitorSettings.MAPPING_FILE.getName(), BED_MM9_SORTED);
        pars.put(FluxCapacitorSettings.ANNOTATION_MAPPING.getName(), AnnotationMapping.PAIRED);
        pars.put(FluxCapacitorSettings.READ_DESCRIPTOR.getName(), UniversalReadDescriptor.DESCRIPTORID_SIMULATOR);

        File parFile = FluxCapacitorRunner.createTestDir(currentTestDirectory,pars);

        MappingStats stats = FluxCapacitorRunner.runCapacitor(parFile, null);

        assertNotNull(stats);
        assertEquals(1, stats.getSingleTxLoci());
        assertEquals(2, stats.getLociExp());
        assertEquals(4, stats.getTxsExp());
        assertEquals(0, stats.getEventsExp());
        //assertEquals(283, stats.getReadsSingleTxLoci());
        assertEquals(283, stats.getMappingsSingleTxLoci());
        assertEquals(586, stats.getMappingPairsSingleTxLoci());
        assertEquals(8005, stats.getMappingsTotal());
        assertEquals(8184, stats.getMappingsMapped());
        assertEquals(0, stats.getMappingPairsNoTx());
        assertEquals(208, stats.getPairsWrongOrientation());
        assertEquals(0, stats.getMappingsWrongStrand());

    }

    @Test
    public void testOutputProfiles() throws Exception {
        File proFile = new File(currentTestDirectory, FileHelper.append(FluxCapacitorRunner.DEFAULT_OUTPUT_FILE.toString(), ".profiles", true, null));

        Map pars = new HashMap();
        pars.put(FluxCapacitorSettings.ANNOTATION_FILE.getName(), GTF_MM9_SORTED);
        pars.put(FluxCapacitorSettings.MAPPING_FILE.getName(), BED_MM9_SORTED);
        pars.put(FluxCapacitorSettings.ANNOTATION_MAPPING.getName(), AnnotationMapping.PAIRED);
        pars.put(FluxCapacitorSettings.READ_DESCRIPTOR.getName(), UniversalReadDescriptor.DESCRIPTORID_SIMULATOR);
        pars.put(FluxCapacitorSettings.PROFILE_FILE.getName(), proFile);

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
        pars.put(FluxCapacitorSettings.ANNOTATION_FILE.getName(), ENSEMBLE_SORTED);
        pars.put(FluxCapacitorSettings.MAPPING_FILE.getName(), BED_ENSEMBLE);
        pars.put(FluxCapacitorSettings.ANNOTATION_MAPPING.getName(), AnnotationMapping.PAIRED);
        pars.put(FluxCapacitorSettings.READ_DESCRIPTOR.getName(), UniversalReadDescriptor.DESCRIPTORID_PAIRED);
        pars.put(FluxCapacitorSettings.PROFILE_FILE.getName(), proFile);
        pars.put(FluxCapacitorSettings.COVERAGE_FILE.getName(), coveragefile.getAbsolutePath());

        File parFile = FluxCapacitorRunner.createTestDir(currentTestDirectory,pars);
        String[] params = {"--profile", "-p", parFile.getAbsolutePath()};
        FluxCapacitorRunner.runCapacitor(parFile, params);
        assertEquals(0, coveragefile.length());

    }
    @Test
    public void testCreatingCoverageFile() throws Exception {
        // TEST for BARNA-269
        File proFile = new File(currentTestDirectory, FileHelper.append(FluxCapacitorRunner.DEFAULT_OUTPUT_FILE.toString(), "_profiles", true, "txt"));
        File coveragefile = File.createTempFile("coverage", ".test");
        coveragefile.deleteOnExit();
        Map pars = new HashMap();
        pars.put(FluxCapacitorSettings.ANNOTATION_FILE.getName(), GTF_MM9_SORTED);
        pars.put(FluxCapacitorSettings.MAPPING_FILE.getName(), BED_MM9_SORTED);
        pars.put(FluxCapacitorSettings.ANNOTATION_MAPPING.getName(), AnnotationMapping.PAIRED);
        pars.put(FluxCapacitorSettings.READ_DESCRIPTOR.getName(), UniversalReadDescriptor.DESCRIPTORID_PAIRED);
        pars.put(FluxCapacitorSettings.PROFILE_FILE.getName(), proFile);
        pars.put(FluxCapacitorSettings.COVERAGE_FILE.getName(), coveragefile.getAbsolutePath());

        File parFile = FluxCapacitorRunner.createTestDir(currentTestDirectory,pars);
        String[] params = {"--profile", "-p", parFile.getAbsolutePath()};
        FluxCapacitorRunner.runCapacitor(parFile, params);
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
        pars.put(FluxCapacitorSettings.ANNOTATION_FILE.getName(), GTF_HG_SORTED);
        pars.put(FluxCapacitorSettings.MAPPING_FILE.getName(), BED_HG_SORTED);
        pars.put(FluxCapacitorSettings.PROFILE_FILE.getName(), BED_HG_PROFILE);
        pars.put(FluxCapacitorSettings.ANNOTATION_MAPPING.getName(), AnnotationMapping.PAIRED);
        pars.put(FluxCapacitorSettings.READ_DESCRIPTOR.getName(), UniversalReadDescriptor.DESCRIPTORID_PAIRED);

        File parFile = FluxCapacitorRunner.createTestDir(currentTestDirectory, pars);

        FluxCapacitorRunner.runCapacitor(parFile, null);
    }

    @Test
    public void testProfileFile() throws Exception {
        File proFile = new File(currentTestDirectory, FileHelper.append(FluxCapacitorRunner.DEFAULT_OUTPUT_FILE.toString(), ".profile", true, null));

        Map pars = new HashMap();
        pars.put(FluxCapacitorSettings.ANNOTATION_FILE.getName(), GTF_MM9_SORTED);
        pars.put(FluxCapacitorSettings.MAPPING_FILE.getName(), BED_MM9_SORTED);
        pars.put(FluxCapacitorSettings.ANNOTATION_MAPPING.getName(), AnnotationMapping.PAIRED);
        pars.put(FluxCapacitorSettings.READ_DESCRIPTOR.getName(), UniversalReadDescriptor.DESCRIPTORID_SIMULATOR);
        pars.put(FluxCapacitorSettings.PROFILE_FILE.getName(), proFile);

        File parFile = FluxCapacitorRunner.createTestDir(currentTestDirectory,pars);
        String[] params = {"--profile", "-p", parFile.getAbsolutePath()};

        FluxCapacitorRunner.runCapacitor(parFile, params);

        Profile runProfile = BiasProfiler.readProfile(proFile, true);
        Profile refProfile = BiasProfiler.readProfile(BED_MM9_PROFILE, true);

        assertEquals(runProfile, refProfile);
    }

    @Test
    public void testFluxGtf() throws Exception {
        File proFile = new File(currentTestDirectory, FileHelper.append(FluxCapacitorRunner.DEFAULT_OUTPUT_FILE.toString(), ".profile", true, null));

        Map pars = new HashMap();
        pars.put(FluxCapacitorSettings.ANNOTATION_FILE.getName(), GTF_MM9_SORTED);
        pars.put(FluxCapacitorSettings.MAPPING_FILE.getName(), BED_MM9_SORTED);
        pars.put(FluxCapacitorSettings.PROFILE_FILE.getName(), proFile);
        pars.put(FluxCapacitorSettings.READ_DESCRIPTOR.getName(), UniversalReadDescriptor.DESCRIPTORID_SIMULATOR);
        pars.put(FluxCapacitorSettings.ANNOTATION_MAPPING.getName(), AnnotationMapping.PAIRED);
        pars.put(FluxCapacitorSettings.COUNT_ELEMENTS.getName(),"[]");

        File parFile = FluxCapacitorRunner.createTestDir(currentTestDirectory, pars);
        String[] params = {"--profile", "-p", parFile.getAbsolutePath()};
        FluxCapacitorRunner.runCapacitor(parFile, params);
        File f= (File) pars.get(FluxCapacitorSettings.STDOUT_FILE.getName());
        f.delete(); // to avoid confirmation check
        FluxCapacitorRunner.runCapacitor(parFile,null);

        // check
        assertTrue(GTF_MM9_SORTED.exists());
        assertTrue(BED_MM9_SORTED.exists());
        File output = new File(currentTestDirectory, FluxCapacitorRunner.DEFAULT_OUTPUT_FILE.toString());
        assertTrue(output.exists());
        assertTrue(new File(currentTestDirectory, FluxCapacitorRunner.DEFAULT_PARAMETER_FILE.toString()).exists());

        BufferedReader refGtf = new BufferedReader(new FileReader(getClass().getResource("/mm9_chr1_chrX_flux.gtf").getFile()));
        BufferedReader runGtf = new BufferedReader(new FileReader(output));

        List<String> runLines = new ArrayList<String>();
        List<String> refLines = new ArrayList<String>();

        String line;
        while ((line = runGtf.readLine()) != null) {
            runLines.add(line);
        }
        runGtf.close();
        while ((line = refGtf.readLine()) != null) {
            refLines.add(line);
        }
        refGtf.close();

        assertEquals(refLines.size(),runLines.size());

        for (int i = 0; i < refLines.size(); i++) {
            assertEquals(refLines.get(i),runLines.get(i));
        }
    }

    @Test
    public void testJunctionsStranded() throws Exception {
        Map pars = new HashMap();
        pars.put(FluxCapacitorSettings.ANNOTATION_FILE.getName(), GTF_HG_JUNCTION);
        pars.put(FluxCapacitorSettings.MAPPING_FILE.getName(), BED_HG_JUNCTION);
        pars.put(FluxCapacitorSettings.READ_DESCRIPTOR.getName(), UniversalReadDescriptor.DESCRIPTORID_MATE2_SENSE);
        pars.put(FluxCapacitorSettings.ANNOTATION_MAPPING.getName(), AnnotationMapping.PAIRED_STRANDED);
        pars.put(FluxCapacitorSettings.DECONVOLUTE.getName(), false);
        pars.put(FluxCapacitorSettings.COUNT_ELEMENTS.getName(),FluxCapacitorSettings.CountElements.SPLICE_JUNCTIONS);

        File parFile = FluxCapacitorRunner.createTestDir(currentTestDirectory, pars);
        FluxCapacitorRunner.runCapacitor(parFile,null);

        assertTrue(GTF_HG_JUNCTION.exists());
        assertTrue(BED_HG_JUNCTION.exists());
        File output = new File(currentTestDirectory, FluxCapacitorRunner.DEFAULT_OUTPUT_FILE.toString());
        assertTrue(output.exists());
        assertTrue(new File(currentTestDirectory, FluxCapacitorRunner.DEFAULT_PARAMETER_FILE.toString()).exists());

        BufferedReader runGtf = new BufferedReader(new FileReader(output));

        List<String> runLines = new ArrayList<String>();

        String line;
        while ((line = runGtf.readLine()) != null) {
            runLines.add(line);
        }
        runGtf.close();

        assertEquals(runLines.size(),2);
        assertEquals("1\tflux\tjunction\t320653\t320881\t.\t+\t.\tgene_id \"ENSG00000237094.3\"; locus_id \"1:320162-328453W\"; reads 4.000000", runLines.get(0));
        assertEquals("1\tflux\tjunction\t320938\t321032\t.\t+\t.\tgene_id \"ENSG00000237094.3\"; locus_id \"1:320162-328453W\"; reads 2.000000",runLines.get(1));
    }

    @Test
    public void testMultiMapsPrimaryOnly() throws Exception {

        Map pars = new HashMap();
        pars.put(FluxCapacitorSettings.ANNOTATION_FILE.getName(), GTF_HG_MULTI);
        pars.put(FluxCapacitorSettings.MAPPING_FILE.getName(), BAM_HG_MULTI);
        pars.put(FluxCapacitorSettings.ANNOTATION_MAPPING.getName(), AnnotationMapping.PAIRED);
        pars.put(FluxCapacitorSettings.SAM_PRIMARY_ONLY.getName(), true);

        File parFile = FluxCapacitorRunner.createTestDir(currentTestDirectory,pars);

        MappingStats stats = FluxCapacitorRunner.runCapacitor(parFile, null);

        assertNotNull(stats);
        assertEquals(1, stats.getSingleTxLoci());
        //assertEquals(1, stats.getReadsSingleTxLoci());
        assertEquals(1, stats.getMappingsSingleTxLoci());
        assertEquals(2, stats.getMappingPairsSingleTxLoci());
        assertEquals(6, stats.getMappingsTotal());
        assertEquals(2, stats.getMappingsMapped());
        assertEquals(0, stats.getMappingPairsNoTx());
        assertEquals(0, stats.getPairsWrongOrientation());
        assertEquals(0, stats.getMappingsWrongStrand());

    }

    @Test
    public void testMultiMapsMatesOnly() throws Exception {

        Map pars = new HashMap();
        pars.put(FluxCapacitorSettings.ANNOTATION_FILE.getName(), GTF_HG_MULTI);
        pars.put(FluxCapacitorSettings.MAPPING_FILE.getName(), BAM_HG_MULTI);
        pars.put(FluxCapacitorSettings.ANNOTATION_MAPPING.getName(), AnnotationMapping.PAIRED);
        pars.put(FluxCapacitorSettings.SAM_MATES_ONLY.getName(), true);
        pars.put(FluxCapacitorSettings.SAM_PRIMARY_ONLY.getName(), false);

        File parFile = FluxCapacitorRunner.createTestDir(currentTestDirectory,pars);

        MappingStats stats = FluxCapacitorRunner.runCapacitor(parFile, null);

        assertNotNull(stats);
        assertEquals(1, stats.getSingleTxLoci());
        assertEquals(3, stats.getReadsSingleTxLoci());
        assertEquals(1, stats.getMappingsSingleTxLoci());
        assertEquals(2, stats.getMappingPairsSingleTxLoci());
        assertEquals(6, stats.getMappingsTotal());
        assertEquals(6, stats.getMappingsMapped());
        assertEquals(0, stats.getMappingPairsNoTx());
        assertEquals(0, stats.getPairsWrongOrientation());
        assertEquals(0, stats.getMappingsWrongStrand());

    }

    @Test
    public void testMultiMapsWeighted() throws Exception {

        Map pars = new HashMap();
        pars.put(FluxCapacitorSettings.ANNOTATION_FILE.getName(), GTF_HG_MULTI);
        pars.put(FluxCapacitorSettings.MAPPING_FILE.getName(), BAM_HG_MULTI);
        pars.put(FluxCapacitorSettings.ANNOTATION_MAPPING.getName(), AnnotationMapping.PAIRED);
        pars.put(FluxCapacitorSettings.WEIGHTED_COUNT.getName(), true);
        pars.put(FluxCapacitorSettings.SAM_PRIMARY_ONLY.getName(), false);

        File parFile = FluxCapacitorRunner.createTestDir(currentTestDirectory,pars);

        MappingStats stats = FluxCapacitorRunner.runCapacitor(parFile, null);

        assertNotNull(stats);
        assertEquals(1, stats.getSingleTxLoci());
        assertEquals(3, stats.getReadsSingleTxLoci());
        assertEquals(0, stats.getMappingsSingleTxLoci());
        assertEquals(1, stats.getMappingPairsSingleTxLoci());
        assertEquals(6, stats.getMappingsTotal());
        assertEquals(2, stats.getMappingsMapped());
        assertEquals(0, stats.getMappingPairsNoTx());
        assertEquals(0, stats.getPairsWrongOrientation());
        assertEquals(0, stats.getMappingsWrongStrand());

    }

    @Test
    public void testMultiMapsMatesOnlyWeighted() throws Exception {

        Map pars = new HashMap();
        pars.put(FluxCapacitorSettings.ANNOTATION_FILE.getName(), GTF_HG_MULTI);
        pars.put(FluxCapacitorSettings.MAPPING_FILE.getName(), BAM_HG_MULTI);
        pars.put(FluxCapacitorSettings.ANNOTATION_MAPPING.getName(), AnnotationMapping.PAIRED);
        pars.put(FluxCapacitorSettings.SAM_MATES_ONLY.getName(), true);
        pars.put(FluxCapacitorSettings.WEIGHTED_COUNT.getName(),true);
        pars.put(FluxCapacitorSettings.SAM_PRIMARY_ONLY.getName(), false);

        File parFile = FluxCapacitorRunner.createTestDir(currentTestDirectory,pars);

        MappingStats stats = FluxCapacitorRunner.runCapacitor(parFile, null);

        assertNotNull(stats);
        assertEquals(1, stats.getSingleTxLoci());
        assertEquals(3, stats.getReadsSingleTxLoci());
        assertEquals(0, stats.getMappingsSingleTxLoci());   // 0.3 rounded down
        assertEquals(1, stats.getMappingPairsSingleTxLoci());
        assertEquals(6, stats.getMappingsTotal());
        assertEquals(2, stats.getMappingsMapped());
        assertEquals(0, stats.getMappingPairsNoTx());
        assertEquals(0, stats.getPairsWrongOrientation());
        assertEquals(0, stats.getMappingsWrongStrand());

    }

    @Test
    public void testMultiMapsUniqueOnly() throws Exception {

        Map pars = new HashMap();
        pars.put(FluxCapacitorSettings.ANNOTATION_FILE.getName(), GTF_HG_MULTI);
        pars.put(FluxCapacitorSettings.MAPPING_FILE.getName(), BAM_HG_MULTI);
        pars.put(FluxCapacitorSettings.ANNOTATION_MAPPING.getName(), AnnotationMapping.PAIRED);
        pars.put(FluxCapacitorSettings.SAM_UNIQUE_ONLY.getName(), true);
        pars.put(FluxCapacitorSettings.SAM_MATES_ONLY.getName(), false);    // 2 unique mappings, but not in a sam pair
        pars.put(FluxCapacitorSettings.SAM_PRIMARY_ONLY.getName(), false);

        File parFile = FluxCapacitorRunner.createTestDir(currentTestDirectory,pars);

        MappingStats stats = FluxCapacitorRunner.runCapacitor(parFile, null);

        assertNotNull(stats);
        assertEquals(1, stats.getSingleTxLoci());
        assertEquals(1, stats.getReadsSingleTxLoci());
        assertEquals(1, stats.getMappingsSingleTxLoci());
        assertEquals(0, stats.getMappingPairsSingleTxLoci());
        assertEquals(6, stats.getMappingsTotal());
        assertEquals(2, stats.getMappingsMapped());
        assertEquals(0, stats.getMappingPairsNoTx());
        assertEquals(0, stats.getPairsWrongOrientation());
        assertEquals(0, stats.getMappingsWrongStrand());

    }

    @Test
    public void testMultiMapsMatesOnlyUniqueOnly() throws Exception {

        Map pars = new HashMap();
        pars.put(FluxCapacitorSettings.ANNOTATION_FILE.getName(), GTF_HG_MULTI);
        pars.put(FluxCapacitorSettings.MAPPING_FILE.getName(), BAM_HG_MULTI);
        pars.put(FluxCapacitorSettings.ANNOTATION_MAPPING.getName(), AnnotationMapping.PAIRED);
        pars.put(FluxCapacitorSettings.SAM_MATES_ONLY.getName(), true);
        pars.put(FluxCapacitorSettings.SAM_UNIQUE_ONLY.getName(), true);
        //pars.put("WEIGHTED_COUNT", true);

        File parFile = FluxCapacitorRunner.createTestDir(currentTestDirectory,pars);

        MappingStats stats = FluxCapacitorRunner.runCapacitor(parFile, null);

        assertNotNull(stats);
        assertEquals(1, stats.getSingleTxLoci());
        assertEquals(1, stats.getReadsSingleTxLoci());
        assertEquals(1, stats.getMappingsSingleTxLoci());
        assertEquals(0, stats.getMappingPairsSingleTxLoci());  // profiling only allows sam-pairing
        assertEquals(6, stats.getMappingsTotal());
        assertEquals(0, stats.getMappingsMapped());
        assertEquals(0, stats.getMappingPairsNoTx());
        assertEquals(0, stats.getPairsWrongOrientation());
        assertEquals(0, stats.getMappingsWrongStrand());

    }
}
