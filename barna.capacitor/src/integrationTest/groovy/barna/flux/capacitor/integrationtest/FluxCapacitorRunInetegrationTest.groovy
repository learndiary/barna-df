	
package barna.flux.capacitor.integrationtest

import barna.commons.Execute
import barna.commons.system.OSChecker
import barna.flux.capacitor.reconstruction.FluxCapacitorSettings
import barna.flux.capacitor.reconstruction.FluxCapacitorSettings.AnnotationMapping
import barna.io.FileHelper
import barna.io.Sorter
import barna.io.rna.UniversalReadDescriptor
import org.junit.AfterClass
import org.junit.BeforeClass
import org.junit.Test

import java.util.concurrent.Future
import java.util.zip.GZIPOutputStream
import java.util.zip.ZipEntry
import java.util.zip.ZipOutputStream

import static junit.framework.Assert.assertTrue
import static org.junit.Assert.fail
import org.junit.Before
import org.junit.After

/**
 * 
 * @author Thasso Griebel (thasso.griebel@gmail.com)
 */

class FluxCapacitorRunInetegrationTest {

	static final int SORTED= -1, UNSORT_GTF= 2, UNSORT_BED= 9;
	final File GTF_SORTED= new File(getClass().getResource("/mm9_chr1_chrX.gtf").getFile());
	final File BED_SORTED= new File(getClass().getResource("/chr1_chrX.bed").getFile());

    static String executable
    static File tmpDir

    @BeforeClass
    public static void setUp(){
        executable = System.getProperty("dist.exe")
        if(executable == null){
            fail("No capacitor executable specified")
        }
        Execute.initialize(2);

    }

    @AfterClass
    public static void shutdownExecuter() {
        Execute.shutdown();
    }

	protected String runCapacitor() throws Exception{
        File parFile = FluxCapacitorRunner.getParFile();
        tmpDir = FluxCapacitorRunner.getTmpDir();
        def pb = new ProcessBuilder()
        pb.environment().put("FLUX_MEM", "1G")
        if (tmpDir != null){
            pb.environment().put("JAVA_OPTS", "-Dflux.io.deny.tmpdir=yes")
        }
        def cmd = [executable, "-p", parFile.getAbsolutePath()]
        if (OSChecker.isWindows()) {
            cmd = ["cmd", "/c", executable, "-p", parFile.getAbsolutePath()]
        }
        def process = pb.directory(tmpDir != null ? tmpDir : parFile.getParentFile())
                .redirectErrorStream(true)
                .command(cmd)
                .start()
        String output = process.inputStream.text
        process.waitFor()
		return output;
		
		
	}

	static final String[] STDERR_MAPPED= ["8009","8192"]
	static final String[] STDERR_ACCESS_DENIED= ["access denied"]
	
	void assertFiles(int nrFilesInGTF, int nrFilesInBED, String stderr, String[] occurrences) {

        System.out.println("Error output : " + stderr);
		
		// current stderr output
		//		[INFO] 	7897 reads, 8009 mappings: R-factor 1.0141826
		//		[INFO] 	6145 entire, 1864 split mappings (2.3273816%)
		//
		//		1 single transcript loci
		//		8009 mappings in file
		//		566 mappings fall in single transcript loci
		//		283 mappings map to annotation
		//		1172 mappings in annotation-mapped pairs
		//		40,75 min/max read length
		//
		//		8009 mappings read from file
		//		8044 mapping pairs map to annotation
		//		0 pairs without tx evidence
		//		208 pairs in wrong orientation
		//		0 single mappings forced
		//
		//		4 transcripts, 4 detected
				
		//		println "EXIT: ${process.exitValue()}"
		//		println "STDOUT: ${process.in.text}"
		//		// alternative to combine wait and stdout
		//		println "ls ${parameter}".execute().text
		for (int i = 0; i < occurrences.length; i++) {
			assertTrue(stderr.contains(occurrences[i]))
		}

		if (occurrences!= STDERR_ACCESS_DENIED)
			assertTrue(FluxCapacitorRunner.getOutFile().exists());

	}

    File currentTestDirectory = null
    @Before
    public void setUpTest(){
        currentTestDirectory = FileHelper.createTempDir("FluxCapacitorIntegration", "", null)
    }

    @After
    public void cleanup(){
        if(currentTestDirectory != null){
            FileHelper.rmDir(currentTestDirectory)
        }
    }

    @Test
    public void testIOflatSortedWritableGTFflatSortedWritableBEDnoKeep() {
        UniversalReadDescriptor descriptor = new UniversalReadDescriptor();
        descriptor.init(UniversalReadDescriptor.getDescriptor("SIMULATOR"));
        FluxCapacitorRunner.createTestDir(currentTestDirectory, [
               "ANNOTATION_FILE" : GTF_SORTED,
               "MAPPING_FILE" : BED_SORTED,
               "ANNOTATION_MAPPING" : AnnotationMapping.PAIRED,
               "READ_DESCRIPTOR" : descriptor,
               "KEEP_SORTED_FILES" : false
        ])
        String stderr= runCapacitor();
        assertFiles(2, 1, stderr, STDERR_MAPPED);
    }

	@Test
	public void testIOgzippedSortedWritableGTFflatSortedWritableBEDnoKeep() {

        UniversalReadDescriptor descriptor = new UniversalReadDescriptor();
        descriptor.init(UniversalReadDescriptor.getDescriptor("SIMULATOR"));
        FluxCapacitorRunner.createTestDir(currentTestDirectory, [
                "ANNOTATION_FILE" : GTF_SORTED,
                "MAPPING_FILE" : BED_SORTED,
                "ANNOTATION_MAPPING" : AnnotationMapping.PAIRED,
                "READ_DESCRIPTOR" : descriptor,
                "KEEP_SORTED_FILES" : false
        ])
        String stderr= runCapacitor();
        assertFiles(2, 1, stderr, STDERR_MAPPED);
	}
	
    @Test
	public void testIOzippedSortedWritableGTFflatSortedWritableBEDnoKeep() {

        UniversalReadDescriptor descriptor = new UniversalReadDescriptor();
        descriptor.init(UniversalReadDescriptor.getDescriptor("SIMULATOR"));
        FluxCapacitorRunner.createTestDir(currentTestDirectory, [
                "ANNOTATION_FILE" : GTF_SORTED,
                "MAPPING_FILE" : BED_SORTED,
                "ANNOTATION_MAPPING" : AnnotationMapping.PAIRED,
                "READ_DESCRIPTOR" : descriptor,
                "KEEP_SORTED_FILES" : false
        ])
        String stderr= runCapacitor();
        assertFiles(2, 1, stderr, STDERR_MAPPED);
	
	}
	
	@Test
	public void testIOflatSortedWritableGTFgzippedSortedWritableBEDnoKeep() {

        UniversalReadDescriptor descriptor = new UniversalReadDescriptor();
        descriptor.init(UniversalReadDescriptor.getDescriptor("SIMULATOR"));
        FluxCapacitorRunner.createTestDir(currentTestDirectory, [
                "ANNOTATION_FILE" : GTF_SORTED,
                "MAPPING_FILE" : BED_SORTED,
                "ANNOTATION_MAPPING" : AnnotationMapping.PAIRED,
                "READ_DESCRIPTOR" : descriptor,
                "KEEP_SORTED_FILES" : false
        ])
        String stderr= runCapacitor();
        assertFiles(2, 1, stderr, STDERR_MAPPED);

	}

	@Test
	public void testIOflatSortedWritableGTFzippedSortedWritableBEDnoKeep() {

        UniversalReadDescriptor descriptor = new UniversalReadDescriptor();
        descriptor.init(UniversalReadDescriptor.getDescriptor("SIMULATOR"));
        FluxCapacitorRunner.createTestDir(currentTestDirectory, [
                "ANNOTATION_FILE" : GTF_SORTED,
                "MAPPING_FILE" : BED_SORTED,
                "ANNOTATION_MAPPING" : AnnotationMapping.PAIRED,
                "READ_DESCRIPTOR" : descriptor,
                "KEEP_SORTED_FILES" : false
        ])
        String stderr= runCapacitor();
        assertFiles(2, 1, stderr, STDERR_MAPPED);

	}

	@Test
	public void testIOgzippedSortedWritableGTFzippedSortedWritableBEDnoKeep() {

        UniversalReadDescriptor descriptor = new UniversalReadDescriptor();
        descriptor.init(UniversalReadDescriptor.getDescriptor("SIMULATOR"));
        FluxCapacitorRunner.createTestDir(currentTestDirectory, [
                "ANNOTATION_FILE" : GTF_SORTED,
                "MAPPING_FILE" : BED_SORTED,
                "ANNOTATION_MAPPING" : AnnotationMapping.PAIRED,
                "READ_DESCRIPTOR" : descriptor,
                "KEEP_SORTED_FILES" : false
        ])
        String stderr= runCapacitor();
        assertFiles(2, 1, stderr, STDERR_MAPPED);

	}

	@Test
	public void testIOzippedSortedWritableGTFgzippedSortedWritableBEDnoKeep() {

        UniversalReadDescriptor descriptor = new UniversalReadDescriptor();
        descriptor.init(UniversalReadDescriptor.getDescriptor("SIMULATOR"));
        FluxCapacitorRunner.createTestDir(currentTestDirectory, [
                "ANNOTATION_FILE" : GTF_SORTED,
                "MAPPING_FILE" : BED_SORTED,
                "ANNOTATION_MAPPING" : AnnotationMapping.PAIRED,
                "READ_DESCRIPTOR" : descriptor,
                "KEEP_SORTED_FILES" : false
        ])
        String stderr= runCapacitor();
        assertFiles(2, 1, stderr, STDERR_MAPPED);

	}

	@Test
	public void testIOflatUnsortedWritableGTFflatSortedWritableBEDnoKeep() {

        UniversalReadDescriptor descriptor = new UniversalReadDescriptor();
        descriptor.init(UniversalReadDescriptor.getDescriptor("SIMULATOR"));
        FluxCapacitorRunner.createTestDir(currentTestDirectory, [
                "ANNOTATION_FILE" : GTF_SORTED,
                "MAPPING_FILE" : BED_SORTED,
                "ANNOTATION_MAPPING" : AnnotationMapping.PAIRED,
                "READ_DESCRIPTOR" : descriptor,
                "KEEP_SORTED_FILES" : false
        ])
        String stderr= runCapacitor();
        assertFiles(2, 1, stderr, STDERR_MAPPED);

	}

	@Test
	public void testIOflatSortedWritableGTFflatUnsortedWritableBEDnoKeep() {

        UniversalReadDescriptor descriptor = new UniversalReadDescriptor();
        descriptor.init(UniversalReadDescriptor.getDescriptor("SIMULATOR"));
        FluxCapacitorRunner.createTestDir(currentTestDirectory, [
                "ANNOTATION_FILE" : GTF_SORTED,
                "MAPPING_FILE" : BED_SORTED,
                "ANNOTATION_MAPPING" : AnnotationMapping.PAIRED,
                "READ_DESCRIPTOR" : descriptor,
                "KEEP_SORTED_FILES" : false
        ])
        String stderr= runCapacitor();
        assertFiles(2, 1, stderr, STDERR_MAPPED);
	}

	@Test
	public void testIOflatUnSortedWritableGTFflatUnsortedWritableBEDkeep() {

        UniversalReadDescriptor descriptor = new UniversalReadDescriptor();
        descriptor.init(UniversalReadDescriptor.getDescriptor("SIMULATOR"));
        FluxCapacitorRunner.createTestDir(currentTestDirectory, [
                "ANNOTATION_FILE" : GTF_SORTED,
                "MAPPING_FILE" : BED_SORTED,
                "ANNOTATION_MAPPING" : AnnotationMapping.PAIRED,
                "READ_DESCRIPTOR" : descriptor,
                "KEEP_SORTED_FILES" : false
        ])
        String stderr= runCapacitor();
        assertFiles(2, 1, stderr, STDERR_MAPPED);

	}

	@Test
	public void testIOflatUnSortedReadOnlyGTFflatUnsortedReadOnlyBEDkeep() {

        UniversalReadDescriptor descriptor = new UniversalReadDescriptor();
        descriptor.init(UniversalReadDescriptor.getDescriptor("SIMULATOR"));
        FluxCapacitorRunner.createTestDir(currentTestDirectory, [
                "ANNOTATION_FILE" : GTF_SORTED,
                "MAPPING_FILE" : BED_SORTED,
                "ANNOTATION_MAPPING" : AnnotationMapping.PAIRED,
                "READ_DESCRIPTOR" : descriptor,
                "KEEP_SORTED_FILES" : false
        ])
        String stderr= runCapacitor();
        assertFiles(2, 1, stderr, STDERR_MAPPED);

	}
	
	@Test
	public void testIOzippedUnSortedReadOnlyGTFgzippedUnsortedReadOnlyBEDnoKeep() {

        UniversalReadDescriptor descriptor = new UniversalReadDescriptor();
        descriptor.init(UniversalReadDescriptor.getDescriptor("SIMULATOR"));
        FluxCapacitorRunner.createTestDir(currentTestDirectory, [
                "ANNOTATION_FILE" : GTF_SORTED,
                "MAPPING_FILE" : BED_SORTED,
                "ANNOTATION_MAPPING" : AnnotationMapping.PAIRED,
                "READ_DESCRIPTOR" : descriptor,
                "KEEP_SORTED_FILES" : false
        ])
        String stderr= runCapacitor();
        assertFiles(2, 1, stderr, STDERR_MAPPED);

	}

	@Test
	public void testTmpDir() {

        UniversalReadDescriptor descriptor = new UniversalReadDescriptor();
        descriptor.init(UniversalReadDescriptor.getDescriptor("SIMULATOR"));
        FluxCapacitorRunner.createTestDir(currentTestDirectory, [
                "ANNOTATION_FILE" : GTF_SORTED,
                "MAPPING_FILE" : BED_SORTED,
                "ANNOTATION_MAPPING" : AnnotationMapping.PAIRED,
                "READ_DESCRIPTOR" : descriptor,
                "KEEP_SORTED_FILES" : false,
                "TMP_DIR" : new File(System.getProperty("java.io.tmpdir"))
        ])
        String stderr= runCapacitor();
        assertFiles(2, 1, stderr, STDERR_ACCESS_DENIED);

	}

}