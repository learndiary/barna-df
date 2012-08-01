	
package barna.flux.capacitor.integrationtest

import barna.commons.Execute
import barna.commons.system.OSChecker

import barna.flux.capacitor.reconstruction.FluxCapacitorSettings.AnnotationMapping
import barna.io.FileHelper

import barna.io.rna.UniversalReadDescriptor

import org.junit.*

import static junit.framework.Assert.assertTrue
import static org.junit.Assert.fail
import static junit.framework.Assert.fail

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
	
	void assertStdErr(String stderr, String[] occurrences) {

		for (int i = 0; i < occurrences.length; i++) {
			assertTrue(stderr.contains(occurrences[i]))
		}

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

        File parFile = FluxCapacitorRunner.createTestDir(currentTestDirectory, [
               "ANNOTATION_FILE" : FluxCapacitorRunner.testData['gtf/mm9_chr1_chrX_sorted.gtf'],
               "MAPPING_FILE" : FluxCapacitorRunner.testData['bed/chr1_chrX.bed'],
               "ANNOTATION_MAPPING" : AnnotationMapping.PAIRED,
               "READ_DESCRIPTOR" : "SIMULATOR",
        ])

        String stderr= FluxCapacitorRunner.runCapacitor(currentTestDirectory,parFile);

        assertStdErr(stderr, STDERR_MAPPED);
        assertDir(currentTestDirectory, [
                FluxCapacitorRunner.DEFAULT_PARAMETER_FILE,
                FluxCapacitorRunner.DEFAULT_OUTPUT_FILE,
                "/tmp/asdasd"
        ])

        assertDir(currentTestDirectory, [
                "parameters.par" : {File file -> return file.exists()},
                "output/result.gtf" : [
                        "lines":20,
                        "contains": ["transcript", "exon"],
                        "equals": "expected-result.gtf",
                        "md5": "34masd314",
                        "eachLine" : {line-> line.endsWith("abc")}
                ]
        ])

    }

    void assertDir(File cwd, Map value){

        def allFile = []
        cwd.eachFileRecurse{ allFile << it.getAbsolutePath()}


        for (Map.Entry e  : value.entrySet()) {
            if(e.key.equals("eachLine")){
                if(e.value instanceof Closure){
                    def v = e.value
                    for (String line  : new File(file).readLines()) {
                        if(!v(line)) fail("LIne comparison failed : ")
                    }
                }else{

                }
            }
        }
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
        assertStdErr(2, 1, stderr, STDERR_MAPPED);
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
        assertStdErr(2, 1, stderr, STDERR_MAPPED);
	
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
        assertStdErr(2, 1, stderr, STDERR_MAPPED);

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
        assertStdErr(2, 1, stderr, STDERR_MAPPED);

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
        assertStdErr(2, 1, stderr, STDERR_MAPPED);

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
        assertStdErr(2, 1, stderr, STDERR_MAPPED);

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
        assertStdErr(2, 1, stderr, STDERR_MAPPED);

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
        assertStdErr(2, 1, stderr, STDERR_MAPPED);
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
        assertStdErr(2, 1, stderr, STDERR_MAPPED);

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
        assertStdErr(2, 1, stderr, STDERR_MAPPED);

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
        assertStdErr(2, 1, stderr, STDERR_MAPPED);

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
        assertStdErr(2, 1, stderr, STDERR_ACCESS_DENIED);

	}

}