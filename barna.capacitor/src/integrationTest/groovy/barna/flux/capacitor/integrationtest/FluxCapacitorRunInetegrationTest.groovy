	
package barna.flux.capacitor.integrationtest

import barna.commons.Execute
import barna.commons.system.OSChecker
import barna.flux.capacitor.reconstruction.FluxCapacitorSettings.AnnotationMapping
import barna.io.FileHelper
import barna.io.rna.UniversalReadDescriptor
import groovy.io.FileType
import org.junit.*

import static junit.framework.Assert.assertTrue
import static junit.framework.Assert.fail

/**
 * 
 * @author Thasso Griebel (thasso.griebel@gmail.com)
 */

class FluxCapacitorRunInetegrationTest {

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
	
	void assertStdErr(String stderr, String[] occurrences, Boolean debug = false) {

        if (debug) System.err.println(stderr)

		for (int i = 0; i < occurrences.length; i++) {
			assertTrue(stderr.contains(occurrences[i]))
		}

	}

    void assertDir(File cwd, List value){

        def allFile = []
        cwd.eachFileRecurse(FileType.FILES, { allFile << it.getAbsolutePath() })

        if (allFile.size()!=value.size()) fail("""
            Differemt number of files.
            Expected: ${value.each {println it}}
            Found: ${allFile.each {println it}}
        """)

//        for (Map.Entry e  : value.entrySet()) {
//            if(e.key.toString().contains(File.separator)){
//                if(e.value instanceof Closure){
//                    def v = e.value
//                    for (String line  : new File(file).readLines()) {
//                        if(!v(line)) fail("LIne comparison failed : ")
//                    }
//                }else{
//
//                }
//            }
//        }
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
    public void testIOflatSortedGTFflatSortedBEDnoKeep() {
        File parFile = FluxCapacitorRunner.createTestDir(currentTestDirectory, [
               "ANNOTATION_FILE" : FluxCapacitorRunner.testData['gtf/mm9_chr1_chrX_sorted.gtf'],
               "MAPPING_FILE" : FluxCapacitorRunner.testData['bed/mm9_chr1_chrX_sorted.bed'],
               "ANNOTATION_MAPPING" : AnnotationMapping.PAIRED,
               "READ_DESCRIPTOR" : "SIMULATOR",
        ])

        String stderr= FluxCapacitorRunner.runCapacitor(currentTestDirectory,parFile);

        assertStdErr(stderr, STDERR_MAPPED);
        assertDir(currentTestDirectory, [
                FluxCapacitorRunner.DEFAULT_PARAMETER_FILE,
                FluxCapacitorRunner.DEFAULT_OUTPUT_FILE,
        ])
//
//        assertDir(currentTestDirectory, [
//                "parameters.par" : {File file -> return file.exists()},
//                "output/result.gtf" : [
//                        "lines":20,
//                        "contains": ["transcript", "exon"],
//                        "equals": "expected-result.gtf",
//                        "md5": "34masd314",
//                        "eachLine" : {line-> line.endsWith("abc")}
//                ]
//        ])

    }

	@Test
	public void testIOgzippedSortedGTFflatSortedBEDnoKeep() {
        File parFile = FluxCapacitorRunner.createTestDir(currentTestDirectory, [
                "ANNOTATION_FILE" : FluxCapacitorRunner.testData['gtf/mm9_chr1_chrX_sorted.gtf.gz'],
                "MAPPING_FILE" : FluxCapacitorRunner.testData['bed/mm9_chr1_chrX_sorted.bed'],
                "ANNOTATION_MAPPING" : AnnotationMapping.PAIRED,
                "READ_DESCRIPTOR" : "SIMULATOR",
        ])

        String stderr= FluxCapacitorRunner.runCapacitor(currentTestDirectory,parFile);

        assertStdErr(stderr, STDERR_MAPPED);
	}

	@Test
	public void testIOflatSortedGTFgzippedSortedBEDnoKeep() {
        File parFile = FluxCapacitorRunner.createTestDir(currentTestDirectory, [
                "ANNOTATION_FILE" : FluxCapacitorRunner.testData['gtf/mm9_chr1_chrX_sorted.gtf'],
                "MAPPING_FILE" : FluxCapacitorRunner.testData['bed/mm9_chr1_chrX_sorted.bed.gz'],
                "ANNOTATION_MAPPING" : AnnotationMapping.PAIRED,
                "READ_DESCRIPTOR" : "SIMULATOR",
        ])
        String stderr= FluxCapacitorRunner.runCapacitor(currentTestDirectory,parFile);

        assertStdErr(stderr, STDERR_MAPPED);
	}

	@Test
	public void testIOgzippedSortedGTFgzippedSortedBEDnoKeep() {
        File parFile = FluxCapacitorRunner.createTestDir(currentTestDirectory, [
                "ANNOTATION_FILE" : FluxCapacitorRunner.testData['gtf/mm9_chr1_chrX_sorted.gtf.gz'],
                "MAPPING_FILE" : FluxCapacitorRunner.testData['bed/mm9_chr1_chrX_sorted.bed.gz'],
                "ANNOTATION_MAPPING" : AnnotationMapping.PAIRED,
                "READ_DESCRIPTOR" : "SIMULATOR",
        ])

        String stderr= FluxCapacitorRunner.runCapacitor(currentTestDirectory,parFile);

        assertStdErr(stderr, STDERR_MAPPED);
	}

	@Test
	public void testIOflatUnsortedGTFflatSortedBEDnoKeep() {
        File parFile = FluxCapacitorRunner.createTestDir(currentTestDirectory, [
                "ANNOTATION_FILE" : FluxCapacitorRunner.testData['gtf/mm9_chr1_chrX.gtf'],
                "MAPPING_FILE" : FluxCapacitorRunner.testData['bed/mm9_chr1_chrX_sorted.bed'],
                "ANNOTATION_MAPPING" : AnnotationMapping.PAIRED,
                "READ_DESCRIPTOR" : "SIMULATOR",
        ])

        String stderr= FluxCapacitorRunner.runCapacitor(currentTestDirectory,parFile);

        assertStdErr(stderr, STDERR_MAPPED);
	}

	@Test
	public void testIOflatSortedGTFflatUnsortedBEDnoKeep() {
        File parFile = FluxCapacitorRunner.createTestDir(currentTestDirectory, [
                "ANNOTATION_FILE" : FluxCapacitorRunner.testData['gtf/mm9_chr1_chrX_sorted.gtf'],
                "MAPPING_FILE" : FluxCapacitorRunner.testData['bed/mm9_chr1_chrX.bed'],
                "ANNOTATION_MAPPING" : AnnotationMapping.PAIRED,
                "READ_DESCRIPTOR" : "SIMULATOR",
        ])

        String stderr= FluxCapacitorRunner.runCapacitor(currentTestDirectory,parFile);

        assertStdErr(stderr, STDERR_MAPPED);
	}

	@Test
	public void testIOflatUnSortedGTFflatUnsortedBEDkeep() {
        File parFile = FluxCapacitorRunner.createTestDir(currentTestDirectory, [
                "ANNOTATION_FILE" : FluxCapacitorRunner.testData['gtf/mm9_chr1_chrX.gtf'],
                "MAPPING_FILE" : FluxCapacitorRunner.testData['bed/mm9_chr1_chrX.bed'],
                "ANNOTATION_MAPPING" : AnnotationMapping.PAIRED,
                "READ_DESCRIPTOR" : "SIMULATOR",
                "KEEP_SORTED_FILES" : "true",
        ])

        String stderr= FluxCapacitorRunner.runCapacitor(currentTestDirectory,parFile);

        assertStdErr(stderr, STDERR_MAPPED);
	}

	@Test
	public void testTmpDir() {
        UniversalReadDescriptor descriptor = new UniversalReadDescriptor();
        descriptor.init(UniversalReadDescriptor.getDescriptor("SIMULATOR"));
        File parFile = FluxCapacitorRunner.createTestDir(currentTestDirectory, [
                "ANNOTATION_FILE" : FluxCapacitorRunner.testData['gtf/mm9_chr1_chrX_sorted.gtf'],
                "MAPPING_FILE" : FluxCapacitorRunner.testData['bed/mm9_chr1_chrX_sorted.bed'],
                "ANNOTATION_MAPPING" : AnnotationMapping.PAIRED,
                "READ_DESCRIPTOR" : "SIMULATOR",
                "TMP_DIR" : new File(System.getProperty("java.io.tmpdir"))
        ])

        String stderr= FluxCapacitorRunner.runCapacitor(currentTestDirectory,parFile,true);

        assertStdErr(stderr, STDERR_ACCESS_DENIED);
	}

}