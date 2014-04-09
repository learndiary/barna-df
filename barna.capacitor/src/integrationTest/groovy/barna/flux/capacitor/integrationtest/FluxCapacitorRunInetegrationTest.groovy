	
package barna.flux.capacitor.integrationtest

import barna.commons.Execute
import barna.commons.system.OSChecker
import barna.flux.capacitor.reconstruction.FluxCapacitorSettings.AnnotationMapping
import barna.io.FileHelper
import barna.model.rna.UniversalReadDescriptor
import org.junit.*

import static junit.framework.Assert.assertTrue
import static org.junit.Assert.fail

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
            def resource = FluxCapacitorRunInetegrationTest.class.getResource("/integration-build.properties")
            if(!resource){
                fail("No capacitor executable specified")
            }else{
                Properties p = new Properties()
                p.load(resource.openStream())
                executable = p.getProperty("dist.exe", null)
                if(executable == null){
                    fail("No capacitor executable specified and integration tests properties do not contain a path")
                }
            }
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

	static final String[] STDERR_MAPPED= ["8005","8184"]
	static final String[] STDERR_ACCESS_DENIED= ["access denied"]
	
	void assertStdErr(String stderr, String[] occurrences, Boolean debug = false) {
        if (debug) System.err.println(stderr)

		for (int i = 0; i < occurrences.length; i++) {
			assertTrue(stderr.contains(occurrences[i]))
		}
	}

    static void assertFileExist(File cwd, Map files){

        for (Map.Entry e  : files.entrySet()) {
            String fileName = e.key.toString();
            File file = new File(fileName)
            boolean isRelative = true;
            for (File root : File.listRoots()) {
                if(fileName.startsWith(root.getAbsolutePath())){
                    isRelative = false
                    break;
                }
            }
            if(isRelative) {
                file = new File(cwd, fileName)
            }
            if(e.value instanceof Closure){
                Closure v = e.value
                if(!v(file)) fail("""
            File does not exsits.
            Expected: ${fileName}
        """)
            }
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
    public void testIOflatSortedGTFflatSortedBEDnoKeep() {
        File parFile = FluxCapacitorRunner.createTestDir(currentTestDirectory, [
               "ANNOTATION_FILE" : FluxCapacitorRunner.testData['gtf/mm9_chr1_chrX_sorted.gtf'],
               "MAPPING_FILE" : FluxCapacitorRunner.testData['bed/mm9_chr1_chrX_sorted.bed'],
               "ANNOTATION_MAPPING" : AnnotationMapping.PAIRED,
               "READ_DESCRIPTOR" : "SIMULATOR",
        ])

        String stderr= FluxCapacitorRunner.runCapacitor(currentTestDirectory,parFile);

        assertStdErr(stderr, STDERR_MAPPED);
        assertFileExist(currentTestDirectory, [
                (FluxCapacitorRunner.testData['gtf/mm9_chr1_chrX_sorted.gtf']) : {File file -> return file.exists()},
                (FluxCapacitorRunner.testData['bed/mm9_chr1_chrX_sorted.bed']) : {File file -> return file.exists()},
                (FluxCapacitorRunner.DEFAULT_PARAMETER_FILE) : {File file -> return file.exists()},
                (FluxCapacitorRunner.DEFAULT_OUTPUT_FILE) : {File file -> return file.exists()},
        ])

//        Possible use case>
//
//        assertFileExist(currentTestDirectory, [
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
        assertFileExist(currentTestDirectory, [
                (FluxCapacitorRunner.DEFAULT_PARAMETER_FILE) : {File file -> return file.exists()},
                (FluxCapacitorRunner.DEFAULT_OUTPUT_FILE) : {File file -> return file.exists()},
        ])
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
        assertFileExist(currentTestDirectory, [
                (FluxCapacitorRunner.DEFAULT_PARAMETER_FILE) : {File file -> return file.exists()},
                (FluxCapacitorRunner.DEFAULT_OUTPUT_FILE) : {File file -> return file.exists()},
        ])
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
        assertFileExist(currentTestDirectory, [
                (FluxCapacitorRunner.DEFAULT_PARAMETER_FILE) : {File file -> return file.exists()},
                (FluxCapacitorRunner.DEFAULT_OUTPUT_FILE) : {File file -> return file.exists()},
        ])
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
        assertFileExist(currentTestDirectory, [
                (FluxCapacitorRunner.DEFAULT_PARAMETER_FILE) : {File file -> return file.exists()},
                (FluxCapacitorRunner.DEFAULT_OUTPUT_FILE) : {File file -> return file.exists()},
        ])
	}

    @Test
    public void testIOflatUnsortedGTFflatSortedBEDKeep() {
        File parFile = FluxCapacitorRunner.createTestDir(currentTestDirectory, [
                "ANNOTATION_FILE" : FluxCapacitorRunner.testData['gtf/mm9_chr1_chrX.gtf'],
                "MAPPING_FILE" : FluxCapacitorRunner.testData['bed/mm9_chr1_chrX_sorted.bed'],
                "ANNOTATION_MAPPING" : AnnotationMapping.PAIRED,
                "READ_DESCRIPTOR" : "SIMULATOR",
                "KEEP_SORTED" : "tmp_sorted",
        ])

        String stderr= FluxCapacitorRunner.runCapacitor(currentTestDirectory,parFile);

        assertStdErr(stderr, STDERR_MAPPED);
        assertFileExist(currentTestDirectory, [
                (FluxCapacitorRunner.DEFAULT_PARAMETER_FILE) : {File file -> return file.exists()},
                (FluxCapacitorRunner.DEFAULT_OUTPUT_FILE) : {File file -> return file.exists()},
                "tmp_sorted/mm9_chr1_chrX_sorted.gtf" : {File file -> return file.exists()},
        ])
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
        assertFileExist(currentTestDirectory, [
                (FluxCapacitorRunner.DEFAULT_PARAMETER_FILE) : {File file -> return file.exists()},
                (FluxCapacitorRunner.DEFAULT_OUTPUT_FILE) : {File file -> return file.exists()},
        ])
	}

    @Test
    public void testIOflatSortedGTFflatUnsortedBEDKeep() {
        File parFile = FluxCapacitorRunner.createTestDir(currentTestDirectory, [
                "ANNOTATION_FILE" : FluxCapacitorRunner.testData['gtf/mm9_chr1_chrX_sorted.gtf'],
                "MAPPING_FILE" : FluxCapacitorRunner.testData['bed/mm9_chr1_chrX.bed'],
                "ANNOTATION_MAPPING" : AnnotationMapping.PAIRED,
                "READ_DESCRIPTOR" : "SIMULATOR",
                "KEEP_SORTED" : "tmp_sorted",
        ])

        String stderr= FluxCapacitorRunner.runCapacitor(currentTestDirectory,parFile);

        assertStdErr(stderr, STDERR_MAPPED);
        assertFileExist(currentTestDirectory, [
                (FluxCapacitorRunner.DEFAULT_PARAMETER_FILE) : {File file -> return file.exists()},
                (FluxCapacitorRunner.DEFAULT_OUTPUT_FILE) : {File file -> return file.exists()},
                "tmp_sorted/mm9_chr1_chrX_sorted.bed" : {File file -> return file.exists()},
        ])
    }

    @Test
    public void testIOflatSortedGTFflatUnsortedBEDKeepAbsolute() {
        File parFile = FluxCapacitorRunner.createTestDir(currentTestDirectory, [
                "ANNOTATION_FILE" : FluxCapacitorRunner.testData['gtf/mm9_chr1_chrX_sorted.gtf'],
                "MAPPING_FILE" : FluxCapacitorRunner.testData['bed/mm9_chr1_chrX.bed'],
                "ANNOTATION_MAPPING" : AnnotationMapping.PAIRED,
                "READ_DESCRIPTOR" : "SIMULATOR",
                "KEEP_SORTED" : new File(System.getProperty("java.io.tmpdir")).getAbsolutePath(),
        ])

        String stderr= FluxCapacitorRunner.runCapacitor(currentTestDirectory,parFile);

        try {
        assertStdErr(stderr, STDERR_MAPPED);
        assertFileExist(currentTestDirectory, [
                (FluxCapacitorRunner.DEFAULT_PARAMETER_FILE) : {File file -> return file.exists()},
                (FluxCapacitorRunner.DEFAULT_OUTPUT_FILE) : {File file -> return file.exists()},
                    (System.getProperty("java.io.tmpdir")+File.separator+"mm9_chr1_chrX_sorted.bed") : {File file -> return file.exists()},
        ])
        } catch (Exception e) {
        }
        finally {
            File f = new File(System.getProperty("java.io.tmpdir")+File.separator+"mm9_chr1_chrX_sorted.bed");
            if (f.exists())
                f.delete();
        }
    }

	@Test
	public void testIOflatUnSortedGTFflatUnsortedBEDkeep() {
        File parFile = FluxCapacitorRunner.createTestDir(currentTestDirectory, [
                "ANNOTATION_FILE" : FluxCapacitorRunner.testData['gtf/mm9_chr1_chrX.gtf'],
                "MAPPING_FILE" : FluxCapacitorRunner.testData['bed/mm9_chr1_chrX.bed'],
                "ANNOTATION_MAPPING" : AnnotationMapping.PAIRED,
                "READ_DESCRIPTOR" : "SIMULATOR",
                "KEEP_SORTED" : "tmp_sorted",
        ])

        String stderr= FluxCapacitorRunner.runCapacitor(currentTestDirectory,parFile);

        assertStdErr(stderr, STDERR_MAPPED);
        assertFileExist(currentTestDirectory, [
                (FluxCapacitorRunner.DEFAULT_PARAMETER_FILE) : {File file -> return file.exists()},
                (FluxCapacitorRunner.DEFAULT_OUTPUT_FILE) : {File file -> return file.exists()},
                "tmp_sorted/mm9_chr1_chrX_sorted.gtf" : {File file -> return file.exists()},
                "tmp_sorted/mm9_chr1_chrX_sorted.bed" : {File file -> return file.exists()},
        ])
	}

	@Test
	public void testTmpDir() {
        UniversalReadDescriptor descriptor = new UniversalReadDescriptor(UniversalReadDescriptor.getDescriptor("SIMULATOR"));
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