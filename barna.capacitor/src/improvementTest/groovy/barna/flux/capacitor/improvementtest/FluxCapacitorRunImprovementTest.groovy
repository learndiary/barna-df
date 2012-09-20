package barna.flux.capacitor.improvementtest

import barna.commons.Execute

import java.util.concurrent.Executors

import org.junit.*

/**
 *
 * @author Thasso Griebel (thasso.griebel@gmail.com)
 */
class FluxCapacitorRunImprovementTest {

    static String executable
    static File tmpDir
    static final THREADS = 5
    static final SAMPLES = 35

    @BeforeClass
    public static void setUp(){
        executable = System.getProperty("dist.exe")
        if(executable == null){
            def resource = FluxCapacitorRunImprovementTest.class.getResource("/integration-build.properties")
            if(!resource){
                Assert.fail("No capacitor executable specified")
            }else{
                Properties p = new Properties()
                p.load(resource.openStream())
                executable = p.getProperty("dist.exe", null)
                if(executable == null){
                    Assert.fail("No capacitor executable specified and integration tests properties do not contain a path")
                }
            }
        }
        Execute.initialize(2);

    }

    @AfterClass
    public static void shutdownExecuter() {
        Execute.shutdown();
    }

	static final String[] STDERR_MAPPED= ["8005","8184"]
	static final String[] STDERR_ACCESS_DENIED= ["access denied"]

	void assertStdErr(String stderr, String[] occurrences, Boolean debug = false) {
        if (debug) System.err.println(stderr)

		for (int i = 0; i < occurrences.length; i++) {
			Assert.assertTrue(stderr.contains(occurrences[i]))
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
                if(!v(file)) Assert.fail("""
            File does not exsits.
            Expected: ${fileName}
        """)
            }
        }
    }

    File currentTestDirectory = null
    @Before
    public void setUpTest(){
        //currentTestDirectory = FileHelper.createTempDir("FluxCapacitorImprovement", "", null)
        currentTestDirectory = new File("/tmp/FluxCapacitorImprovement")
    }

    @After
    public void cleanup(){
        /*if(currentTestDirectory != null){
            FileHelper.rmDir(currentTestDirectory)
        } */
    }

    @Test
    public void testDeconvolution() {
        println "======Running deconvolution======"
        def pool = Executors.newFixedThreadPool(THREADS)
        println "Setting up ${SAMPLES/THREADS} pools with $THREADS slots"
        def bams = new File(FluxCapacitorRunner.testData['bam'])
        def c = 0
        bams.eachFileMatch(~/^.*\.bam$/) {
            def name = it.toString().split("/").last().replace(".bam","")
            //def fileDir = FileHelper.createTempDir(name, "",currentTestDirectory)
            /*File parFile = barna.flux.capacitor.improvementtest.FluxCapacitorRunner.createTestDir(fileDir, [
                    "ANNOTATION_FILE" : barna.flux.capacitor.improvementtest.FluxCapacitorRunner.testData['gtf/gencode_v12.gtf'],
                    "MAPPING_FILE" : it,
                    "ANNOTATION_MAPPING" : FluxCapacitorSettings.AnnotationMapping.PAIRED,
                    "READ_DESCRIPTOR" : "PAIRED",
            ])*/

            def i = c++%THREADS
            /*def t = new Thread() {
                public void run() {
                    barna.flux.capacitor.improvementtest.FluxCapacitorRunner.runCapacitor(fileDir,parFile)
                }
            } */

            //pool.submit(t);
            println "$name submitted to pool"
            if (i == (THREADS-1)) {
                //pool.shutdown()
                print "Executing threads in pool ${c/THREADS}..."
                //while (!pool.isTerminated()) {}
                println "done"
                pool = Executors.newFixedThreadPool(THREADS)
            }
        }
        println "======Deconvolution completed======"
        println "======Building transcripts RPKM table======"
        def i = 0
        def file = []
        file[0] = "TargetID\tGene_Symbol\tChr\tCoord"
        currentTestDirectory.eachDir {
            def sample = it.toString().split("/").last().substring(0,26)
            file[0]= file[0].toString().concat("\t$sample")
            File gtf = new File(it.getAbsolutePath()+"/output/results.gtf")
            def j = 1
            gtf.eachLine {
                def tokens = it.split("\t")
                def chr = tokens[0]
                def coord = tokens[3]
                def options = tokens.last().split(";")
                def optionMap = [:]
                options.each {
                    def option = it.trim().split(" ")
                    optionMap.put(option[0],option[1].replaceAll("\"",""))
                }
                def transcript = optionMap["transcript_id"]
                def gene = optionMap["gene_id"]
                def rpkm = optionMap["RPKM"]
                if (i == 0)
                    file[j]="$transcript\t$gene\t$chr\t$coord"
                file[j] = file[j].toString().concat("\t$rpkm")
                i++
                j++
            }
        }
        File table = new File(currentTestDirectory.absolutePath+"/transcriptsRPKM.txt")
        file.each {
            table.append(it)
            table.append(barna.commons.system.OSChecker.NEW_LINE)
        }
        println "======RPKM table completed======"
        println "======Building correlation matrix of samples======"
        println "======Correlation matrix completed======"
        println "======SVD of the corelation matrix======"
        println "======SVD completed======"
        println "======Evaluating dataset======"
        println "======Evaluation completed======"

        /*File parFile = barna.flux.capacitor.improvementtest.FluxCapacitorRunner.createTestDir(currentTestDirectory, [
               "ANNOTATION_FILE" : barna.flux.capacitor.improvementtest.FluxCapacitorRunner.testData['gtf/mm9_chr1_chrX_sorted.gtf'],
               "MAPPING_FILE" : barna.flux.capacitor.improvementtest.FluxCapacitorRunner.testData['bed/mm9_chr1_chrX_sorted.bed'],
               "ANNOTATION_MAPPING" : FluxCapacitorSettings.AnnotationMapping.PAIRED,
               "READ_DESCRIPTOR" : "SIMULATOR",
        ])

        String stderr= barna.flux.capacitor.improvementtest.FluxCapacitorRunner.runCapacitor(currentTestDirectory,parFile);

        assertStdErr(stderr, STDERR_MAPPED);
        assertFileExist(currentTestDirectory, [
                (barna.flux.capacitor.improvementtest.FluxCapacitorRunner.testData['gtf/mm9_chr1_chrX_sorted.gtf']) : {File file -> return file.exists()},
                (barna.flux.capacitor.improvementtest.FluxCapacitorRunner.testData['bed/mm9_chr1_chrX_sorted.bed']) : {File file -> return file.exists()},
                (barna.flux.capacitor.improvementtest.FluxCapacitorRunner.DEFAULT_PARAMETER_FILE) : {File file -> return file.exists()},
                (barna.flux.capacitor.improvementtest.FluxCapacitorRunner.DEFAULT_OUTPUT_FILE) : {File file -> return file.exists()},
        ]) */

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

}
