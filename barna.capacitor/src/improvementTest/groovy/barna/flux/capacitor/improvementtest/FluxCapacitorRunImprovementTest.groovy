package barna.flux.capacitor.improvementtest

import barna.commons.Execute
import org.rosuda.JRI.Rengine

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
        }*/
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
            }*/

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
        println "======Building transcripts RPKM table======"
        def file = []
        def i = 0
        file[0] = "TargetID\tGene_Symbol\tChr\tCoord"
        currentTestDirectory.eachDir {
            def sample = it.toString().split("/").last().substring(0,26)
            file[0]= file[0].toString().concat("\t$sample")
            File gtf = new File(it.getAbsolutePath()+"/output/results.gtf")
            def j = 1
            gtf.eachLine {
                def tokens = it.split("\t")
                def chr = tokens[0].replace("chr","")
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
                j++
            }
            i++
        }
        def table = new File(currentTestDirectory.absolutePath+"/transcriptsRPKM.txt")
        if (table.exists())
            table.delete()
        file.each {
            table.append(it)
            table.append(barna.commons.system.OSChecker.NEW_LINE)
        }
        println "======Building correlation matrix of samples======"
        def corMat = new File(currentTestDirectory.absolutePath+"/corMat.txt")
        def pb = new ProcessBuilder()
        def cmd = ["Rscript", FluxCapacitorRunner.testData["R/matrixCalculateCorOps.R"], table.absolutePath, corMat.absolutePath,"1",SAMPLES.toString()]
        /*def process = pb.directory(new File(FluxCapacitorRunner.testData["R"]))
                .redirectErrorStream(true)
                .command(cmd)
                .start()
        String output = process.inputStream.text
        process.waitFor()*/
        //Runtime.getRuntime().exec("Rscript " + FluxCapacitorRunner.testData["R/matrixCalculateCorOps.R"] + " " + table.absoluteFile + " " + corMat.absoluteFile + " 1 " + SAMPLES, null, new File(FluxCapacitorRunner.testData["R"]))
        println "======Performing SVD of the corelation matrix======"
        def args = new String[1]
        args[0] = "--vanilla"
        Rengine re = new Rengine(args,false,null)
        re.eval("mat_raw_trrpkm <- read.table(\"" + corMat.absolutePath + "\", header=T, row.names=1, sep=\"\t\")")
        re.eval("data <- list(mat_raw_trrpkm)")
        re.eval("names <- c(\"TR_RPKM by Lab - chr1\",\"TR_RPKM by Population\")")
        re.eval("ind <- rownames(mat_raw_trrpkm)")
        re.eval("qc <- read.table(\"" + FluxCapacitorRunner.testData["stats/mRNA.stats"] + "\", header=T, row.names=1, sep=\"\t\")")
        re.eval("qcind <- rownames(qc)")
        re.eval("qc <- qc[qcind,]")
        re.eval("data <- lapply(data, function(x) {x[qcind,qcind]})")
        re.eval("colvec <- as.character(qc[,\"LabColor\"])")
        re.eval("m <- cmdscale(as.dist(1-data[[1]]))")
        def svd = re.eval("m")
//        re.eval("X11()")
//        re.eval("plot(m[,1],-m[,2], main=names[1], xlab=\"1st dimension\", ylab=\"2nd dimension\", cex.axis=1.2, cex.lab=1.2, col=colvec)")
//        re.eval("colvec1 <- as.character(qc[,\"PopColor\"])")
//        re.eval("X11()")
//        re.eval("plot(m[,1], -m[,2], main=names[2], xlab=\"1st dimension\", ylab=\"2nd dimension\", cex.axis=1.2, cex.lab=1.2, col=colvec1)")

        println "======Evaluating dataset======"
        def dist = 0;
        svd.asDoubleMatrix().each {
            dist += Math.sqrt(Math.pow(it[0],  2) + Math.pow(it[1], 2))
        }
        println "Overall distance: $dist"

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
