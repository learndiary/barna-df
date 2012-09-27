package barna.flux.capacitor.improvementtest

import barna.commons.Execute
import barna.flux.capacitor.reconstruction.FluxCapacitorSettings
import barna.io.FileHelper
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
        currentTestDirectory = FileHelper.createTempDir("FluxCapacitorImprovement", "", null)
    }

    @After
    public void cleanup(){
        if(currentTestDirectory != null){
            FileHelper.rmDir(currentTestDirectory)
        }
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
            def fileDir = FileHelper.createTempDir(name, "",currentTestDirectory)
            File parFile = barna.flux.capacitor.improvementtest.FluxCapacitorRunner.createTestDir(fileDir, [
                    "ANNOTATION_FILE" : barna.flux.capacitor.improvementtest.FluxCapacitorRunner.testData['gtf/gencode_v12.gtf'],
                    "MAPPING_FILE" : it,
                    "ANNOTATION_MAPPING" : FluxCapacitorSettings.AnnotationMapping.PAIRED,
                    "READ_DESCRIPTOR" : "PAIRED",
            ])

            def i = c++%THREADS
            def t = new Thread() {
                public void run() {
                    barna.flux.capacitor.improvementtest.FluxCapacitorRunner.runCapacitor(fileDir,parFile)
                }
            }

            pool.submit(t);
            println "$name submitted to pool"
            if (i == (THREADS-1)) {
                pool.shutdown()
                print "Executing threads in pool ${c/THREADS}..."
                while (!pool.isTerminated()) {}
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
        def cmd = ["/usr/bin/Rscript", FluxCapacitorRunner.testData["R/matrixCalculateCorOps.R"], table.absolutePath, corMat.absolutePath,"1",SAMPLES.toString()]
        def process = pb.directory(new File(FluxCapacitorRunner.testData["R"]))
                .redirectErrorStream(true)
                .command(cmd)
                .start()
        String output = process.inputStream.text
        process.waitFor()
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
        re.eval("m <- cmdscale(as.dist(1-data[[1]]))")
        //Plot for debugging purposes
/*
        re.eval("colvec <- as.character(qc[,\"LabColor\"])")
        re.eval("X11()")
        re.eval("plot(m[,1],-m[,2], main=names[1], xlab=\"1st dimension\", ylab=\"2nd dimension\", cex.axis=1.2, cex.lab=1.2, col=colvec)")
        re.eval("colvec1 <- as.character(qc[,\"PopColor\"])")
        re.eval("X11()")
        re.eval("plot(m[,1], -m[,2], main=names[2], xlab=\"1st dimension\", ylab=\"2nd dimension\", cex.axis=1.2, cex.lab=1.2, col=colvec1)")
*/

        println "======Evaluating dataset======"
        def svd = re.eval("m")
        def overDist = computeDist(svd.asDoubleMatrix(),[0,0])
        println "Overall distance from origin: ${String.format('%1$2.6s',overDist)}\n"


        def labvec = re.eval("as.character(qc[,\"LabColor\"])")
        def labcolors = [] as Set
        def dist = 0;
        labvec.asStringArray().each {labcolors.add(it)}
        labcolors.eachWithIndex {color,n ->
            def ind = labvec.asStringArray().findIndexValues {it == color}
            dist = computeDist(svd.asDoubleMatrix()[ind],computeCentroid(svd.asDoubleMatrix()[ind]))
            println "Lab ${n+1} - overall distance from centroid (normalized): ${String.format('%1$2.6s',dist/overDist)}"
        }

        println ""

        def popvec = re.eval("as.character(qc[,\"PopColor\"])")
        def popcolors = [] as Set
        dist = 0;
        popvec.asStringArray().each {popcolors.add(it)}
        popcolors.eachWithIndex {color,n ->
            def ind = popvec.asStringArray().findIndexValues {it == color}
            dist = computeDist(svd.asDoubleMatrix()[ind],computeCentroid(svd.asDoubleMatrix()[ind]))
            println "Population ${n+1} - overall distance from centroid (normalized): ${String.format('%1$2.6s',dist/overDist)}"
        }
    }

    private double computeDist(matrix, centroid) {
        def dist = 0
        def diff = []
        matrix.each {
            dist += Math.sqrt(Math.pow(it[0]-centroid[0],  2) + Math.pow(it[1]-centroid[1], 2))
        }
        return dist
    }

    private double[] computeCentroid(matrix) {
        def x=0,y=0,n=0
        matrix.each {
            x+=it[0]
            y+=it[1]
            n++
        }
        return [x/n,y/n]
    }

}
