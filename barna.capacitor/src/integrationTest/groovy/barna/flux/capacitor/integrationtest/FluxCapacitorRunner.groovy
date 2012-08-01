package barna.flux.capacitor.integrationtest

import barna.commons.system.OSChecker
import groovy.json.JsonSlurper

import java.security.MessageDigest
import java.util.zip.ZipInputStream

import com.martiansoftware.jsap.RequiredParameterMissingException
import barna.io.FileHelper
import barna.flux.capacitor.reconstruction.FluxCapacitorSettings
import barna.io.rna.UniversalReadDescriptor

/**
 *
 * @author  Emilio Palumbo (emiliopalumbo@gmail.com)
 */
class FluxCapacitorRunner {
    
    /**
     * Map from the relative filename to the absolute filename within
     * the test data directory
     */
    private static Map testData
    /**
     * The default URL used to download Capacitor test data
     */
    private static String artifactoryUrl;

    /**
     * The repository to use
     */
    private static String repository

    /**
     * The test data artifact
     */
    private static String artifact

    /**
     * The test data target directory
     */
    private static String targetDirectory

    static {
        // initialize defaults
        artifactoryUrl = System.getProperty("testdata.artifactory", "http://sammeth.net/artifactory")
        repository = System.getProperty("testdata.repository", "repo")
        artifact = System.getProperty("testdata.artifact", "barna/test_data-1.0.zip")
        targetDirectory = System.getProperty("testdata.target", new File("").getAbsolutePath())
    }

    /**
     * The path to the capacitor executable
     */
    private static String executable

    /**
     * Synchronization lock
     */
    private static final Object lock = new Object()

    /**
     * Creates a temporary directory structure for the current run of the Capacitor
     *
     * @param cwd current working directory
     * @param parameters Map of Capacitor parameters
     * @throws RequiredParameterMissingException if param name is not a known parameter name
     * @return the current paramter file
     */
    static File createTestDir(File cwd, Map parameters) throws RequiredParameterMissingException {

        //checking for mandatory parameters
        if (!parameters.containsKey("ANNOTATION_FILE"))
            throw new RequiredParameterMissingException("The parameter for annotation file is missing")
        if (!parameters.containsKey("MAPPING_FILE"))
            throw new RequiredParameterMissingException("The parameter for mapping file is missing")

        //getting instance for the read descriptor
        if (parameters.containsKey("READ_DESCRIPTOR")) {
            UniversalReadDescriptor descriptor = new UniversalReadDescriptor();
            descriptor.init(UniversalReadDescriptor.getDescriptor("SIMULATOR"));
        }

        //setting up the output file
        File outDir = new File(cwd, "output");
        outDir.mkdir();
        if (!parameters.containsKey("STDOUT_FILE")) {
            File outFile = new File(outDir,"FluxCapacitor.gtf");
            parameters.put("STDOUT_FILE", outFile);
            outFile.delete();
        }

        //writing the parameter file
        File parFile= new File(cwd,"parameters.par");
        FluxCapacitorSettings settings= new FluxCapacitorSettings();
        parameters.each{k,v->
            settings.set(k,v)
        }
        settings.validate();
        parameters.each{k,v->
            parFile.append("${k}\t${v}\n")
        }

        return parFile;
    }

    /**
     * Run the capacitor as external process in the given working directory
     * with the specified parameter file.
     *
     * @param cwd working directory
     * @param parameterFile the parameter file
     * @param denyTempDir if true, access to system tmp dir will be forbidden
     * @return output the output of the capacitor
     */
    public static String runCapacitor(File cwd, File parameterFile, boolean denyTempDir = false){
        synchronized (lock){
            if(executable == null){
                executable = System.getProperty("dist.exe")
                if(executable == null){
                    throw new RuntimeException("No capacitor executable specified")
                }
            }
        }

        def pb = new ProcessBuilder()
        pb.environment().put("FLUX_MEM", "1G")

        if (denyTempDir){
            pb.environment().put("JAVA_OPTS", "-Dflux.io.deny.tmpdir=yes")
        }

        def cmd = [executable, "-p", parameterFile.getAbsolutePath()]
        if (OSChecker.isWindows()) {
            cmd = ["cmd", "/c", executable, "-p", parameterFile.getAbsolutePath()]
        }
        def process = pb.directory(cwd)
                .redirectErrorStream(true)
                .command(cmd)
                .start()
        String output = process.inputStream.text
        process.waitFor()
        return output;
    }

    /**
     * Access the test data map. The map is initialized lazily. If test data
     * is not available locally, the data file is downloaded from the repository.
     * <p>
     *     The data map contains relative paths to the test data pointing to the absolut
     *     files and directories. For example, to access a gtf file
     *     {@code FluxCapacitorRunner.getTestData()['gtf/mydata.gtf']} returns the absolute path
     *     to 'mydata.gtf'
     *
     *
     *
     *
     * @return
     */
    static synchronized Map getTestData() {
        if (this.testData == null) {
            this.testData = [:]
            File testDir = prepareTestData()
            testDir.eachFileRecurse {file ->
                def relPath = file.getAbsolutePath().replace(testDir.getAbsolutePath(), "")
                relPath = relPath.replaceAll("\\\\", "/")
                if(relPath.startsWith("/")) relPath = relPath.substring(1)
                this.testData[relPath] = file.getAbsolutePath()
            }
        }
        return this.testData
    }

    /**
     * Checks for valid local test data or downloads test data from the repository
     *
     * @return test_data the unzipped test_data directory
     */
    private static File prepareTestData() {
        println("Checking for test data")
        def slurper = new JsonSlurper()
        def dataFQN = new URL("${artifactoryUrl}/api/storage/${repository}/${artifact}")
        def fileName = dataFQN.getFile().split("/").last()


        def targetFile = new File(targetDirectory, fileName)
        File test_data_dir = new File(targetDirectory, "test_data")
        def md5File = new File(targetDirectory, "${fileName}.md5")

        def metaData = slurper.parseText(dataFQN.openStream().text)
        def md5sum = metaData.checksums['md5']

        boolean invalid_file = false
        if (!targetFile.exists() || !md5File.exists() || md5File.readLines()[0].trim() != md5sum) {
            invalid_file = true
        }

        if (invalid_file) {
            println("Downloading test data from ${dataFQN.toExternalForm()}")
            if (targetFile.exists()) {
                targetFile.delete()
            }
            if (test_data_dir.exists()) {
                FileHelper.rmDir(new File(targetDirectory, "test_data"))
            }
            if (md5File.exists()) {
                md5File.delete()
            }
            // download
            def out = new BufferedOutputStream(new FileOutputStream(targetFile))
            out << new URL(metaData['downloadUri']).openStream()
            out.close()

            // compare md5
            if (md5sum != generateMD5(targetFile)) {
                throw new RuntimeException("Test Data downlaoded but md5 sums do not match !")
            }

            // write md5 file
            md5File.write(md5sum)
        }


        if (!test_data_dir.exists()) {
            println "Unzipping test data"
            unzip(targetFile, targetDirectory)
        }
        println("Test data available")
        return test_data_dir
    }

    /**
     * Helper method to generate the md5 sum for a given file
     *
     * @param file the file
     * @return md5sum the hex representation of the md5sum
     */
    private static String generateMD5(final file) {
        MessageDigest digest = MessageDigest.getInstance("MD5")
        file.withInputStream() {is ->
            byte[] buffer = new byte[8192]
            int read = 0
            while ((read = is.read(buffer)) > 0) {
                digest.update(buffer, 0, read);
            }
        }
        byte[] md5sum = digest.digest()
        BigInteger bigInt = new BigInteger(1, md5sum)
        return bigInt.toString(16).padLeft(32, '0')
    }
    /**
     * Helper to unzip unzip a file
     */
    private static unzip = { File file, String dest ->
        def result = new ZipInputStream(new FileInputStream(file))
        def destFile = new File(dest)
        if (!destFile.exists()) {
            destFile.mkdir();
        }
        result.withStream {
            def entry
            while (entry = result.nextEntry) {
                if (!entry.isDirectory()) {
                    new File(dest + File.separator + entry.name).parentFile?.mkdirs()
                    def output = new FileOutputStream(dest + File.separator
                            + entry.name)
                    output.withStream {
                        int len = 0;
                        byte[] buffer = new byte[4096]
                        while ((len = result.read(buffer)) > 0) {
                            output.write(buffer, 0, len);
                        }
                    }
                } else {
                    new File(dest + File.separator + entry.name).mkdir()
                }
            }
        }
    }

    public static void main(String[] args) {
        //println getTestData()
        FluxCapacitorSettings s =  new FluxCapacitorSettings();
        s.set("ANNOTATION_FILE", new File(""));
        s.set("ANNOTATION_FILE", "ABF");
    }
}
