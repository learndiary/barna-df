package barna.flux.capacitor.improvementtest

import barna.commons.system.OSChecker
import barna.flux.capacitor.reconstruction.FluxCapacitorSettings
import barna.io.rna.UniversalReadDescriptor
import com.martiansoftware.jsap.RequiredParameterMissingException
import org.kamranzafar.jtar.TarInputStream

import java.util.zip.GZIPInputStream

/**
 *
 * @author  Emilio Palumbo (emiliopalumbo@gmail.com)
 */
class FluxCapacitorRunner {

    /**
     * Default common files
     */
    public static DEFAULT_PARAMETER_FILE = "params.par"
    public static DEFAULT_OUTPUT_FILE = "output"+File.separator+"results.gtf"

    /**
     * Map from the relative filename to the absolute filename within
     * the test data directory
     */
    private static Map testData

    /**
     * The test data target directory
     */
    private static String testDir

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

        //check for mandatory parameters
        if (!parameters.containsKey("ANNOTATION_FILE"))
            throw new RequiredParameterMissingException("The parameter for annotation file is missing")
        if (!parameters.containsKey("MAPPING_FILE"))
            throw new RequiredParameterMissingException("The parameter for mapping file is missing")

        //get instance for the read descriptor
        if (parameters.containsKey("READ_DESCRIPTOR")) {
            UniversalReadDescriptor descriptor = new UniversalReadDescriptor();
            descriptor.init(UniversalReadDescriptor.getDescriptor("SIMULATOR"));
        }

        //check if sorted files should be kept and set up directory
        if (parameters.containsKey("KEEP_SORTED")) {
            String sortedPath = parameters['KEEP_SORTED']
            if (!sortedPath.startsWith(File.separator)) {
                File sortDir = new File(cwd, parameters["KEEP_SORTED"]);
                sortDir.mkdir();
                parameters.put("KEEP_SORTED", sortDir);
            }
        }

        //set up the output file
        File outDir = new File(cwd, "output");
        outDir.mkdir();
        if (!parameters.containsKey("STDOUT_FILE")) {
            File outFile = new File(cwd, DEFAULT_OUTPUT_FILE);
            parameters.put("STDOUT_FILE", outFile);
            outFile.delete();
        }

        //write the parameter file
        File parFile= new File(cwd,DEFAULT_PARAMETER_FILE);
        FluxCapacitorSettings settings= new FluxCapacitorSettings();
        parameters.each{k,v->
            settings.set(k,v)
        }
        settings.validate();
        parameters.each{k,v->
            parFile.append("${k}\t${v}"+OSChecker.NEW_LINE)
        }

        return parFile;
    }

    /**
     * Run the capacitor as external process in the given working directory
     * with the specified parameter file.
     *
     * @param cwd working directory
     * @param parameterFile the parameter file
     * @return output the output of the capacitor
     */
    public static InputStream runCapacitor(File cwd, File parameterFile){
        synchronized (lock){
            if(executable == null){
                executable = System.getProperty("dist.exe")
                if(executable == null){
                    if(executable == null){
                        // try to load from properties
                        def resource = FluxCapacitorRunner.class.getResource("/integration-build.properties")
                        if(!resource){
                            throw new RuntimeException("No capacitor executable specified")
                        }else{
                            Properties p = new Properties()
                            p.load(resource.openStream())
                            executable = p.getProperty("dist.exe", null)
                            if(executable == null){
                                throw new RuntimeException("No capacitor executable specified and integration tests properties do not contain a path")
                            }
                        }
                    }


                }
            }
        }

        def pb = new ProcessBuilder()
        pb.environment().put("FLUX_MEM", "2G")

        def cmd = [executable, "-p", parameterFile.getAbsolutePath()]
        if (OSChecker.isWindows()) {
            cmd = ["cmd", "/c", executable, "-p", parameterFile.getAbsolutePath()]
        }
        def process = pb.directory(cwd)
                .redirectErrorStream(true)
                .command(cmd)
                .start()
        InputStream output = process.inputStream
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
            testDir = System.getProperty("improvedata.path")
            if (testDir==null) {
                throw new RuntimeException("No path for test data specified")
            }
            File testDir = new File(testDir)
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
     * Helper to untar a file
     */
    private static untar = { File file, String dest ->
        def result = new TarInputStream(new BufferedInputStream(new GZIPInputStream(new FileInputStream(file))))
        def destFile = new File(dest)
        if (!destFile.exists()) {
            destFile.mkdir();
        }
        result.withStream {
            def entry
            while ((entry = result.getNextEntry()) != null) {
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
        println getTestData()
    }
}
