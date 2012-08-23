package barna.flux.simulator.integration

import barna.commons.system.OSChecker
import barna.flux.simulator.FluxSimulatorSettings
import barna.io.FileHelper
import org.junit.BeforeClass
import org.junit.Test

import static org.junit.Assert.*

/**
 * Created with IntelliJ IDEA.
 * User: thasso
 * Date: 7/31/12
 * Time: 10:06 PM
 * To change this template use File | Settings | File Templates.
 */
class FluxSimulatorIntegrationTest {
    final File SPIKE_SEQUENCES= new File(getClass().getResource("/spike_sequences.gtf").getFile());
    final File GENOME_SPIKE= new File(getClass().getResource("/genome_Spikes").getFile());

    static String executable

    @BeforeClass
    public static void setUp(){
        executable = System.getProperty("dist.exe")
        if(executable == null){
            fail("No simulator executable specified")
        }
    }

    protected FluxSimulatorSettings createSettings(File tmpDir) {
        FluxSimulatorSettings settings= new FluxSimulatorSettings()
        settings.set(FluxSimulatorSettings.REF_FILE,
                new File(getClass().getResource("/spike_sequences.gtf").getFile()))
        settings.set(FluxSimulatorSettings.GEN_DIR,
                new File(getClass().getResource("/genome_Spikes").getFile()))
        settings.set(FluxSimulatorSettings.NB_MOLECULES,
                10000)
        settings.set(FluxSimulatorSettings.READ_NUMBER,
                700176)
        settings.set(FluxSimulatorSettings.READ_LENGTH,
                36)
        settings.set(FluxSimulatorSettings.PAIRED_END,
                true)
        if (tmpDir!= null)
            settings.set(FluxSimulatorSettings.TMP_DIR,
                    tmpDir)

        return settings
    }

    protected void writeParFile(File parFile, FluxSimulatorSettings settings) {
        BufferedWriter buffy= new BufferedWriter(new FileWriter(parFile));
        buffy.write(settings.toString());
        buffy.close();
    }

    public Process runSimulator(File directory, File parameterFile, boolean tmpDirDeny = false){
        def pb = new ProcessBuilder()
        def out = new HashMap<String,String>()
        pb.environment().put("FLUX_MEM", "1G")
        if (tmpDirDeny) {
            pb.environment().put("JAVA_OPTS", "-Dflux.io.deny.tmpdir=yes")
        }
        def cmd = [executable, "-p", parameterFile.getAbsolutePath()]
        if (OSChecker.isWindows()) {
            cmd = ["cmd", "/c", executable, "-p", parameterFile.getAbsolutePath()]
        }
        def process = pb.directory(directory)
                .redirectErrorStream(true)
                .command(cmd)
                .start()
        process.waitFor()
        return process
    }

    @Test
    public void testSimpleRun(){
        def dir = FileHelper.createTempDir("Tests-Barna-166", "", null)
        def targetParams = new File(dir, "params.par")
        try {
            targetParams.append("""
            REF_FILE_NAME   ${SPIKE_SEQUENCES.getAbsolutePath()}
            GEN_DIR         ${GENOME_SPIKE.getAbsolutePath()}

            NB_MOLECULES    1000

            # Fragmentation
            FRAG_SUBSTRATE    RNA
            FRAG_METHOD    UR
            FRAG_UR_ETA     350

            # Reverse Transcription
            RTRANSCRIPTION    NO

            # Amplification
            PCR_DISTRIBUTION default
            GC_MEAN      NaN
            PCR_PROBABILITY  0.05

            # Size Filtering
            FILTERING     NO

            ## Sequencing

            READ_NUMBER    150
            READ_LENGTH    76
            PAIRED_END    NO
            FASTA    YES
            ERR_FILE    76
            """.split(barna.commons.system.OSChecker.NEW_LINE).collect {it.trim()}.join(barna.commons.system.OSChecker.NEW_LINE)) // get rid of the spaces at the beginning
            assertEquals(0, runSimulator(dir, targetParams).exitValue());

        } catch (Exception e) {
            e.printStackTrace()
            fail(e.getMessage())
        }finally{
            if (dir!=null)
                FileHelper.rmDir(dir)
        }
    }

    @Test
    public void testTmpDirFail() {

        File tmpDir= new File(System.getProperty("java.io.tmpdir")); //FileHelper.createTempDir(getClass().getSimpleName(), "_tmpdir", null)
        FluxSimulatorSettings settings= createSettings(tmpDir)
        File parFile= FileHelper.createTempFile(getClass().getSimpleName(), ".par")
        writeParFile(parFile, settings)
        try{
            String stdOut= runSimulator(tmpDir, parFile, true).inputStream.readLines().join(OSChecker.NEW_LINE)
            String[] denied= ["access denied"]
            assertTrue(stdOut,stdOut.contains(denied[0]));
        }catch (Exception e){
            e.printStackTrace()
            fail()
        }finally{
            if (parFile!=null)
                FileHelper.rmDir(parFile)
        }

    }

}

