package test;


import barna.genome.io.FileHelper
import barna.genome.sequencing.rnaseq.simulation.FluxSimulatorSettings
import org.junit.Test
import static org.junit.Assert.assertTrue

class FluxSimulatorTest {

	final File MINIMAL_PAR= new File(getClass().getResource("/minimal.par").getFile());
	final File SPIKE_SEQUENCES= new File(getClass().getResource("/spike_sequences.gtf").getFile());
	final File GENOME_SPIKE= new File(getClass().getResource("/genome_spikes").getFile());
	
	@Test
	public void testTmpDirFail() {
	
		File tmpDir= new File(System.getProperty("java.io.tmpdir")); //FileHelper.createTempDir(getClass().getSimpleName(), "_tmpdir", null)
		FluxSimulatorSettings settings= createSettings(tmpDir)
		File parFile= FileHelper.createTempFile(getClass().getSimpleName(), ".par")
		writeParFile(parFile, settings)
		
		String stdErr= runSimulator(parFile, tmpDir)
		String[] denied= ["access denied"]
		assertResult(stdErr, denied)
	}
	
	@Test
	public void testTmpDirSuccess() {
	
		File tmpDir= FileHelper.createTempDir(getClass().getSimpleName(), "_tmpdir", null)
		FluxSimulatorSettings settings= createSettings(tmpDir)
		File parFile= FileHelper.createTempFile(getClass().getSimpleName(), ".par")
		writeParFile(parFile, settings)
		
		String stdErr= runSimulator(parFile, tmpDir)
		String[] end= ["[END]"]
		assertResult(stdErr, end)
	}
	
	protected void assertResult(String stream, String[] contains) {
		for (int i = 0; i < contains.length; i++) {
			assertTrue(stream.contains(contains[i]))
		}
	}
	
	protected FluxSimulatorSettings createSettings(File tmpDir) {
		FluxSimulatorSettings settings= new FluxSimulatorSettings()
//		## File locations
//REF_FILE_NAME   spike_sequences.gtf
//GEN_DIR         genome_Spikes
//
//## Expression
//NB_MOLECULES	10000
//
//
//## Sequencing
//READ_NUMBER	700176
//READ_LENGTH	36
//PAIRED_END	YES
		settings.set(FluxSimulatorSettings.REF_FILE,
				new File(getClass().getResource("/spike_sequences.gtf").getFile()))
		settings.set(FluxSimulatorSettings.GEN_DIR,
				new File(getClass().getResource("/genome_spikes").getFile()))
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
	
	protected String runSimulator(File parFile, File tmpDir) {
		String cmd= "java -cp "+System.getProperty("java.class.path")
		if (tmpDir!= null)
			cmd+= " -Dflux.io.deny.tmpdir=yes"
		cmd+= " -Xmx1G barna.commons.launcher.Flux -t simulator -p "+parFile.getAbsolutePath()

		Process process= cmd.execute()
		process.waitFor()

		String stderr= "STDERR: ${process.err.text}"
		return stderr;

	}
}
