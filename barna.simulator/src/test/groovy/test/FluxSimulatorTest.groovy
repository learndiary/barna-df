/*
 * Copyright (c) 2010, Micha Sammeth
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * The names of its contributors may be not used to endorse or promote
 *       products derived from this software without specific prior written
 *       permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL MICHA SAMMETH BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

package test;


import barna.flux.simulator.FluxSimulatorSettings
import barna.io.FileHelper
import org.junit.Test

import static org.junit.Assert.assertTrue

class FluxSimulatorTest {

	final File MINIMAL_PAR= new File(getClass().getResource("/minimal.par").getFile());
	final File SPIKE_SEQUENCES= new File(getClass().getResource("/spike_sequences.gtf").getFile());
	final File GENOME_SPIKE= new File(getClass().getResource("/genome_Spikes").getFile());
	
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
		println (stdErr)
		//assertResult(stdErr, end)
	}
	
	protected void assertResult(String stream, String[] contains) {
		for (int i = 0; i < contains.length; i++) {		
			assertTrue("Stream ${stream} does not contain ${contains[i]}", stream.contains(contains[i]))
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
	
	protected String runSimulator(File parFile, File tmpDir) {
		String cmd= "java -cp ";
        String cp= System.getProperty("java.class.path");
        String[] cpp= cp.split(":");
        String cp2="";
        for(int i= 0; i< cpp.length; ++i) {
            if (cpp[i].indexOf(" ")< 0)
                cp2+= ":"+ cpp[i];
        }
        cmd+= cp2.substring(1);
		if (tmpDir!= null)
			cmd+= " -Dflux.io.deny.tmpdir=yes"
		cmd+= " -Xmx1G barna.commons.launcher.Flux -t simulator -p "+parFile.getAbsolutePath()

		Process process= cmd.execute()
		process.waitFor()

		String stderr= "STDERR: ${process.err.text}"
		return stderr;

	}
}
