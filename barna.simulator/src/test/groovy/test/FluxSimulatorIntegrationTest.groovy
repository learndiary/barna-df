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


import barna.flux.simulator.SimulationPipeline
import barna.io.FileHelper
import org.junit.Test

import static org.junit.Assert.fail

class FluxSimulatorIntegrationTest {

	final File MINIMAL_PAR= new File(getClass().getResource("/minimal.par").getFile());
	final File SPIKE_SEQUENCES= new File(getClass().getResource("/spike_sequences.gtf").getFile());
	final File GENOME_SPIKE= new File(getClass().getResource("/genome_Spikes").getFile());
	
    @Test
	public void testBarna166ErrorModelLoading(){
        File params = new File(getClass().getResource("/BARNA-166.par").getFile())
        def dir = FileHelper.createTempDir("Tests-Barna-166", "", null)
        def targetParams = new File(dir, "params.par")

        try {
            FileHelper.copy(params, targetParams)
            FileHelper.copy(SPIKE_SEQUENCES, new File(dir, "id2.gtf"))
            def simulator = new SimulationPipeline()
            simulator.setFile(targetParams)
            simulator.call()
        } catch (Exception e) {
            fail()
        }finally{
            dir.eachFile {it.delete()}
            dir.delete()
        }

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
