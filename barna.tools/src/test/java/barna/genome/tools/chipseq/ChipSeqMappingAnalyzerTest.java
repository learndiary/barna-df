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

package barna.genome.tools.chipseq;

import barna.commons.Execute;
import barna.commons.launcher.Flux;
import barna.io.FileHelper;
import barna.io.rna.UniversalReadDescriptor;
import com.martiansoftware.jsap.JSAP;
import com.martiansoftware.jsap.JSAPException;
import com.martiansoftware.jsap.JSAPResult;
import com.martiansoftware.jsap.Parameter;
import org.junit.Test;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Future;

import static org.junit.Assert.fail;

public class ChipSeqMappingAnalyzerTest {

    public static void setUp() {
        Execute.initialize(2);
    }

    public static void shutdown() {
        Execute.shutdown();
    }

    @Test
    public void testPaired() {
        setUp();
        // copy input file
        File f = new File(ChipSeqMappingAnalyzerTest.class.getResource("/Paired_sorted_chrY.bed").getFile());
        long start = System.currentTimeMillis();
        // instantiate and run
        ChipSeqMappingAnalyzer myRun = new ChipSeqMappingAnalyzer();
        myRun.fileInput= f; 
        myRun.descriptor = new UniversalReadDescriptor();
        myRun.descriptor.init(UniversalReadDescriptor.getDescriptor(UniversalReadDescriptor.DESCRIPTORID_PAIRED));

        Future<int[]> captain = Execute.getExecutor().submit(myRun);
        int[] distr = null;
        try {
            distr = captain.get();
        } catch (InterruptedException e) {
            e.printStackTrace();
        } catch (ExecutionException e) {
            e.printStackTrace();
        }

        if (distr != null)
            System.err.println("peak " + distr[1] + ", bounds=[" + distr[0] + "," + distr[2] + "]");
        System.out.println((System.currentTimeMillis()-start) / 1000+"s");
        shutdown();
    }

    @Test
    public void testPairedFlux() {

        // copy input file
        File f = new File(ChipSeqMappingAnalyzerTest.class.getResource("/Paired_sorted_chrY.bed").getFile());
        File tmpF = null;
        try {
            tmpF = FileHelper.createTempFile(FileHelper.stripExtension(f.getName()), FileHelper.getExtension(f));
        } catch (IOException e) {
            e.printStackTrace();
        }
        FileHelper.copy(f, tmpF);


        // write par
        File parF = null;
        BufferedWriter writer = null;
        try {
            parF = FileHelper.createTempFile(ChipSeqMappingAnalyzerTest.class.getSimpleName(), ".par");
            writer = new BufferedWriter(new FileWriter(parF));
            writer.write(ChipSeqSettings.FILE_INPUT.getName() + " " + tmpF.getAbsolutePath() + "\n");
            writer.write(ChipSeqSettings.FILE_OUTPUT.getName() + " " + tmpF.getAbsolutePath()+".out" + "\n");
            writer.write(ChipSeqSettings.READ_DESCRIPTOR.getName() + " " + UniversalReadDescriptor.DESCRIPTORID_PAIRED + "\n");
        } catch (Exception e) {
            e.printStackTrace();
        } finally {
            if (writer != null)
                try {
                    writer.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
        }


        // instantiate and run
        String[] args = {"-t", "mapdist", "-p", parF.getAbsolutePath()};
        long start = System.currentTimeMillis();
        Flux.main(args);

        System.out.println((System.currentTimeMillis()-start) / 1000+"s");

    }
    @Test
    public void testThatTheFileOutputParameterMustBeSpecified() {

        // copy input file
        File f = new File(ChipSeqMappingAnalyzerTest.class.getResource("/Paired_sorted_chrY.bed").getFile());
        File tmpF = null;
        try {
            tmpF = FileHelper.createTempFile(FileHelper.stripExtension(f.getName()), FileHelper.getExtension(f));
        } catch (IOException e) {
            e.printStackTrace();
        }
        FileHelper.copy(f, tmpF);

        ChipSeqMappingAnalyzer analyser = new ChipSeqMappingAnalyzer(f, null, null);
        JSAP jsap = new JSAP();
        for (Parameter parameter : analyser.getParameter()) {
            try {
                jsap.registerParameter(parameter);
            } catch (JSAPException e) {
                e.printStackTrace();
                fail();
            }
        }

        File parF = null;
        BufferedWriter writer = null;
        try {
            parF = FileHelper.createTempFile(ChipSeqMappingAnalyzerTest.class.getSimpleName(), ".par");
            writer = new BufferedWriter(new FileWriter(parF));
            writer.write(ChipSeqSettings.FILE_INPUT.getName() + " " + tmpF.getAbsolutePath() + "\n");

            writer.write(ChipSeqSettings.READ_DESCRIPTOR.getName() + " " + UniversalReadDescriptor.DESCRIPTORID_PAIRED + "\n");
        } catch (Exception e) {
            e.printStackTrace();
        } finally {
            if (writer != null)
                try {
                    writer.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
        }

        JSAPResult result = jsap.parse("-t mapdist -p "+parF.getAbsolutePath());
        try {
            analyser.validateParameter(result);
            fail();
        } catch (Exception e) {

        }
    }

}
