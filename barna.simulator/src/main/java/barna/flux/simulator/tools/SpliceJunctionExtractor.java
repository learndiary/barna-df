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

package barna.flux.simulator.tools;

import barna.commons.launcher.Options;
import barna.commons.log.Log;
import barna.io.SpliceGraphIO;
import barna.model.IntronModel;

import java.io.File;

/**
 * Extract splice junctions from genome
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
//@Cli(name = "extract", description = "Splice junction extraction", restrict = false)
public class SpliceJunctionExtractor { //implements Tool<Void> {
    /**
     * Default 5' flank
     */
    private static final int FLANK5_LEN = 50;
    /**
     * Default 3' flank
     */
    private static final int FLANK3_LEN = 50;
    /**
     * The Genome File
     */
    private File gffFile;
    /**
     * The Flanks
     */
    private int[] eFlanks;
    /**
     * The model file
     */
    private File modelFile;
    /**
     * The options
     */
    private Options options;

    public Void call() throws Exception {
        if (gffFile == null || !gffFile.exists()) {
            throw new RuntimeException("No valid GFF input file specified!");
        } else {
            // todo : refactor this to not use a static variable
            barna.model.Graph.overrideSequenceDirPath = gffFile.getAbsolutePath();
        }

        int[] flanks = getEFlanks();
        flanks[0] = flanks[0] < 0 ? FLANK5_LEN : flanks[0];
        flanks[1] = flanks[1] < 0 ? FLANK3_LEN : flanks[1];

        IntronModel iModel = new IntronModel();
        if (modelFile != null) {
            Log.info("Reading Model file " + modelFile.getAbsolutePath());
            iModel.read(modelFile);
        }
        Log.info("Extracting splice junctions, 5'sequence " + flanks[0]
                + ", 3'sequence " + eFlanks[1] + ", intron model " + iModel.getName());
        SpliceGraphIO.extractSpliceJunctions(flanks[0], flanks[1], iModel, gffFile, null);
        return null;
    }


    /**
     * Get the current flanks. Returns a two element array, index 0 is 5' flank,  index 1 is the 3' flank
     *
     * @return flanks the flanks
     */
    public int[] getEFlanks() {
        if (eFlanks == null) {
            eFlanks = new int[2];
            eFlanks[0] = -1;
            eFlanks[1] = -1;
        }

        return eFlanks;
    }

    /**
     * Set the 5' flank
     *
     * @param length flank length  @code{&gt; 0}
     */

//    @Option(name = "5", longName = "5flank", description = "exonic flank 5' of intron")
    public void set5flank(int length) {
        getEFlanks()[0] = -1;
        if (length <= 0) {
            throw new IllegalArgumentException("Not a valid length for 5' exon flank: " + length);
        }
        getEFlanks()[0] = length;
    }

    /**
     * Set the 3' flank
     *
     * @param length flank length @code{&gt; 0}
     */
//    @Option(name = "3", longName = "3flank", description = "exonic flank 3' of intron")
    public void set3flank(int length) {
        getEFlanks()[1] = -1;
        if (length <= 0) {
            throw new IllegalArgumentException("Not a valid length for 3' exon flank: " + length);
        }
        getEFlanks()[1] = length;
    }

    /**
     * Get the GFF input file
     *
     * @return gffFile the gff input file
     */
    public File getGffFile() {
        return gffFile;
    }

    /**
     * Set the GFF input file
     *
     * @param gffFile the input file
     */
//    @Option(name = "i", longName = "input", description = "GFF input file")
    public void setGffFile(File gffFile) {
        this.gffFile = gffFile;
    }

    /**
     * Get the model file
     *
     * @return model the model file
     */
    public File getModelFile() {
        return modelFile;
    }

    /**
     * Set the model file
     *
     * @param modelFile model file
     */
//    @Option(name = "m", longName = "model", description = "specify the intron model")
    public void setModelFile(File modelFile) {
        this.modelFile = modelFile;
    }

//    public boolean validateParameters(HelpPrinter printer, ArgumentProcessor toolArguments) {
//        return true;
//    }
}
