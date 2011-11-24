/*
 * This file is part of the Flux Library.
 *
 * The code of the Flux Library may be freely distributed and modified under the terms of the
 * European Union Public Licence (EUPL) published on the web site <http://www.osor.eu/eupl/european-union-public-licence-eupl-v.1.1>.
 * Copyright for the code is held jointly by the individual authors, who should be listed
 * in @author doc comments. According to Article 5 and Article 11 of the EUPL, publications that
 * include results produced by the Flux Library are liable to reference the Work,
 * see the Flux Library homepage <http://flux.sammeth.net> for more information.
 */

package fbi.genome.sequencing.rnaseq.simulation.tools;

import fbi.commons.Log;
import fbi.commons.flux.FluxTool;
import fbi.commons.options.HelpPrinter;
import fbi.commons.options.Options;
import fbi.genome.io.SpliceGraphIO;
import fbi.genome.model.IntronModel;
import org.cyclopsgroup.jcli.ArgumentProcessor;
import org.cyclopsgroup.jcli.annotation.Cli;
import org.cyclopsgroup.jcli.annotation.Option;

import java.io.File;

/**
 * Extract splice junctions from genome
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
//@Cli(name = "extract", description = "Splice junction extraction", restrict = false)
public class SpliceJunctionExtractor { //implements FluxTool<Void> {
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
            fbi.genome.model.Graph.overrideSequenceDirPath = gffFile.getAbsolutePath();
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
     * @param length flank length  @code{> 0}
     */

    @Option(name = "5", longName = "5flank", description = "exonic flank 5' of intron")
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
     * @param length flank length @code{> 0}
     */
    @Option(name = "3", longName = "3flank", description = "exonic flank 3' of intron")
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
    @Option(name = "i", longName = "input", description = "GFF input file")
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
    @Option(name = "m", longName = "model", description = "specify the intron model")
    public void setModelFile(File modelFile) {
        this.modelFile = modelFile;
    }

    public boolean validateParameters(HelpPrinter printer, ArgumentProcessor toolArguments) {
        return true;
    }
}
