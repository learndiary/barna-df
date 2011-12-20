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

package fbi.genome.sequencing.rnaseq.simulation.error;

import fbi.commons.ByteArrayCharSequence;
import fbi.commons.tools.Qualities;

import java.util.Random;

/**
 * The model pool applies error model qualities and mutates the sequence according to
 * the crosstalk table
 */
public class ModelPool {
    /**
     * Random generator fot mutations
     */
    private Random rndMutator = new Random();
    /**
     * True if a fastq file should be generated
     */
    private boolean fastaOutput = true;
    /**
     * The error model
     */
    private QualityErrorModel errorModel;
    /**
     * The nucleotides written
     */
    private long writtenNucleotides;
    /**
     * Sum of the qualities
     */
    private long sumQualities;
    /**
     * Sum of mutations
     */
    private long sumMutations;

    /**
     * Create a new model pool
     *
     * @param fastOutput create fastq output
     * @param errorModel the error model
     */
    public ModelPool(final boolean fastOutput, final QualityErrorModel errorModel) {
        if(errorModel == null) throw new NullPointerException("Null error model not permitted");
        this.fastaOutput = fastOutput;
        this.errorModel = errorModel;
    }

    /**
     * Apply qualities and mutations to the sequence and
     * return the averages for [mutations,quality]
     *
     * @param cs the current line
     * @param seqStart the sequence start within the line
     */
    public void apply(ByteArrayCharSequence cs, int seqStart) {
        // generate qualities (FASTQ)
        int seqEnd = cs.end;
        int len = seqEnd - seqStart;

        // prepare sequence for fastQ
        if (fastaOutput && errorModel != null) {
            cs.append("\n+\n");
            seqEnd += 3;
            cs.ensureLength(cs.end, len);    // for qualities
        }


        byte[] a = cs.chars;
        int quality = -1;

        if(errorModel != null){
            for (int i = 0; i < len; i++) {
                // iterate over the sequence and generate qualities
                char character = (char) a[seqStart+i];
                double r = rndMutator.nextDouble();
                quality = errorModel.getQualityModel().getQuality(i, quality, r);
                sumQualities +=quality;
                writtenNucleotides++;
                // check if we have to mutate
                double pe = Qualities.getPropability(quality);
                double random = rndMutator.nextDouble();
                if ( random <= pe) {
                    // mutate using crosstalk
                    a[seqStart + i] = (byte) errorModel.getCrossTalk().getTransition(quality, (char) a[seqStart + i], random);
                    // SIMULATOR-29 make sure we count only for "real" mutations
                    if (a[seqStart + i] != character){
                        sumMutations++;
                    }
                }

                // if fastq, write quality value
                if(fastaOutput){
                    a[cs.end++] = Qualities.ascii(quality);
                }
            }
        }else{
            writtenNucleotides += len;
        }
    }

    /**
     * Returns the overall average mutation rate
     *
     * @return mutation mutation rate
     */
    public double getAverageMutations() {
        return (double)sumMutations/(double)writtenNucleotides;

    }

    /**
     * Return the overall average quality value
     *
     * @return quality overall average quality
     */
    public double getAverageQuality() {
        return (double)sumQualities/(double)writtenNucleotides;
    }

    /**
     * Returns true if an error model exists
     *
     * @return model the error model
     */
    public boolean hasErrorModel() {
        return errorModel != null;
    }
}