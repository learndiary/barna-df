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

package barna.flux.simulator.error;

import barna.commons.ByteArrayCharSequence;
import barna.commons.RandomFactory;
import barna.commons.system.OSChecker;
import barna.model.Qualities;

import java.util.Random;

/**
 * The model pool applies error model qualities and mutates the sequence according to
 * the crosstalk table
 */
public class ModelPool {
    /**
     * Random generator fot mutations
     */
    private Random rndMutator = RandomFactory.get();
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
     * Length of the simulated reads
     */
    private int simulatedReadLength;

    /**
     * Create a new model pool
     *
     * @param fastOutput create fastq output
     * @param errorModel the error model
     * @param simulatedReadLength length of the simulated reads (used for scaling)
     */
    public ModelPool(final boolean fastOutput, final QualityErrorModel errorModel, int simulatedReadLength) {
        this.simulatedReadLength = simulatedReadLength;
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
            cs.append(OSChecker.NEW_LINE+"+"+ OSChecker.NEW_LINE);
            seqEnd += 3;
            cs.ensureLength(cs.end, len);    // for qualities
        }


        byte[] a = cs.chars;
        int quality = -1;

        if(errorModel != null){
            for (int i = 0; i < len; i++) {
                // iterate over the sequence and generate qualities
                int characterPosition = seqStart + i;
                char character = (char) a[characterPosition];
                double r = rndMutator.nextDouble();
                // implement BARNA-167 and scale the error model up/down
                int scaledPosition = scalePosition(i, errorModel.getReadLength(), simulatedReadLength);
                quality = errorModel.getQualityModel().getQuality(scaledPosition, quality, r);
                sumQualities +=quality;
                writtenNucleotides++;
                // check if we have to mutate
                double pe = Qualities.getPropability(quality);
                double random = rndMutator.nextDouble();
                if ( random <= pe) {
                    // mutate using crosstalk
                    a[characterPosition] = (byte) errorModel.getCrossTalk().getTransition(quality, (char) a[characterPosition], random);
                    // SIMULATOR-29 make sure we count only for "real" mutations
                    if (a[characterPosition] != character){
                        sumMutations++;
                        // BARNA-184 make sure mutations are lowe case
                        a[characterPosition] = (byte) Character.toLowerCase(a[characterPosition]);
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

    /**
     * Scale the position to fit within the error model range
     * Also scales down
     *
     * @param position source position
     * @param modelLength supported length of the model
     * @param readLength length of the simulated reads
     * @return scaled scaled position
     */
    protected static int scalePosition(int position, int modelLength, int readLength){
        if(position == 0) return 0;
        if(modelLength == readLength) return position;
        position = (int) (position * (((modelLength-1)/(double)readLength)));
        position = Math.min(modelLength-1, position);
        position = Math.max(0, position);
        return position;
    }
}