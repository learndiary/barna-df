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

/**
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class QualityTransitions {

    /**
     * Transition matrix
     */
    private long[][][] transitions;
    /**
     * Number of reads
     */
    private int[][] reads;

    /**
     * The distribution for position 0
     */
    private int[] initialDistribution;

    /**
     * Number of quality states
     */
    private int size;

    /**
     * The read length
     */
    private int readLength;

    /**
     * The number of reads added
     */
    private int numReads;

    QualityTransitions(int size, int readLength) {
        this.size = size;
        this.readLength = readLength;
        transitions = new long[readLength][size][size];
        reads = new int[readLength][size];
        initialDistribution = new int[size];
    }

    void addRead(Read r) {
        numReads++;
        int[] readQualities = r.getQualities();

        // position 0 distribution
        initialDistribution[readQualities[0]]++;

        // positions 1..n transitions
        for (int i = 1; i < r.getLength(); i++) {
            int q0 = readQualities[i - 1];
            int q1 = readQualities[i];

            transitions[i][q0][q1]++;
            reads[i][q0]++;
        }
    }

    /**
     * Supported read length
     *
     * @return read length the supported read length
     */
    public int getLength(){
        return transitions.length;
    }

    public long[][][] getTransitions() {
        return transitions;
    }

    public int getQuality(int position, int lastQualityValue, double random) {
        // if start, sample from initial position 0 distribution
        if (position == 0) return getInitialQuality(random);

        // otherwise sample from transitions
        double sum = 0;
        for (int i = 0; i < size; i++) {

            double numberOfReadsFromQ0 = (double) reads[position][lastQualityValue];
            if (numberOfReadsFromQ0 == 0) {
                numberOfReadsFromQ0 = 1d / size; // equally distributed ?
            }
            sum += transitions[position][lastQualityValue][i] / numberOfReadsFromQ0;
            if (sum >= random) return i;
        }
        // catch case where there is no such transition
        return getInitialQuality(random);
    }

    private int getInitialQuality(double r) {
        double sum = 0;
        for (int i = 0; i < size; i++) {
            sum += initialDistribution[i] / (double) numReads;
            if (sum >= r) return i;
        }
        return 38; // return a rather good quality by default
    }
}
