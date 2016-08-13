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

import barna.commons.log.Log;
import barna.commons.utils.TableFormatter;

import java.util.Arrays;
import java.util.List;

/**
 * Cross talk table contains the probability for a genomic character that mutated in a read.
 * Say we observe an 'A' in the genomic sequence and want to add an error to the simulated read.
 * We can then use the transition table to check to which character we should mutate the 'A' to.
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class CrossTalkModel {
    /**
     * Supported characters
     */
    static final char[] SYMBOLS = new char[]{'A', 'C', 'G', 'N', 'T'};

    /**
     * Global flag whether (further) unknown symbols should be reported.
     */
    static boolean reportUnknownSymbol= true;

    /**
     * The transitions per quality in the form [quality][from][to]
     */
    private long[][][] transitions;
    /**
     * Observed counts
     */
    private long[][] counts;
    /**
     * Use qualities or positions
     */
    private boolean qualities;

    /**
     * The number of scanned nucleotides
     */
    private long numberOfNucleotides;

    /**
     * The number of mutations
     */
    private long numberOfMutations;


    public CrossTalkModel(int states, boolean qualities) {
        super();
        this.qualities = qualities;
        this.transitions = new long[states][SYMBOLS.length][SYMBOLS.length];
        this.counts = new long[states][SYMBOLS.length];
    }

    public void addRead(Read read) {
        List<Read.Mapping> mappings = read.getMappings();
        if (mappings == null || mappings.size() == 0) {
            return; // no mapping - skip this read
        }

        numberOfNucleotides += read.getLength();
        // add missmatches
        for (Read.Mapping mapping : mappings) {
            if (mapping.getMissmatches() != null) {
                for (Read.Missmatch missmatch : mapping.getMissmatches()) {
                    int i = missmatch.getPosition() - 1;
                    char readChar = Character.toUpperCase(read.getSequence().charAt(i));
                    char genomicCharacter = Character.toUpperCase(missmatch.getGenomicCharacter());
                    int p0 = Arrays.binarySearch(SYMBOLS, genomicCharacter);
                    int p1 = Arrays.binarySearch(SYMBOLS, readChar);
                    numberOfMutations++;

                    int state = transitions.length == 1 ? 0 : qualities ? read.getQualities()[i] : i;
                    transitions[state][p0][p1]++;
                    counts[state][p0]++;
                }
            }
        }
    }

    /**
     * Retrieves a falsified base for the original character.
     * @param state a state
     * @param from original character
     * @param random a randomly generated number
     * @return a falsified base
     */
    public char getTransition(int state, char from, double random) {

        // make sure its upper case
        from = Character.toUpperCase(from);
        int p0 = Arrays.binarySearch(SYMBOLS, from);
        if (p0< 0|| p0>= SYMBOLS.length) {
            if (reportUnknownSymbol) {
                Log.warn("Unknown symbol "+ from+ " has been mapped to N!\n" +
                    "Forthcoming unknown symbols will not anymore be reported.");
            }
            reportUnknownSymbol= false;
            p0= 3;  // map unknown symbols to 'N'
        }

        state = transitions.length == 1 ? 0 : state;

        if (p0 >= 0) {
            double sum = 0;
            for (int i = 0; i < SYMBOLS.length; i++) {
                double numberOfReadsFromQ0 = counts[state][p0];
                if (numberOfReadsFromQ0 == 0) {
                    numberOfReadsFromQ0 = 1d / transitions.length; // equally distributed ?
                }

                sum += transitions[state][p0][i] / numberOfReadsFromQ0;
                if (sum >= random) {
                    return SYMBOLS[i];
                }
            }
        }
        return from;
    }

    private double p(int state, int from, int to) {
        double numberOfReadsFromQ0 = counts[state][from];
        if (numberOfReadsFromQ0 == 0) return 0;
        return transitions[state][from][to] / numberOfReadsFromQ0;
    }

    public String toString() {
        StringBuffer b = new StringBuffer();

        for (int s = 0; s < transitions.length; s++) {
            b.append("State " + s).append(barna.commons.system.OSChecker.NEW_LINE);
            b.append(toTable(s));
        }
        b.append(barna.commons.system.OSChecker.NEW_LINE);
        return b.toString();
    }

    private String toTable(int state) {
        TableFormatter table = new TableFormatter(SYMBOLS.length + 1);
        String[] header = new String[SYMBOLS.length + 1];
        header[0] = "";
        for (int i = 0; i < SYMBOLS.length; i++) {
            header[i + 1] = Character.toString(SYMBOLS[i]);
        }
        table.addRow(header);

        for (int i = 0; i < SYMBOLS.length; i++) {
            String[] row = new String[SYMBOLS.length + 1];
            row[0] = Character.toString(SYMBOLS[i]);
            for (int j = 0; j < SYMBOLS.length; j++) {
                row[j + 1] = "" + p(state, i, j);
            }
            table.addRow(row);
        }
        return table.toString();
    }

    public String printStateTable() {
        StringBuffer b = new StringBuffer();
        for (int i = 0; i < transitions.length; i++) {
            b.append(i).append("\t");

            for (int k = 0; k < SYMBOLS.length; k++) {
                b.append(p(i, k, k)).append("\t");
            }

            b.append(barna.commons.system.OSChecker.NEW_LINE);
        }
        return b.toString();
    }

    public double[][] getDistribution(char from, char to) {
        double[][] d = new double[2][transitions.length];
        int p0 = Arrays.binarySearch(SYMBOLS, from);
        int p1 = Arrays.binarySearch(SYMBOLS, to);
        for (int k = 0; k < transitions.length; k++) {
            d[0][k] = k;
            d[1][k] = p(k, p0, p1);
        }
        return d;
    }

    /**
     * Get the average number of mutations
     *
     * @return mutations average number of mutations
     */
    public double getAverageMutations() {
        return (double) numberOfMutations / (double) numberOfNucleotides;
    }
}
