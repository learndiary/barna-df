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

package barna.genome.sequencing.rnaseq.simulation.error;

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

    public char getTransition(int state, char from, double random) {

        // make sure its upper case
        from = Character.toUpperCase(from);
        int p0 = Arrays.binarySearch(SYMBOLS, from);
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
            b.append("State " + s).append("\n");
            b.append(toTable(s));
        }
        b.append("\n");
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

            b.append("\n");
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
