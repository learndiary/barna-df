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

/**
 * Distribution of all quality values, independent of the position
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
class QualityDistribution extends Distribution {
    /**
     * If true, we trim the end of reads to cut away all qualities scores == 2
     */
    private static final boolean TRIM_ILLUMINA_BAD_ENDS = false;
    /**
     * The number of nucleotides
     */
    private long numberOfNucleotides;
    /**
     * The sum of qualities
     */
    private long sumQualities;

    /**
     * Create a quality distribution
     *
     * @param size number of quality values
     */
    public QualityDistribution(int size) {
        super(size);
    }

    public void addRead(Read read) {
        // trim ends

        int finalPositon = read.getLength() - 1;
        int[] q = read.getQualities();
        if (TRIM_ILLUMINA_BAD_ENDS) {
            for (; finalPositon > 0; finalPositon--) {
                if (q[finalPositon] != 2) break;
            }
        }

        reads += finalPositon + 1;

        //for (int i = 0; i < q.length; i++) {
        for (int i = 0; i < finalPositon; i++) {
            numberOfNucleotides++;
            sumQualities += q[i];
            values[q[i]]++;
        }
    }

    /**
     * Returns the average quality
     *
     * @return average quality
     */
    public double getAverageQuality() {
        return (double) sumQualities / (double) numberOfNucleotides;
    }
}
