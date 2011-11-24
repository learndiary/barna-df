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

/**
 * Distribution of a Character to its quality values
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
class CharacterQualityDistribution extends Distribution {
    /**
     * If true, we trim the end of reads to cu away all qualities scores == 2
     */
    private static final boolean TRIM_ILLUMINA_BAD_ENDS = false;
    /**
     * The character
     */
    private char character;

    /**
     * Create a new character distribution
     *
     * @param character the character
     * @param size      the size of the distribution (i.e. number of quality values)
     */
    public CharacterQualityDistribution(char character, int size) {
        super(size);
        this.character = Character.toUpperCase(character);
    }

    /**
     * Add a read
     *
     * @param read the read
     */
    public void addRead(Read read) {
        // trim ends

        int finalPositon = read.getLength() - 1;
        int[] q = read.getQualities();
        if (TRIM_ILLUMINA_BAD_ENDS) {
            for (; finalPositon > 0; finalPositon--) {
                if (q[finalPositon] != 2) break;
            }
        }

        for (int i = 0; i <= finalPositon; i++) {
            char readChar = Character.toUpperCase(read.getSequence().charAt(i));
            if (readChar == character) {
                reads++;
                values[q[i]]++;
            }
        }
    }
}
