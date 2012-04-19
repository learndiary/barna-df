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
