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

package barna.model;

/**
 * Phred/Solexa/Illumina quality translator
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class Qualities {
    /**
     * General minimum ASCII value
     */
    private static final int MIN = 33;

    /**
     * General maximum ASCII value
     */
    private static final int MAX = 126;

    public static final int[] PHRED_RANGE = {0, 93};
    private static final int[] PHRED_RANGE_ASCII = {33, 126};

    public static final int[] ILLUMINA_18_RANGE = PHRED_RANGE;
    private static final int[] ILLUMINA_18_RANGE_ASCII = PHRED_RANGE_ASCII;

    public static final int[] SOLEXA_RANGE = {-5, 62};
    private static final int[] SOLEXA_RANGE_ASCII = {59, 126};

    public static final int[] ILLUMINA_13_RANGE = {0, 62};
    private static final int[] ILLUMINA_13_RANGE_ASCII = {64, 126};


    public static enum Technology {Phred, Solexa, Illumina18, Illumina13}


    /**
     * Returns the Phred quality for the given character. Solexa qualities are converted
     *
     * @param tech  the technology
     * @param value the value
     * @return quality the phred quality
     */
    public static int quality(Technology tech, char value) {
        switch (tech) {
            case Phred:
                return (value - PHRED_RANGE_ASCII[0]);
            case Illumina13:
                return (value - ILLUMINA_13_RANGE_ASCII[0]);
            case Illumina18:
                return (value - ILLUMINA_13_RANGE_ASCII[0]);
            case Solexa:
                return solexa2phredQuality((value - SOLEXA_RANGE_ASCII[0]));
        }
        return -1;
    }

    public static byte ascii(int value) {
        return ascii(Technology.Illumina18, value);
    }

    /**
     * Returns the Phred quality for the given character. Solexa qualities are converted
     *
     * @param tech  the technology
     * @param value the value
     * @return quality the phred quality
     */
    public static byte ascii(Technology tech, int value) {
        switch (tech) {
            case Phred:
                return (byte) (value + PHRED_RANGE_ASCII[0]);
            case Illumina13:
                return (byte) (value + ILLUMINA_13_RANGE_ASCII[0]);
            case Illumina18:
                return (byte) (value + ILLUMINA_18_RANGE_ASCII[0]);
            case Solexa:
                return (byte) (value + SOLEXA_RANGE_ASCII[0]); // todo check solexa quality
        }
        return '-';
    }

    /**
     * Get the probability for a given phred quality value
     *
     * @param q the quality value
     * @return p the probability
     */
    public static double getPropability(int q) {
        return Math.pow(10.0, ((double) -q / 10.0));
    }

    public static int solexa2phredQuality(int qualSol) {
        // 10*log(10**(solexa_quality/10.0) + 1, 10)

        int q = (int) (10 *
                Math.log(1 + Math.pow(10, (qualSol / 10d)))
                / Math.log(10));
        assert (q >= 0);
        return q;
    }

    public static int phred2solexaQuality(int qualPhred) {
        //10*log(10**(phred_quality/10.0) - 1, 10)
        int q = (int) (10 *
                Math.log(Math.pow(10, (qualPhred / 10d)) - 1)
                / Math.log(10));
        return q;
    }


}
