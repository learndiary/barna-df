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

package fbi.commons.tools;

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