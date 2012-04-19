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

package barna.flux.simulator;

import barna.model.Graph;
import org.junit.Test;

import static junit.framework.Assert.assertEquals;

/**
 * Test sequence read extensions in graph genome thingy
 *
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class GraphTest {

    @Test
    public void testSimpleRead() throws Exception {
        Graph.overrideSequenceDirPath = getClass().getResource("/genome_Spikes").getFile();

        String first100 = Graph.readSequence(
                null,
                "Lambdaclone23-2",
                true,
                1,
                100);

        assertEquals(100, first100.length());
        assertEquals(
                "CCACGTAAGCGAAACAAAAACGGGGTTTACCTTACCGAAATCGGTACGGATACCGCGAAAGAGCAGATTTATAACCGCTTCACACTGACGCCGGAAGGGG",
                first100);

        // check negative start
        String first200 = Graph.readSequence(
                null,
                "Lambdaclone23-2",
                true,
                -100,
                100);

        assertEquals(100, first200.length());
        assertEquals(
                "CCACGTAAGCGAAACAAAAACGGGGTTTACCTTACCGAAATCGGTACGGATACCGCGAAAGAGCAGATTTATAACCGCTTCACACTGACGCCGGAAGGGG",
                first200);

        String last100 = Graph.readSequence(
                null,
                "Lambdaclone23-2",
                true,
                9687,
                10000,
                false);

        assertEquals(100, last100.length());
        assertEquals(
                "GAAACAGTCCAGCGTGAAGGTGTCTGCGGGCGATCGTCAGGAAGACAGTGCTCATGCTGCCCTGCTGACGCTTCAGGCAGAACTCCGGACGCTGGAGAAG",
                last100);



    }
}
