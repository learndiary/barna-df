package barna.genome.sequencing.rnaseq.simulation;

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
