package barna.flux.capacitor.diffexp;

import barna.flux.capacitor.diffexp.math.FishersExactTest;
import org.junit.Test;

import static junit.framework.Assert.assertEquals;

public class FishersExactTestTest {

    @Test
    public void testFisherRComparison() throws Exception {
        double[] doubles = FishersExactTest.fishersExactTest(756, 9014419, 1308, 11319059);
        assertEquals(1.4241909876202e-12, doubles[0], 0.00000000001);
    }
}
