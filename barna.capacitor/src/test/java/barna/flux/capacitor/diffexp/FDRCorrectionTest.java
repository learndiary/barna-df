package barna.flux.capacitor.diffexp;

import org.junit.Test;

import java.util.ArrayList;

import static junit.framework.Assert.assertEquals;

public class FDRCorrectionTest {
    @Test
    public void testFDR() throws Exception {
        ArrayList<DifferentialExpression> list = new ArrayList<DifferentialExpression>();

        list.add(new DifferentialExpression(null, null, 0.004, 0,0));
        list.add(new DifferentialExpression(null, null, 0.05, 0,0));
        list.add(new DifferentialExpression(null, null, 0.06, 0,0));

        Corrections.fdr(list);
        assertEquals(0.012, list.get(0).getFdrP(), 0.00000001);
        assertEquals(0.06, list.get(1).getFdrP(), 0.00000001);
        assertEquals(0.06, list.get(2).getFdrP(), 0.00000001);
    }
}
