package barna.flux.capacitor.diffexp;

import org.junit.Test;

import static junit.framework.Assert.assertEquals;
import static junit.framework.Assert.assertNotNull;
import static junit.framework.Assert.assertTrue;

public class TranscriptComparatorTest {
    @Test
    public void testBothSpecified() throws Exception {
        Transcript source = new Transcript("chr1", "TEST", "transcript", 1, 100, (short) 0, '+', '.', "reads 756.000000; RPKM 22.246664; transcript_id \"a\"");
        Transcript target = new Transcript("chr1", "TEST", "transcript", 1, 100, (short) 0, '+', '.', "reads 1307.999878; RPKM 32.432495; transcript_id \"a\"");
        double sourceReads = 9015174.868426;
        double targetReads = 11320367.100461;

        TranscriptComparator cmp = new TranscriptComparator(source, target, sourceReads, targetReads);
        DifferentialExpression df = cmp.call();
        assertNotNull(df);
        assertEquals(1.4241909876202e-12, df.getP(), 0.000000001);
        assertEquals(10.185831, df.getDifference(), 0.000000001);
        assertEquals(0.54385100031018, df.getFoldChange(), 0.000000001);

    }
    @Test
    public void testOneNotExpressed() throws Exception {
        Transcript source = new Transcript("chr1", "TEST", "transcript", 1, 100, (short) 0, '+', '.', "reads 756.000000; RPKM 22.246664; transcript_id \"a\"");
        Transcript target = new Transcript("chr1", "TEST", "transcript", 1, 100, (short) 0, '+', '.', "reads 0; RPKM 0; transcript_id \"a\"");
        double sourceReads = 9015174.868426;
        double targetReads = 11320367.100461;

        TranscriptComparator cmp = new TranscriptComparator(source, target, sourceReads, targetReads);
        DifferentialExpression df = cmp.call();
        assertNotNull(df);
        assertEquals(8.155749649539851E-268, df.getP(), 0.0000000001);
        assertEquals(-22.246664, df.getDifference(), 0.00000001);
        assertEquals(Double.NEGATIVE_INFINITY, df.getFoldChange());
    }
}
