package barna.flux.capacitor.diffexp;

import org.junit.Test;

import static junit.framework.Assert.assertEquals;
import static junit.framework.Assert.assertNotNull;

public class GFFEntryTest {
    @Test
    public void testBasicParsing() throws Exception {
        GFFEntry e = GFFEntry.parse("chr1\tHAVANA\ttranscript\t12010\t13670\t.\t+\t.\ttranscript_id \"ENST00000450305.2\"; locus_id \"chr1:11869-14412W\"; gene_id \"ENSG00000223972.4\"; reads 86.362259; length 632; RPKM 5.362027\n");
        assertNotNull(e);
        assertEquals("chr1", e.getChromosome());
        assertEquals("HAVANA", e.getSource());
        assertEquals("transcript", e.getFeature());
        assertEquals(12010, e.getStart());
        assertEquals(13670, e.getEnd());
        assertEquals(-1, e.getScore());
        assertEquals('+', e.getStrand());
        assertEquals('.', e.getCodingFrame());
        assertEquals("transcript_id \"ENST00000450305.2\"; locus_id \"chr1:11869-14412W\"; gene_id \"ENSG00000223972.4\"; reads 86.362259; length 632; RPKM 5.362027", e.getGroup());
    }
    @Test
    public void testTranscriptParsing() throws Exception {
        GFFEntry gffEntry = GFFEntry.parse("chr1\tHAVANA\ttranscript\t12010\t13670\t.\t+\t.\ttranscript_id \"ENST00000450305.2\"; locus_id \"chr1:11869-14412W\"; gene_id \"ENSG00000223972.4\"; reads 86.362259; length 632; RPKM 5.362027\n");
        QuantificationEntry e = new QuantificationEntry("transcript_id", gffEntry);
        assertEquals(5.362027, e.getRpkm(), 0.000001);
        assertEquals(86.362259, e.getReadCount(), 0.000001);
    }
}
