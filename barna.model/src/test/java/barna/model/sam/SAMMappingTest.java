package barna.model.sam;

import net.sf.samtools.SAMRecord;
import org.junit.Before;
import org.junit.Test;

import static org.junit.Assert.*;

/**
 * @author Emilio Palumbo (emiliopalumbo@gmail.com)
 */
public class SAMMappingTest {
    private static SAMRecord r = null;
    SAMMapping mapping = null;

    @Before
    public void setUp() throws Exception {
        r = new SAMRecord(null);
        r.setReadName("r005");
        r.setReferenceName("ref");
        r.setAlignmentStart(16);
        r.setReadString("ATAGCTTCAGT");
        r.setCigarString("6M14N5M");
        r.setMappingQuality(30);

        mapping = new SAMMapping(r);
    }

    @Test
    public void testName() throws Exception {
        assertEquals("r005",mapping.getName());
    }

    @Test
    public void testChromosome() throws Exception {
        assertEquals("ref",mapping.getChromosome());
    }

    @Test
    public void testPositions() throws Exception {
        assertEquals(16,mapping.getStart());
        assertEquals(40,mapping.getEnd());
    }

    @Test
    public void testLength() throws Exception {
        assertEquals(11,mapping.getLength());
    }

    @Test
    public void testBlocks() throws Exception {
        assertEquals(2, mapping.getBlockCount());
        assertEquals(6,mapping.getNextBlockSize());
        assertEquals(36,mapping.getNextBlockStart());
        assertEquals(5,mapping.getNextBlockSize());
        assertEquals(-1,mapping.getNextBlockStart());
    }

    @Test
    public void testQuality() throws Exception {
        assertEquals(30, mapping.getScore());
    }

    @Test
    public void testStrand() throws Exception {
        assertEquals(1, mapping.getStrand());

        r.setReadNegativeStrandFlag(true);
        mapping = new SAMMapping(r);
        assertEquals(-1, mapping.getStrand());
    }
}
