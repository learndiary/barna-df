package barna.model.sam;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMSequenceRecord;
import org.junit.Before;
import org.junit.Test;

import static org.junit.Assert.assertEquals;

/**
 * @author Emilio Palumbo (emiliopalumbo@gmail.com)
 */
public class SAMMappingTest {
    private static SAMRecord r = null;
    SAMMapping mapping = null;

    @Before
    public void setUp() throws Exception {
        SAMFileHeader h = new SAMFileHeader();
        SAMSequenceRecord rec = new SAMSequenceRecord("ref", 100);
        h.addSequence(rec);

        r = new SAMRecord(null);
        r.setHeader(h);
        r.setReadName("r005");
        r.setReferenceName("ref");
        r.setAlignmentStart(16);
        r.setReadString("ATAGCTTCAGT");
        r.setCigarString("6M14N5M");
        r.setMappingQuality(30);
        r.setAttribute("NH", 10);

        mapping = new SAMMapping(r);
    }

    @Test
    public void testName() throws Exception {
        assertEquals("r005",mapping.getName(false));
        assertEquals("r005/1",mapping.getName(true));
    }

    @Test
    public void testChromosome() throws Exception {
        assertEquals("ref",mapping.getChromosome());
    }

    @Test
    public void testPositions() throws Exception {
        assertEquals(15,mapping.getStart());
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
        assertEquals(20,mapping.getNextBlockStart());
        assertEquals(5,mapping.getNextBlockSize());
        assertEquals(0,mapping.getNextBlockStart());
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

    @Test
    public void testSequence() throws Exception {
        assertEquals("ATAGCTTCAGT", mapping.getSequence());
    }

    @Test
    public void testCigar() throws Exception {
        assertEquals("6M14N5M", mapping.getCigar());
    }

    @Test
    public void testCount() throws Exception {
        assertEquals(1.0/10, mapping.getCount(true),0.0);
        assertEquals(1.0, mapping.getCount(false),0.0);
    }

    @Test
    public void testUniqueXT() throws Exception {
        r.setAttribute("XT", 'U');
        assertEquals(true, new SAMMapping(r).isUnique());
        assertEquals(false, mapping.isUnique());
    }

    @Test
    public void testUniqueNH() throws Exception {
        r.setAttribute("NH", 1);
        assertEquals(true, new SAMMapping(r).isUnique());
        assertEquals(false, mapping.isUnique());
    }
}
