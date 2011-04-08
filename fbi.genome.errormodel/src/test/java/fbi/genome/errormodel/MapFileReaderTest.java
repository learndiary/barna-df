package fbi.genome.errormodel;

import fbi.commons.tools.Qualities;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.File;
import java.io.IOException;

import static junit.framework.Assert.*;

/**
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class MapFileReaderTest {

    private static File testfile;

    @BeforeClass
    public static void setUp(){
        testfile = new File(MapFileReaderTest.class.getResource("/test.map").getFile());
    }
    @Test
    public void testReadOne(){
        MapFileReader reader = new MapFileReader(testfile, Qualities.Technology.Phred);
        try {
            Read read1 = reader.parseNext();
            assertNotNull(read1);

            assertEquals("HWUSI-EAS627_1:2:1:4:1299/1", read1.getName());
            assertEquals("ATTTTNNTCAAAAACTTTGTCTTTTTTTCTTTCCTCCCCTAAATTTTCCCCAATTTAAATTTTTCCCCCAGGGGTC", read1.getSequence());
            assertEquals(76, read1.getLength());
            // only one read
            assertNull(reader.parseNext());
        } catch (IOException e) {
            e.printStackTrace();
            fail();
        }
    }
}
