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
            Read read1 = reader.parseNext(false);
            assertNotNull(read1);

            assertEquals("HWUSI-EAS627_1:2:1:4:1299/1", read1.getName());
            assertEquals("ATTTTNNTCAAAAACTTTGTCTTTTTTTCTTTCCTCCCCTAAATTTTCCCCAATTTAAATTTTTCCCCCAGGGGTC", read1.getSequence());
            assertEquals(76, read1.getLength());
        } catch (IOException e) {
            e.printStackTrace();
            fail();
        }finally {
            reader.close();
        }
    }
    @Test
    public void testReadMappings(){
        MapFileReader reader = new MapFileReader(testfile, Qualities.Technology.Phred);
        try {
            Read read1 = reader.parseNext(false);
            assertNotNull(read1);

            read1 = reader.parseNext(false);
            // second read
            // should have one perfect mapping without missmatches
            assertNotNull(read1.getMappings());
            assertEquals(1, read1.getMappings().size());
            assertEquals("chr8:R11218885", read1.getMappings().get(0).getName());

            assertNull(read1.getMappings().get(0).getMissmatches());

            read1 = reader.parseNext(false);
            // third read
            // should have one mapping with 8 missmatches
            assertNotNull(read1.getMappings());
            assertEquals(3, read1.getMappings().size());
            assertEquals(8, read1.getMappings().get(0).getMissmatches().size());
            assertEquals("chr12:F10906850", read1.getMappings().get(0).getName());


        } catch (IOException e) {
            e.printStackTrace();
            fail();
        }finally {
            reader.close();
        }
    }
}
