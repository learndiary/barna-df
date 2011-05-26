package fbi.genome.io.bed;

import fbi.commons.Execute;
import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.File;

import static junit.framework.Assert.assertEquals;

/**
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class BEDwrapperTest {

    private static File testfile;

    @BeforeClass
    public static void setUp(){
        testfile = new File(BEDwrapperTest.class.getResource("/test.bed").getFile());
        Execute.initialize(4);
    }

    @AfterClass
    public static void tearDown() throws Exception {
        Execute.shutdown();
    }


    @Test
    public void testScanFile(){
        BEDwrapper wrapper = new BEDwrapper(testfile.getAbsolutePath());
        wrapper.scanFile();

        //scanFileReadLines= 0;
        //countAll= 0; countEntire= 0; countSplit= 0; countReads= 0;
        assertEquals(17, wrapper.nrUniqueLinesRead );
        assertEquals(17, wrapper.countAll );
        assertEquals(5, wrapper.countSplit );
        assertEquals(12, wrapper.countEntire );
        assertEquals(17, wrapper.countReads );
    }
}
