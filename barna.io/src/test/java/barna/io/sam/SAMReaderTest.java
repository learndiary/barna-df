package barna.io.sam;

import barna.commons.Execute;
import barna.io.MSIterator;
import barna.model.sam.SAMMapping;
import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.File;

import static org.junit.Assert.assertEquals;

public class SAMReaderTest {
	
	private File testfile = new File(getClass().getResource("/test.bam").getFile());

	@BeforeClass
	public static void setUp() {
        Execute.initialize(4);
	}

	@AfterClass
	public static void tearDown() throws Exception {
		Execute.shutdown();
	}

	@Test
	public void testRead() {
        SAMReader reader = new SAMReader(testfile, true, true);
        MSIterator iter = reader.read("chr22", 24030323, 24041363);
        SAMMapping mapping = (SAMMapping)iter.next();

        int c = 1;
        while (iter.hasNext()) {
            ++c;
            mapping = (SAMMapping)iter.next();
        }


        assertEquals(181, c);
        assertEquals("chr22", mapping.getChromosome());
        assertEquals(24034554, mapping.getStart());
        assertEquals(24034629, mapping.getEnd());
	}

}
