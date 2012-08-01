package barna.io.sam;

import java.io.File;

import barna.io.MSIterator;
import barna.model.sam.SAMMapping;
import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;

import barna.commons.Execute;

public class SAMReaderTest {
	
	private static File testfile;

	@BeforeClass
	public static void setUp() {
		testfile = new File("/home/emilio/fromMicha/test.bam");
        Execute.initialize(4);
	}

	@AfterClass
	public static void tearDown() throws Exception {
		Execute.shutdown();
	}

	@Test
	public void testRead() {
        SAMReader reader = new SAMReader(testfile, true);
        MSIterator iter = reader.read("chr22", 24030323, 24041363);
        SAMMapping mapping = (SAMMapping)iter.next();

        int c = 1;
        while (iter.hasNext()) {
            ++c;
            mapping = (SAMMapping)iter.next();
        }


        assertEquals(181, c);
        assertEquals("chr22", mapping.getChromosome());
        assertEquals(24236618, mapping.getStart());
        assertEquals(24236692, mapping.getEnd());
	}

}
