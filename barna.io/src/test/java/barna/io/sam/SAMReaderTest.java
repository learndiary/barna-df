package barna.io.sam;

import java.io.File;

import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;

import barna.commons.Execute;

public class SAMReaderTest {
	
	private static File testfile;

	@BeforeClass
	public static void setUp() {
		//testfile = new File("/home/emilio/test.sam");
        Execute.initialize(4);
	}

	@AfterClass
	public static void tearDown() throws Exception {
		Execute.shutdown();
	}

	@Test
	public void testRead() {
		//SAMReader wrapper = new SAMReader(testfile);
		//wrapper.read();
		//assertNotNull(wrapper.beds);
		//assertTrue(wrapper.beds.length > 0);
		//assertEquals("chrM",wrapper.beds[0].getChrom().toString());
	}

}
