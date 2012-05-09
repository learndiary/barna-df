package barna.io.sam;

import static org.junit.Assert.*;

import java.io.File;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import barna.commons.Execute;
import barna.io.bed.BEDwrapperTest;
import barna.io.sam.SAMWrapper;

public class SAMWrapperTest {
	
	private static File testfile;

	@BeforeClass
	public static void setUp() {
		testfile = new File("/home/emilio/test.sam");
        Execute.initialize(4);
	}

	@AfterClass
	public static void tearDown() throws Exception {
		Execute.shutdown();
	}

	@Test
	public void testRead() {
		SAMWrapper wrapper = new SAMWrapper(testfile);
		wrapper.read();
		assertNotNull(wrapper.beds);
		assertTrue(wrapper.beds.length > 0);		
		//assertEquals("chrM",wrapper.beds[0].getChrom().toString());
	}

}
