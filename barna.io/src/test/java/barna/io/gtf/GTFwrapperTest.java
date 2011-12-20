package barna.io.gtf;

import barna.genome.io.gtf.GTFwrapper;
import org.junit.Test;

import java.io.File;

import static junit.framework.Assert.*;

public class GTFwrapperTest {

    @Test
    public void testThatTheWrapperWorksWithGzippedFiles() throws Exception {
        File gzippedGtf = new File(getClass().getResource("/gzipped-gtf.gtf.gz").getFile());
        assertTrue(gzippedGtf.exists());

        GTFwrapper wrapperFile = new GTFwrapper(gzippedGtf);
        assertFalse(wrapperFile.isApplicable());

        File sorted = wrapperFile.sort();
        assertNotNull(sorted);
        assertTrue(new GTFwrapper(sorted).isApplicable());
        sorted.delete();
    }
}
