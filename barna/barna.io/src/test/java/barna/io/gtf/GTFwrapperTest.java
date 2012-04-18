package barna.io.gtf;

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
    @Test
    public void testThatNormaltabsForFirstFieldsAreApplied() throws Exception {
        File gzippedGtf = new File(getClass().getResource("/valid-tab-space.gtf").getFile());
        assertTrue(gzippedGtf.exists());
        GTFwrapper wrapperFile = new GTFwrapper(gzippedGtf);
        File sorted = wrapperFile.sort();
        assertNotNull(sorted);
        assertTrue(new GTFwrapper(sorted).isApplicable());
        sorted.delete();
    }
    @Test
    public void testThatOnlySpaceSeparatedFields() throws Exception {
        File gzippedGtf = new File(getClass().getResource("/only-space.gtf").getFile());
        assertTrue(gzippedGtf.exists());
        GTFwrapper wrapperFile = new GTFwrapper(gzippedGtf);
        File sorted = wrapperFile.sort();
        assertNotNull(sorted);
        assertTrue(new GTFwrapper(sorted).isApplicable());
        sorted.delete();
    }
    @Test
    public void testMixedFileWithTabsAndSpaces() throws Exception {
        File gzippedGtf = new File(getClass().getResource("/mixed-space.gtf").getFile());
        assertTrue(gzippedGtf.exists());
        GTFwrapper wrapperFile = new GTFwrapper(gzippedGtf);
        File sorted = wrapperFile.sort();
        assertNotNull(sorted);
        assertTrue(new GTFwrapper(sorted).isApplicable());
        sorted.delete();
    }
    @Test
    public void testWithTabsInTheGTFFieldsAtTheEnd() throws Exception {
        File gzippedGtf = new File(getClass().getResource("/tab-fields.gtf").getFile());
        assertTrue(gzippedGtf.exists());
        GTFwrapper wrapperFile = new GTFwrapper(gzippedGtf);
        File sorted = wrapperFile.sort();
        assertNotNull(sorted);
        assertTrue(new GTFwrapper(sorted).isApplicable());
        sorted.delete();
    }
    @Test
    public void testThatInvalidGTFFailsWithReport() throws Exception {
        File gzippedGtf = new File(getClass().getResource("/invalid.gtf").getFile());
        assertTrue(gzippedGtf.exists());
        GTFwrapper wrapperFile = new GTFwrapper(gzippedGtf);
        try{
            File sorted = wrapperFile.sort();
            fail();
        }catch (Exception error){
        }
    }
}
