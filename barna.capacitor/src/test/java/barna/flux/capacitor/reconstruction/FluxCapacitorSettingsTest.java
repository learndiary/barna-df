package barna.flux.capacitor.reconstruction;

import org.junit.Test;

import java.io.ByteArrayInputStream;

import static junit.framework.Assert.assertEquals;
import static junit.framework.Assert.fail;

public class FluxCapacitorSettingsTest {
    @Test
    public void testThatSpacesAreAcceptedInReadDescriptor() throws Exception {
        FluxCapacitorSettings settings = new FluxCapacitorSettings();
        ByteArrayInputStream input = new ByteArrayInputStream(("READ_DESCRIPTOR\t{ID} {MATE}" + barna.commons.system.OSChecker.NEW_LINE +
                "MAPPING_FILE "+getClass().getResource("/mm9_chr1_chrX_sorted.bed").getFile()+barna.commons.system.OSChecker.NEW_LINE +
                "ANNOTATION_FILE "+getClass().getResource("/mm9_chr1_chrX_sorted.gtf").getFile()+barna.commons.system.OSChecker.NEW_LINE).getBytes());

        try{
            settings.parse(input);
            settings.validate();
            assertEquals("{ID} {MATE}[1,2]", settings.get(FluxCapacitorSettings.READ_DESCRIPTOR).toString());
        }catch (Exception err){
            err.printStackTrace();
            fail();
        }
    }
}
