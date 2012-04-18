package barna.flux.capacitor.reconstruction;

import org.junit.Test;

import java.io.ByteArrayInputStream;
import java.io.StringReader;

import static junit.framework.Assert.assertEquals;
import static junit.framework.Assert.fail;

public class FluxCapacitorSettingsTest {
    @Test
    public void testThatSpacesAreAcceptedInReadDescriptor() throws Exception {
        FluxCapacitorSettings settings = new FluxCapacitorSettings();
        ByteArrayInputStream input = new ByteArrayInputStream(("READ_DESCRIPTOR\t{ID} {MATE}\n" +
                "MAPPING_FILE "+getClass().getResource("/chr1_chrX.bed").getFile()+"\n" +
                "ANNOTATION_FILE "+getClass().getResource("/mm9_chr1_chrX.gtf").getFile()+"\n").getBytes());

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
