package barna.flux.capacitor;

import barna.flux.capacitor.reconstruction.FluxCapacitor;
import barna.flux.capacitor.reconstruction.FluxCapacitorSettings;
import org.junit.Test;

import java.io.File;

/**
 * Test writing sam/bam output
 * @author Thasso Griebel <thasso.griebel@gmail.com>
 */
public class WriteSamTests {

    @Test
    public void testWritingPairedEndSam() throws Exception {
        FluxCapacitorSettings settings = new FluxCapacitorSettings();
        //settings.set(FluxCapacitorSettings.ANNOTATION_FILE, new File("/Users/thasso/data/capacitor-test-data/chr21.gtf"));
        settings.set(FluxCapacitorSettings.ANNOTATION_FILE, new File(getClass().getResource("/single_multimap.gtf").getFile()));
        settings.set(FluxCapacitorSettings.MAPPING_FILE, new File(getClass().getResource("/single_multimap.bam").getFile()));
        settings.set(FluxCapacitorSettings.SORT_IN_RAM, true);
        FluxCapacitor capacitor = new FluxCapacitor(settings);
        capacitor.call();
        System.out.println("Done");
    }
}
