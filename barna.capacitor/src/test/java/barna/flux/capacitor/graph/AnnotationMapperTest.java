package barna.flux.capacitor.graph;

import barna.flux.capacitor.reconstruction.FluxCapacitorSettings;
import barna.io.BufferedIterator;
import barna.io.BufferedIteratorRAM;
import barna.io.bed.BEDwrapper;
import barna.io.gtf.GTFwrapper;
import barna.io.rna.UniversalReadDescriptor;
import barna.model.Gene;
import junit.framework.Assert;
import junit.framework.TestCase;
import org.junit.Test;

import java.io.*;
import java.util.zip.ZipOutputStream;

/**
 * Created with IntelliJ IDEA.
 * User: Emilio Palumbo
 * Date: 5/31/12
 * Time: 12:39 PM
 */
public class AnnotationMapperTest extends TestCase {

    private final File gtfFile = new File(getClass().getResource("/mm9_chr1_chrX.gtf").getFile());
    private final File bedFile = new File(getClass().getResource("/chr1_chrX.bed").getFile());

    @Test
    public void testWriteSJReads() throws Exception {
        UniversalReadDescriptor descriptor= new UniversalReadDescriptor();
        descriptor.init(UniversalReadDescriptor.getDescriptor("SIMULATOR"));
        FluxCapacitorSettings settings= new FluxCapacitorSettings();
        settings.set(FluxCapacitorSettings.ANNOTATION_FILE,
                new File(gtfFile.getAbsolutePath()));
        settings.set(FluxCapacitorSettings.MAPPING_FILE,
                new File(bedFile.getAbsolutePath()));
        settings.set(FluxCapacitorSettings.READ_DESCRIPTOR,
                descriptor);
        settings.set(FluxCapacitorSettings.SORT_IN_RAM,
                false);
        settings.set(FluxCapacitorSettings.KEEP_SORTED_FILES,
                false);
        settings.set(FluxCapacitorSettings.ANNOTATION_MAPPING,
                FluxCapacitorSettings.AnnotationMapping.PAIRED);
        settings.set(FluxCapacitorSettings.STDOUT_FILE,
                null);
        settings.set(FluxCapacitorSettings.STATS_FILE,
                null);
        GTFwrapper gtf = new GTFwrapper(gtfFile);
        BEDwrapper bed = new BEDwrapper(bedFile);
        try {
            gtf.read();
            Gene g = gtf.getGenes()[0];
            AnnotationMapper a = new AnnotationMapper(g);
            BufferedIterator iter = bed.readBedFile(g, g.getStart(), g.getEnd(), true, settings.get(FluxCapacitorSettings.READ_DESCRIPTOR),null);
            a.map(iter, settings);
            ZipOutputStream out = new ZipOutputStream(new BufferedOutputStream(new FileOutputStream("/home/emilio/test.zip")));
            a.writeSJReads(out);
            out.close();
        } catch (Exception e) {
            e.printStackTrace();
            Assert.fail();
        }
    }
}
