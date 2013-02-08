package barna.flux.capacitor.reconstruction;

import barna.commons.parameters.ParameterException;
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

    @Test
    public void testNullReadDescriptor() throws Exception {
        FluxCapacitorSettings setting = new FluxCapacitorSettings();

        setting.set(FluxCapacitorSettings.MAPPING_FILE.getName(),getClass().getResource("/mm9_chr1_chrX_sorted.bed").getFile());
        setting.set(FluxCapacitorSettings.ANNOTATION_FILE.getName(), getClass().getResource("/mm9_chr1_chrX_sorted.gtf").getFile());
        setting.set(FluxCapacitorSettings.READ_STRAND.getName(), "NONE");
        setting.set(FluxCapacitorSettings.ANNOTATION_MAPPING.getName(),"PAIRED");

        try {
            setting.validate();
        } catch (ParameterException ex) {
            fail(ex.getMessage());
        }
    }

    @Test
    public void testNonNullReadDescriptorStranded() throws Exception {
        FluxCapacitorSettings setting = new FluxCapacitorSettings();

        setting.set(FluxCapacitorSettings.MAPPING_FILE.getName(),getClass().getResource("/mm9_chr1_chrX_sorted.bed").getFile());
        setting.set(FluxCapacitorSettings.ANNOTATION_FILE.getName(), getClass().getResource("/mm9_chr1_chrX_sorted.gtf").getFile());
        setting.set(FluxCapacitorSettings.READ_DESCRIPTOR.getName(),"MATE2_SENSE");
        setting.set(FluxCapacitorSettings.READ_STRAND.getName(), "NONE");
        setting.set(FluxCapacitorSettings.ANNOTATION_MAPPING.getName(),"PAIRED_STRANDED");

        try {
            setting.validate();
        } catch (ParameterException ex) {
            fail(ex.getMessage());
        }
    }

    @Test
    public void testNonNullReadDescriptorStrandInfo() throws Exception {
        FluxCapacitorSettings setting = new FluxCapacitorSettings();

        setting.set(FluxCapacitorSettings.MAPPING_FILE.getName(),getClass().getResource("/mm9_chr1_chrX_sorted.bed").getFile());
        setting.set(FluxCapacitorSettings.ANNOTATION_FILE.getName(), getClass().getResource("/mm9_chr1_chrX_sorted.gtf").getFile());
        setting.set(FluxCapacitorSettings.READ_DESCRIPTOR.getName(),"PAIRED");
        setting.set(FluxCapacitorSettings.READ_STRAND.getName(), "MATE2_SENSE");
        setting.set(FluxCapacitorSettings.ANNOTATION_MAPPING.getName(),"PAIRED_STRANDED");

        try {
            setting.validate();
        } catch (ParameterException ex) {
            fail(ex.getMessage());
        }
    }

    @Test
    public void testAnnotationMappingStrandedNoStrand() throws Exception {
        FluxCapacitorSettings setting = new FluxCapacitorSettings();

        setting.set(FluxCapacitorSettings.MAPPING_FILE.getName(),getClass().getResource("/mm9_chr1_chrX_sorted.bed").getFile());
        setting.set(FluxCapacitorSettings.ANNOTATION_FILE.getName(), getClass().getResource("/mm9_chr1_chrX_sorted.gtf").getFile());
        setting.set(FluxCapacitorSettings.READ_STRAND.getName(), "NONE");


        String[] ams = {"SINGLE_STRANDED", "PAIRED_STRANDED"};

        for (String am : ams) {
            try {
                setting.set(FluxCapacitorSettings.ANNOTATION_MAPPING.getName(),am);
                setting.validate();
                fail("Exception not throw!!!");
            } catch(ParameterException ex) {
                assertEquals("Annotation mapping " + am + " requires strand information.", ex.getMessage());
            }
        }
    }

    @Test
    public void testAnnotationMappingPairedNoPaired() throws Exception {
        FluxCapacitorSettings setting = new FluxCapacitorSettings();

        setting.set(FluxCapacitorSettings.MAPPING_FILE.getName(),getClass().getResource("/mm9_chr1_chrX_sorted.bed").getFile());
        setting.set(FluxCapacitorSettings.ANNOTATION_FILE.getName(), getClass().getResource("/mm9_chr1_chrX_sorted.gtf").getFile());
        setting.set(FluxCapacitorSettings.READ_STRAND.getName(), "SENSE");


        String[] ams = {"PAIRED", "PAIRED_STRANDED"};

        for (String am : ams) {
            try {
                setting.set(FluxCapacitorSettings.ANNOTATION_MAPPING.getName(),am);
                setting.validate();
                fail("Exception not throw!!!");
            } catch(ParameterException ex) {
                assertEquals("Annotation mapping " + am + " requires paired reads.", ex.getMessage());
            }
        }
    }
}
