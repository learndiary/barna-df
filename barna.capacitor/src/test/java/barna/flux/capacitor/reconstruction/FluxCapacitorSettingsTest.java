package barna.flux.capacitor.reconstruction;

import barna.commons.parameters.Parameter;
import barna.commons.parameters.ParameterSchema;
import org.junit.Test;

import java.io.*;
import java.lang.reflect.Field;
import java.util.EnumSet;

import static junit.framework.Assert.*;

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

    @Test
    public void testReflectionForCountElements() throws Exception {
        FluxCapacitorSettings settings = new FluxCapacitorSettings();
        settings.set(FluxCapacitorSettings.COUNT_ELEMENTS,EnumSet.of(FluxCapacitorSettings.CountElements.INTRONS));
        BufferedWriter wbuffy = new BufferedWriter(new FileWriter("/home/emilio/test.par"));
        wbuffy.write(settings.toString());
        wbuffy.close();

        FluxCapacitorSettings settings2 =  new FluxCapacitorSettings();
        Field field = settings2.getClass().getField("COUNT_ELEMENTS");
        Parameter p = (Parameter)field.get(null);
        assertTrue(p.getValuesString().length() >0);
        assertNotNull(p);
        BufferedInputStream rbuffy = new BufferedInputStream(new FileInputStream("/home/emilio/test.par"));
        settings2.parse(rbuffy);
        //assertNotNull(settings2.get(FluxCapacitorSettings.COUNT_ELEMENTS));
    }
}
