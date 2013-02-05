package barna.astalavista;

import barna.commons.Execute;
import barna.io.FileHelper;
import barna.io.gtf.GTFwrapper;
import barna.model.ASEvent;
import barna.model.Gene;
import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.EnumSet;
import java.util.Vector;

/**
 * Created with IntelliJ IDEA.
 * User: micha
 * Date: 11/23/12
 * Time: 7:07 PM
 * To change this template use File | Settings | File Templates.
 */
public class EventExtractorTest {

    @BeforeClass
    public static void initExecuter() {
        Execute.initialize(2);
    }

    @AfterClass
    public static void shutdownExecuter() {
        Execute.shutdown();
    }

    private static Gene[] getGene(String[] gtf) throws Exception {
        if (gtf== null)
            return null;

        // write file
        final File f= FileHelper.createTempFile(EventExtractorTest.class.getName(), ".gtf");
        BufferedWriter buffy= new BufferedWriter(new FileWriter(f));
        for (int i = 0; i < gtf.length; i++)
            buffy.write(gtf[i]+ "\n");
        buffy.close();

        // read file
        GTFwrapper myReader= new GTFwrapper(f);
        myReader.read();
        return myReader.getGenes();
    }

    @Test
    public void testAlternativeTranscriptStart() throws Exception {
        final String[] gtf= new String[]
                {"XXX\tfoo\texon\t10\t15\t.\t+\t.\ttranscript_id=\"aTranscript\";",
                "XXX\tfoo\texon\t5\t15\t.\t+\t.\ttranscript_id=\"anotherTranscript\";"};
        Gene[] ge= getGene(gtf);

        AStaSettings settings= new AStaSettings();
        settings.set(AStaSettings.EVENTS, EnumSet.of(
                AStaSettings.EventTypes.ASI,
                AStaSettings.EventTypes.ASE,
                AStaSettings.EventTypes.VST)
        );
        EventExtractor extractor= new EventExtractor(ge[0], settings);

        extractor.constructGraph();
        extractor.contractGraph(2);
        extractor.getEventsByPartitions(2);
        Vector<ASEvent> v= extractor.getEventV();

        assertEquals(1, v.size());
        assertEquals("1[,2[", v.elementAt(0).toString());
    }

    @Test
    public void testAlternativeCleavage() throws Exception {
        final String[] gtf= new String[]
                {"XXX\tfoo\texon\t10\t15\t.\t+\t.\ttranscript_id=\"aTranscript\";",
                        "XXX\tfoo\texon\t10\t20\t.\t+\t.\ttranscript_id=\"anotherTranscript\";"};
        Gene[] ge= getGene(gtf);
        AStaSettings settings= new AStaSettings();
        settings.set(AStaSettings.EVENTS, EnumSet.of(
                AStaSettings.EventTypes.ASI,
                AStaSettings.EventTypes.ASE,
                AStaSettings.EventTypes.VST)
        );
        EventExtractor extractor= new EventExtractor(ge[0], settings);

        extractor.constructGraph();
        extractor.contractGraph(2);

        extractor.getEventsByPartitions(2);
        Vector<ASEvent> v= extractor.getEventV();

        assertEquals(1, v.size());
        assertEquals("1],2]", v.elementAt(0).toString());
    }


}
