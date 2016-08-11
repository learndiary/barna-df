package barna.astafunk;

import barna.astalavista.AStaSettings;
import barna.astalavista.EventExtractor;
import barna.astafunk.utils.OrderingEventComparator;
import barna.io.gtf.GTFwrapper;
import barna.model.ASEvent;
import barna.model.Gene;
import barna.model.Graph;
import org.junit.Test;

import java.io.File;
import java.util.Collections;
import java.util.EnumSet;
import java.util.List;

/**
 * @author vitorcoelho
 * @version 2
 * @since 07/12/15
 */
public class GroupSplitTest {

    String gtf = "";
    Gene[] genes = null;
    String genome = "";

    private static GTFwrapper wrapper = null;

    public void parseGTF() {

        File f = new File(gtf);

        wrapper = new GTFwrapper(f);

        if (!wrapper.isApplicable()) {
            f = wrapper.sort();
            wrapper = new GTFwrapper(f);
        }

        //Read all genes
        wrapper.setReadAll(true);

        //wrapper.setChromosomeWise(false);
        //wrapper.setGeneWise(false);
        wrapper.setReadAheadTranscripts(-1);

        //Read the file
        wrapper.read();

        //Get genes
        genes = wrapper.getGenes();

        //Define path to fastas
        Graph.overrideSequenceDirPath = genome;
    }

    protected List<ASEvent> extractAllEvents(Gene gene) {

        AStaSettings astaSettings = new AStaSettings();

        /**
         * ASI:
         * ASE:
         */
        astaSettings.set(AStaSettings.EVENTS, EnumSet.of(AStaSettings.EventTypes.ASI,
                AStaSettings.EventTypes.ASE));

        astaSettings.set(AStaSettings.EVENTS_DIMENSION, -1);

        EventExtractor extractor = new EventExtractor(gene, astaSettings);

        //Constructing the graph
        extractor.constructGraph();

        //Contract the graph
        extractor.contractGraph(0);

        //Get events by partitions
        extractor.getEventsByPartitions(-1);
        List<ASEvent> l = extractor.getEventV();
        Collections.sort(l, new OrderingEventComparator());

        return l;
    }

    @Test
    public void testSplitVariant(){

        // TODO micha was here
        if (1== 1)
            return;

            // Parse GTF

        parseGTF();

        // Extract events



        // Group Events

        // Split Variants

    }
}
