package barna.flux.simulator.error;

import barna.model.Qualities;
import com.thoughtworks.xstream.XStream;
import org.junit.Test;

import java.io.File;
import java.io.IOException;

import static org.junit.Assert.*;

public class MarkovErrorModelTest {

    @Test
    public void testThatErrorModelLoadingWorksAlsoForOldModels() {
        // models before refactoring to barna package
        // (BARNA-86)
        try {
            QualityErrorModel model_35 = MarkovErrorModel.loadErrorModel(
                    new File(getClass().getResource("/before_barna_35.error.model.gz").getFile())
            );
            assertNotNull(model_35);
        } catch (IOException e) {
            e.printStackTrace();
            fail();
        }
        try {
            QualityErrorModel model_76 = MarkovErrorModel.loadErrorModel(
                    new File(getClass().getResource("/before_barna_76.error.model.gz").getFile())
            );
            assertNotNull(model_76);
        } catch (IOException e) {
            e.printStackTrace();
            fail();
        }
    }
    @Test
    public void testThatDistributedModelsCanBeLoaded() {
        // models before refactoring to barna package
        // (BARNA-86)
        try {
            QualityErrorModel model_35 = MarkovErrorModel.loadErrorModel(
                    new File(getClass().getResource("/35_error.model").getFile())
            );
            assertNotNull(model_35);
        } catch (IOException e) {
            e.printStackTrace();
            fail();
        }
        try {
            QualityErrorModel model_76 = MarkovErrorModel.loadErrorModel(
                    new File(getClass().getResource("/76_error.model").getFile())
            );
            assertNotNull(model_76);
        } catch (IOException e) {
            e.printStackTrace();
            fail();
        }
    }

    @Test
    public void testThatTheErrorModelWriterUsesOnlyTheSimpleClassName(){
        QualityErrorModel model = new QualityErrorModel(Qualities.Technology.Illumina13, 36, new QualityTransitions(10, 36), new CrossTalkModel(4, true));
        XStream streamer = MarkovErrorModel.createXStream();
        String string = streamer.toXML(model);
        assertTrue(string.startsWith("<QualityErrorModel>"));

        // test reloading
        Object loaded_model = streamer.fromXML(string);
        assertNotNull(loaded_model);

    }
}
