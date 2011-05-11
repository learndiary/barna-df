package fbi.genome.sequencing.rnaseq.simulation;

import org.junit.Test;

import static junit.framework.Assert.assertEquals;
import static junit.framework.Assert.fail;

/**
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class ProfilerTest {

    @Test
    public void testConstructor() throws Exception {
        try{
            Profiler profiler = new Profiler(null);
            fail();
        }catch (NullPointerException e){
            // expected
            assertEquals("You have to specify settings! NULL not permitted.", e.getMessage());
        }
    }
}
