package barna.genome.sequencing.rnaseq.simulation.fragmentation;

import barna.genome.sequencing.rnaseq.simulation.distributions.EmpiricalDistribution;
import org.junit.Test;

import java.io.IOException;
import java.net.URL;

import static junit.framework.Assert.assertNotNull;
import static junit.framework.Assert.fail;

public class FragmenterTest {


    @Test
    public void testThatTheExpAllSizeDistributionSumIsValid(){
            // load the default
            URL url = Fragmenter.class.getResource("/expAll.isizes");
            try {
                // get the attributes ?
                EmpiricalDistribution dist = EmpiricalDistribution.create(url.openStream(), 209208, 1, 297, 100, EmpiricalDistribution.LINE_PARSER_SINGLE_VALUE);
                assertNotNull(dist);
            } catch (IOException e) {
                e.printStackTrace();
                fail();
            }

    }
}
