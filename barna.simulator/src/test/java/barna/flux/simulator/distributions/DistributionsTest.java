/*
 * This file is part of the Flux Library.
 *
 * The code of the Flux Library may be freely distributed and modified under the terms of the
 * European Union Public Licence (EUPL) published on the web site <http://www.osor.eu/eupl/european-union-public-licence-eupl-v.1.1>.
 * Copyright for the code is held jointly by the individual authors, who should be listed
 * in @author doc comments. According to Article 5 and Article 11 of the EUPL, publications that
 * include results produced by the Flux Library are liable to reference the Work,
 * see the Flux Library homepage <http://flux.sammeth.net> for more information.
 */

package barna.flux.simulator.distributions;

import barna.genome.sequencing.rnaseq.simulation.distributions.AbstractDistribution;
import barna.genome.sequencing.rnaseq.simulation.distributions.Distributions;
import barna.genome.sequencing.rnaseq.simulation.distributions.NormalDistribution;
import org.junit.Test;

import static junit.framework.Assert.assertEquals;
import static junit.framework.Assert.assertTrue;

/**
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class DistributionsTest {

    @Test
    public void testParser() throws Exception {
        AbstractDistribution d = Distributions.parseDistribution("N(1.0,1.0)");
        assertEquals(1.0, d.getMean(), 0.00001);
        //assertEquals(1.0, d.get(), 0.00001);
        assertTrue(d instanceof NormalDistribution);
        assertEquals(1.0, ((NormalDistribution)d).getSd(), 0.00001);


        d = Distributions.parseDistribution("N(0.1,0.1)");
        assertEquals(0.1, d.getMean(), 0.00001);
        //assertEquals(1.0, d.get(), 0.00001);
        assertTrue(d instanceof NormalDistribution);
        assertEquals(0.1, ((NormalDistribution)d).getSd(), 0.00001);

    }


}
