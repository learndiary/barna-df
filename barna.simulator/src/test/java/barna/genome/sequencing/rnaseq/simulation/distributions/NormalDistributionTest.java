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

package barna.genome.sequencing.rnaseq.simulation.distributions;

import org.junit.Test;

import static junit.framework.Assert.assertEquals;

/**
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class NormalDistributionTest {

    @Test
    public void testName() throws Exception {
        NormalDistribution d = new NormalDistribution(0, 1, false);
        EmpiricalDistribution e = new EmpiricalDistribution(d, 1000, 4.0);
        assertEquals(0.0, e.getP(0), 0.000001);
        assertEquals(0.0, e.getP(1), 0.000001);
        assertEquals(0.0031, e.getP(0.5), 0.0001);
        assertEquals(0.999, e.getSum(), 0.001);
    }
}
