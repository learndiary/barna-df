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

package fbi.genome.sequencing.rnaseq.simulation.distributions;

import org.junit.Test;

import static junit.framework.Assert.assertEquals;

/**
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class EmpiricalDistributionTest {

    @Test
    public void testPropabilityAndFrequency() throws Exception {
        double[] dist = new double[]{1.0,2.0,2.0,3.0};

        double[] edist = new double[3];
        for (int i = 0; i < dist.length; i++) {
            double v = dist[i];
            EmpiricalDistribution.addToBin(v, edist, 1.0, 3.0);
        }
        EmpiricalDistribution distribution = new EmpiricalDistribution(edist, 1.0, 3.0, 2.0);


        assertEquals(0.5, distribution.getP(2.0), 0.00001);
        assertEquals(1.0, distribution.getRelFreq(2.0), 0.00001);
    }
}
