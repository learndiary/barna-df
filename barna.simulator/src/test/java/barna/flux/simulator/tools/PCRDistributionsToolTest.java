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

package barna.flux.simulator.tools;

import barna.genome.sequencing.rnaseq.simulation.distributions.GCPCRDistribution;
import barna.genome.sequencing.rnaseq.simulation.tools.PCRDistributionsTool;
import org.junit.Test;

import java.io.InputStream;

import static junit.framework.Assert.assertEquals;
import static junit.framework.Assert.assertNotNull;

/**
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class PCRDistributionsToolTest {

    @Test
    public void testLoadDefault() throws Exception {
        InputStream inputStream = getClass().getResource("/pcr_15_20.dat").openStream();
        GCPCRDistribution def = PCRDistributionsTool.load(inputStream);
        assertNotNull(def);
        assertEquals(15, def.getGenerations());
    }
}
