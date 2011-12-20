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

package barna.genome.sequencing.rnaseq.simulation.error;

/**
 * Distribution of the average read quality
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
class ReadQualityDistribution extends Distribution {

    public ReadQualityDistribution(int size) {
        super(size);
    }

    public void addRead(Read read) {
        reads++;
        int[] q = read.getQualities();
        int s = 0;
        for (int i = 0; i < q.length; i++) {
            s += q[i];
        }
        int avg = s / read.getLength();
        values[avg]++;
    }
}
