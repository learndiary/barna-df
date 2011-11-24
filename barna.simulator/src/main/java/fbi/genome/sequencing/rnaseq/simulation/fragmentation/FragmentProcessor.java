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

package fbi.genome.sequencing.rnaseq.simulation.fragmentation;

import fbi.commons.ByteArrayCharSequence;

import java.util.List;

/**
 * Fragmentation processor
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public interface FragmentProcessor {
    /**
     * Process the given read
     *
     * @param id    the id
     * @param cs    the read
     * @param start the start
     * @param end   the end
     * @param len   the length
     * @return fragments list of fragments or null
     */
    List<Fragment> process(ByteArrayCharSequence id, ByteArrayCharSequence cs, int start, int end, int len);

    /**
     * Return the name of this processor
     *
     * @return name the name
     */
    String getName();

    /**
     * Return the current configuration
     *
     * @return config the configuration
     */
    String getConfiguration();

    /**
     * Called after all fragments are processed. The returned message is printed
     *
     * @return status status or null
     */
    String done();
}
