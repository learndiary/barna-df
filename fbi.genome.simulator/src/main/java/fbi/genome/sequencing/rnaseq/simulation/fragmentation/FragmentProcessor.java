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
     * @param id the id
     * @param cs the read
     * @param start the start
     * @param end the end
     * @param len the length
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
     * Return status string or null when processing is done
     *
     * @return status status or null
     */
    String done();
}
