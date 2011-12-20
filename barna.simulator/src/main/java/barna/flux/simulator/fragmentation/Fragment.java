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

package barna.flux.simulator.fragmentation;

/**
 * Model fragments based on the global id ({@code <chromosome>:<position>@<transcriptID>})
 * start and end coordinates.
 */
public class Fragment {
    /**
     * The global id
     */
    private CharSequence id;
    /**
     * The start
     */
    private int start;

    /**
     * The end
     */
    private int end;

    /**
     * Number of duplications created during amplification
     */
    private int duplicates;

    /**
     * Create a new Fragment
     *
     * @param id    the id
     * @param start the start
     * @param end   the end
     */
    public Fragment(final CharSequence id, final int start, final int end) {
        if (id == null) {
            throw new NullPointerException("Fragment with NULL id not permitted");
        }
        if (id.length() == 0) {
            throw new IllegalArgumentException("Fragment with empty ID not permitted");
        }
        if (start > end) {
            throw new IllegalArgumentException("Fragment start > end! " + start + "->" + end);
        }
        this.id = id;
        this.start = start;
        this.end = end;
    }

    /**
     * Get the ID
     *
     * @return id the global fragment id
     */
    public CharSequence getId() {
        return id;
    }

    /**
     * Get the start
     *
     * @return start the start
     */
    public int getStart() {
        return start;
    }

    /**
     * Get the end
     *
     * @return end the end
     */
    public int getEnd() {
        return end;
    }

    /**
     * The length
     *
     * @return length the fragment length
     */
    public int length() {
        return end - start + 1;
    }

    /**
     * Get the number of duplicates created during amplification
     *
     * @return duplicates the number of duplicates
     */
    public int getDuplicates() {
        return duplicates;
    }

    /**
     * Set the number of duplicates
     *
     * @param duplicates number of duplicates
     */
    public void setDuplicates(final int duplicates) {
        if(duplicates < 0) throw new IllegalArgumentException("Number of Duplicates must be >= 0");
        this.duplicates = duplicates;
    }

    @Override
    public String toString() {
        return start + "\t" + end + "\t" + id + "\t" + duplicates;
    }
}
