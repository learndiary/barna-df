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

package fbi.genome.errormodel;

import java.util.ArrayList;
import java.util.List;

/**
 * Read information. Mutable class so we do not have to create a new object for each read.
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
class Read {
    /**
     * The read name
     */
    private CharSequence name;
    /**
     * The sequence
     */
    private CharSequence sequence;
    /**
     * read length
     */
    private int length;
    /**
     * Qualities per position
     */
    private int[] qualities;
    /**
     * Mappings
     */
    private List<Mapping> mappings;


    /**
     * Get the read ID
     *
     * @return ID the read id
     */
    public CharSequence getName() {
        return name;
    }

    /**
     * Set the read ID
     *
     * @param name the read ID
     */
    public void setName(CharSequence name) {
        this.name = name;
    }

    /**
     * Get the sequence
     *
     * @return seq the sequnece
     */
    public CharSequence getSequence() {
        return sequence;
    }

    /**
     * Set the sequence
     *
     * @param sequence the sequence
     */
    public void setSequence(CharSequence sequence) {
        this.sequence = sequence;
        setLength(sequence.length());

        // init the quality array
        if (qualities == null || qualities.length < getLength()) qualities = new int[getLength()];
    }

    /**
     * Get the read length
     *
     * @return length the read length
     */
    public int getLength() {
        return length;
    }

    /**
     * Set the read length
     *
     * @param length the length
     */
    public void setLength(int length) {
        this.length = length;
    }

    /**
     * Get the qualities
     *
     * @return qualities the qualities per position
     */
    public int[] getQualities() {
        return qualities;
    }

    /**
     * Set the qualities
     *
     * @param qualities the qualities
     */
    public void setQualities(int[] qualities) {
        this.qualities = qualities;
    }

    /**
     * Add a new empty mapping and return it
     *
     * @return mapping the created mapping
     */
    public Mapping addMapping() {
        if (mappings == null) mappings = new ArrayList<Mapping>();
        Mapping mapping = new Mapping();
        mappings.add(mapping);
        return mapping;
    }

    /**
     * Get the list of mappings
     *
     * @return mappings the list of mappings or null
     */
    public List<Mapping> getMappings() {
        return mappings;
    }

    /**
     * Reset this read. Clears ID, sequence, length and the mappings
     */
    public void reset() {
        if (mappings != null) mappings.clear();
        length = 0;
        name = null;
        sequence = null;
    }

    /**
     * Mapping
     */
    class Mapping {
        /**
         * List of missmatches
         */
        private List<Missmatch> missmatches;
        /**
         * The quality value
         */
        private int quality;
        /**
         * The name
         */
        private String name;

        /**
         * Add a missmatch
         *
         * @param position         the position
         * @param genomicCharacter the genomic character
         */
        public void addMissmatch(int position, char genomicCharacter) {
            if (missmatches == null) missmatches = new ArrayList<Missmatch>();
            missmatches.add(new Missmatch(position, genomicCharacter, getSequence().charAt(position - 1)));
        }

        /**
         * Get a list of missmatches
         *
         * @return missmatches the list of missmatches
         */
        public List<Missmatch> getMissmatches() {
            return missmatches;
        }

        /**
         * Get the genomic characters for missmatches at the given position
         *
         * @param position the position
         * @return missmatches the genomic characters for missmachtes at the given position
         */
        public List<Character> getMissmatches(int position) {
            if (missmatches == null || missmatches.size() == 0) return null;
            ArrayList<Character> cc = new ArrayList<Character>();
            for (Missmatch missmatch : missmatches) {
                if (missmatch.getPosition() == position) cc.add(missmatch.getGenomicCharacter());
            }
            return cc.size() == 0 ? null : cc;
        }

        /**
         * Set the quality
         *
         * @param quality the quality
         */
        public void setQuality(int quality) {
            this.quality = quality;
        }

        /**
         * Get the quality
         *
         * @return quality the quality
         */
        public int getQuality() {
            return quality;
        }

        /**
         * Set the name
         *
         * @param name the name
         */
        public void setName(String name) {
            this.name = name;
        }

        /**
         * Get the name
         *
         * @return name the name
         */
        public String getName() {
            return name;
        }
    }

    /**
     * Missmatch
     */
    class Missmatch {
        /**
         * The position
         */
        private int position;
        /**
         * The genomic character
         */
        private char genomicCharacter;
        /**
         * The read character
         */
        private char readCharacter;

        /**
         * Create a new missmatch
         *
         * @param position         the position
         * @param genomicCharacter the genomic character
         * @param readCharacter    the read character
         */
        Missmatch(int position, char genomicCharacter, char readCharacter) {
            this.position = position;
            this.genomicCharacter = genomicCharacter;
            this.readCharacter = readCharacter;
        }

        /**
         * Get the positions
         *
         * @return position the position
         */
        public int getPosition() {
            return position;
        }

        /**
         * Set the position
         *
         * @param position the position
         */
        public void setPosition(int position) {
            this.position = position;
        }

        /**
         * Get the genomic character
         *
         * @return genomic the genomic character
         */
        public char getGenomicCharacter() {
            return genomicCharacter;
        }

        /**
         * Set the genomic character
         *
         * @param genomicCharacter the genomic character
         */
        public void setGenomicCharacter(char genomicCharacter) {
            this.genomicCharacter = genomicCharacter;
        }

        /**
         * Get the read character
         *
         * @return readCharacter the read character
         */
        public char getReadCharacter() {
            return readCharacter;
        }

        /**
         * Set the read character
         *
         * @param readCharacter the read character
         */
        public void setReadCharacter(char readCharacter) {
            this.readCharacter = readCharacter;
        }
    }
}
