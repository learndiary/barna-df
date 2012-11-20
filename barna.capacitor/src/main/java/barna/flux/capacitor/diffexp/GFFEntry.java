/*
 * Copyright (c) 2012, Micha Sammeth, Thasso Griebel, Emilio Palumbo
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *      * Redistributions of source code must retain the above copyright
 *        notice, this list of conditions and the following disclaimer.
 *      * Redistributions in binary form must reproduce the above copyright
 *        notice, this list of conditions and the following disclaimer in the
 *        documentation and/or other materials provided with the distribution.
 *      * The names of its contributors may be not used to endorse or promote
 *        products derived from this software without specific prior written
 *        permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 *  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 *  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 *  DISCLAIMED. IN NO EVENT SHALL MICHA SAMMETH BE LIABLE FOR ANY
 *  DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 *  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 *  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 *  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 *  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 *  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

package barna.flux.capacitor.diffexp;

import java.io.IOException;
import java.util.*;

/**
 * Basic GTF Entry with the default fields as specified at
 *
 * <a href="http://genome.ucsc.edu/FAQ/FAQformat.html#format3"> UCSC</a>
 * <p>
 * The class provides a {@link #parse(String)} method to create
 * an instance from a GFF line. Note that the line MUST consist of 9
 * tab separated fields.
 * </p>
 *
 * @author Thasso Griebel <thasso.griebel@gmail.com>
 */
class GFFEntry {
    /**
     * !SORTED! Valid strand entries
     */
    private static final char[] VALID_STRANDS = new char[]{'+', '-', '.'};
    /**
     * !SORTED! Valid coding frame entries
     */
    private static final char[] VALID_CODING_FRAMES = new char[]{'.','0', '1', '2'};

    /**
     * The Chromosome
     */
    private final String chromosome;
    /**
     * The entry source
     */
    private final String source;

    /**
     * The feature name
     */
    private final String feature;

    /**
     * Start position (1-based)
     */
    private final long start;

    /**
     * End position (inclusive) (1-based)
     */
    private final long end;

    /**
     * Score between 0 and 1000 or -1 for no score
     */
    private final short score;

    /**
     * +/-/. strand information
     */
    private final char strand;

    /**
     * Coding frame: -1 for no coding frame
     * otherwise a value of 0-2 for the reading frame
     */
    private final char codingFrame;

    /**
     * GTF group
     */
    private final String group;
    /**
     * Additional attributes
     */
    private Map<String, String> attributes;

    /**
     * Create a new GFFEntry
     *
     * @param chromosome the chromosome
     * @param source the source name
     * @param start start position 1-based
     * @param end end position 1-based inclusive
     * @param score the score (between 0 and 1000) or -1 for no score
     * @param strand the strand (+/- or . for unknown)
     * @param codingFrame the coding frame 0,1,2 or . for unknown
     * @param group optional group entry, can be null or or empty string
     */
    public GFFEntry(final String chromosome, final String source, final String feature, final long start, final long end, final short score,
             final char strand, final char codingFrame, final String group) {
        if(chromosome == null || chromosome.isEmpty()) throw new IllegalArgumentException("No valid Chromosome specified: " + chromosome);
        if(source == null || source.isEmpty()) throw new IllegalArgumentException("No valid Source specified: " + source);
        if(feature == null || feature.isEmpty()) throw new IllegalArgumentException("No feature specified: " + feature);
        if(start <= 0) throw new IllegalArgumentException("Start has to be > 0: " + start);
        if(end <= 0 || end < start) throw new IllegalArgumentException("End has to be 0 < start <= end: " + end);
        if(score < -1 || score > 1000) throw new IllegalArgumentException("Score has to be -1 <= score <= 1000: " + score);
        if(Arrays.binarySearch(VALID_STRANDS, strand) < 0) throw new IllegalArgumentException("Invalid strand: " + strand);
        if(Arrays.binarySearch(VALID_CODING_FRAMES, codingFrame) < 0) throw new IllegalArgumentException("Invalid coding frame: " + codingFrame);

        this.chromosome = chromosome;
        this.source = source;
        this.feature = feature;
        this.start = start;
        this.end = end;
        this.score = score;
        this.strand = strand;
        this.codingFrame = codingFrame;
        this.group = group;
    }

    /**
     * Get the chromosome name of the entry
     *
     * @return chr the chromosome name
     */
    public String getChromosome() {
        return chromosome;
    }

    /**
     * Get the source value
     * @return source the source value
     */
    public String getSource() {
        return source;
    }

    /**
     * Get the feature name
     * @return feature the feature name
     */
    public String getFeature() {
        return feature;
    }

    /**
     * Get the 1-based start position
     * @return start the 1-based start
     */
    public long getStart() {
        return start;
    }

    /**
     * Get the 1-based end position (inclusive)
     * @return end 1-based end (inclusive)
     */
    public long getEnd() {
        return end;
    }

    /**
     * Get the score between 0 and 1000
     * @return score the score between 0 and 1000
     */
    public short getScore() {
        return score;
    }

    /**
     * Get the strand as +/-/.
     * @return strand one of +/-/.
     */
    public char getStrand() {
        return strand;
    }

    /**
     * The coding frame as 0/1/2/.
     *
     * @return frame the coding frame as 0/1/2/.
     */
    public char getCodingFrame() {
        return codingFrame;
    }

    /**
     * Optional group entry, might be null or empty
     *
     * @return group the group entry
     */
    public String getGroup() {
        return group;
    }

    /**
     * Get additional attributes from group field
     *
     * @return attributes additional attributes from group field
     */
    public Map<String, String> getAttributes() {
        if(attributes == null){
            attributes = parseAttributes();
        }
        return attributes;
    }

    /**
     * Length of the entry as end - start + 1
     *
     * @return length the length of the entry
     */
    public long getLength(){
        if(end <=0 || start <= 0) throw new IllegalArgumentException("Start or end are <= 0 ! Transcript is not initialized!");
        return end - start + 1; // +1 because end is inclusive
    }

    /**
     * Parse the attributes in the group field and
     * return the attribute map with all key in lower case.
     * <p>
     *     This parser does not assume an attribute map in the group field. Therefore,
     *     if the group can not be splitted properly into key value pairs, the
     *     attribute is silently ignored and not added to the map.
     * </p>
     * @return attributes attribute map
     */
    private Map<String, String> parseAttributes(){
        if(group == null || group.isEmpty()) return Collections.EMPTY_MAP;
        String[] split = group.split("; ");
        HashMap<String, String> atrs = new HashMap<String, String>();
        for (String s : split) {
            String[] key_value = s.split("\\s+", 2);
            if(key_value.length == 2){
                if(key_value[1].startsWith("\"") && key_value[1].endsWith("\"")){
                    key_value[1] = key_value[1].substring(1, key_value[1].length()-1);
                }
                atrs.put(key_value[0].toLowerCase(), key_value[1]);
            }
        }
        return atrs;
    }

    /**
     * Parse a GFF line
     *
     * @param line the line
     * @return gff the entry
     */
    public static GFFEntry parse(String line) throws IOException{
        line = line.trim();
        String[] split = line.split("\t");
        if(split.length != 9) throw new IllegalArgumentException("The line does not contain 9 tab separated fields !");

        String chromosome = split[0];
        String source = split[1];
        String feature = split[2];
        long start = parseLong(split[3], "start");
        long end = parseLong(split[4], "end");
        short score = parseStrand(split[5], "score");
        char strand = parseChar(split[6], "strand");
        char codingFrame = parseChar(split[7], "frame");
        String group = split[8];
        return new GFFEntry(chromosome, source, feature, start, end, score, strand, codingFrame, group);
    }

    /**
     * Helper to parse long and catch exception. Takes
     * the field name to create a better exception message
     *
     * @param value the value
     * @param name the field name
     * @return result parsed result
     */
    private static long parseLong(String value, String name){
        try {
            return Long.parseLong(value);
        } catch (NumberFormatException e) {
            throw new IllegalArgumentException("Error while parsing " + name +" gtf field from value: " + value);
        }
    }

    /**
     * Helper to parse strand, get the '.' and catch exception. Takes
     * the field name to create a better exception message
     *
     * @param value the value
     * @param name the field name
     * @return result parsed result
     */
    private static short parseStrand(String value, String name){
        if(value.equals(".")) return -1;
        try {
            return Short.parseShort(value);
        } catch (NumberFormatException e) {
            throw new IllegalArgumentException("Error while parsing " + name +" gtf field from value: " + value);
        }
    }

    /**
     * Helper to parse char and catch exception. Takes
     * the field name to create a better exception message
     *
     * @param value the value
     * @param name the field name
     * @return result parsed result
     */
    private static char parseChar(String value, String name){
        try {
            return value.charAt(0);
        } catch (Exception e) {
            throw new IllegalArgumentException("Error while parsing " + name +" gtf field from value: " + value);
        }
    }

}
