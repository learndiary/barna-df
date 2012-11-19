package barna.flux.capacitor.diffexp;

import java.util.Map;

/**
 * GFF extension that wraps around a transcript entry. The constructor
 * tries to read readCount and rpkm value from the gtf attributes
 */
class Transcript extends GFFEntry{

    /**
     * The transcript id
     */
    private final String id;

    /**
     * The transcripts read count
     */
    private double readCount;

    /**
     * The transcripts RPKM value
     */
    private double rpkm;

    /**
     * Create a new Transcript. This assumes an existing "transcript_id" attribute
     * and throws an Exception if that attribute is not found!
     *
     * @param chromosome  the chromosome
     * @param source      the source name
     * @param start       start position 1-based
     * @param end         end position 1-based inclusive
     * @param score       the score (between 0 and 1000)
     * @param strand      the strand (+/- or . for unknown)
     * @param codingFrame the coding frame 0,1,2 or . for unknown
     * @param group       optional group entry, can be null or or empty string
     */
    public Transcript(final String chromosome, final String source, final String feature, final long start, final long end, final short score, final char strand, final char codingFrame, final String group) {
        super(chromosome, source, feature, start, end, score, strand, codingFrame, group);
        this.id = getAttributes().get("transcript_id");
        if(id == null || id.isEmpty()) throw new IllegalArgumentException("No valid transcript id specified : " + id);

        String reads = getAttributes().get("reads");
        String rpkm = getAttributes().get("rpkm");
        if(reads != null){
            this.readCount = Double.parseDouble(reads);
        }
        if(rpkm != null){
            this.rpkm = Double.parseDouble(rpkm);
        }
    }

    /**
     * Get the transcript id
     *
     * @return id the transcript id
     */
    public String getId() {
        return id;
    }

    /**
     * Get the read count for this transcript
     *
     * @return readCount the read count
     */
    public double getReadCount() {
        return readCount;
    }

    /**
     * Get the rpkm value for this transcript
     *
     * @return rpkm the rpkm value
     */
    public double getRpkm() {
        return rpkm;
    }
}
