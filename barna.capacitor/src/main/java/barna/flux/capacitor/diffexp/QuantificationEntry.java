package barna.flux.capacitor.diffexp;

import java.util.Map;

/**
 * GFF extension that wraps around a transcript entry. The constructor
 * tries to read readCount and rpkm value from the gtf attributes
 *
 * @author Thasso Griebel <thasso.griebel@gmail.com>
 */
class QuantificationEntry {
    /**
     * Grouping name
     */
    private final String name;

    /**
     * Unique key to identify this entry
     */
    private final String key;
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
     *
     */
    public QuantificationEntry(String id_key, GFFEntry gffEntry) {
        super();
        this.id = gffEntry.getAttributes().get(id_key);
        if(id == null || id.isEmpty()) throw new IllegalArgumentException("No valid transcript id specified : " + id);
        String reads = gffEntry.getAttributes().get("reads");
        String rpkm = gffEntry.getAttributes().get("rpkm");
        this.key = this.id + ":::" + gffEntry.getAttributes().get("locus_id");
        this.name = gffEntry.getChromosome();
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

    /**
     * Get the key that identifies this entry uniquely
     *
     * @return key the key that identifies this entry uniquely
     */
    public String getKey() {
        return key;
    }

    /**
     * Get the grouping name
     *
     * @return name the grouping name (i.e. the chromosome)
     */
    public String getName() {
        return name;
    }
}
