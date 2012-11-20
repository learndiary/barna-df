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
