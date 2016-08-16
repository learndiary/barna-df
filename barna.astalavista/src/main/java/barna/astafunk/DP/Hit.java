package barna.astafunk.DP;

/**
 * @version 2.0
 * @autor Vitor Coelho
 * @since 2014-07-17
 */

/**
 * Hit describes a predicted hit of a reference transcript (A transcript representing a set of transcript
 * with the same exon/intron structure between a source and sink).
 */
public class Hit {

    /**
     * Chromossome ID (see GTF format).
     */
    private String chr;

    /**
     * Gene ID. It is a string with concatenated IDs of genes/transcripts with AS.
     */
    private String geneID;

    /**
     * Name of family/domain/repeat/etc of Pfam.
     */
    private String family;

    /**
     * Accession number ACC of Pfam.
     */
    private String acc;

    /**
     * Description of profile HMM.
     */
    private String desc;

    /**
     * Bit score of the hit
     */
    private double score;

    /**
     * Start position of alignment in the sequence.
     */
    private int startAlignment;

    /**
     * End position of alignment in the sequence.
     */
    private int endAlignment;

    /**
     * Genomic start position of alignment.
     */

    private int genomicStartAlignment;

    /**
     * Genomic end position of alignment.
     */

    private int genomicEndAlignment;

    /**
     * Event source.
     */

    private int source;

    /**
     * Event sink.
     */
    private int sink;

    /**
     * Start position of alignment in model.
     */
    private int startModel;

    /**
     * End position of alignment in model.
     */
    private int endModel;

    /**
     * Number of states of the profile HMM or length.
     */
    private int lengthModel;

    /**
     * ID of transcript from where the sequence was translated to search a profile HMM.
     */
    private String transcriptID;

    /**
     * Group ID: a string with all transcripts ID that have the same exon/intron structure
     * between a source and sink.
     */

    private String groupID;

    /**
     * Constructor of Hit object with genomic information.
     * @param geneID GTF annotation of gene id
     * @param family Family or domain name
     * @param acc Accession number from Pfam
     * @param desc Description of family/domain
     * @param score Bit score.
     * @param startAlignment Start position of hit.
     * @param endAlignment End position of hit.
     * @param genomicStartAlignment Genomic start position of hit.
     * @param genomicEndAlignment Genomic end position of hit.
     * @param startModel Start position of model.
     * @param endModel End position of model.
     * @param lengthModel Number of states of the model.
     * @param transcriptID Transcript ID.
     */
    public Hit(String geneID, String family, String acc, String desc, double score,
               int startAlignment, int endAlignment, int genomicStartAlignment,
               int genomicEndAlignment, int startModel, int endModel, int lengthModel,
               String transcriptID) {
        this.geneID = geneID;
        this.family = family;
        this.acc = acc;
        this.desc = desc;
        this.score = score;
        this.startAlignment = startAlignment;
        this.endAlignment = endAlignment;
        this.genomicStartAlignment = genomicStartAlignment;
        this.genomicEndAlignment = genomicEndAlignment;
        this.startModel = startModel;
        this.endModel = endModel;
        this.lengthModel = lengthModel;
        this.transcriptID = transcriptID;
    }

    /**
     * Constructor of Hit object without genomic information.
     * @param family Family or domain name
     * @param acc Accession number from Pfam
     * @param desc Description of family/domain
     * @param score Bit score
     * @param startAlignment Start position of alignment
     * @param endAlignment End position of alignment
     * @param startModel Start position of model
     * @param endModel End position of model
     * @param lengthModel Number of states of the model
     * @param transcriptID Transcript ID
     */
    public Hit(String family, String acc, String desc, double score,
               int startAlignment, int endAlignment, int startModel, int endModel,
               int lengthModel, String transcriptID) {

        this.family = family;
        this.acc = acc;
        this.desc = desc;
        this.score = score;
        this.startAlignment = startAlignment;
        this.endAlignment = endAlignment;
        this.startModel = startModel;
        this.endModel = endModel;
        this.lengthModel = lengthModel;
        this.transcriptID = transcriptID;
    }

    /**
     * @return Family/Domain name.
     */
    public String getFamily() {
        return family;
    }

    /**
     * @return ACC number of profile HMM.
     */
    public String getAcc() {
        return acc;
    }

    /**
     * @param acc Accession number of profile HMM to set.
     */
    public void setAcc(String acc) {
        this.acc = acc;
    }

    /**
     * @return Description of profile HMM.
     */
    public String getDesc() {
        return desc;
    }

    /**
     * @return Bit score of the hit.
     */
    public double getScore() {
        return score;
    }

    /**
     * @param score Bit score of the hit to set.
     */
    public void setScore(double score) {
        this.score = score;
    }

    /**.
     * @return Start position of the alignment (Dynamic Matrix column number or
     *                       first alignment position in the sequence) to set.
     */
    public int getStartAlignment() {
        return startAlignment;
    }

    /**
     * @return End position of alignment (Dynamic Matrix column number
     * or last alignment position in the sequence).
     */
    public int getEndAlignment() {
        return endAlignment;
    }

    /**
     * @return Start position of the alignment in genomic coordinates.
     */
    public int getGenomicStartAlignment() {
        return genomicStartAlignment;
    }

    /**
     * @return End position of the alignment in genomic coordinates.
     */
    public int getGenomicEndAlignment() {
        return genomicEndAlignment;
    }

    /**
     * @return Number of the first state of alignment in HMM (first Dynamic Matrix row number in the
     * alignment).
     */
    public int getStartModel() {
        return startModel;
    }

    /**
     * @return Number of the last  state of alignment in the HMM ( last Dynamic Matrix row number
     * in the alignment)
     */
    public int getEndModel() {
        return endModel;
    }

    /**
     * Method to return length model (number of states).
     * @return Length model (number of states).
     */
    public int getLengthModel() {
        return lengthModel;
    }

    /**
     * Method to set length model (number of states).
     * @param lengthModel Length model (number of states).
     */
    public void setLengthModel(int lengthModel) {
        this.lengthModel = lengthModel;
    }

    /**
     * Method to return transcript ID.
     * @return Transcript ID.
     */
    public String getTranscriptID() {
        return transcriptID;
    }

    /**
     * @return Source genomic coordinate of hit.
     */
    public int getSource() {
        return source;
    }

    public void setSource(int source) {
        this.source = source;
    }

    public void setSink(int sink) {
        this.sink = sink;
    }

    /**
     * Method to set genomic information to the hit. These information are obtained from annotation.
     * @param chr Chromosome ID.
     * @param geneID Gene ID.
     * @param genomicStartAlignment Start position of the alignment in genomic coordinates.
     * @param genomicEndAlignment End position of the alignment in genomic coordinates.
     */
    public void setGenomicInfo(String chr, String geneID,
                               int genomicStartAlignment, int genomicEndAlignment){
        this.chr = chr;
        this.geneID =geneID;
        this.genomicStartAlignment = genomicStartAlignment;
        this.genomicEndAlignment = genomicEndAlignment;

    }

    /**
     * @return Chromosome ID.
     */
    public String getChr() {
        return chr;
    }

    /**
     * @param chr Chromosome ID to set.
     */
    public void setChr(String chr) {
        this.chr = chr;
    }

    /**
     * @param groupID Group ID (a string with all transcripts ID that have the
     *                same exon/intron structure) to set.
     */
    public void setGroupID(String groupID) {
        this.groupID = groupID;
    }

    /**
     * @return Default output of the hit.
     */
    @Override
    public String toString() {
        String output = chr + '\t' +
                geneID + '\t' +
                family + '\t' +
                acc + '\t' +
                desc + '\t' +
                score + '\t' +
                startAlignment + '\t' +
                endAlignment + '\t' +
                genomicStartAlignment + '\t' +
                genomicEndAlignment + '\t' +
                source + '\t'+
                sink + '\t' +
                startModel + '\t' +
                endModel + '\t' +
                lengthModel + '\t' +
                transcriptID + '\t' +
                groupID;
        return output;
    }
}