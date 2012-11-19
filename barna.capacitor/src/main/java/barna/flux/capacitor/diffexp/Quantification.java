package barna.flux.capacitor.diffexp;

import java.io.*;
import java.util.HashMap;
import java.util.Map;
import java.util.zip.GZIPInputStream;

/**
 * Quantification container for differential expression. You can create a new instance
 * using the {@link #read(java.io.File)} method.
 *
 */
class Quantification {
    /**
     * Absolute read count
     */
    private double totalReads;

    /**
     * Total RPKM count
     */

    private double totalRpkm;

    /**
     * Map from transcript identifier to
     * the transcript
     */
    private Map<String, Transcript> transcripts;


    /**
     * Private constructor
     */
    private Quantification() {
        transcripts = new HashMap<String, Transcript>();
    }

    /**
     * Add a transcript to the quantification
     *
     * @param transcript the transcript
     */
    void addTranscript(Transcript transcript) {
        if(transcripts.containsKey(transcript.getId())){
            throw new RuntimeException("Duplicated transcript id " + transcript.getId());
        }else{
            transcripts.put(transcript.getId(), transcript);
            this.totalReads += transcript.getReadCount();
            this.totalRpkm  += transcript.getRpkm();
        }
    }

    /**
     * Iterate over the transcripts
     *
     * @return transcripts iterate over the transcripts
     */
    public Iterable<Transcript> transcripts(){
        return transcripts.values();
    }

    /**
     * Get transcript with the given id or null
     * @param id the transcript id
     * @return transcript the transcript
     */
    public Transcript getTranscript(String id) {
        if(transcripts == null) return null;
        return transcripts.get(id);
    }

    /**
     * Get the total read count
     *
     * @return totalReads the total read count
     */
    public double getTotalReads() {
        return totalReads;
    }

    /**
     * Create a quantification from a given file
     *
     * @param file the input file
     * @return quantification quantification created from the given file
     * @throws IOException in case the file can not be read
     */
    public static Quantification read(File file) throws IOException{
        if(file == null) throw new NullPointerException();
        if(!file.exists()) throw new FileNotFoundException("Quantification input " + file.getAbsolutePath() + " not found");
        Reader reader = null;
        if(file.getName().endsWith(".gz")){
            reader = new InputStreamReader(new GZIPInputStream(new FileInputStream(file)));
        }else{
            reader = new FileReader(file);
        }
        try {
            return read(reader);
        } catch (IOException e) {
            try {reader.close();} catch (IOException ignore) {}
            throw e;
        }
    }

    /**
     * Create a quantification from a given input reader
     *
     * @param reader the input
     * @return quantification quantification created from the given file
     * @throws IOException in case the file can not be read
     */
    private static Quantification read(Reader reader) throws IOException {
        BufferedReader bb = new BufferedReader(reader);
        String line = null;
        Quantification q = new Quantification();
        while((line = bb.readLine()) != null){
            if(!line.isEmpty()){
                GFFEntry entry = GFFEntry.parse(line);
                if(entry.getFeature().equals("transcript")){
                    Transcript transcript = (Transcript) entry;
                    q.addTranscript(transcript);
                }
            }
        }
        return q;
    }

}
