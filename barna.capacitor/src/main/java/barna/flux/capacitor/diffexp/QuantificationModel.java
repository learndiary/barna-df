package barna.flux.capacitor.diffexp;

import java.io.*;
import java.util.HashMap;
import java.util.Map;
import java.util.zip.GZIPInputStream;

/**
 * Quantification Model is a container created from a base annotation used for
 * the quantification. The model can be queried for transcripts and genes.
 * <p>
 * An instance can be created using the {@link #read(java.io.File)} method.
 * </p>
 */
class QuantificationModel {

    /**
     * Map from transcript identifier to
     * the transcript
     */
    private Map<String, Transcript> transcripts;

    /**
     * Map from transcript identifier to
     * the transcript
     */
    private Map<String, GFFEntry> genes;


    /**
     * Private constructor
     */
    private QuantificationModel() {
        transcripts = new HashMap<String, Transcript>();
        genes = new HashMap<String, GFFEntry>();
    }

    /**
     * Add a transcript to the quantification
     *
     * @param entry the gff entry
     */
    void addEntry(GFFEntry entry) {
        if(entry.getFeature().equals("transcript")){
            Transcript transcript = (Transcript) entry;
            if(transcripts.containsKey(transcript.getId())){
                throw new RuntimeException("Duplicated transcript id " + transcript.getId());
            }else{
                transcripts.put(transcript.getId(), transcript);
            }
        }else if (entry.getFeature().equals("gene")){
            String gene_id = entry.getAttributes().get("gene_id");
            if(genes.containsKey(gene_id)){
                throw new RuntimeException("Duplicated gene id " + gene_id);
            }else{
                genes.put(gene_id, entry);
            }
        }
    }

    /**
     * Create a quantification from a given file
     *
     * @param file the input file
     * @return quantification quantification created from the given file
     * @throws java.io.IOException in case the file can not be read
     */
    public static QuantificationModel read(File file) throws IOException{
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
     * @throws java.io.IOException in case the file can not be read
     */
    private static QuantificationModel read(Reader reader) throws IOException {
        BufferedReader bb = new BufferedReader(reader);
        String line = null;
        QuantificationModel q = new QuantificationModel();
        while((line = bb.readLine()) != null){
            if(!line.isEmpty()){
                q.addEntry(GFFEntry.parse(line));
            }
        }
        return q;
    }

}
