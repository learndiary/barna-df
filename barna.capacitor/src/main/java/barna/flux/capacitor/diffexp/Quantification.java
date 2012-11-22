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

import java.io.*;
import java.util.HashMap;
import java.util.Map;
import java.util.zip.GZIPInputStream;

/**
 * Quantification container for differential expression. You can create a new instance
 * using the {@link #read(java.io.File)} method. The quantification should be created
 * from a gff file produced by the capacitor.
 *
 * @author Thasso Griebel <thasso.griebel@gmail.com>
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
    private Map<String, QuantificationEntry> transcripts;


    /**
     * Private constructor
     */
    private Quantification() {
        transcripts = new HashMap<String, QuantificationEntry>();
    }

    /**
     * Add a transcript to the quantification
     *
     * @param transcript the transcript
     */
    void addTranscript(QuantificationEntry transcript) {
        if(transcripts.containsKey(transcript.getId())){
            throw new RuntimeException("Duplicated transcript id " + transcript.getId());
        }else{
            transcripts.put(transcript.getKey(), transcript);
            this.totalReads += transcript.getReadCount();
            this.totalRpkm  += transcript.getRpkm();
        }
    }

    /**
     * Iterate over the transcripts
     *
     * @return transcripts iterate over the transcripts
     */
    public Iterable<QuantificationEntry> transcripts(){
        return transcripts.values();
    }

    /**
     * Get transcript with the given id or null
     * @param id the transcript id
     * @return transcript the transcript
     */
    public QuantificationEntry getTranscript(String id) {
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
                    QuantificationEntry transcript = new QuantificationEntry("transcript_id", entry);
                    q.addTranscript(transcript);
                }
            }
        }
        return q;
    }

}
