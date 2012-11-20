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
 * Quantification Model is a container created from a base annotation used for
 * the quantification. The model can be queried for transcripts and genes.
 * <p>
 * An instance can be created using the {@link #read(java.io.File)} method.
 * </p>
 *
 * @author Thasso Griebel <thasso.griebel@gmail.com>
 */
class QuantificationModel {

    /**
     * Map from transcript identifier to
     * the transcript
     */
    private Map<String, GFFEntry> transcripts;

    /**
     * Map from transcript identifier to
     * the transcript
     */
    private Map<String, GFFEntry> genes;


    /**
     * Private constructor
     */
    private QuantificationModel() {
        transcripts = new HashMap<String, GFFEntry>();
        genes = new HashMap<String, GFFEntry>();
    }

    /**
     * Add a transcript to the quantification
     *
     * @param entry the gff entry
     */
    void addEntry(GFFEntry entry) {
        if(entry.getFeature().equals("transcript")){
            String transcript_id = entry.getAttributes().get("transcript_id");
            if(transcripts.containsKey(transcript_id)){
                throw new RuntimeException("Duplicated transcript id " + transcript_id);
            }else{
                transcripts.put(transcript_id, entry);
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
