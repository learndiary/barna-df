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
import java.util.*;
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
     * Feature map that maps from feature type (transcript, gene) to
     * a map with feature_id mapping the the feature
     */
    private Map<String, Map<String, Feature>> features;

    /**
     * Private constructor
     */
    private QuantificationModel() {
        features = new HashMap<String, Map<String, Feature>>();
    }

    /**
     * Add a a feature from a GFF entry. If an id_attribute is
     * specified, the entries attributes map is queried for the value and
     * it is used as id. If no id_attribute is specified, a unique one is
     * created as {@literal <chr>:<start>:<end>:<strand>}.
     *
     * @param entry the gff entry
     * @param id_attribute the key of the attribute used as id
     */
    void addGFFEntry(GFFEntry entry, String id_attribute) {
        String type = entry.getFeature();
        String id = null;
        if(id_attribute == null){
            id = entry.getChromosome() + ":"+entry.getStart() +":" + entry.getEnd() +":" +entry.getEnd();
        }else{
            id = entry.getAttributes().get(id_attribute);
        }
        getFeatures(type).put(id, new Feature(id, entry.getAttributes()));
    }

    /**
     * Returns the map that stores the feature entries for
     * the specified feature type
     *
     * @param feature the feature type
     * @return features map of features
     */
    Map<String, Feature> getFeatures(String feature){
        Map<String, Feature> ff = features.get(feature);
        if(ff == null){
            ff = new HashMap<String, Feature>();
            features.put(feature, ff);
        }
        return ff;
    }

    /**
     * Returns the value for the attribute from the feature with given type and id
     *
     * @param type the feature type
     * @param id the id
     * @param attribute the attribute key
     * @return value the value or null
     */
    String get(String type, String id, String attribute){
        Map<String, Feature> ff = features.get(type);
        if(ff == null) return null;
        Feature feature = ff.get(id);
        if(feature == null) return null;
        return feature.get(attribute);
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
            if(file.getName().toLowerCase().endsWith("gtf") || file.getName().toLowerCase().endsWith("gff") ||
                    file.getName().toLowerCase().endsWith("gtf.gz") || file.getName().toLowerCase().endsWith("gff.gz")){
                return readGFF(reader);
            }else{
                // reads tab separated file
                return readTable(reader);
            }
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
    public static QuantificationModel readGFF(Reader reader) throws IOException {
        BufferedReader bb = new BufferedReader(reader);
        String line = null;
        QuantificationModel q = new QuantificationModel();
        while((line = bb.readLine()) != null){
            if(!line.isEmpty()){
                GFFEntry gff = GFFEntry.parse(line);
                if(gff.getFeature().equals("transcript")){ // todo : remove this when we open up for genes
                    String id_attribute = null;
                    if(gff.getFeature().equals("transcript")){
                        id_attribute = "transcript_id";
                    }else if (gff.getFeature().equals("gene")){
                        id_attribute = "gene_id";
                    }
                    q.addGFFEntry(gff, id_attribute);
                }
            }
        }
        return q;
    }

    /**
     * Create a quantification from a tab separated table file.
     * <p>
     *  Lines starting with {@literal #} are ignored. The first non comment
     *  line in the file must be a header line that defines the attribute names. First
     *  column is always the ID, second column is always the feature type (gene, transcript etc)
     * </p>
     *
     * @param reader the input
     * @return quantification quantification created from the given file
     * @throws java.io.IOException in case the file can not be read
     */
    public static QuantificationModel readTable(Reader reader) throws IOException {
        BufferedReader bb = new BufferedReader(reader);
        String line = null;
        QuantificationModel q = new QuantificationModel();
        List<String> header = null;
        long lineCount = 0;
        while((line = bb.readLine()) != null){
            lineCount++;
            line = line.trim();
            if(!line.isEmpty() && line.charAt(0) != '#'){
                String[] split = line.split("\t");
                if(split.length < 2){
                    throw new RuntimeException("Error in line "+lineCount + ": You have to specify at least 2 columns (ID and type)");
                }

                if(header == null){
                    header = new ArrayList<String>();
                    for (String s : split) {
                        header.add(s.trim());
                    }
                    continue;
                }

                String id = split[0].trim();
                String type = split[1].trim();
                if(type.isEmpty()){
                    throw new RuntimeException("Error in line "+lineCount + ": No type specified!");
                }
                if(id.isEmpty()){
                    throw new RuntimeException("Error in line "+lineCount + ": No id specified!");
                }
                HashMap<String, String> attrs = new HashMap<String, String>();
                for (int i = 2; i < split.length; i++) {
                    if(header.size() <= i){
                        throw new RuntimeException("Error in line "+lineCount + ": No enough fields specified in your header!");
                    }
                    String value = split[i].trim();
                    if(!value.isEmpty()){
                        attrs.put(header.get(i), value);
                    }
                }
                q.getFeatures(type).put(id, new Feature(id, attrs));
            }
        }
        return q;
    }

}
