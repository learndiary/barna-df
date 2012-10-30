/*
 * Copyright (c) 2010, Micha Sammeth
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * The names of its contributors may be not used to endorse or promote
 *       products derived from this software without specific prior written
 *       permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL MICHA SAMMETH BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

package barna.flux.capacitor.profile;

import barna.commons.log.Log;
import com.google.gson.Gson;
import com.google.gson.GsonBuilder;

import java.io.*;
import java.nio.channels.FileChannel;
import java.nio.channels.FileLock;

/**
 * Capacitor stats wrapper
 *
 * @author Thasso Griebel (thasso.griebel@gmail.com)
 */
public class MappingStats {
    /*
    public final static String
            JSON_LOCI_SINGLE= "LOCI_NO-AS",
            JSON_MAPPINGS_SL= "MAPPINGS_NO-AS",
            JSON_MAPPINGS_SL_PAIRS= "MAPPING-PAIRS_NO-AS",
            JSON_MAPPINGS_SL_PAIRS_MAPPED= "MAPPING-PAIRS_NO-AS_MAPPED",

    JSON_MAPPINGS_TOTAL= "MAPPINGS_TOTAL",
            JSON_MAPPINGS_MAPPED= "MAPPINGS_MAPPED",
            JSON_MAPPINGS_PAIRS_NA= "MAPPINGS_PAIRS-WO",
            JSON_MAPPINGS_PAIRS_WO= "MAPPINGS_PAIRS-WO",
            JSON_MAPPINGS_NOTSENSE= "MAPPINGS_NOTSENSE",
            JSON_LOCI_EXP= "LOCI_EXP",
            JSON_TX_EXP= "TX_EXP",
            JSON_EVENTS_EXP= "EVENTS_EXP";
    */

    //Total
    private long readsTotal;
    private long mappingsTotal;

    //Learning
    private long readsSingleTxLoci;
    private long mappingsSingleTxLoci;
    private long mappingPairsSingleTxLoci;
    private long mappingsSingleTxLociNoAnn;
    private long mappingPairs;
    private long mappingsWrongStrand;

    //Learning & Annotation Mapping
    private long pairsWrongOrientation;

    //Annotation Mapping
    private long singleTxLoci;
    private long readsLoci;
    private long mappingsMapped;
    private long mappingPairsNoTx;

    //Deconvolution
    private long lociExp;
    private long txsExp;
    private long eventsExp;
    private long lociUnsolved;


    public long getSingleTxLoci() {
        return singleTxLoci;
    }

    public void setSingleTxLoci(long singleTxLoci) {
        this.singleTxLoci = singleTxLoci;
    }

    public void incrSingleTxLoci() {
        ++this.singleTxLoci;
    }

    public long getLociExp() {
        return this.lociExp;
    }

    public void setLociExp(long lociExp) {
        this.lociExp = lociExp;
    }

    public void incrLociExp() {
        ++this.lociExp;
    }

    public long getTxsExp() {
        return txsExp;
    }

    public void setTxsExp(long txsExp) {
        this.txsExp = txsExp;
    }

    public void incrTxsExp() {
        ++this.txsExp;
    }

    public long getEventsExp() {
        return eventsExp;
    }

    public void setEventsExp(long eventsExp) {
        this.eventsExp = eventsExp;
    }

    public void incrEventsExp() {
        ++this.eventsExp;
    }

    public long getReadsSingleTxLoci() {
        return readsSingleTxLoci;
    }

    public void setReadsSingleTxLoci(long readsSingleTxLoci) {
        this.readsSingleTxLoci = readsSingleTxLoci;
    }

    public void incrReadsSingleTxLoci() {
        ++this.readsSingleTxLoci;
    }

    public long getMappingsSingleTxLoci() {
        return mappingsSingleTxLoci;
    }

    public void setMappingsSingleTxLoci(long mappingsSingleTxLoci) {
        this.mappingsSingleTxLoci = mappingsSingleTxLoci;
    }

    public void incrMappingsSingleTxLoci() {
        ++this.mappingsSingleTxLoci;
    }

    public long getMappingPairs() {
        return mappingPairs;
    }

    public void setMappingPairs(long mappingPairs) {
        this.mappingPairs = mappingPairs;
    }

    public void incrMappingPairs() {
        ++this.mappingPairs;
    }

    public long getMappingsTotal() {
        return mappingsTotal;
    }

    public void setMappingsTotal(long mappingsTotal) {
        this.mappingsTotal = mappingsTotal;
    }

    public void incrMappingsTotal() {
        ++this.mappingsTotal;
    }

    public long getMappingsMapped() {
        return mappingsMapped;
    }

    public void setMappingsMapped(long mappingsMapped) {
        this.mappingsMapped = mappingsMapped;
    }

    public void incrMappingsMapped() {
        ++this.mappingsMapped;
    }

    public long getMappingPairsNoTx() {
        return mappingPairsNoTx;
    }

    public void setMappingPairsNoTx(long mappingPairsNoTx) {
        this.mappingPairsNoTx = mappingPairsNoTx;
    }

    public void incrMappingPairsNoTx() {
        ++this.mappingPairsNoTx;
    }

    public long getPairsWrongOrientation() {
        return pairsWrongOrientation;
    }

    public void setPairsWrongOrientation(long mappingsPairsWa) {
        this.pairsWrongOrientation = mappingsPairsWa;
    }

    public void incrPairsWrongOrientation() {
        this.pairsWrongOrientation+=2;
    }

    public long getMappingsWrongStrand() {
        return mappingsWrongStrand;
    }

    public void setMappingsWrongStrand(long mappingsWrongStrand) {
        this.mappingsWrongStrand = mappingsWrongStrand;
    }

    public void incrMappingsWrongStrand() {
        ++this.mappingsWrongStrand;
    }

    public long getReadsTotal() {
        return readsTotal;
    }

    public void setReadsTotal(long readsTotal) {
        this.readsTotal = readsTotal;
    }

    public void incrReadsTotal() {
        ++this.readsTotal;
    }

    public long getMappingsSingleTxLociNoAnn() {
        return mappingsSingleTxLociNoAnn;
    }

    public void setMappingsSingleTxLociNoAnn(long mappingsSingleTxLociNoAnn) {
        this.mappingsSingleTxLociNoAnn = mappingsSingleTxLociNoAnn;
    }

    public void incrMappingsSingleTxLociNoAnn() {
        ++this.mappingsSingleTxLociNoAnn;
    }

    public long getReadsLoci() {
        return readsLoci;
    }

    public void setReadsLoci(long readsLoci) {
        this.readsLoci = readsLoci;
    }

    public void incrReadsLoci() {
        ++this.readsLoci;
    }

    public long getLociUnsolved() {
        return lociUnsolved;
    }

    public void setLociUnsolved(long lociUnsolved) {
        this.lociUnsolved = lociUnsolved;
    }

    public void incrLociUnsolved() {
        ++this.lociUnsolved;
    }

    public long getMappingPairsSingleTxLoci() {
        return mappingPairsSingleTxLoci;
    }

    public void setMappingPairsSingleTxLoci(long mappingPairsSingleTxLoci) {
        this.mappingPairsSingleTxLoci = mappingPairsSingleTxLoci;
    }

    public void incrMappingPairsSingleTxLoci() {
        this.mappingPairsSingleTxLoci+=2;
    }

    /**
     * Add the numbers from the given stats to this instance.
     *
     * @param other the other stats
     */
    public void add(MappingStats other){
        if(other == null) return;
        this.singleTxLoci               += other.singleTxLoci;
        this.lociExp                    += other.lociExp;
        this.txsExp                     += other.txsExp;
        this.eventsExp                  += other.eventsExp;
        this.readsSingleTxLoci          += other.readsSingleTxLoci;
        this.mappingsSingleTxLoci       += other.mappingsSingleTxLoci;
        this.mappingPairs               += other.mappingPairs;
        this.mappingsTotal              += other.mappingsTotal;
        this.mappingsMapped             += other.mappingsMapped;
        this.mappingPairsNoTx           += other.mappingPairsNoTx;
        this.pairsWrongOrientation      += other.pairsWrongOrientation;
        this.mappingsWrongStrand        += other.mappingsWrongStrand;
        this.readsLoci                  += other.readsLoci;
        this.readsTotal                 += other.readsTotal;
        this.mappingsSingleTxLociNoAnn  += other.mappingsSingleTxLociNoAnn;
        this.lociUnsolved               += other.lociUnsolved;
    }

    /**
     * Write the stats to file in JSON format
     *
     * @param statsFile the file to write to
     * @param append whether append the current stats to an existing file
     * @throws Exception
     */
    public void writeStats(File statsFile, boolean append) throws Exception {// BARNA-103 : write stats to file
        if (statsFile != null) {
            MappingStats statsToWrite = this;
            BufferedWriter writer = null;
            BufferedReader reader = null;


            File lockFile = new File(statsFile.getAbsolutePath() + ".lock");
    //            if(!lockFile.exists()) lockFile.createNewFile();
            FileChannel channel = new RandomAccessFile(lockFile, "rw").getChannel();
            FileLock lock = channel.lock();


            try {
                Gson gson = new GsonBuilder().setPrettyPrinting().create();
                if (statsFile.exists() && append) {
                    // read stats file and append
                    reader = new BufferedReader(new FileReader(statsFile));
                    MappingStats existingStats = gson.fromJson(reader, MappingStats.class);
                    reader.close();
                    existingStats.add(this);
                    statsToWrite = existingStats;
                }
                Log.info((append ? "Appending stats to " : "Writing stats to ") + statsFile.getAbsolutePath());
                writer = new BufferedWriter(new FileWriter(statsFile));
                gson.toJson(statsToWrite, writer);
                writer.close();
            } catch (Exception e) {
                Log.error("Unable to " + (append ? "append stats to " : "write stats to ") + statsFile.getAbsolutePath() + " : " + e.getMessage(), e);
            } finally {
                if (reader != null) reader.close();
                if (writer != null) writer.close();
                // release the lock
                try {
                    lock.release();
                } catch (IOException e) {
                    Log.error("Unable to release lock");
                }
                channel.close();
            }
        }
    }

    public void readStats(File statsFile) {
        final String MSG_READING_STATS = "reading mapping stats";

        Log.progressStart(MSG_READING_STATS);
        String status = "OK";

        try {
            BufferedReader buffy = new BufferedReader(new FileReader(statsFile));
            BiasProfile profile = new BiasProfile();
            Gson gson = new GsonBuilder().serializeSpecialFloatingPointValues().create();
            this.add(gson.fromJson(buffy,MappingStats.class));
        } catch (Exception e) {
            Log.error("Cannot read stats from file: " + statsFile.getAbsolutePath());
            status = "KO";
        } finally {
            Log.progressFinish(status, false);
        }
    }
}
