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
    private double mappingsSingleTxLoci;
    private double mappingPairsSingleTxLoci;
    private long mappingsSingleTxLociNoAnn;
    private long mappingsWrongStrand;
    private int readLenMin;
    private int readLenMax;


    //Learning & Annotation Mapping
    private long pairsWrongOrientation;

    //Annotation Mapping
    private long singleTxLoci;
    private long readsLoci;
    private double mappingsMapped;
    private long mappingPairsNoTx;

    // new

    public double getCtrHits() {
        return ctrHits;
    }

    public void setCtrHits(double ctrHits) {
        this.ctrHits = ctrHits;
    }

    public long getCtrHitsNone() {
        return ctrHitsNone;
    }

    public void setCtrHitsNone(long ctrHitsNone) {
        this.ctrHitsNone = ctrHitsNone;
    }

    public long getCtrHitsMultiLocus() {
        return ctrHitsMultiLocus;
    }

    public void setCtrHitsMultiLocus(long ctrHitsMultiLocus) {
        this.ctrHitsMultiLocus = ctrHitsMultiLocus;
    }

    public long getCtrHitsMultiGenome() {
        return ctrHitsMultiGenome;
    }

    public void setCtrHitsMultiGenome(long ctrHitsMultiGenome) {
        this.ctrHitsMultiGenome = ctrHitsMultiGenome;
    }

    public long getCtrHitsMultiLocusAndGenome() {
        return ctrHitsMultiLocusAndGenome;
    }

    public void setCtrHitsMultiLocusAndGenome(long ctrHitsMultiLocusAndGenome) {
        this.ctrHitsMultiLocusAndGenome = ctrHitsMultiLocusAndGenome;
    }

    public long getCtrHitsNoneMultiGenome() {
        return ctrHitsNoneMultiGenome;
    }

    public void setCtrHitsNoneMultiGenome(long ctrHitsNoneMultiGenome) {
        this.ctrHitsNoneMultiGenome = ctrHitsNoneMultiGenome;
    }

    /**
     * The number of reads (single end) or read pairs (paired-end) that map to the locus, with at least one valid
     * mapping. The number contains the count of <code>ctrHitsMultiLocus</code>, but not any other counters as it
     * focuses on the reads of which we are really sure they stem from this locus. Hits are possibly weighted.
     */
    private double ctrHits= 0;

    /**
     * The number of reads (single end) or read pairs (paired-end) without any valid mapping within the locus. The
     * redundant distribution of annotation mapping errors is stored in the corresponding mapping counter.
     */
    private long ctrHitsNone= 0;

    /**
     * The number of reads (single end) or read pairs (paired-end) that map multiple times only within the locus,
     * with at least one valid mapping.
     */
    private long ctrHitsMultiLocus= 0;

    /**
     * The number of reads (single end) or read pairs (paired-end) that map multiple times within the genome but
     * not within the locus and that have a valid mapping in the locus.
     */
    private long ctrHitsMultiGenome= 0;

    /**
     * The number of reads (single end) or read pairs (paired-end) that map multiple times within the locus and
     * also in the genome.
     */
    private long ctrHitsMultiLocusAndGenome= 0;

    /**
     * The number of reads (single end) or read pairs (paired-end) that do not map within the locus, but multiple
     * times within the genome.
     */
    private long ctrHitsNoneMultiGenome= 0;


    //Deconvolution
    private long lociExp;
    private long txsExp;
    private long eventsExp;
    private long lociUnsolved;
    
    //Annotation
    private long txs;
    private long loci;
    private long events;


    public long getSingleTxLoci() {
        return singleTxLoci;
    }

    public void setSingleTxLoci(long singleTxLoci) {
        this.singleTxLoci = singleTxLoci;
    }

    public void incrSingleTxLoci(long value) {
        this.singleTxLoci+=value;
    }

    public long getLoci() {
        return this.loci;
    }

    public void setLoci(long loci) {
        this.loci = loci;
    }

    public void incrLoci(long value) {
        this.loci+=value;
    }

    public long getTxs() {
        return this.txs;
    }

    public void setTxs(long txs) {
        this.txs = txs;
    }

    public void incrTxs(long value) {
        this.txs+=value;
    }

    public long getEvents() {
        return this.events;
    }

    public void setEvents(long events) {
        this.events = events;
    }

    public void incrEvents(long value) {
        this.events+=value;
    }

    public long getLociExp() {
        return this.lociExp;
    }

    public void setLociExp(long lociExp) {
        this.lociExp = lociExp;
    }

    public void incrLociExp(long value) {
        this.lociExp+=value;
    }

    public long getTxsExp() {
        return txsExp;
    }

    public void setTxsExp(long txsExp) {
        this.txsExp = txsExp;
    }

    public void incrTxsExp(long value) {
        this.txsExp+=value;
    }

    public long getEventsExp() {
        return eventsExp;
    }

    public void setEventsExp(long eventsExp) {
        this.eventsExp = eventsExp;
    }

    public void incrEventsExp(long value) {
        this.eventsExp+=value;
    }

    public long getReadsSingleTxLoci() {
        return readsSingleTxLoci;
    }

    public void setReadsSingleTxLoci(long readsSingleTxLoci) {
        this.readsSingleTxLoci = readsSingleTxLoci;
    }

    public void incrReadsSingleTxLoci(long value) {
        this.readsSingleTxLoci+=value;
    }

    public long getMappingsSingleTxLoci() {
        if (!Double.isInfinite(mappingsSingleTxLoci) && mappingsSingleTxLoci == Math.floor(mappingsSingleTxLoci))
            return (long)mappingsSingleTxLoci;
        return Math.round(mappingsSingleTxLoci);
    }

    public void setMappingsSingleTxLoci(long mappingsSingleTxLoci) {
        this.mappingsSingleTxLoci = mappingsSingleTxLoci;
    }

    public void incrMappingsSingleTxLoci(double value) {
        this.mappingsSingleTxLoci+=value;
    }

    public long getMappingsTotal() {
        return mappingsTotal;
    }

    public void setMappingsTotal(long mappingsTotal) {
        this.mappingsTotal = mappingsTotal;
    }

    public void incrMappingsTotal(long value) {
        this.mappingsTotal+=value;
    }

    public long getMappingsMapped() {
        if (!Double.isInfinite(mappingsMapped) && mappingsMapped == Math.floor(mappingsMapped))
            return (long)mappingsMapped;
        return Math.round(mappingsMapped);
    }

    public void setMappingsMapped(double mappingsMapped) {
        this.mappingsMapped = mappingsMapped;
    }

    public void incrMappingsMapped(long value) {
        this.mappingsMapped+=value;
    }

    public long getMappingPairsNoTx() {
        return mappingPairsNoTx;
    }

    /**
     * @deprecated not supported anymore, see <code>mappingsNotMapped</code>
     * @param mappingPairsNoTx
     */
    public void setMappingPairsNoTx(long mappingPairsNoTx) {
        this.mappingPairsNoTx = mappingPairsNoTx;
    }

    public void incrMappingPairsNoTx(long value) {
        this.mappingPairsNoTx+=value;
    }

    public long getPairsWrongOrientation() {
        return pairsWrongOrientation;
    }

    /**
     * @deprecated not supported anymore, see <code>mappingsWrongStrand</code>
     * @param mappingsPairsWa
     */
    public void setPairsWrongOrientation(long mappingsPairsWa) {
        this.pairsWrongOrientation = mappingsPairsWa;
    }

    public void incrPairsWrongOrientation(long value) {
        this.pairsWrongOrientation+=value;
    }

    public long getMappingsWrongStrand() {
        return mappingsWrongStrand;
    }

    public void setMappingsWrongStrand(long mappingsWrongStrand) {
        this.mappingsWrongStrand = mappingsWrongStrand;
    }

    public void incrMappingsWrongStrand(long value) {
        this.mappingsWrongStrand+=value;
    }

    public long getReadsTotal() {
        return readsTotal;
    }

    public void setReadsTotal(long readsTotal) {
        this.readsTotal = readsTotal;
    }

    public void incrReadsTotal(long value) {
        this.readsTotal+=value;
    }

    public long getMappingsSingleTxLociNoAnn() {
        return mappingsSingleTxLociNoAnn;
    }

    public void setMappingsSingleTxLociNoAnn(long mappingsSingleTxLociNoAnn) {
        this.mappingsSingleTxLociNoAnn = mappingsSingleTxLociNoAnn;
    }

    public void incrMappingsSingleTxLociNoAnn(long value) {
        this.mappingsSingleTxLociNoAnn+=value;
    }

    public long getReadsLoci() {
        return readsLoci;
    }

    public void setReadsLoci(long readsLoci) {
        this.readsLoci = readsLoci;
    }

    public void incrReadsLoci(long value) {
        this.readsLoci+=value;
    }

    public long getLociUnsolved() {
        return lociUnsolved;
    }

    public void setLociUnsolved(long lociUnsolved) {
        this.lociUnsolved = lociUnsolved;
    }

    public void incrLociUnsolved(long value) {
        this.lociUnsolved+=value;
    }

    public long getMappingPairsSingleTxLoci() {
        if (!Double.isInfinite(mappingPairsSingleTxLoci) && mappingPairsSingleTxLoci == Math.floor(mappingPairsSingleTxLoci))
            return (long)mappingPairsSingleTxLoci;
        return Math.round(mappingPairsSingleTxLoci);
    }

    public void incrMappingPairsSingleTxLoci(double value) {
        this.mappingPairsSingleTxLoci+=value;
    }

    public int getReadLenMin() {
        return readLenMin;
    }

    public void setReadLenMin(int readLenMin) {
        this.readLenMin = readLenMin;
    }

    public int getReadLenMax() {
        return readLenMax;
    }

    public void setReadLenMax(int readLenMax) {
        this.readLenMax = readLenMax;
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
        this.mappingsTotal              += other.mappingsTotal;
        this.mappingsMapped             += other.mappingsMapped;
        this.mappingPairsNoTx           += other.mappingPairsNoTx;
        this.pairsWrongOrientation      += other.pairsWrongOrientation;
        this.mappingsWrongStrand        += other.mappingsWrongStrand;
        this.readsLoci                  += other.readsLoci;
        this.readsTotal                 += other.readsTotal;
        this.mappingsSingleTxLociNoAnn  += other.mappingsSingleTxLociNoAnn;
        this.lociUnsolved               += other.lociUnsolved;
        this.mappingPairsSingleTxLoci   += other.mappingPairsSingleTxLoci;
    }

    /**
     * Add the numbers from the given stats referring to a single locus to this instance.
     *
     * @param other the other stats
     */
    public void addLocus(MappingStats other){
        if(other == null) return;

        // TODO reflective solution

        this.singleTxLoci               += other.singleTxLoci;
        this.lociExp                    += other.lociExp;
        this.txsExp                     += other.txsExp;
        this.eventsExp                  += other.eventsExp;
        this.mappingsMapped             += other.mappingsMapped;
        this.mappingPairsNoTx           += other.mappingPairsNoTx;
        this.pairsWrongOrientation      += other.pairsWrongOrientation;
        this.readsLoci                  += other.readsLoci;
        this.lociUnsolved               += other.lociUnsolved;
        this.txs                        += other.txs;
        this.loci                       += other.loci;
        this.events                     += other.events;

        this.ctrHits                    += other.ctrHits;
        this.ctrHitsMultiGenome         += other.ctrHitsMultiGenome;
        this.ctrHitsMultiLocus          += other.ctrHitsMultiLocus;
        this.ctrHitsMultiLocusAndGenome += other.ctrHitsMultiLocusAndGenome;
        this.ctrHitsNone                += other.ctrHitsNone;
        this.ctrHitsNoneMultiGenome     += other.ctrHitsNoneMultiGenome;
    }

    /**
     * Write the stats to file in JSON format
     *
     * @param statsFile the file to write to
     * @throws Exception
     */
    public void writeStats(File statsFile) throws Exception {// BARNA-103 : write stats to file
        if (statsFile != null) {
            MappingStats statsToWrite = this;
            BufferedWriter writer = null;
            try {
                Gson gson = new GsonBuilder().setPrettyPrinting().create();
                Log.info("Writing stats to " + statsFile.getAbsolutePath());
                writer = new BufferedWriter(new FileWriter(statsFile));
                gson.toJson(statsToWrite, writer);
                writer.close();
            } catch (Exception e) {
                Log.error("Unable to write stats to " + statsFile.getAbsolutePath() + " : " + e.getMessage(), e);
            } finally {
                if (writer != null) writer.close();
            }
        }
    }

    public void readStats(File statsFile) {
        final String MSG_READING_STATS = "reading mapping stats";

        Log.progressStart(MSG_READING_STATS);
        String status = "OK";

        try {
            BufferedReader buffy = new BufferedReader(new FileReader(statsFile));
            Gson gson = new GsonBuilder().create();
            MappingStats other = gson.fromJson(buffy,MappingStats.class);
            if(other == null) return;
            this.singleTxLoci               = other.singleTxLoci;
            this.lociExp                    = other.lociExp;
            this.txsExp                     = other.txsExp;
            this.eventsExp                  = other.eventsExp;
            this.readsSingleTxLoci          = other.readsSingleTxLoci;
            this.mappingsSingleTxLoci       = other.mappingsSingleTxLoci;
            this.mappingsTotal              = other.mappingsTotal;
            this.mappingsMapped             = other.mappingsMapped;
            this.mappingPairsNoTx           = other.mappingPairsNoTx;
            this.pairsWrongOrientation      = other.pairsWrongOrientation;
            this.mappingsWrongStrand        = other.mappingsWrongStrand;
            this.readsLoci                  = other.readsLoci;
            this.readsTotal                 = other.readsTotal;
            this.mappingsSingleTxLociNoAnn  = other.mappingsSingleTxLociNoAnn;
            this.lociUnsolved               = other.lociUnsolved;
            this.mappingPairsSingleTxLoci   = other.mappingPairsSingleTxLoci;
            this.readLenMin                 = other.readLenMin;
            this.readLenMax                 = other.readLenMax;
        } catch (Exception e) {
            Log.error("Cannot read stats from file: " + statsFile.getAbsolutePath());
            status = "KO";
        } finally {
            Log.progressFinish(status, false);
        }
    }

    @Override
    public boolean equals(Object other) {
        if (!(other instanceof MappingStats))
            return false;
        MappingStats otherStats = (MappingStats)other;
        if (this.singleTxLoci != otherStats.singleTxLoci)
            return false;
        if (this.eventsExp != otherStats.eventsExp)
            return false;
        if (this.readsSingleTxLoci != otherStats.readsSingleTxLoci)
            return false;
        if (this.mappingsSingleTxLoci != otherStats.mappingsSingleTxLoci)
            return false;
        if (this.mappingsTotal != otherStats.mappingsTotal)
            return false;
        if (this.txsExp != otherStats.txsExp)
            return false;
        if (this.mappingsMapped != otherStats.mappingsMapped)
            return false;
        if (this.mappingPairsNoTx != otherStats.mappingPairsNoTx)
            return false;
        if (this.pairsWrongOrientation != otherStats.pairsWrongOrientation)
            return false;
        if (this.mappingsWrongStrand != otherStats.mappingsWrongStrand)
            return false;
        if (this.readsLoci != otherStats.readsLoci)
            return false;
        if (this.readsTotal != otherStats.readsTotal)
            return false;
        if (this.mappingsSingleTxLociNoAnn != otherStats.mappingsSingleTxLociNoAnn)
            return false;
        if (this.lociUnsolved != otherStats.lociUnsolved)
            return false;
        if (this.mappingPairsSingleTxLoci != otherStats.mappingPairsSingleTxLoci)
            return false;
        if (this.readLenMin != otherStats.readLenMin)
            return false;
        if (this.readLenMax != otherStats.readLenMax)
            return false;
        if (this.lociExp != otherStats.lociExp)
            return false;
        return true;
    }

    public void reset() {
        //Annotation Mapping
        readsLoci = 0;
        mappingsMapped = 0;
        mappingPairsNoTx = 0;
        pairsWrongOrientation = 0;
        singleTxLoci = 0;

        //After deconvolution
        loci = 0;
        lociExp = 0;
        lociUnsolved = 0;
        txs = 0;
        txsExp = 0;
        events = 0;
        eventsExp = 0;
    }
}
