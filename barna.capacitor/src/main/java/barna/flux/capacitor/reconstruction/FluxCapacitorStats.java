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

package barna.flux.capacitor.reconstruction;

/**
 * Capacitor stats wrapper
 *
 * @author Thasso Griebel (thasso.griebel@gmail.com)
 */
public class FluxCapacitorStats {
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
    private long lociSingle;
    private long lociExp;
    private long txExp;
    private long eventsExp;
    private long mappingsSingle;
    private long mappingsSinglePairs;
    private long mappingsSinglePairsMapped;
    private long mappingsTotal;
    private long mappingsMapped;
    private long mappingsPairsNa;
    private long mappingsPairsWo;
    private long mappingsNotSens;

    public long getLociSingle() {
        return lociSingle;
    }

    public void setLociSingle(long lociSingle) {
        this.lociSingle = lociSingle;
    }

    public long getLociExp() {
        return lociExp;
    }

    public void setLociExp(long lociExp) {
        this.lociExp = lociExp;
    }

    public long getTxExp() {
        return txExp;
    }

    public void setTxExp(long txExp) {
        this.txExp = txExp;
    }

    public long getEventsExp() {
        return eventsExp;
    }

    public void setEventsExp(long eventsExp) {
        this.eventsExp = eventsExp;
    }

    public long getMappingsSingle() {
        return mappingsSingle;
    }

    public void setMappingsSingle(long mappingsSingle) {
        this.mappingsSingle = mappingsSingle;
    }

    public long getMappingsSinglePairs() {
        return mappingsSinglePairs;
    }

    public void setMappingsSinglePairs(long mappingsSinglePairs) {
        this.mappingsSinglePairs = mappingsSinglePairs;
    }

    public long getMappingsSinglePairsMapped() {
        return mappingsSinglePairsMapped;
    }

    public void setMappingsSinglePairsMapped(long mappingsSinglePairsMapped) {
        this.mappingsSinglePairsMapped = mappingsSinglePairsMapped;
    }

    public long getMappingsTotal() {
        return mappingsTotal;
    }

    public void setMappingsTotal(long mappingsTotal) {
        this.mappingsTotal = mappingsTotal;
    }

    public long getMappingsMapped() {
        return mappingsMapped;
    }

    public void setMappingsMapped(long mappingsMapped) {
        this.mappingsMapped = mappingsMapped;
    }

    public long getMappingsPairsNa() {
        return mappingsPairsNa;
    }

    public void setMappingsPairsNa(long mappingsPairsNa) {
        this.mappingsPairsNa = mappingsPairsNa;
    }

    public long getMappingsPairsWo() {
        return mappingsPairsWo;
    }

    public void setMappingsPairsWo(long mappingsPairsWa) {
        this.mappingsPairsWo = mappingsPairsWa;
    }

    public long getMappingsNotSens() {
        return mappingsNotSens;
    }

    public void setMappingsNotSens(long mappingsNotSens) {
        this.mappingsNotSens = mappingsNotSens;
    }

    /**
     * Add the numbers from the given stats to this instance.
     *
     * @param other the other stats
     */
    public void add(FluxCapacitorStats other){
        if(other == null) return;
        this.lociSingle                 += other.lociSingle               ;
        this.lociExp                    += other.lociExp                  ;
        this.txExp                      += other.txExp                    ;
        this.eventsExp                  += other.eventsExp                ;
        this.mappingsSingle             += other.mappingsSingle           ;
        this.mappingsSinglePairs        += other.mappingsSinglePairs      ;
        this.mappingsSinglePairsMapped  += other.mappingsSinglePairsMapped;
        this.mappingsTotal              += other.mappingsTotal            ;
        this.mappingsMapped             += other.mappingsMapped           ;
        this.mappingsPairsNa            += other.mappingsPairsNa          ;
        this.mappingsPairsWo            += other.mappingsPairsWo          ;
        this.mappingsNotSens            += other.mappingsNotSens          ;
    }
}
