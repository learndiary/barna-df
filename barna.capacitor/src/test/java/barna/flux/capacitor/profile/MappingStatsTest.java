package barna.flux.capacitor.profile;

import barna.flux.capacitor.reconstruction.FluxCapacitor;
import org.junit.Test;

import static junit.framework.Assert.assertEquals;

public class MappingStatsTest {

    static {
        FluxCapacitor.DEBUG= false;}

    @Test
    public void testAppendingToStats() throws Exception {
        MappingStats other = new MappingStats();
        other.setEventsExp(1);
        other.setLociExp(1);
        other.setSingleTxLoci(1);
        other.setMappingsMapped(1);
        other.setMappingsWrongStrand(1);
        other.setMappingPairsNoTx(1);
        other.setPairsWrongOrientation(1);
        other.setReadsSingleTxLoci(1);
        other.setMappingsSingleTxLoci(1);
//        other.setMappingPairs(1);
        other.setMappingsTotal(1);
        other.setTxsExp(1);

        MappingStats stats = new MappingStats();
        stats.add(other);
        assertEquals(1, stats.getEventsExp());
        assertEquals(1, stats.getLociExp());
        assertEquals(1, stats.getSingleTxLoci());
        assertEquals(1, stats.getMappingsMapped());
        assertEquals(1, stats.getMappingsWrongStrand());
        assertEquals(1, stats.getMappingPairsNoTx());
        assertEquals(1, stats.getPairsWrongOrientation());
        assertEquals(1, stats.getReadsSingleTxLoci());
        assertEquals(1, stats.getMappingsSingleTxLoci());
//        assertEquals(1, stats.getMappingPairs());
        assertEquals(1, stats.getMappingsTotal());
        assertEquals(1, stats.getTxsExp());

        stats.add(other);
        assertEquals(2, stats.getEventsExp());
        assertEquals(2, stats.getLociExp());
        assertEquals(2, stats.getSingleTxLoci());
        assertEquals(2, stats.getMappingsMapped());
        assertEquals(2, stats.getMappingsWrongStrand());
        assertEquals(2, stats.getMappingPairsNoTx());
        assertEquals(2, stats.getPairsWrongOrientation());
        assertEquals(2, stats.getReadsSingleTxLoci());
        assertEquals(2, stats.getMappingsSingleTxLoci());
//        assertEquals(2, stats.getMappingPairs());
        assertEquals(2, stats.getMappingsTotal());
        assertEquals(2, stats.getTxsExp());
    }
}
