package barna.flux.capacitor.reconstruction;

import barna.flux.capacitor.profile.MappingStats;
import org.junit.Test;

import static junit.framework.Assert.assertEquals;

public class MappingStatsTest {

    @Test
    public void testAppendingToStats() throws Exception {
        MappingStats other = new MappingStats();
        other.setEventsExp(1);
        other.setLociExp(1);
        other.setLociSingle(1);
        other.setMappingsMapped(1);
        other.setMappingsNotSens(1);
        other.setMappingsPairsNa(1);
        other.setMappingsPairsWo(1);
        other.setMappingsSingle(1);
        other.setMappingsSinglePairs(1);
        other.setMappingsSinglePairsMapped(1);
        other.setMappingsTotal(1);
        other.setTxExp(1);

        MappingStats stats = new MappingStats();
        stats.add(other);
        assertEquals(1, stats.getEventsExp());
        assertEquals(1, stats.getLociExp());
        assertEquals(1, stats.getLociSingle());
        assertEquals(1, stats.getMappingsMapped());
        assertEquals(1, stats.getMappingsNotSens());
        assertEquals(1, stats.getMappingsPairsNa());
        assertEquals(1, stats.getMappingsPairsWo());
        assertEquals(1, stats.getMappingsSingle());
        assertEquals(1, stats.getMappingsSinglePairs());
        assertEquals(1, stats.getMappingsSinglePairsMapped());
        assertEquals(1, stats.getMappingsTotal());
        assertEquals(1, stats.getTxExp());

        stats.add(other);
        assertEquals(2, stats.getEventsExp());
        assertEquals(2, stats.getLociExp());
        assertEquals(2, stats.getLociSingle());
        assertEquals(2, stats.getMappingsMapped());
        assertEquals(2, stats.getMappingsNotSens());
        assertEquals(2, stats.getMappingsPairsNa());
        assertEquals(2, stats.getMappingsPairsWo());
        assertEquals(2, stats.getMappingsSingle());
        assertEquals(2, stats.getMappingsSinglePairs());
        assertEquals(2, stats.getMappingsSinglePairsMapped());
        assertEquals(2, stats.getMappingsTotal());
        assertEquals(2, stats.getTxExp());
    }
}
