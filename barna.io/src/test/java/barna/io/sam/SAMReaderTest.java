package barna.io.sam;

import barna.commons.Execute;
import barna.io.MSIterator;
import barna.io.rna.UniversalReadDescriptor;
import barna.model.Mapping;
import barna.model.sam.SAMMapping;
import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.File;
import java.util.Iterator;

import static org.junit.Assert.assertEquals;

public class SAMReaderTest {
	
	private File testfile = new File(getClass().getResource("/test.bam").getFile());
    private File testMultiMaps = new File(getClass().getResource("/single_multimap.bam").getFile());

	@BeforeClass
	public static void setUp() {
        Execute.initialize(4);
	}

	@AfterClass
	public static void tearDown() throws Exception {
		Execute.shutdown();
	}

	@Test
	public void testRead() {
        SAMReader reader = new SAMReader(testfile, true, true);
        MSIterator iter = reader.read("chr22", 24030323, 24041363);
        SAMMapping mapping = (SAMMapping)iter.next();

        int c = 1;
        while (iter.hasNext()) {
            ++c;
            mapping = (SAMMapping)iter.next();
        }


        assertEquals(181, c);
        assertEquals("chr22", mapping.getChromosome());
        assertEquals(24034554, mapping.getStart());
        assertEquals(24034629, mapping.getEnd());
	}

    @Test
    public void testPrimaryMaps() {
        SAMReader reader = new SAMReader(testMultiMaps, true, false, true, false, false);
        MSIterator iter = reader.read("chr21", 34924516, 34924516+1000);
        SAMMapping mapping;

        int c = 0;
        while (iter.hasNext()) {
            ++c;
            mapping = (SAMMapping)iter.next();
        }

        assertEquals(c, 2);
    }

    @Test
    public void testReadMultiMapsSortedMates() {
        SAMReader reader = new SAMReader(testMultiMaps, true, false, false, true, false);
        MSIterator<Mapping> iter = reader.read("chr21", 34924516, 34924516+1000);
        SAMMapping mapping;

        UniversalReadDescriptor desc = new UniversalReadDescriptor();
        desc.init(UniversalReadDescriptor.DESCRIPTORID_PAIRED);

        int c = 0;
        while (iter.hasNext()) {
            ++c;
            mapping = (SAMMapping)iter.next();
            if (mapping.getName().endsWith("1")) {
                int d = 0;
                Iterator<Mapping> mates = iter.getMates(mapping, desc);
                Mapping m;
                while (mates.hasNext()) {
                    m = mates.next();
                    ++d;
                }
                assertEquals(1,d);
            }
        }

        assertEquals(3,c);
    }

    @Test
    public void testReadMultiMapsSortedPrimary() {
        SAMReader reader = new SAMReader(testMultiMaps, true, false, true, false, false);
        MSIterator<Mapping> iter = reader.read("chr21", 34924516, 34924516+1000);
        SAMMapping mapping;

        UniversalReadDescriptor desc = new UniversalReadDescriptor();
        desc.init(UniversalReadDescriptor.DESCRIPTORID_PAIRED);

        int c = 0;
        while (iter.hasNext()) {
            ++c;
            mapping = (SAMMapping)iter.next();
            if (mapping.getName().endsWith("1")) {
                int d = 0;
                Iterator<Mapping> mates = iter.getMates(mapping, desc);
                Mapping m;
                while (mates.hasNext()) {
                    m = mates.next();
                    ++d;
                }
                assertEquals(1,d);
            }
        }

        assertEquals(2,c);
    }
    @Test
    public void testReadMultiMapsSortedUnique() {
        SAMReader reader = new SAMReader(testMultiMaps, true, false, false, false, true);
        MSIterator<Mapping> iter = reader.read("chr21", 34924516, 34924516+1000);
        SAMMapping mapping;

        UniversalReadDescriptor desc = new UniversalReadDescriptor();
        desc.init(UniversalReadDescriptor.DESCRIPTORID_PAIRED);

        int c = 0;
        while (iter.hasNext()) {
            ++c;
            mapping = (SAMMapping)iter.next();
            if (mapping.getName().endsWith("1")) {
                int d = 0;
                Iterator<Mapping> mates = iter.getMates(mapping, desc);
                Mapping m;
                while (mates.hasNext()) {
                    m = mates.next();
                    ++d;
                }
                assertEquals(1,d);
            }
        }

        assertEquals(2,c);
    }

    @Test
    public void testReadMultiMapsSortedMatesUnique() {
        SAMReader reader = new SAMReader(testMultiMaps, true, false, false, true, true);
        MSIterator<Mapping> iter = reader.read("chr21", 34924516, 34924516+1000);
        SAMMapping mapping;

        UniversalReadDescriptor desc = new UniversalReadDescriptor();
        desc.init(UniversalReadDescriptor.DESCRIPTORID_PAIRED);

        int c = 0;
        while (iter.hasNext()) {
            ++c;
            mapping = (SAMMapping)iter.next();
            if (mapping.getName().endsWith("1")) {
                int d = 0;
                Iterator<Mapping> mates = iter.getMates(mapping, desc);
                Mapping m;
                while (mates.hasNext()) {
                    m = mates.next();
                    ++d;
                }
                assertEquals(0,d);
            }
        }

        assertEquals(2,c);
    }


    @Test
    public void testReadMultiMapsMates() {
        SAMReader reader = new SAMReader(testMultiMaps, true, true, false, true, false);
        MSIterator<Mapping> iter = reader.read("chr21", 34924516, 34924516+1000);
        SAMMapping mapping;

        UniversalReadDescriptor desc = new UniversalReadDescriptor();
        desc.init(UniversalReadDescriptor.DESCRIPTORID_PAIRED);

        int c = 0;
        while (iter.hasNext()) {
            ++c;
            mapping = (SAMMapping)iter.next();
            if (mapping.getName().endsWith("1")) {
                int d = 0;
                Iterator<Mapping> mates = iter.getMates(mapping, desc);
                Mapping m;
                while (mates.hasNext()) {
                    m = mates.next();
                    ++d;
                }
                assertEquals(1,d);
            }
        }

        assertEquals(3,c);  // 3 pairs when considering sam pairing
    }

    @Test
    public void testReadMultiMapsPrimary() {
        SAMReader reader = new SAMReader(testMultiMaps, true, true, true, false, false);
        MSIterator<Mapping> iter = reader.read("chr21", 34924516, 34924516+1000);
        SAMMapping mapping;

        UniversalReadDescriptor desc = new UniversalReadDescriptor();
        desc.init(UniversalReadDescriptor.DESCRIPTORID_PAIRED);

        int c = 0;
        while (iter.hasNext()) {
            ++c;
            mapping = (SAMMapping)iter.next();
            if (mapping.getName().endsWith("1")) {
                int d = 0;
                Iterator<Mapping> mates = iter.getMates(mapping, desc);
                Mapping m;
                while (mates.hasNext()) {
                    m = mates.next();
                    ++d;
                }
                assertEquals(1,d);
            }
        }

        assertEquals(2,c);
    }

    @Test
    public void testReadMultiMapsUnique() {
        SAMReader reader = new SAMReader(testMultiMaps, true, true, false, false, true);
        MSIterator<Mapping> iter = reader.read("chr21", 34924516, 34924516+1000);
        SAMMapping mapping;

        UniversalReadDescriptor desc = new UniversalReadDescriptor();
        desc.init(UniversalReadDescriptor.DESCRIPTORID_PAIRED);

        int c = 0;
        while (iter.hasNext()) {
            ++c;
            mapping = (SAMMapping)iter.next();
            if (mapping.getName().endsWith("1")) {
                int d = 0;
                Iterator<Mapping> mates = iter.getMates(mapping, desc);
                Mapping m;
                while (mates.hasNext()) {
                    m = mates.next();
                    ++d;
                }
                assertEquals(1,d);
            }
        }

        assertEquals(2,c);
    }

    @Test
    public void testReadMultiMapsMatesUnique() {
        SAMReader reader = new SAMReader(testMultiMaps, true, true, false, true, true);
        MSIterator<Mapping> iter = reader.read("chr21", 34924516, 34924516+1000);
        SAMMapping mapping;

        UniversalReadDescriptor desc = new UniversalReadDescriptor();
        desc.init(UniversalReadDescriptor.DESCRIPTORID_PAIRED);

        int c = 0;
        while (iter.hasNext()) {
            ++c;
            mapping = (SAMMapping)iter.next();
            if (mapping.getName().endsWith("1")) {
                int d = 0;
                Iterator<Mapping> mates = iter.getMates(mapping, desc);
                Mapping m;
                while (mates.hasNext()) {
                    m = mates.next();
                    ++d;
                }
                assertEquals(0,d);
            }
        }

        assertEquals(2,c);
    }
}
