package barna.io.sam;

import barna.commons.Execute;
import barna.io.MSIterator;
import barna.model.Mapping;
import barna.model.sam.SAMMapping;
import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.*;
import java.util.Iterator;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

public class SAMReaderTest {

    // file 'test.bam' is sorted by position
	private File testfile = new File(getClass().getResource("/test.bam").getFile());

    // file 'testID.bam' is sorted by position
    private File testIDfile = new File(getClass().getResource("/testID.bam").getFile());

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

    /**
     * Sorts a lexicographically presorted file by position, checks reference ID and alignment position
     * for correct sorting.
     * @throws Exception
     */
    @Test
    public void testSort() throws Exception {
        SAMReader reader = new SAMReader(testIDfile, true, true);
        PipedInputStream pipi= new PipedInputStream();
        PipedOutputStream pipo= new PipedOutputStream(pipi);
        final BufferedReader buffy= new BufferedReader(new InputStreamReader(pipi));
        Thread tt= new Thread() {
            @Override
            public void run() {
                try {
                    for(String line= null, last= null; (line= buffy.readLine())!= null; ) {
                        if (line.startsWith("@"))
                            continue;
                        if (last!= null) {
                            String[] s1= last.split("\\s");
                            String[] s2= line.split("\\s");
                            assertTrue(s1[2].compareTo(s2[2])<= 0);
                            assertTrue(Integer.parseInt(s1[3])<= Integer.parseInt(s2[3]));
                        }
                        last= line;
                    }
                } catch (IOException e) {
                    e.printStackTrace();
                    fail();
                }
            }
        };

        tt.start();
        reader.sort(pipo);
        tt.join();

    }

    /**
     * Sorts alphabetically by read IDs a file that has been presorted by position.
     * Lexicographical ordering of read IDs is checked.
     * @throws Exception
     */
    @Test
    public void testSortByReadID() throws Exception {
        SAMReader reader = new SAMReader(testfile, true, true);
        PipedInputStream pipi= new PipedInputStream();
        PipedOutputStream pipo= new PipedOutputStream(pipi);
        final BufferedReader buffy= new BufferedReader(new InputStreamReader(pipi));
        Thread tt= new Thread() {
            @Override
            public void run() {
                try {
                    for(String line= null, last= null; (line= buffy.readLine())!= null; ) {
                        if (line.startsWith("@"))
                            continue;
                        if (last!= null) {
                            String[] s1= last.split("\\s");
                            String[] s2= line.split("\\s");
                            assertTrue(s1[0].compareTo(s2[0])<= 0);
                        }
                        last= line;
                    }
                } catch (IOException e) {
                    e.printStackTrace();
                    fail();
                }
            }
        };

        tt.start();
        reader.sortByReadID(pipo);
        tt.join();

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

        int c = 0;
        while (iter.hasNext()) {
            ++c;
            mapping = (SAMMapping)iter.next();
            if (mapping.getName(true).endsWith("1")) {
                int d = 0;
                Iterator<Mapping> mates = iter.getMates(mapping);
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

        int c = 0;
        while (iter.hasNext()) {
            ++c;
            mapping = (SAMMapping)iter.next();
            if (mapping.getName(true).endsWith("1")) {
                int d = 0;
                Iterator<Mapping> mates = iter.getMates(mapping);
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

        int c = 0;
        while (iter.hasNext()) {
            ++c;
            mapping = (SAMMapping)iter.next();
            if (mapping.getName(true).endsWith("1")) {
                int d = 0;
                Iterator<Mapping> mates = iter.getMates(mapping);
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

        int c = 0;
        while (iter.hasNext()) {
            ++c;
            mapping = (SAMMapping)iter.next();
            if (mapping.getName(true).endsWith("1")) {
                int d = 0;
                Iterator<Mapping> mates = iter.getMates(mapping);
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

        int c = 0;
        while (iter.hasNext()) {
            ++c;
            mapping = (SAMMapping)iter.next();
            if (mapping.getName(true).endsWith("1")) {
                int d = 0;
                Iterator<Mapping> mates = iter.getMates(mapping);
                Mapping m;
                while (mates.hasNext()) {
                    m = mates.next();
                    ++d;
                }
                assertEquals(1,d);
            }
        }

        assertEquals(6,c);  // 3 pairs when considering sam pairing, but mates iterated redundantly
    }

    @Test
    public void testReadMultiMapsPrimary() {
        SAMReader reader = new SAMReader(testMultiMaps, true, true, true, false, false);
        MSIterator<Mapping> iter = reader.read("chr21", 34924516, 34924516+1000);
        SAMMapping mapping;

        int c = 0;
        while (iter.hasNext()) {
            ++c;
            mapping = (SAMMapping)iter.next();
            if (mapping.getName(true).endsWith("1")) {
                int d = 0;
                Iterator<Mapping> mates = iter.getMates(mapping);
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

        int c = 0;
        while (iter.hasNext()) {
            ++c;
            mapping = (SAMMapping)iter.next();
            if (mapping.getName(true).endsWith("1")) {
                int d = 0;
                Iterator<Mapping> mates = iter.getMates(mapping);
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

        int c = 0;
        while (iter.hasNext()) {
            ++c;
            mapping = (SAMMapping)iter.next();
            if (mapping.getName(true).endsWith("1")) {
                int d = 0;
                Iterator<Mapping> mates = iter.getMates(mapping);
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
