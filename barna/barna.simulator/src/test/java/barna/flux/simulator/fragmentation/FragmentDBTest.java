package barna.flux.simulator.fragmentation;

import barna.commons.ByteArrayCharSequence;
import barna.io.FileHelper;
import org.junit.Test;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.List;

import static org.junit.Assert.assertEquals;

public class FragmentDBTest {

    @Test
    public void testThatIndexCreationWorks() throws Exception {
        File libFile = new File(getClass().getResource("/fragmentdb_test.txt").getFile());
        long lines = FileHelper.countLines(libFile);

        FragmentDB fragmentDB = new FragmentDB(libFile);
        fragmentDB.createIndex();
        assertEquals(lines, fragmentDB.getNumberOfLines());
        assertEquals(9060, fragmentDB.getNumberOfLines());
        assertEquals(3751, fragmentDB.getNumberOfEntries());
        assertEquals(14072, fragmentDB.getNumberOfFragments());
    }

    @Test
    public void testThatTheIDIteratorHasTheRightSize() throws Exception {
        File libFile = new File(getClass().getResource("/fragmentdb_test.txt").getFile());

        FragmentDB fragmentDB = new FragmentDB(libFile);
        fragmentDB.createIndex();
        long c = 0;
        for (String s : fragmentDB.fragmentIDs()) {
            c++;
        }

        assertEquals(3751, c);

    }

    @Test
    public void testThatAllIteratedEntriesHaveTheSameIDAndTheNumebrOfFragmentsmatches() throws Exception {
        File libFile = new File(getClass().getResource("/fragmentdb_test.txt").getFile());

        FragmentDB fragmentDB = new FragmentDB(libFile);
        fragmentDB.createIndex();
        long entryCounter = 0;
        long fragmentCounter = 0;
        for (String id : fragmentDB.fragmentIDs()) {
            entryCounter++;
            Iterable<ByteArrayCharSequence> entries = fragmentDB.getEntries(id);
            for (ByteArrayCharSequence fragment : entries) {
                assertEquals(id, fragment.getToken(2).toString());
                fragmentCounter++;
            }
        }

        assertEquals(3751, entryCounter);
        assertEquals(9060, fragmentCounter);

    }
    @Test
    public void test_that_the_fragmentsets_are_complete() throws Exception {
        File libFile = new File(getClass().getResource("/fragmentdb_test.txt").getFile());

        FragmentDB fragmentDB = new FragmentDB(libFile);
        fragmentDB.createIndex();

        BufferedReader libFileReader = null;
        try {
            // read the sorted file and put it in a zip form
            libFileReader = new BufferedReader(new FileReader(libFile));
            ByteArrayCharSequence cs = new ByteArrayCharSequence(300);
            ByteArrayCharSequence nextID = null;
            ByteArrayCharSequence currentID = null;
            String line = null;
            List<String> frags = new ArrayList<String>();
            while ((line = libFileReader.readLine()) != null) {
                cs.set(line);
                cs.resetFind();
                nextID = cs.getToken(2);

                // initialize current
                if(currentID == null){
                    currentID = nextID.cloneCurrentSeq();
                }

                if (!nextID.equals(currentID)) {
                    Iterable<ByteArrayCharSequence> entries = fragmentDB.getEntries(currentID.toString());
                    int c = 0;
                    for (ByteArrayCharSequence entry : entries) {
                        assertEquals(entry.toString(), frags.get(c++));
                    }
                    assertEquals(c, frags.size());
                    frags.clear();
                    currentID = nextID.cloneCurrentSeq();
                }
                frags.add(cs.toString());

            }
        } finally {
            libFileReader.close();
        }

    }
}
