package fbi.commons.tools;

import fbi.commons.Execute;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.util.Comparator;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Future;

import static junit.framework.Assert.assertEquals;
import static junit.framework.Assert.fail;

/**
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class SorterTest {
    private static final String SIMPLE_INPUT =
            "C\tY\t2\n" +
            "B\tZ\t1\n" +
            "A\tX\t3\n";

    private static final String SIMPLE_INPUT_WITH_NEWLINES =
            "\n" +
            "C\tY\t2\n" +
            "\n" +
            "B\tZ\t1\n" +
            "\n" +
            "A\tX\t3\n";


    @Before
    public void setUp() throws Exception {
        Execute.initialize(4);
    }

    @After
    public void tearDown() throws Exception {
        Execute.shutdown();
    }

    @Test
    public void testSmallSort(){
        ByteArrayInputStream in = new ByteArrayInputStream(SIMPLE_INPUT.getBytes());
        ByteArrayOutputStream out = new ByteArrayOutputStream();
        try {
            Sorter.create(in, out, true).field(0, false).sort();

            String outString = new String(out.toByteArray());
            assertEquals(
                    "A\tX\t3\n"+
                    "B\tZ\t1\n"+
                    "C\tY\t2\n"
                    ,outString);
        } catch (IOException e) {
            e.printStackTrace();
            fail();
        }

    }
    @Test
    public void testSmallSortNewLIneLines(){
        ByteArrayInputStream in = new ByteArrayInputStream(SIMPLE_INPUT_WITH_NEWLINES.getBytes());
        ByteArrayOutputStream out = new ByteArrayOutputStream();

        try {
            Sorter.create(in, out, true).field(0, false).sort();
            String outString = new String(out.toByteArray());
            assertEquals(
                    "\n"+
                    "\n"+
                    "\n"+
                    "A\tX\t3\n"+
                    "B\tZ\t1\n"+
                    "C\tY\t2\n"
                    ,outString);
        } catch (IOException e) {
            e.printStackTrace();
            fail();
        }
    }

    @Test
    public void testSmallSortNumbers(){
        ByteArrayInputStream in = new ByteArrayInputStream(SIMPLE_INPUT.getBytes());
        ByteArrayOutputStream out = new ByteArrayOutputStream();

        try {
            Sorter.create(in, out, true).field(2, true).sort();

            String outString = new String(out.toByteArray());
            assertEquals(
                    "B\tZ\t1\n"+
                    "C\tY\t2\n"+
                    "A\tX\t3\n"
                    ,outString);
        } catch (IOException e) {
            e.printStackTrace();
            fail();
        }
    }

    @Test
    public void testSmallSortThreaded(){
        ByteArrayInputStream in = new ByteArrayInputStream(SIMPLE_INPUT.getBytes());
        ByteArrayOutputStream out = new ByteArrayOutputStream();
        Future future = Sorter.create(in, out, true).field(0, false).sortInBackground();
        try {
            future.get();
        } catch (InterruptedException e) {
            e.printStackTrace();
        } catch (ExecutionException e) {
            e.printStackTrace();
            fail();
            return;
        }

        String outString = new String(out.toByteArray());
        assertEquals(
                "A\tX\t3\n"+
                "B\tZ\t1\n"+
                "C\tY\t2\n"
                ,outString);

    }

    @Test
    public void testSmallSortCustomComparator(){
        ByteArrayInputStream in = new ByteArrayInputStream(SIMPLE_INPUT_WITH_NEWLINES.getBytes());
        ByteArrayOutputStream out = new ByteArrayOutputStream();

        try {
            Sorter.create(in, out, true).field(new Comparator<CharSequence>() {
                @Override
                public int compare(final CharSequence o1, final CharSequence o2) {
                    return o1.toString().compareTo(o2.toString());
                }
            }).sort();
            String outString = new String(out.toByteArray());
            assertEquals(
                    "\n"+
                    "\n"+
                    "\n"+
                    "A\tX\t3\n"+
                    "B\tZ\t1\n"+
                    "C\tY\t2\n"
                    ,outString);
        } catch (IOException e) {
            e.printStackTrace();
            fail();
        }
    }



}
