package fbi.commons.tools;

import org.junit.Test;

import java.io.*;

import static junit.framework.Assert.assertEquals;
import static junit.framework.Assert.fail;

/**
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class UnixSortTest {
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


    @Test
    public void testSmallSort(){
        ByteArrayInputStream in = new ByteArrayInputStream(SIMPLE_INPUT.getBytes());
        ByteArrayOutputStream out = new ByteArrayOutputStream();

        StreamSorter sorter = new UnixSort();
        try {
            sorter.sort(in, out, 0, false, "\t");

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

        StreamSorter sorter = new UnixSort();
        try {
            sorter.sort(in, out, 0, false, "\t");

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

        StreamSorter sorter = new UnixSort();
        try {
            sorter.sort(in, out, 2, true, "\t");

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
    public void testrealWorld() throws FileNotFoundException {
        FileInputStream in = new FileInputStream("/home/thasso/data/big.dat");
        FileOutputStream out = new FileOutputStream("/home/thasso/data/big-java.dat");

        StreamSorter sorter = new UnixSort(false);
        try {
            sorter.sort(in, out, -1, false, "\t");
        } catch (IOException e) {
            e.printStackTrace();
            fail();
        }
    }


    public static void main(String[] args) throws FileNotFoundException {

        try {
            Thread.sleep(20000);
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
        FileInputStream in = new FileInputStream("/home/thasso/data/big.dat");
        FileOutputStream out = new FileOutputStream("/home/thasso/data/big-java.dat");

        StreamSorter sorter = new UnixSort(false);
        try {
            sorter.sort(in, out, -1, false, "\t");
        } catch (IOException e) {
            e.printStackTrace();
            fail();
        }

    }
}
