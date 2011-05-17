package fbi.commons.tools;

import org.junit.Test;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;

import static junit.framework.Assert.assertEquals;
import static junit.framework.Assert.fail;

/**
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class UnixStreamSorterTest {
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

        StreamSorter sorter = new UnixStreamSorter(0, false, "\t");
        try {
            sorter.sort(in, out);

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
    public void testSmallSortSmallMem(){
        byte[] bytes = SIMPLE_INPUT.getBytes();
        System.out.println("Len: " + bytes.length);
        ByteArrayInputStream in = new ByteArrayInputStream(bytes);
        ByteArrayOutputStream out = new ByteArrayOutputStream();

        StreamSorter sorter = new UnixStreamSorter(10, 0, false, "\t");
        try {
            sorter.sort(in, out);

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

        StreamSorter sorter = new UnixStreamSorter(0, false, "\t");
        try {
            sorter.sort(in, out);

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

        StreamSorter sorter = new UnixStreamSorter(2, true, "\t");
        try {
            sorter.sort(in, out);

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
    public void testRandomNumbers(){
        Random r = new Random();
        StringBuffer bb = new StringBuffer();
        List<Double> data = new ArrayList<Double>();
        for (int i = 0; i < 1000; i++) {
            double d = r.nextDouble();
            bb.append(d).append("\n");
            data.add(d);
        }

        Collections.sort(data);

        ByteArrayInputStream in = new ByteArrayInputStream(bb.toString().getBytes());
        ByteArrayOutputStream out = new ByteArrayOutputStream();

        StreamSorter sorter = new UnixStreamSorter(0, true, "\t");
        try {
            sorter.sort(in, out);

            String outString = new String(out.toByteArray());
            String[] lines = outString.split("\\n");
            for (int i = 0; i < lines.length; i++) {
                String line = lines[i];
                assertEquals(data.get(i).toString(), line);
            }
        } catch (IOException e) {
            e.printStackTrace();
            fail();
        }
    }

}
