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

package barna.io;

import barna.commons.Execute;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

import java.io.*;
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


    @Before
    public void before(){
        Execute.initialize(4);
    }
    @After
    public void after(){
        Execute.shutdown();
    }

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

    public static void main(String[] args) throws Exception {
        Execute.initialize(16);
        UnixStreamSorter sorter = new UnixStreamSorter(200*1024*1024, false, -1, false, "\t");
//        UnixStreamSorter sorter = new UnixStreamSorter(200*1024*1024, new Comparator<String>(){
//            @Override
//            public int compare(String o1, String o2) {
//                return o1.compareTo(o2);
//            }
//        });
        String input = "/home/thasso/data/simulator/genomes/mm9/mm9.fasta";
        String output = "/home/thasso/data/simulator/genomes/mm9/mm9.fasta.sorted";
        long t = System.currentTimeMillis();
        sorter.sort(new BufferedInputStream(new FileInputStream(input)), new BufferedOutputStream(new FileOutputStream(output)));
        System.out.println((System.currentTimeMillis()-t)/1000+"s");
        Execute.shutdown();
    }


}
