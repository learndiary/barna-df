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
import barna.commons.system.OSChecker;
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
            "C\tY\t2"+ OSChecker.NEW_LINE +
            "B\tZ\t1"+ OSChecker.NEW_LINE +
            "A\tX\t3"+ OSChecker.NEW_LINE;

    private static final String SIMPLE_INPUT_WITH_NEWLINES =
            barna.commons.system.OSChecker.NEW_LINE +
            "C\tY\t2"+ OSChecker.NEW_LINE +
            barna.commons.system.OSChecker.NEW_LINE +
            "B\tZ\t1"+ OSChecker.NEW_LINE +
            barna.commons.system.OSChecker.NEW_LINE +
            "A\tX\t3"+ OSChecker.NEW_LINE;


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
            Sorter.create(in, out, true, null).field(0, false).sort();

            String outString = new String(out.toByteArray());
            assertEquals(
                    "A\tX\t3"+ OSChecker.NEW_LINE+
                    "B\tZ\t1"+ OSChecker.NEW_LINE+
                    "C\tY\t2"+ OSChecker.NEW_LINE
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
            Sorter.create(in, out, true, null).field(0, false).sort();
            String outString = new String(out.toByteArray());
            assertEquals(
                    barna.commons.system.OSChecker.NEW_LINE+
                    barna.commons.system.OSChecker.NEW_LINE+
                    barna.commons.system.OSChecker.NEW_LINE+
                    "A\tX\t3"+ OSChecker.NEW_LINE+
                    "B\tZ\t1"+ OSChecker.NEW_LINE+
                    "C\tY\t2"+ OSChecker.NEW_LINE
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
            Sorter.create(in, out, true, null).field(2, true).sort();

            String outString = new String(out.toByteArray());
            assertEquals(
                    "B\tZ\t1"+ OSChecker.NEW_LINE+
                    "C\tY\t2"+ OSChecker.NEW_LINE+
                    "A\tX\t3"+ OSChecker.NEW_LINE
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
        Future future = Sorter.create(in, out, true, null).field(0, false).sortInBackground();
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
                "A\tX\t3"+ OSChecker.NEW_LINE+
                "B\tZ\t1"+ OSChecker.NEW_LINE+
                "C\tY\t2"+ OSChecker.NEW_LINE
                ,outString);

    }

    @Test
    public void testSmallSortCustomComparator(){
        ByteArrayInputStream in = new ByteArrayInputStream(SIMPLE_INPUT_WITH_NEWLINES.getBytes());
        ByteArrayOutputStream out = new ByteArrayOutputStream();

        try {
            Sorter.create(in, out, true, null).field(new Comparator<CharSequence>() {
                @Override
                public int compare(final CharSequence o1, final CharSequence o2) {
                    return o1.toString().compareTo(o2.toString());
                }
            }).sort();
            String outString = new String(out.toByteArray());
            assertEquals(
                    barna.commons.system.OSChecker.NEW_LINE+
                    barna.commons.system.OSChecker.NEW_LINE+
                    barna.commons.system.OSChecker.NEW_LINE+
                    "A\tX\t3"+ OSChecker.NEW_LINE+
                    "B\tZ\t1"+ OSChecker.NEW_LINE+
                    "C\tY\t2"+ OSChecker.NEW_LINE
                    ,outString);
        } catch (IOException e) {
            e.printStackTrace();
            fail();
        }
    }



}
