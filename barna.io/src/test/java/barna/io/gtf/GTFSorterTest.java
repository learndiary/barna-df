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

package barna.io.gtf;

import barna.commons.Execute;
import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

import static junit.framework.Assert.assertEquals;

/**
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class GTFSorterTest {

    private static File unsorted;
    private static File sorted;

    @BeforeClass
    public static void setUp(){
        unsorted = new File(GTFSorterTest.class.getResource("/testGtfSort.gtf").getFile());
        sorted = new File(GTFSorterTest.class.getResource("/testGtfSort.gtf.sorted").getFile());
        Execute.initialize(4);
    }

    @AfterClass
    public static void tearDown() throws Exception {
        Execute.shutdown();
    }

    @Test
    public void testGFFSort() throws IOException {
    	GTFwrapper wrapper= new GTFwrapper(unsorted);
        File sortedFile = wrapper.sort();	//sort(unsorted);

        BufferedReader o = new BufferedReader(new FileReader(sorted));
        BufferedReader s = new BufferedReader(new FileReader(sortedFile));
        String l1, l2;
        try{
            int c = 0;
            while( (l1 = o.readLine()) != null && (l2 = s.readLine()) != null){
                assertEquals("Missmatch in line " + c,l1, l2);
                c++;
            }
        }finally {
            o.close();
            s.close();
        }
    }
    @Test
    public void testGFFSort4() throws IOException {
    	GTFwrapper wrapper= new GTFwrapper(new File(GTFSorterTest.class.getResource("/testGtf4.gtf").getFile()).getAbsolutePath());
        File sortedFile = wrapper.sort();
        BufferedReader s = new BufferedReader(new FileReader(sortedFile));
        String l1;
        try{
            int c = 0;
            while( (l1 = s.readLine()) != null){
                System.out.println(l1);
            }
        }finally {
            s.close();
        }
    }

    @Test
    public void testSpikeSorting() throws IOException{
        GTFwrapper wrapper= new GTFwrapper(new File(GTFSorterTest.class.getResource("/spikes.gtf").getFile()).getAbsolutePath());
        File sortedFile = wrapper.sort();
        BufferedReader s = new BufferedReader(new FileReader(sortedFile));
        String l1;
        try{
            int c = 0;
            while( (l1 = s.readLine()) != null){
                System.out.println(l1);
            }
        }finally {
            s.close();
        }

    }
}
