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

package barna.flux.simulator;

import barna.commons.ByteArrayCharSequence;
import barna.commons.Execute;
import barna.io.gtf.GTFwrapper;
import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;

import static junit.framework.Assert.assertEquals;
import static junit.framework.Assert.fail;

/**
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class GtfSorterTest {

    @BeforeClass
    public static void setUp() throws Exception {
        Execute.initialize(2);
    }

    @AfterClass
    public static void tearDown() throws Exception {
        Execute.shutdown();
    }

    /**
     * Test #56
     * @throws Exception
     */
    @Test
    public void testSortErrorNPE() throws Exception {
        File out = File.createTempFile("sortertest", ".out");
        try{
            File in = new File(getClass().getResource("/sort-err56.gtf").getFile());
            GTFwrapper w = new GTFwrapper(in);
            w.sort(out);

            // check the result
            BufferedReader reader = new BufferedReader(new FileReader(out));
            StringBuffer bb = new StringBuffer();
            String l = null;
            while((l = reader.readLine()) != null){
                bb.append(l).append("\n");
                System.out.println(l);
            }
        }catch (Exception e){
            e.printStackTrace();
            fail();
        }finally {
            out.delete();
        }


    }
    //@Test
    public void testSortBigNPE() throws Exception {
        File out = File.createTempFile("sortertest", ".out");
        try{
            File in = new File(getClass().getResource("/mutations.gtf").getFile());
            if(!in.exists()) return;

            GTFwrapper w = new GTFwrapper(in);
            w.sort(out);


            // check the result
            BufferedReader reader = new BufferedReader(new FileReader(out));
            StringBuffer bb = new StringBuffer();
            String l = null;
            while((l = reader.readLine()) != null){
                bb.append(l).append("\n");
                System.out.println(l);
            }
        }catch (Exception e){
            e.printStackTrace();
            fail();
        }finally {
            out.delete();
        }


    }

    @Test
    public void testFinds() throws Exception {
        String s = "chrY\tGenomeSimulator\texon\t9343455\t9343664\t.\t+\t.\tgene_id \"ENSG00000228927\";\ttranscript_id \"ENST00000451548\";\tgene_name \"9272328-9343664\";";
        ByteArrayCharSequence ss = new ByteArrayCharSequence(s);
        assertEquals("chrY",ss.getToken(0).toString());
        assertEquals("GenomeSimulator",ss.getToken(1).toString());
        assertEquals("exon",ss.getToken(2).toString());
        assertEquals("9343455",ss.getToken(3).toString());
        assertEquals("9343664",ss.getToken(4).toString());
        assertEquals(".",ss.getToken(5).toString());
        assertEquals("+",ss.getToken(6).toString());
        assertEquals(".",ss.getToken(7).toString());
        assertEquals("gene_id \"ENSG00000228927\";",ss.getToken(8).toString());
        assertEquals("transcript_id \"ENST00000451548\";",ss.getToken(9).toString());
        assertEquals("gene_name \"9272328-9343664\";",ss.getToken(10).toString());

    }
}
