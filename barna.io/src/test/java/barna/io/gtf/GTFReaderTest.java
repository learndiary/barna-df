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

import barna.commons.ByteArrayCharSequence;
import barna.commons.Execute;
import barna.model.Gene;
import barna.model.Transcript;
import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.File;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.fail;

/**
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class GTFReaderTest {


    @BeforeClass
    public static void setUp(){
        Execute.initialize(4);
    }

    @AfterClass
    public static void tearDown() throws Exception {
        Execute.shutdown();
    }

    @Test
    public void testReadCorrectTranscriptStart(){
        // file
        File currentRefFile = new File(getClass().getResource("/sacCer2_sorted.gtf").getFile());
        GTFwrapper reader = new GTFwrapper(currentRefFile.getAbsolutePath());
        // make sure the gtf is valid and sorted
        if (!reader.isApplicable()) {
            fail();
        }
        reader.setSilent(true);
        reader.setStars(false);

        try {
            reader.read();
        } catch (Exception e) {
            e.printStackTrace();
            fail();
        }

        int counter = 0;

        try {
            for (Gene[] g; (g = reader.getGenes()) != null; reader.read()) {
                for (int i = 0; i < g.length; i++) {
                    Gene gene = g[i];
                    for (int j = 0; j < gene.getTranscripts().length; j++) {
                        counter++;
                        Transcript transcript = gene.getTranscripts()[j];
                        ByteArrayCharSequence id = new ByteArrayCharSequence(transcript.getTranscriptID());
                        ByteArrayCharSequence locName = new ByteArrayCharSequence(gene.getGeneID());
                        int exonicLength = transcript.getExonicLength();

                        if(id.toString().equals("YAL068W-A")){
                            assertEquals(255, exonicLength);
                            assertEquals("chrI:335-792W",locName.toString());
                        }

                    }
                }
            }
        } catch (Exception e) {
            e.printStackTrace();
            fail();
        }
        System.out.println(counter);


    }
    @Test
    public void testReadCorrectTranscriptStartSimple(){
        // file
        File currentRefFile = new File(getClass().getResource("/testGtf1.gtf").getFile());

        GTFwrapper reader = new GTFwrapper(currentRefFile.getAbsolutePath());
        reader.setSilent(true);
        reader.setStars(false);

        // make sure the gtf is valid and sorted
        if (!reader.isApplicable()) {
            fail();
        }

        try {
            reader.read();
        } catch (Exception e) {
            e.printStackTrace();
            fail();
        }

        int counter = 0;

        try {
            for (Gene[] g; (g = reader.getGenes()) != null; reader.read()) {
                for (int i = 0; i < g.length; i++) {
                    Gene gene = g[i];
                    for (int j = 0; j < gene.getTranscripts().length; j++) {
                        counter++;
                        Transcript transcript = gene.getTranscripts()[j];
                        ByteArrayCharSequence id = new ByteArrayCharSequence(transcript.getTranscriptID());
                        ByteArrayCharSequence locName = new ByteArrayCharSequence(gene.getGeneID());
                        int exonicLength = transcript.getExonicLength();

                        if(id.toString().equals("YAL068W-A")){
                            assertEquals(255, exonicLength);
                            assertEquals("chrI:538-792W",locName.toString());
                        }

                    }
                }
            }
        } catch (Exception e) {
            e.printStackTrace();
            fail();
        }
        System.out.println(counter);


    }
    @Test
    public void testReadCorrectTranscriptStartSimple2(){
        // file
        File currentRefFile = new File(getClass().getResource("/testGtf2.gtf").getFile());

        GTFwrapper reader = new GTFwrapper(currentRefFile.getAbsolutePath());
        reader.setSilent(true);
        reader.setStars(false);

        // make sure the gtf is valid and sorted
        if (!reader.isApplicable()) {
            fail();
        }

        try {
            reader.read();
        } catch (Exception e) {
            e.printStackTrace();
            fail();
        }

        int counter = 0;

        try {
            for (Gene[] g; (g = reader.getGenes()) != null; reader.read()) {
                for (int i = 0; i < g.length; i++) {
                    Gene gene = g[i];
                    for (int j = 0; j < gene.getTranscripts().length; j++) {
                        counter++;
                        Transcript transcript = gene.getTranscripts()[j];
                        ByteArrayCharSequence id = new ByteArrayCharSequence(transcript.getTranscriptID());
                        ByteArrayCharSequence locName = new ByteArrayCharSequence(gene.getGeneID());
                        int exonicLength = transcript.getExonicLength();

                        if(id.toString().equals("YAL068W-A")){
                            assertEquals(255, exonicLength);
                            assertEquals("chrI:538-792W",locName.toString());
                        }

                    }
                }
            }
        } catch (Exception e) {
            e.printStackTrace();
            fail();
        }
        System.out.println(counter);


    }


    @Test
    public void testReadCorrectTranscriptStartSimple3(){
        // file
        File currentRefFile = new File(getClass().getResource("/testGtf3.gtf").getFile());

        GTFwrapper reader = new GTFwrapper(currentRefFile.getAbsolutePath());
        reader.setSilent(true);
        reader.setStars(false);
        reader.setClusterGenes(false);

        // make sure the gtf is valid and sorted
        if (!reader.isApplicable()) {
            fail();
        }

        try {
            reader.read();
        } catch (Exception e) {
            e.printStackTrace();
            fail();
        }

        int counter = 0;

        try {
            for (Gene[] g; (g = reader.getGenes()) != null; reader.read()) {
                for (int i = 0; i < g.length; i++) {
                    Gene gene = g[i];
                    for (int j = 0; j < gene.getTranscripts().length; j++) {
                        counter++;
                        Transcript transcript = gene.getTranscripts()[j];
                        ByteArrayCharSequence id = new ByteArrayCharSequence(transcript.getTranscriptID());
                        ByteArrayCharSequence locName = new ByteArrayCharSequence(gene.getGeneID());
                        int exonicLength = transcript.getExonicLength();

                        if(id.toString().equals("YAL068W-A")){
                            assertEquals(255, exonicLength);
                            assertEquals("chrI:538-792W",locName.toString());
                        }

                    }
                }
            }
        } catch (Exception e) {
            e.printStackTrace();
            fail();
        }
        System.out.println(counter);


    }

}
