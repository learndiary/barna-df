/*
 * This file is part of the Flux Library.
 *
 * The code of the Flux Library may be freely distributed and modified under the terms of the
 * European Union Public Licence (EUPL) published on the web site <http://www.osor.eu/eupl/european-union-public-licence-eupl-v.1.1>.
 * Copyright for the code is held jointly by the individual authors, who should be listed
 * in @author doc comments. According to Article 5 and Article 11 of the EUPL, publications that
 * include results produced by the Flux Library are liable to reference the Work,
 * see the Flux Library homepage <http://flux.sammeth.net> for more information.
 */

package fbi.genome.io.gff;

import fbi.commons.ByteArrayCharSequence;
import fbi.commons.Execute;
import fbi.genome.model.Gene;
import fbi.genome.model.Transcript;
import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.File;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.fail;

/**
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class GFFReaderTest {


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
        GFFReader reader = new GFFReader(currentRefFile.getAbsolutePath());
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

        GFFReader reader = new GFFReader(currentRefFile.getAbsolutePath());
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

        GFFReader reader = new GFFReader(currentRefFile.getAbsolutePath());
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

        GFFReader reader = new GFFReader(currentRefFile.getAbsolutePath());
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
