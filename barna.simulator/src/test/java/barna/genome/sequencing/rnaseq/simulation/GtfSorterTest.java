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

package barna.genome.sequencing.rnaseq.simulation;

import barna.commons.ByteArrayCharSequence;
import barna.commons.Execute;
import barna.genome.io.gtf.GTFwrapper;
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
