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

import fbi.commons.Execute;
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
public class GFFSorterTest {

    private static File unsorted;
    private static File sorted;

    @BeforeClass
    public static void setUp(){
        unsorted = new File(GFFSorterTest.class.getResource("/testGtfSort.gtf").getFile());
        sorted = new File(GFFSorterTest.class.getResource("/testGtfSort.gtf.sorted").getFile());
        Execute.initialize(4);
    }

    @AfterClass
    public static void tearDown() throws Exception {
        Execute.shutdown();
    }

    @Test
    public void testGFFSort() throws IOException {
        File sortedFile = GFFSorter.sort(unsorted);

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
}
