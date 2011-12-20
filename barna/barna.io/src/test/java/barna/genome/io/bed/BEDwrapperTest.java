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

package barna.genome.io.bed;

import barna.commons.Execute;
import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.File;

import static junit.framework.Assert.assertEquals;

/**
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class BEDwrapperTest {

    private static File testfile;

    @BeforeClass
    public static void setUp(){
        testfile = new File(BEDwrapperTest.class.getResource("/test.bed").getFile());
        Execute.initialize(4);
    }

    @AfterClass
    public static void tearDown() throws Exception {
        Execute.shutdown();
    }


    @Test
    public void testScanFile(){
        BEDwrapper wrapper = new BEDwrapper(testfile.getAbsolutePath());
        wrapper.scanFile();

        //scanFileReadLines= 0;
        //countAll= 0; countEntire= 0; countSplit= 0; countReads= 0;
        assertEquals(17, wrapper.nrUniqueLinesRead );
        assertEquals(17, wrapper.countAll );
        assertEquals(5, wrapper.countSplit );
        assertEquals(12, wrapper.countEntire );
        assertEquals(17, wrapper.countReads );
    }
}
