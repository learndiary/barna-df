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

import org.junit.Test;

import java.io.File;

import static junit.framework.Assert.*;

/**
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class ProfilerTest {

    @Test
    public void testConstructor() throws Exception {
        try{
            Profiler profiler = new Profiler(null);
            fail();
        }catch (NullPointerException e){
            // expected
            assertEquals("You have to specify settings! NULL not permitted.", e.getMessage());
        }
    }


    @Test
    public void testExpression() throws Exception {
        FluxSimulatorSettings settings = FluxSimulatorSettings.createSettings(new File(getClass().getResource("/minimal.par").getFile()));
        Profiler profiler = new Profiler(settings);
        profiler.resetProfile();
        profiler.call();

        /*
        Lambdaclone23-2	spike	exon	1	9786	.	+	.	transcript_id "Lambdaclone23-2";
        Lambdaclone1-1	spike	exon	1	11934	.	+	.	transcript_id "Lambdaclone1-1";
        Apetala2	spike	exon	1	1405	.	+	.	transcript_id "Apetala2";
        OBF5	spike	exon	1	1429	.	+	.	transcript_id "OBF5";
        EPR-1	spike	exon	1	1451	.	+	.	transcript_id "EPR-1";
        AGP	spike	exon	1	325	.	+	.	transcript_id "AGP";
        VATG	spike	exon	1	376	.	+	.	transcript_id "VATG";

         */
        assertEquals(7, profiler.size());
        assertEquals(9786, profiler.getLength(0));
        assertEquals(11934, profiler.getLength(1));
        assertEquals(1405, profiler.getLength(2));
        assertEquals(1429, profiler.getLength(3));
        assertEquals(1451, profiler.getLength(4));
        assertEquals(325, profiler.getLength(5));
        assertEquals(376, profiler.getLength(6));

        assertEquals(11934, profiler.getMaxMoleculeLength());

        int s = 0;
        for(int i=0; i<profiler.size(); i++){
            s += profiler.getMolecules(i);
        }
        assertTrue(s > 0);

    }
    @Test
    public void testLoadingProfile() throws Exception {
        FluxSimulatorSettings settings = FluxSimulatorSettings.createSettings(new File(getClass().getResource("/minimal.par").getFile()));
        File testProFile = new File(getClass().getResource("/test.pro").getFile());
        settings.set(FluxSimulatorSettings.PRO_FILE, testProFile);
        Profiler profiler = new Profiler(settings);
        //profiler.initializeProfiler(testProFile);

        /*
        chrI:335-792W   YAL069W CDS     315     3.7997796127824584E-4   190
        chrI:335-792W   YAL068W-A       CDS     255     5.5996752188373074E-5   28
        chrI:2480-2707W YAL067W-A       CDS     228     2.9598283299568624E-4   148
        chrI:10092-10400W       YAL066W CDS     309     3.5997912121096976E-5   18
        chrI:12047-12427W       YAL064W-B       CDS     381     1.419917644776603E-4    71
         */
        assertTrue(profiler.initializeProfiler(testProFile));
        assertEquals(5, profiler.size());
        assertEquals(315, profiler.getLength(0));

        // test number of molecules
        assertEquals(190, profiler.getMolecules(0));
        assertEquals(28, profiler.getMolecules(1));
        assertEquals(148, profiler.getMolecules(2));
        assertEquals(18, profiler.getMolecules(3));
        assertEquals(71, profiler.getMolecules(4));

    }

}
