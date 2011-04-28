package fbi.genome.sequencing.rnaseq.simulation;

import org.junit.BeforeClass;
import org.junit.Test;

import java.io.File;

import static junit.framework.Assert.*;

/**
 * Test flus simulator paramter file loading
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class FluxSimulatorSettingsTest {

    private static File par1;

    @BeforeClass
    public static void setUp(){
        par1 = new File(FluxSimulatorTest.class.getResource("/test_sacCer2.par").getFile());
    }

    @Test
    public void testPar1(){
        FluxSimulatorSettings s1 = FluxSimulatorSettings.createSettings(par1);

        /*
        REF_FILE_NAME	sacCer2_sgdGene_fromUCSC110329.gtf
        PRO_FILE_NAME	test_sacCer2.pro
        LIB_FILE_NAME	test_sacCer2.lib
        SEQ_FILE_NAME	test_sacCer2.bed
        ERR_FILE_NAME	error_model_qualities.err
        GEN_DIR	genome_sacCer2
        NB_MOLECULES	1000000
        EXPRESSION_K	-0.6
        EXPRESSION_X0	5.0E7
        EXPRESSION_X1	9500.0
        TSS_MEAN	25.0
        POLYA_SHAPE	5.0
        POLYA_SCALE	50.0
        RT_MIN	500
        RT_MAX	5500
        RT_PRIMER	RANDOM
        FRAGMENTATION	YES
        FRAG_B4_RT	NO
        FRAG_MODE	PHYSICAL
        FRAG_LAMBDA	900.0
        FRAG_SIGMA	0.05
        FRAG_THRESHOLD	0.1
        FILTERING	YES
        LOAD_CODING	YES
        LOAD_NONCODING	YES
        FILT_MIN	100
        FILT_MAX	250
        READ_NUMBER	700176
        READ_LENGTH	36
        PAIRED_END	YES
        FASTQ	YES
         */

        assertEquals("sacCer2_sgdGene_fromUCSC110329.gtf",s1.getRefFile().getName());
        assertEquals("test_sacCer2.pro",s1.getProFile().getName());
        assertEquals("test_sacCer2.lib",s1.getFrgFile().getName());
        assertEquals("test_sacCer2.bed",s1.getSeqFile().getName());
        assertEquals("error_model_qualities.err",s1.getErrFile().getName());
        assertEquals("genome_sacCer2", s1.getGenDir().getName());
        assertEquals(1000000, s1.getNbMolecules());
        assertEquals(-0.6, s1.getExpDistrP1(), 0.000001); // EXPRESSION K
        assertEquals(50000000, s1.getExpDistrP2(), 0.000001); // EXPRESSION X0
        assertEquals(9500.0, s1.getDecDistrP1(), 0.000001); // EXPRESSION X1
        assertEquals(25.0, s1.getTssMean(), 0.000001);
        assertEquals(5.0, s1.getPolyAshape(), 0.000001);
        assertEquals(50.0, s1.getPolyAscale(), 0.000001);
        assertEquals(500, s1.getRTminLen(), 0.000001);
        assertEquals(5500, s1.getRTmaxLen(), 0.000001);
        assertEquals("RANDOM", s1.getRtMode());
        assertTrue(s1.isFragment());
        assertFalse(s1.isFragB4RT());
        assertEquals("PHYSICAL", s1.getFragMode());
        assertEquals(900.0, s1.getFragNBlambda(), 0.000001); // FRAG LAMDA
        assertEquals(0.05, s1.getFragNBsigma(), 0.000001); // FRAG SIGMA
        assertEquals(0.1, s1.getThold(), 0.000001); // FRAG THRESHOLD
        assertTrue(s1.isFilter());
        assertTrue(s1.isLoadCoding());
        assertTrue(s1.isLoadNoncoding());
        assertEquals(100, s1.getFiltMin());
        assertEquals(250, s1.getFiltMax());
        assertEquals(700176, s1.getReadNr());
        assertEquals(36, s1.getReadLength());
        assertTrue(s1.isPairedEnd());
        assertTrue(s1.isFastQ());
        assertEquals(5500, s1.getRTmaxLen());
    }
}
