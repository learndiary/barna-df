package barna.flux.capacitor.reconstruction;

import barna.commons.Execute;
import barna.commons.system.OSChecker;
import barna.io.FileHelper;
import barna.io.gtf.GTFwrapper;
import barna.model.Gene;
import org.junit.*;
import org.junit.Test;

import static junit.framework.Assert.assertEquals;
import static junit.framework.Assert.assertTrue;
import static junit.framework.Assert.fail;
import java.io.File;
import java.util.HashMap;
import java.util.Locale;
import java.util.Random;

/**
 * Created with IntelliJ IDEA.
 * User: micha
 * Date: 11/8/13
 * Time: 3:24 PM
 * To change this template use File | Settings | File Templates.
 */
public class PreprocessorTest {

    //getClass().getResource("/test.bam").getFile()
    // /Volumes/Raptor/annotation/hg19/gencode_v12.gtf
    // BigMac
    // /Volumes/Raptor/annotation/hg19/gencode_v12.gtf
    // LNCC
    // /home/micha/gencode_v12.gtf
    // /home/micha/ERR030892.rnd24M.filtered.sorted.bam
    private File gencode12 = new File("/Volumes/Raptor/annotation/hg19/gencode_v12.gtf");
    // file with mappings presorted by name
    private File mappingsPsort = new File("/Volumes/Raptor/scratch/hg19_gencode_paired_sim01.filtered.bam");
    // file with mappings presorted by position
    private File mappingsQsort = new File("/Volumes/Raptor/scratch/hg19_gencode_paired_sim01.filtered.qsort.bam");



    @BeforeClass
    public static void initExecuter() {
            //Force en-US locale to use "." as the decimal separator in Windows OS
            if (OSChecker.isWindows()) {
            Locale.setDefault(new Locale("en", "US"));
    }
            Execute.initialize(2);
    }

    @AfterClass
    public static void shutdownExecuter() {
            Execute.shutdown();
    }

    File currentTestDirectory = null;
    @Before
    public void setUpTest() throws Exception {
            currentTestDirectory = FileHelper.createTempDir(getClass().getSimpleName(), "", null);
    }

    @After
    public void cleanup(){
            if(currentTestDirectory != null){
                FileHelper.rmDir(currentTestDirectory);
            }
    }

    /**
     * Gene loci derived from the Gencode annotation.
     */
    private Gene[] gencodeGenes= null;

    private Gene[] getGencodeGenes() {
        if (gencodeGenes== null) {
            GTFwrapper wrapper= new GTFwrapper(gencode12);
            wrapper.loadAllGenes();
            gencodeGenes= wrapper.getGenes();
        }
        return gencodeGenes;
    }

    @Test
    public void testCollapse() throws Exception {
        Gene[] oGenes= getGencodeGenes();
        Gene[] cGenes= PreProcessor.collapse(oGenes);
        System.err.println(((oGenes.length- cGenes.length)/ 2)+ " antisense loci.");
        Gene[] dGenes= PreProcessor.collapse(cGenes);
        assertTrue(dGenes.length== cGenes.length);
    }

    @Test
    public void testIndex() throws Exception {

        int readLength= 75;

        Gene[] oGenes= getGencodeGenes();
        Gene[] cGenes= PreProcessor.collapse(oGenes);
        HashMap<String,Gene[]> hGenes= new HashMap<String, Gene[]>();
        PreProcessor.index(cGenes, hGenes);

        Random r= new Random();
        Gene qGene= null;
        int rPos= -1;
        for (int i = 0; i < 100; i++) {
            int p= r.nextInt(cGenes.length);
            int gLen= Math.abs(cGenes[p].getEnd())- Math.abs(cGenes[p].getStart());
            // check contained
            rPos= Math.abs(cGenes[p].getStart())
                    + r.nextInt(gLen);
            // ensure that read is contained at both ends
            qGene= PreProcessor.getGene(hGenes,
                    cGenes[p].getChromosome(),
                    rPos,
                    Math.min(rPos + readLength, Math.abs(cGenes[p].getEnd())),
                    true);
            assertEquals(cGenes[p], qGene);
            // check overlap at start
            rPos= Math.abs(cGenes[p].getStart())+ r.nextInt(readLength- 1);
            int min= (p> 0? Math.abs(cGenes[p- 1].getEnd())+ 1: 0);
            qGene= PreProcessor.getGene(hGenes,
                    cGenes[p].getChromosome(),
                    Math.max(rPos - readLength, min),
                    Math.min(rPos, Math.abs(cGenes[p].getEnd())),
                    false);
            assertEquals(cGenes[p], qGene);
            // check overlap at end
            rPos= Math.abs(cGenes[p].getEnd())+ r.nextInt(readLength- 1);
            int max= (p< cGenes.length- 1? Math.abs(cGenes[p+ 1].getStart())- 1: Integer.MAX_VALUE);
            qGene= PreProcessor.getGene(hGenes,
                    cGenes[p].getChromosome(),
                    Math.max(rPos - readLength, Math.abs(cGenes[p].getStart())),
                    Math.min(rPos, max),
                    false);
            assertEquals(cGenes[p], qGene);
        }
    }

    @Test
    public void testProcess() {
        PreProcessor pp= new PreProcessor(gencode12, mappingsQsort);
        File result= pp.call();
        assertTrue(result!= null);
    }

}
