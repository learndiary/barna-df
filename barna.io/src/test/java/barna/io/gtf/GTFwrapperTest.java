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

import barna.io.FileHelper;
import barna.model.Gene;
import org.junit.Test;

import java.io.*;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.TreeSet;

import static junit.framework.Assert.*;
import static junit.framework.Assert.assertEquals;

public class GTFwrapperTest {

    /**
     * The test data target directory
     */
    private static String targetDirectory;

    public static File getGencodeFile() {
        if (gencodeFile== null)
            prepareTestData();
        return gencodeFile;
    }

    /**
     * Complete Gencode annotation
     */
    private static File gencodeFile= getGencodeFile();

    static {
        // initialize defaults
        //artifactoryUrl = System.getProperty("testdata.artifactory", "http://sammeth.net/artifactory")
        //repository = System.getProperty("testdata.repository", "repo")
        //artifact = System.getProperty("testdata.artifact", "barna/test_data-1.0.zip")
        targetDirectory = System.getProperty("testdata.target", new File("").getAbsolutePath());
    }

    private static File prepareTestData() {

        System.out.println("Checking for test data");
        //JsonSlurper slurper = new JsonSlurper();
        URL dataFQN = null;
        try {
            dataFQN = new URL(
                    // ${artifactoryUrl}/api/storage/${repository}/${artifact}
                    " http://sammeth.net/artifactory/repo/gencode_v12_gtf/gencode_v12_gtf/12/gencode_v12_gtf-12.gz");
        } catch (MalformedURLException e) {
            e.printStackTrace();
        }
        String[] tokens= dataFQN.getFile().split("/");
        String fileName = tokens[tokens.length- 1];


        File targetFile = new File(targetDirectory, fileName);
        File test_data_dir = new File(targetDirectory, "test_data");
        File md5File = new File(targetDirectory, "${fileName}.md5");

        //metaData = slurper.parseText(dataFQN.openStream().text);
        //def md5sum = metaData.checksums['md5']

        boolean invalid_file = false;
        if (!targetFile.exists()) {
            // || !md5File.exists()) || md5File.readLines()[0].trim() != md5sum) {
            invalid_file = true;
        }

        // download
        if (invalid_file) {
            System.out.println("Downloading test data from ${dataFQN.toExternalForm()}");
            if (targetFile.exists()) {
                targetFile.delete();
            }
            if (test_data_dir.exists()) {
                FileHelper.rmDir(new File(targetDirectory, "test_data"));
            }
            if (md5File.exists()) {
                md5File.delete();
            }
            // download
            //out << new URL(metaData['downloadUri']).openStream();
            try {
                InputStream in= dataFQN.openStream();
                OutputStream out = new FileOutputStream(targetFile);
                byte[] buffer = new byte[1024];
                int len;
                while ((len = in.read(buffer)) != -1) {
                    out.write(buffer, 0, len);
                }
                out.close();
            } catch (Exception e) {
                e.printStackTrace();
            }

            // compare md5
            //if (md5sum != generateMD5(targetFile)) {
            //    throw new RuntimeException("Test Data downlaoded but md5 sums do not match !");
            //}

            // write md5 file
            //md5File.write(md5sum)
        }

        gencodeFile= new File(FileHelper.stripExtension(targetFile.getAbsolutePath()));
        if (!gencodeFile.exists()) {
            System.out.println("Unzipping test data");
            //FluxCapacitorRunner.unzip(targetFile, targetDirectory);
            try {
                FileHelper.inflate(targetFile, gencodeFile, FileHelper.getCompression(targetFile));
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
        System.out.println("Test data available");
        return test_data_dir;
    }


    @Test
    public void testIteratorAll() {
        GTFwrapper wrapper= new GTFwrapper(gencodeFile);
        wrapper.setLoadAllGenes();
        assertTrue(wrapper.isApplicable());

        GTFwrapper.GeneIterator iter= wrapper.new GeneIterator();
        int n= 0;
        while (iter.hasNext()) {
            iter.next();
            ++n;
        }

        assertEquals(45013, n);
    }


    @Test
    public void testIteratorProgressive() {
        GTFwrapper wrapper= new GTFwrapper(gencodeFile);
        assertTrue(wrapper.isApplicable());
        wrapper.reset();

        GTFwrapper.GeneIterator iter= wrapper.new GeneIterator();
        int n= 0;
        while (iter.hasNext()) {
            iter.next();
            ++n;
        }

        assertEquals(45013, n);
    }


    @Test
    public void testLoaderProgressive() {
        GTFwrapper wrapper= new GTFwrapper(gencodeFile);
        assertTrue(wrapper.isApplicable());

        GTFwrapper.GeneLoader loader= wrapper.new GeneLoader();
        loader.start();

        int nrGenes= 0;
        while(true) {
            Gene[] g= loader.fetch();
            if (g== null)
                break;
            // else
            nrGenes+= g.length;
            /*for (int i = 0; i < g.length; i++) {
                try {
                    BufferedWriter buffy= new BufferedWriter(new FileWriter("/Volumes/Raptor/scratch/progress.txt", true));
                    buffy.write(g[i].getChromosome()+ "\t"+ g[i].getStart()+ "\t"+ g[i].getEnd()+ "\n");
                    buffy.close();
                } catch (Exception e) {
                    throw new RuntimeException(e);
                }
            }*/
        }
        System.err.println("Genes "+ nrGenes);
        assertEquals(45013, nrGenes);
    }

    @Test
    public void testLoaderAll() {
        GTFwrapper wrapper= new GTFwrapper(gencodeFile);
        wrapper.setLoadAllGenes();
        assertTrue(wrapper.isApplicable());

        GTFwrapper.GeneLoader loader= wrapper.new GeneLoader();
        loader.start();

        int nrGenes= 0;
        while(true) {
            Gene[] g= loader.fetch();
            if (g== null)
                break;
            // else
            nrGenes+= g.length;
            /*for (int i = 0; i < g.length; i++) {
                try {
                    BufferedWriter buffy= new BufferedWriter(new FileWriter("/Volumes/Raptor/scratch/loadall.txt", true));
                    buffy.write(g[i].getChromosome()+ "\t"+ g[i].getStart()+ "\t"+ g[i].getEnd()+ "\n");
                    buffy.close();
                } catch (Exception e) {
                    throw new RuntimeException(e);
                }

            }*/
        }
        System.err.println("Genes "+ nrGenes);
        assertEquals(45013, nrGenes);
    }

    /**
     * Test that loads all genes up to the level of detail necessary to construct splicing graphs.
     */
    @Test
    public void testLoadAll() {

        GTFwrapper wrapper= new GTFwrapper(gencodeFile);
        assertTrue(wrapper.isApplicable());
        long t0= System.currentTimeMillis();
        wrapper.loadAllGenes();
        long t1= System.currentTimeMillis();
        System.err.println((t1- t0)/ 1000+ " sec.");

    }

    @Test
    public void testThatTheWrapperWorksWithGzippedFiles() throws Exception {
        File gzippedGtf = new File(getClass().getResource("/gzipped-gtf.gtf.gz").getFile());
        assertTrue(gzippedGtf.exists());

        GTFwrapper wrapperFile = new GTFwrapper(gzippedGtf);
        assertFalse(wrapperFile.isApplicable());

        File sorted = wrapperFile.sort();
        assertNotNull(sorted);
        assertTrue(new GTFwrapper(sorted).isApplicable());
        sorted.delete();
    }
    @Test
    public void testThatNormaltabsForFirstFieldsAreApplied() throws Exception {
        File gzippedGtf = new File(getClass().getResource("/valid-tab-space.gtf").getFile());
        assertTrue(gzippedGtf.exists());
        GTFwrapper wrapperFile = new GTFwrapper(gzippedGtf);
        File sorted = wrapperFile.sort();
        assertNotNull(sorted);
        assertTrue(new GTFwrapper(sorted).isApplicable());
        sorted.delete();
    }
    @Test
    public void testThatOnlySpaceSeparatedFields() throws Exception {
        File gzippedGtf = new File(getClass().getResource("/only-space.gtf").getFile());
        assertTrue(gzippedGtf.exists());
        GTFwrapper wrapperFile = new GTFwrapper(gzippedGtf);
        File sorted = wrapperFile.sort();
        assertNotNull(sorted);
        assertTrue(new GTFwrapper(sorted).isApplicable());
        sorted.delete();
    }
    @Test
    public void testMixedFileWithTabsAndSpaces() throws Exception {
        File gzippedGtf = new File(getClass().getResource("/mixed-space.gtf").getFile());
        assertTrue(gzippedGtf.exists());
        GTFwrapper wrapperFile = new GTFwrapper(gzippedGtf);
        File sorted = wrapperFile.sort();
        assertNotNull(sorted);
        assertTrue(new GTFwrapper(sorted).isApplicable());
        sorted.delete();
    }
    @Test
    public void testWithTabsInTheGTFFieldsAtTheEnd() throws Exception {
        File gzippedGtf = new File(getClass().getResource("/tab-fields.gtf").getFile());
        assertTrue(gzippedGtf.exists());
        GTFwrapper wrapperFile = new GTFwrapper(gzippedGtf);
        File sorted = wrapperFile.sort();
        assertNotNull(sorted);
        assertTrue(new GTFwrapper(sorted).isApplicable());
        sorted.delete();
    }
    @Test
    public void testThatInvalidGTFFailsWithReport() throws Exception {
        File gzippedGtf = new File(getClass().getResource("/invalid.gtf").getFile());
        assertTrue(gzippedGtf.exists());
        GTFwrapper wrapperFile = new GTFwrapper(gzippedGtf);
        try{
            File sorted = wrapperFile.sort();
            fail();
        }catch (Exception error){
        }
    }
    @Test
    public void testEnsemleGeneLoading() throws Exception {
        File gzippedGtf = new File(getClass().getResource("/BARNA-268-ensemble.gtf").getFile());
        GTFwrapper wrapperFile = new GTFwrapper(gzippedGtf);
        wrapperFile.read();
        Gene[] genes = null;
        TreeSet<String> ids = new TreeSet<String>();
        while((genes = wrapperFile.getGenes()) != null){
            for (Gene gene : genes) {
                ids.add(gene.getGeneID());
            }
            wrapperFile.read();
        }
        assertEquals(3, ids.size());
        assertTrue(ids.contains("ENSG00000223972"));
        assertTrue(ids.contains("ENSG00000227159"));
        assertTrue(ids.contains("ENSG00000233614"));

    }
}
