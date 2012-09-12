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
import barna.flux.simulator.error.MarkovErrorModel;
import barna.flux.simulator.error.ModelPool;
import barna.flux.simulator.error.QualityErrorModel;
import barna.io.FileHelper;
import barna.io.gtf.GTFwrapper;
import barna.model.Exon;
import barna.model.Gene;
import barna.model.Graph;
import barna.model.Transcript;
import barna.model.bed.BEDobject2;
import org.junit.After;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.util.Random;

import static junit.framework.Assert.assertEquals;

public class SequencerTest {

    File currentTestDirectory = null;
    @Before
    public void setUpTest() throws Exception {
        currentTestDirectory = FileHelper.createTempDir("FluxCapacitorUnitTest", "", null);
    }

    @After
    public void cleanup(){
        if(currentTestDirectory != null){
            FileHelper.rmDir(currentTestDirectory);
        }
    }

    @Test
	public void testReadGeneration() throws Exception {

        Sequencer sequencer = new Sequencer(null, null);
        int nrTests= 100; // nr. of tests per transcript
		int maxBoundary= 100; // nt added at beginning/end of tx
		Random rand= new Random();
		
		GTFwrapper reader= new GTFwrapper(new File(
				getClass().getResource("/sacCer2_sgdGene_fromUCSC110329.gtf").getFile()).getAbsolutePath());
		Gene[] genes;
		for (reader.read(); (genes=reader.getGenes())!= null; reader.read()) {
			for (int i = 0; i < genes.length; i++) {
				Gene g= genes[i];
				int n= g.getTranscriptCount();
				for (int j = 0; j < n; j++) {
					Transcript t= g.getTranscripts()[j];
					int tlen= t.getExonicLength();						
					for (int k = 0; k < nrTests; k++) {
						int maxLeft= maxBoundary;
						if (t.getStrand()>= 0)
							maxLeft= Math.min(maxLeft, t.get5PrimeEdge());
						int maxRight= maxBoundary;
						if (t.getStrand()< 0)
							maxRight= Math.min(maxRight, Math.abs(t.get3PrimeEdge()));
						int boundLeft= (maxLeft== 0? 0: rand.nextInt(maxLeft));
						int boundRight= (maxRight== 0? 0: rand.nextInt(maxRight));
						
						byte absDir= t.getStrand();
						int fragStart= rand.nextInt(tlen+ boundLeft+ boundRight)- boundLeft;
						int fragEnd= fragStart+ rand.nextInt(tlen+ boundRight- fragStart);
						int fragLen= fragEnd- fragStart+ 1;
						if (fragLen<= 1)
							continue;
						int readLen= 1+ rand.nextInt(fragLen- 1);
						boolean left= rand.nextBoolean();
						int start= (left? fragStart: fragEnd- readLen+ 1);
						int end= (left? fragStart+ readLen- 1: fragEnd);
						assert(start<= end);
						BEDobject2 obj= new BEDobject2();
						int polyA= sequencer.createRead(obj, start, end, t, k,
                                absDir, fragStart, fragEnd, left, 1);
						
						assertEquals(readLen, obj.getLength());
						int bedLeft= t.getGenomicPosition(start);
						int bedRight= t.getGenomicPosition(end);
						if (absDir< 0) {
							int h= bedLeft;
							bedLeft= -bedRight;
							bedRight= -h;
						}
						--bedLeft;
						
						if (polyA> 0) {
							if (polyA== readLen) {
								assertEquals(0, obj.getStart());
								assertEquals(readLen, obj.getEnd());
							}
						} else {
							assertEquals(bedLeft, obj.getStart());
							assertEquals(bedRight, obj.getEnd());
						}						
						obj.getLength();
					}
				}
			}
		}
	}

    @Test
    public void testSequenceAmbiguities() {
        Sequencer sequencer = new Sequencer(null, null);
        String allChars= "ACGUTWSRYMKBDHV";
        Graph.fileSep = null;

        // write sequence
        File f= null;
        try {
            f= File.createTempFile(this.getClass().getSimpleName(), ".fa",currentTestDirectory);
            BufferedWriter writer= new BufferedWriter(new FileWriter(f));
            writer.write(">"+ FileHelper.stripExtension(f.getName())+ barna.commons.system.OSChecker.NEW_LINE);
            writer.write(allChars+ barna.commons.system.OSChecker.NEW_LINE);
            writer.close();
        } catch (Exception e) {
            e.printStackTrace();
        }

        // generate object
/*        ByteArrayCharSequence cs= new ByteArrayCharSequence(
                FileHelper.stripExtension(f.getName())+ "\t"+
                Integer.toString(1)+ "\t"+ Integer.toString(11)+ "\t"+
                "TestRead+\t"+ Byte.toString(1)+ "\t"+
        );
*/
        BEDobject2 obj= new BEDobject2();
        Graph.overrideSequenceDirPath= f.getParent();
        obj.setChromosome(FileHelper.stripExtension(f.getName()));
        obj.setStart(0);
        obj.setEnd(allChars.length());
        obj.setName("TestRead");
        obj.setStrand((byte) 1);

        // load model
        File eFile= new File(getClass().getResource("/76_error.model").getFile());
        FileInputStream istream= null;
        try {
            istream= new FileInputStream(eFile);
            QualityErrorModel errorModel = MarkovErrorModel.loadErrorModel(eFile.getName(), istream);
            ModelPool babes = new ModelPool(true, errorModel, errorModel.getReadLength());

            // do it
            ByteArrayCharSequence cs= new ByteArrayCharSequence(10);
            sequencer.createQSeq(cs,
                    obj,
                    allChars.length(),  // 3p of tx
                    (byte) 1,  // tx dir
                    allChars.length(), // read length
                    allChars.length(), // fragment length
                    babes);
            String[] scs= cs.toString().split(barna.commons.system.OSChecker.NEW_LINE);
            Assert.assertTrue(allChars.equals(scs[0]));

        } catch (Exception e) {
            e.printStackTrace();
        }


    }

    /**
     * Creates a gene, transcript and exon objects according to the parameters provided.
     * @param chr name of the chromosome
     * @param txID transcript identifier
     * @param eStarts array with exon start coordinates, synchronized with end coordinates
     * @param eEnds array with exon end coordinates, synchronized with start coordinates
     * @return transcript instance with the gene model
     */
    protected static Transcript getModel(String chr, String txID, boolean sense, int[] eStarts, int[] eEnds) {
        // create model
        Gene g= new Gene("myGene");
        g.setChromosome(chr);
        Transcript t= new Transcript(g, txID);
        t.setStrand((byte) (sense? 1: -1));
        t.setSource(".");
        for (int i= 0; i< eStarts.length; ++i) {
            t.addExon(new Exon(t, "exon-"+ (i+1), eStarts[i], eEnds[i]));
        }

        return t;
    }

    @Test
    public void testSequenceRead() {
        Sequencer sequencer = new Sequencer(null, null);
        int nrTests= 100;
        int nrReadTests= 100;
        Random rnd= new Random();
        Graph.fileSep = null;

        for (int j= 0; j< nrTests; ++j) {

            int nrExons= rnd.nextInt(20)+ 1;
            int[] exonStarts= new int[nrExons], exonEnds= new int[nrExons];
            int s= rnd.nextInt(50)+ 1;
            for (int x= 0; x< nrExons; ++x) {
                exonStarts[x]= s;
                s+= 80+ rnd.nextInt(70);    // exon
                exonEnds[x]= s;
                s+= 1000+ rnd.nextInt(2000);    // intron
            }
            boolean sense= rnd.nextBoolean();
            String myTranscriptID= "myTranscript-"+ Integer.toString(j+ 1);

            // init
            Transcript t= getModel("", myTranscriptID, sense, exonStarts, exonEnds);
            int tlen= t.getExonicLength();  // do not add exon int coordinates, model can merge
//            int last3= -1;
//            for (int i= 0; i< t.getExons().length; ++i) {
//                Exon e= t.getExons()[i];
//                if (last3> 0)
//                    System.err.println("intron "+(e.get5PrimeEdge()- last3- 1));
//                System.err.println("exon "+ e.getLength());
//                last3= e.get3PrimeEdge();
//            }

            String chrSeq= getChrSequence(exonEnds[exonEnds.length - 1]);
            File f= writeSequence(chrSeq);
            String chr= FileHelper.stripExtension(f.getName());
            t.getGene().setChromosome(chr); // correct Chromosome to be found when reading sequence
            Graph.overrideSequenceDirPath= f.getParent();

            // test
            for (int i= 0; i< nrReadTests; ++i) {
                int polyA= sense? rnd.nextInt((int) (tlen/ 2d)): 0;    // TODO simulate poly-A only for sense reads
                int truLen= tlen+ polyA;
                int fragStart= rnd.nextInt(truLen)+ 1;
                int fragEnd= fragStart+ (truLen== fragStart? 0: rnd.nextInt(truLen- fragStart));
                int fragLen= fragEnd- fragStart+ 1;
                assert(fragStart<= fragEnd);
                boolean left= rnd.nextBoolean();
                int readLen= (fragStart== fragEnd? 1: rnd.nextInt(fragLen)+ 1);
                int readStart= (left? fragStart: fragEnd- readLen+ 1);
                int readEnd= (left? fragStart+ readLen- 1: fragEnd);
                assert(readStart<= readEnd);
                int molNr= rnd.nextInt(Integer.MAX_VALUE);
                int mate= rnd.nextInt(3);
                byte absDir= t.getStrand();
                if (!left)
                    absDir*= -1;    // anti-sense

                // paste HERE initialization of variables for debug

                // print settings for debugging randomly found fail
                BEDobject2 obj= new BEDobject2();
//                System.err.println("\n=== Simulation ===");
//                System.err.println("sense= "+ sense+ ";\n" +
//                        "chrSeq= getChrSequence("+ chrSeq.length()+");\n" +
//                        "f= writeSequence(chrSeq);\n" +
//                        "chr= FileHelper.stripExtension(f.getName());\n"+
//                        "Gene g= new Gene(\"myGene\");\n" +
//                        "g.setChromosome(chr);\n" +
//                        "t= new Transcript(g, \""+ t.getTranscriptID()+ "\");\n" +
//                        "t.setStrand((byte) (sense? 1: -1));\n" +
//                        "t.setSource(\".\");");
                for (int k = 0; k < t.getExons().length; k++) {
                    Exon ee= t.getExons()[k];
//                    System.err.println("t.addExon(new Exon(t, "+ ee.getExonID()+", "
//                            + Math.abs(ee.getStart())+", "+ Math.abs(ee.getEnd())+"));");

                }
//                System.err.println("tlen= t.getExonicLength();");
//                System.err.println("polyA= "+ polyA+ ";");
//                System.err.println("truLen= tlen+ polyA;");
//                System.err.println("readStart= "+ readStart+ ";");
//                System.err.println("readEnd= "+ readEnd+ ";");
//                System.err.println("readLen= readEnd- readStart+ 1;");
//                System.err.println("molNr= "+ molNr+ ";");
//                System.err.println("absDir= "+ absDir+ ";");
//                System.err.println("fragStart= "+ fragStart+ ";");
//                System.err.println("fragEnd= "+ fragEnd+ ";");
//                System.err.println("fragLen= fragEnd- fragStart+ 1;");
//                System.err.println("left= "+ left+ ";");
//                System.err.println("mate= "+ mate+ ";");
//                System.err.println();

                // create BED
                int polyAcnt= sequencer.createRead(
                        obj,
                        readStart- 1,  // read start in tx (0-based)
                        readEnd- 1, // read end in tx (0-based)
                        t,
                        molNr,       // molecule nr
                        absDir,   // absolute directionality
                        fragStart- 1,  // fragment start (0-based)
                        fragEnd- 1,    // fragment end (0-based)
                        left,   // left
                        mate       // mate ID
                );
                //System.err.println(obj);

                // test BED
                if (readStart> tlen) {
                    Assert.assertEquals(readLen, polyAcnt);
                    Assert.assertEquals("polyA", obj.getChr().toString());
                    Assert.assertEquals(0, obj.getStart());
                    Assert.assertEquals(readLen, obj.getEnd());
                } else {
                    Assert.assertEquals(Math.max(0, readStart- 1+ readLen- tlen), polyAcnt);
                    Assert.assertEquals(chr, obj.getChr().toString());
                    // don't check delimiter coordinates for polyA
                    if (readStart> 0&& readEnd< tlen) {
                        Assert.assertEquals(readStart - 1,
                                t.getExonicPosition(sense ? obj.getStart() + 1 : obj.getEnd()));
                        Assert.assertEquals(readEnd> tlen? readEnd: readEnd- 1,
                                t.getExonicPosition(sense? obj.getEnd(): obj.getStart()+ 1));
                    }
                }
                Assert.assertEquals(absDir, obj.getStrand());

                String[] ids= obj.getName().toString().split(":");
                Assert.assertEquals(8, ids.length);
                Assert.assertEquals(chr, ids[0]);
                Assert.assertEquals(Math.abs(t.getExons()[sense? 0: t.getExons().length- 1].getStart())+ "-"+
                        Math.abs(t.getExons()[sense? t.getExons().length- 1: 0].getEnd())+ (sense? "W": "C"), ids[1]);
                Assert.assertEquals(t.getTranscriptID(), ids[2]);
                Assert.assertEquals(Integer.toString(molNr+ 1), ids[3]);
                Assert.assertEquals(Integer.toString(tlen), ids[4]);
                Assert.assertEquals(Integer.toString(fragStart- 1), ids[5]);
                Assert.assertEquals(Integer.toString(fragEnd- 1), ids[6]);

                if (mate== 0) {
                    Assert.assertEquals(left? "S": "A", ids[7]);
                } else {
                    Assert.assertEquals((left? "S": "A")+ "/"+ Integer.toString(mate), ids[7]);
                }

                // create FASTA
                ByteArrayCharSequence cs= new ByteArrayCharSequence(200);
                sequencer.createQname(obj, cs, null);
                sequencer.createQSeq(cs, obj, t.get3PrimeEdge(), t.getStrand(), readLen, fragLen, null);

                // test FASTA
                //System.err.println(cs);
                String[] fasta= cs.toString().split(barna.commons.system.OSChecker.NEW_LINE);
                Assert.assertEquals(obj.getName().toString(), fasta[0].substring(1));
                Assert.assertEquals(readLen, fasta[1].length());
                if (readStart> tlen) {
                    for (int x= 0; x< readLen; ++x)
                        Assert.assertEquals(Character.toString(left? 'a': 't'), Character.toString(fasta[1].charAt(x)));
                } else {
                    // get transcript sequence
                    if (readEnd> tlen&& !sense)
                        continue;   // TODO cannot test anti-sense partial polyA, BED object misses info
                    String seq= "";
                    Exon[] ee= t.getExons();
                    for (int x= 0; x< ee.length; ++x) {
                        String ss= chrSeq.substring(Math.abs(ee[x].getStart())- 1, Math.abs(ee[x].getEnd()));
                        if (sense)
                            seq+= ss;
                        else
                            seq= ss+ seq;
                    }
                    //System.err.println(chrSeq);
                    //System.err.println(seq);
                    seq= (sense? seq: Graph.reverseSequence(Graph.complementarySequence(seq)));
                    //System.err.println(seq);
                    //System.err.println(readStart);
                    seq= seq.substring(readStart- 1);
                    //System.err.println(seq);
                    int x= 0;
                    for (;x< Math.min(readLen, tlen- readStart+ 1); ++x) {
                        if (left)
                            Assert.assertEquals(Character.toString(seq.charAt(x)),
                                    Character.toString(fasta[1].charAt(x)));
                        else
                            Assert.assertEquals(Character.toString(Graph.complementaryCharacter(seq.charAt(x))),
                                    Character.toString(fasta[1].charAt(readLen- 1- x)));
                    }
                    if (x < readLen) {  // read contains poly-A tail
                        Assert.assertEquals(tlen, readStart+ x- 1);   // check whether last informative position corresponds to the end of the transcript
                    }
                    for (; x< readLen; ++x)
                        if (left)
                            Assert.assertEquals('a', fasta[1].charAt(x));
                        else
                            Assert.assertEquals('t', fasta[1].charAt(readLen- 1- x));
                }

            }
        }
    }
    @Test
    public void testSequenceReadWithUniqueIds() {
        Sequencer sequencer = new Sequencer(null, null);

        BEDobject2 seq = new BEDobject2();
        Gene gene = new Gene("gene1");
        gene.setChromosome("chr1");

        Transcript transcript = new Transcript(gene, "trans-id");
        // single end reads
        sequencer.createRead(seq, 1, 100, transcript, 1, (byte) 1, 1000, 1100,  true, -1);
        assertEquals("chr1:0-0W:trans-id:2:0:1000:1100:S", seq.getName().toString());
        sequencer.createRead(seq, 1, 100, transcript, 1, (byte) -1, 1000, 1100,  false, -1);
        assertEquals("chr1:0-0W:trans-id:2:0:1000:1100:A", seq.getName().toString());
        //paired reads non unique ids
        sequencer.createRead(seq, 1, 100, transcript, 1, (byte) 1, 1000, 1100,  true, 1);
        assertEquals("chr1:0-0W:trans-id:2:0:1000:1100:S/1", seq.getName().toString());
        sequencer.createRead(seq, 1, 100, transcript, 1, (byte) -1, 1000, 1100,  true, 2);
        assertEquals("chr1:0-0W:trans-id:2:0:1000:1100:S/2", seq.getName().toString());
        sequencer.createRead(seq, 1, 100, transcript, 1, (byte) 1, 1000, 1100,  false, 1);
        assertEquals("chr1:0-0W:trans-id:2:0:1000:1100:A/1", seq.getName().toString());
        sequencer.createRead(seq, 1, 100, transcript, 1, (byte) -1, 1000, 1100,  false, 2);
        assertEquals("chr1:0-0W:trans-id:2:0:1000:1100:A/2", seq.getName().toString());
        //paired reads unique ids
        sequencer.setUniqueIds(true);
        sequencer.createRead(seq, 1, 100, transcript, 1, (byte) 1, 1000, 1100,  true, 1);
        assertEquals("chr1:0-0W:trans-id:2:0:1000:1100/1", seq.getName().toString());
        sequencer.createRead(seq, 1, 100, transcript, 1, (byte) -1, 1000, 1100,  true, 2);
        assertEquals("chr1:0-0W:trans-id:2:0:1000:1100/1", seq.getName().toString());
        sequencer.createRead(seq, 1, 100, transcript, 1, (byte) 1, 1000, 1100,  false, 1);
        assertEquals("chr1:0-0W:trans-id:2:0:1000:1100/2", seq.getName().toString());
        sequencer.createRead(seq, 1, 100, transcript, 1, (byte) -1, 1000, 1100,  false, 2);
        assertEquals("chr1:0-0W:trans-id:2:0:1000:1100/2", seq.getName().toString());
    }

     final char[] bases= new char[] {'A', 'C', 'G', 'T'};

    /**
     * Creates a random sequence for a chromosome of a certain length
     * @param maxLen maximum length, highest chromosomal coordinate
     * @return a random sequence for the chromosome
     */
     protected String getChrSequence(int maxLen) {
        char[] seq= new char[maxLen];
        Random rnd= new Random();
        for (int i= 0; i< maxLen; ++i) {
            seq[i]= bases[rnd.nextInt(bases.length)];
        }
        return new String(seq);
    }

    /**
     * Writes chromosome sequence in FASTA format to a temporary file.
     * @param chrSeq the chromosome sequence to be written
     * @return handle of the file to which sequence has been written
     */
     protected File writeSequence(String chrSeq) {
        // write sequence
        File f= null;
        try {
            f= File.createTempFile(SequencerTest.class.getSimpleName(), ".fa", currentTestDirectory);
            BufferedWriter writer= new BufferedWriter(new FileWriter(f));
            writer.write(">"+ FileHelper.stripExtension(f.getName())+ barna.commons.system.OSChecker.NEW_LINE);
            writer.write(chrSeq+ barna.commons.system.OSChecker.NEW_LINE);
            writer.close();
        } catch (Exception e) {
            e.printStackTrace();
        }

        return f;
    }
}
