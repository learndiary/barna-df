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
import com.google.common.io.Resources;
import org.junit.Assert;
import org.junit.Test;

import java.io.*;
import java.util.Random;

import static junit.framework.Assert.*;

public class SequencerTest {

    @Test
	public void testReadGeneration() throws Exception {
		
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
						int polyA= Sequencer.createRead(obj, start, end, t, k,
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

        String allChars= "ACGUTWSRYMKBDHV";

        // write sequence
        File f= null;
        try {
            f= File.createTempFile(this.getClass().getSimpleName(), ".fa");
            BufferedWriter writer= new BufferedWriter(new FileWriter(f));
            writer.write(">"+ FileHelper.stripExtension(f.getName())+ "\n");
            writer.write(allChars+ "\n");
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
            ModelPool babes = new ModelPool(true, errorModel);

            // do it
            ByteArrayCharSequence cs= new ByteArrayCharSequence(10);
            Sequencer.createQSeq(cs,
                    obj,
                    allChars.length(),  // 3p of tx
                    (byte) 1,  // tx dir
                    allChars.length(), // read length
                    allChars.length(), // fragment length
                    babes);
            String[] scs= cs.toString().split("\n");
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

        int nrTests= 100;
        int nrReadTests= 100;
        Random rnd= new Random();

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
            boolean sense= true; // rnd.nextBoolean();  // TODO debug for Crick strand
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
                int truLen= tlen+ rnd.nextInt((int) (tlen/ 2d));
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

                // create BED
                BEDobject2 obj= new BEDobject2();
                if (readStart- 1> tlen- 1)
                    System.currentTimeMillis();
                int polyAcnt= Sequencer.createRead(
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
                    Assert.assertEquals(readStart- 1,
                            t.getExonicPosition(sense? obj.getStart()+ 1: obj.getEnd()));
                    Assert.assertEquals(readEnd> tlen? readEnd: readEnd- 1,
                            t.getExonicPosition(sense? obj.getEnd(): obj.getStart()+ 1));
                }
                Assert.assertEquals(absDir, obj.getStrand());

                String[] ids= obj.getName().toString().split(":");
                Assert.assertEquals(8, ids.length);
                Assert.assertEquals(chr, ids[0]);
                Assert.assertEquals(exonStarts[0]+ "-"+ exonEnds[exonEnds.length- 1]+ (sense? "W": "C"), ids[1]);
                Assert.assertEquals(myTranscriptID, ids[2]);
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
                Sequencer.createQname(obj, cs, null);
                Sequencer.createQSeq(cs, obj, t.get3PrimeEdge(), t.getStrand(), readLen, fragLen, null);

                // test FASTA
                //System.err.println(cs);
                String[] fasta= cs.toString().split("\n");
                Assert.assertEquals(obj.getName().toString(), fasta[0].substring(1));
                Assert.assertEquals(readLen, fasta[1].length());
                if (readStart> tlen) {
                    for (int x= 0; x< readLen; ++x)
                        Assert.assertEquals(left? 'a': 't', fasta[1].charAt(x));
                } else {
                    // get transcript sequence
                    String seq= "";
                    for (int x= 0; x< t.getExons().length; ++x)
                        seq+= chrSeq.substring(Math.abs(t.getExons()[x].getStart())- 1, Math.abs(t.getExons()[x].getEnd()));
                    seq= (sense? seq: Graph.reverseSequence(Graph.complementarySequence(seq)));

                    seq= seq.substring(readStart- 1);
                    //System.err.println(chrSeq);
                    //System.err.println(seq);
                    int x= 0;
                    for (;x< Math.min(readLen, tlen- readStart+ 1); ++x) {
                        //System.err.println(x);
                        if (left)
                            Assert.assertEquals(Character.toString(seq.charAt(x)),
                                    Character.toString(fasta[1].charAt(x)));
                        else
                            Assert.assertEquals(Character.toString(Graph.complementaryCharacter(seq.charAt(x))),
                                    Character.toString(fasta[1].charAt(readLen- 1- x)));
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

    static final char[] bases= new char[] {'A', 'C', 'G', 'T'};

    /**
     * Creates a random sequence for a chromosome of a certain length
     * @param maxLen maximum length, highest chromosomal coordinate
     * @return a random sequence for the chromosome
     */
    static protected String getChrSequence(int maxLen) {
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
    static protected File writeSequence(String chrSeq) {
        // write sequence
        File f= null;
        try {
            f= File.createTempFile(SequencerTest.class.getSimpleName(), ".fa");
            BufferedWriter writer= new BufferedWriter(new FileWriter(f));
            writer.write(">"+ FileHelper.stripExtension(f.getName())+ "\n");
            writer.write(chrSeq+ "\n");
            writer.close();
        } catch (Exception e) {
            e.printStackTrace();
        }

        return f;
    }
}
