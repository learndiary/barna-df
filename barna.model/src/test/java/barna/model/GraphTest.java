package barna.model;

import barna.commons.ByteArrayCharSequence;
import barna.commons.system.OSChecker;
import junit.framework.Assert;
import org.junit.Test;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.Random;

import static org.junit.Assert.fail;


/**
 * Created with IntelliJ IDEA.
 * User: micha
 * Date: 7/16/12
 * Time: 2:55 PM
 * To change this template use File | Settings | File Templates.
 */
public class GraphTest {

    public static String BASIC_SEQUENCE=
            "The two of them keep their eyes on each other...........  "+OSChecker.NEW_LINE+
            "She sits down.  They sit around her.  Nick sits directly  "+OSChecker.NEW_LINE+
            "across from her.  She lights up a cigarette.  They watch  "+OSChecker.NEW_LINE+
            "her.  She is poised, cool, in complete command of herself."+OSChecker.NEW_LINE+
            "------- There is no smoking in this building, Ms. Tramell."+OSChecker.NEW_LINE+
            "------- What are you going to do?  Charge me with smoking?"+OSChecker.NEW_LINE+
            "Ever so casually, she blows her smoke across at Nick......";


    private String writeTmpChromosome(String seq, File tmpFile) throws Exception {
        String chr= tmpFile.getName().substring(0, tmpFile.getName().lastIndexOf('.'));
        BufferedWriter writer= new BufferedWriter(new FileWriter(tmpFile));
        writer.write(">"+ chr+ barna.commons.system.OSChecker.NEW_LINE);
        writer.write(seq+ barna.commons.system.OSChecker.NEW_LINE);
        writer.close();
        return chr;
    }

    @Test
    public void testReadChromosome() {
        File f = null;
        try {
            String seq = BASIC_SEQUENCE;
            f = File.createTempFile(getClass().getSimpleName(), ".fa");
            String chr = writeTmpChromosome(seq, f);
            Graph.overrideSequenceDirPath = f.getParentFile().getAbsolutePath();

            // tests
            Random rnd = new Random();
            ByteArrayCharSequence cs = new ByteArrayCharSequence(seq.length());
            for (int i = 0; i < 100; i++) {
                cs.end = 0;
                testReadWithinLine(rnd, chr, seq, cs);
                cs.end = 0;
                testReadMultiLines(rnd, chr, seq, cs);
                cs.end = 0;
                testReadOverflow(rnd, chr, seq, cs);
            }

        } catch (Exception e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            fail();
        } finally {
            if (f != null && f.exists()) {
                f.delete();
            }
        }


    }

    private void testReadWithinLine(Random rnd, CharSequence chr, String seq, ByteArrayCharSequence cs) {

        int seqLen= 0;
        for (int i = 0; i < seq.length(); i++) {
            if (seq.charAt(i)!= '\n' && seq.charAt(i) != '\r')
                ++seqLen;
        }

        int lineLen= seq.indexOf(OSChecker.NEW_LINE);
        int nrLines= (int) Math.ceil(seqLen/ (lineLen+ 1d));
        int lineNr= rnd.nextInt(nrLines);
        int posStart= rnd.nextInt(lineLen);
        int posEnd= rnd.nextInt(lineLen- posStart);

        int offCR= lineNr* (lineLen+ 1);
        int off= lineNr* lineLen;
        Graph.readSequence(
                null,   // Species spe,
                chr,    // CharSequence chromosome,
                true,   // boolean forwardStrand,
                off + posStart + 1,   // start
                off + posStart + posEnd,  // end
                cs,     // ByteArrayCharSequence cs,
                0,      // int from,
                posEnd     // int to)
        );
        Assert.assertEquals(seq.substring(offCR+ posStart, offCR+ posStart+ posEnd), cs.toString());
    }

    private void testReadMultiLines(Random rnd, CharSequence chr, String seq, ByteArrayCharSequence cs) {

        int seqLen= 0;
        for (int i = 0; i < seq.length(); i++) {
            if (seq.charAt(i)!= '\n' || seq.charAt(i) != '\r')
                ++seqLen;
        }

        int lineLen= seq.indexOf(barna.commons.system.OSChecker.NEW_LINE);
        int nrLines= (int) Math.ceil(seqLen/ (double) (lineLen+ 1));
        int lineNr= rnd.nextInt(nrLines- 1);
        int endLineNr= lineNr+ 1+ rnd.nextInt(nrLines- lineNr- 1);
        int posStart= rnd.nextInt(lineLen);
        int posEnd= rnd.nextInt(lineLen);

        int offCR= lineNr* (lineLen+ 1);
        int off= lineNr* lineLen;
        int off2CR= endLineNr* (lineLen+ 1);
        int off2= endLineNr* lineLen;

//        System.err.println("seq.length()= "+ seq.length()+ ", lineLen= "+ lineLen+", nrLines= "+ nrLines+ ", lineNr= "+ lineNr+ ", endLineNr= "+ endLineNr+ ", posStart= "+ posStart+ ", posEnd= "+posEnd+
//                ", off= "+ off+ ", off2= "+ off2);
        Graph.readSequence(
                null,   // Species spe,
                chr,    // CharSequence chromosome,
                true,   // boolean forwardStrand,
                off + posStart + 1,   // start
                off2 + posEnd,  // end
                cs,     // ByteArrayCharSequence cs,
                0,      // int from,
                (lineLen - posStart) + posEnd + (Math.max(0, (endLineNr - lineNr + 1) - 2) * lineLen)     // int to)
        );
//        System.err.println(cs.toString());
//        System.err.println(seq.substring(offCR+ posStart, offCR+ posStart+ posEnd));
        StringBuffer sb= new StringBuffer(seq.substring(offCR + posStart, off2CR + posEnd));
        for (int i = 0; i < sb.length(); i++) {
            if (sb.charAt(i)== '\n' || sb.charAt(i)== '\r')
                sb.deleteCharAt(i--);
        }

        Assert.assertEquals(sb.toString(), cs.toString());

    }

    private void testReadOverflow(Random rnd, CharSequence chr, String seq, ByteArrayCharSequence cs) {

        int seqLen= 0;
        for (int i = 0; i < seq.length(); i++) {
             if (seq.charAt(i)!= '\n')
                 ++seqLen;
        }

        int lineLen= seq.indexOf(barna.commons.system.OSChecker.NEW_LINE);
        int nrLines= (int) Math.ceil(seq.length()/ (double) (lineLen+ 1));
        int endLineNr= 1+ rnd.nextInt(nrLines- 1);
        int lineNr= rnd.nextInt(endLineNr);
        int posStart= rnd.nextInt(lineLen);
        int posEnd= rnd.nextInt(lineLen);

        int offCR= lineNr* (lineLen+ 1);
        int off= lineNr* lineLen;
        int off2CR= endLineNr* (lineLen+ 1);
        int off2= endLineNr* lineLen;

        //System.err.println("seq.length()= "+ seq.length()+ ", lineLen= "+ lineLen+", nrLines= "+ nrLines+ ", lineNr= "+ lineNr+ ", endLineNr= "+ endLineNr+ ", posStart= "+ posStart+ ", posEnd= "+posEnd+
        //        ", off= "+ off+ ", off2= "+ off2);
        Graph.readSequence(
                null,   // Species spe,
                chr,    // CharSequence chromosome,
                true,   // boolean forwardStrand,
                off2 + posEnd+ 1,   // start
                seqLen + (off + posStart),  // end
                cs,     // ByteArrayCharSequence cs,
                0,      // int from,
                seqLen- ((lineLen- posStart)+ posEnd+ (Math.max(0, (endLineNr- lineNr+ 1)- 2)* lineLen))     // int to)
        );
//        System.err.println(cs.toString());
//        System.err.println(seq.substring(offCR+ posStart, offCR+ posStart+ posEnd));
        StringBuffer sb= new StringBuffer(seq.substring(off2CR+ posEnd)+ seq.substring(0, offCR + posStart));
        for (int i = 0; i < sb.length(); i++) {
            if (sb.charAt(i)== '\n')
                sb.deleteCharAt(i--);
        }
        //System.err.println(sb.toString().length()+ " - "+cs.toString().length());
        Assert.assertEquals(sb.toString(), cs.toString());

    }
}