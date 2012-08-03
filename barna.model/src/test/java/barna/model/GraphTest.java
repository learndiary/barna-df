package barna.model;

import barna.commons.ByteArrayCharSequence;
import barna.commons.system.OSChecker;
import junit.framework.Assert;
import org.junit.Test;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.Random;

import static org.junit.Assert.assertEquals;
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


    private String writeTmpChromosome(String seq, File tmpFile, String newline) throws Exception {
        String chr= tmpFile.getName().substring(0, tmpFile.getName().lastIndexOf('.'));
        BufferedWriter writer= new BufferedWriter(new FileWriter(tmpFile));
        writer.write(">"+ chr+ newline);
        writer.write(seq+ newline);
        writer.close();
        return chr;
    }

    @Test
    public void testReadChromosome() {
        File f = null;
        try {
            String seq = BASIC_SEQUENCE;
            f = File.createTempFile(getClass().getSimpleName(), ".fa");
            String chr = writeTmpChromosome(seq, f, OSChecker.NEW_LINE);
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


    @Test
    public void testReadWithinLineSetup() throws Exception {

        String line_unix = "ABCDE\nFGHIJ\nKLMNO";
        String line_win = "ABCDE\r\nFGHIJ\r\nKLMNO";

        // check full length
        int seqLen_u = getSeqLength(line_unix);
        int seqLen_w = getSeqLength(line_win);
        assertEquals(15, seqLen_u);
        assertEquals(15, seqLen_w);

        // check single line length
        int lineLen_u = line_unix.indexOf("\n");
        int lineLen_w = line_win.indexOf("\r\n");
        assertEquals(5, lineLen_u);
        assertEquals(5, lineLen_w);

        // check number of lines
        int nrLines_u = (int) Math.ceil(seqLen_u / (lineLen_u + 1d));
        int nrLines_w = (int) Math.ceil(seqLen_w / (lineLen_w + 1d));
        assertEquals(3, nrLines_u);
        assertEquals(3, nrLines_w);

        int lineNr= 1;
        int posStart = 2;
        int posEnd = 6;

        int offCR_u = lineNr * (lineLen_u+1);
        int off_u= lineNr* lineLen_u;

        File chr_file_unix = File.createTempFile(getClass().getSimpleName(), ".fa");
        String chr_unix = writeTmpChromosome(line_unix, chr_file_unix, "\n");
        File chr_file_win = File.createTempFile(getClass().getSimpleName(), ".fa");
        String chr_win = writeTmpChromosome(line_win, chr_file_win, "\r\n");

        Graph.overrideSequenceDirPath = chr_file_unix.getParentFile().getAbsolutePath();
        Graph.fileSep = "\n";
        ByteArrayCharSequence result_unix = new ByteArrayCharSequence(10);
        Graph.readSequence(
                null,   // Species spe,
                chr_unix,    // CharSequence chromosome,
                true,   // boolean forwardStrand,
                off_u + posStart + 1,   // start
                off_u + posStart + posEnd,  // end
                result_unix,     // ByteArrayCharSequence cs,
                0,      // int from,
                posEnd     // int to)
        );
        assertEquals("HIJKLM", result_unix.toString());

        Graph.overrideSequenceDirPath = chr_file_win.getParentFile().getAbsolutePath();
        Graph.fileSep = "\r\n";
        ByteArrayCharSequence result_win = new ByteArrayCharSequence(10);
        Graph.readSequence(
                null,   // Species spe,
                chr_win,    // CharSequence chromosome,
                true,   // boolean forwardStrand,
                off_u + posStart + 1,   // start
                off_u + posStart + posEnd,  // end
                result_win,     // ByteArrayCharSequence cs,
                0,      // int from,
                posEnd     // int to)
        );
        assertEquals("HIJKLM", result_win.toString());



    }

    private void testReadWithinLine(Random rnd, CharSequence chr, String seq, ByteArrayCharSequence cs) {

        int seqLen = getSeqLength(seq);

        int lineLen= seq.indexOf(OSChecker.NEW_LINE);
        int nrLines= (int) Math.ceil(seqLen/ (lineLen+ OSChecker.NEW_LINE.length()));
        int lineNr= Math.max(1,rnd.nextInt(nrLines));
        int posStart= rnd.nextInt(lineLen);
        int posEnd= rnd.nextInt(lineLen- posStart);

        int offCR= lineNr * (lineLen+ OSChecker.NEW_LINE.length());
        int off= lineNr * lineLen;
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


    /**
     * Get the number of characters without newline characters
     *
     * @param seq the sequence
     * @return length length if the sequence without newline characters
     */
    private int getSeqLength(String seq) {
        int seqLen= 0;
        for (int i = 0; i < seq.length(); i++) {
            if (seq.charAt(i)!= '\n' && seq.charAt(i) != '\r')
                ++seqLen;
        }
        return seqLen;
    }

    private void testReadMultiLines(Random rnd, CharSequence chr, String seq, ByteArrayCharSequence cs) {

        int seqLen= 0;
        for (int i = 0; i < seq.length(); i++) {
            if (seq.charAt(i)!= '\n' || seq.charAt(i) != '\r')
                ++seqLen;
        }

        int lineLen= seq.indexOf(barna.commons.system.OSChecker.NEW_LINE);
        int nrLines= (int) Math.ceil(seqLen/ (double) (lineLen+ OSChecker.NEW_LINE.length()));
        int lineNr= rnd.nextInt(nrLines- 1);
        int endLineNr= lineNr+ 1+ rnd.nextInt(nrLines- lineNr- 1);
        int posStart= rnd.nextInt(lineLen);
        int posEnd= rnd.nextInt(lineLen);

        int offCR= lineNr* (lineLen+ OSChecker.NEW_LINE.length());
        int off= lineNr* lineLen;
        int off2CR= endLineNr* (lineLen+ OSChecker.NEW_LINE.length());
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

        int seqLen= getSeqLength(seq);

        int lineLen= seq.indexOf(barna.commons.system.OSChecker.NEW_LINE);
        int nrLines= (int) Math.ceil(seq.length()/ (double) (lineLen+ OSChecker.NEW_LINE.length()));
        int endLineNr= 1+ rnd.nextInt(nrLines- 1);
        int lineNr= rnd.nextInt(endLineNr);
        int posStart= rnd.nextInt(lineLen);
        int posEnd= rnd.nextInt(lineLen);

        int offCR= lineNr* (lineLen+ OSChecker.NEW_LINE.length());
        int off= lineNr* lineLen;
        int off2CR= endLineNr* (lineLen+ OSChecker.NEW_LINE.length());
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
            if (sb.charAt(i)== '\n' || sb.charAt(i)== '\r')
                sb.deleteCharAt(i--);
        }
        //System.err.println(sb.toString().length()+ " - "+cs.toString().length());
        Assert.assertEquals(sb.toString(), cs.toString());

    }
}