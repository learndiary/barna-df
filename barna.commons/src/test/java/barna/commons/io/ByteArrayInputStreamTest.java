package barna.commons.io;

import barna.commons.ByteArrayCharSequence;
import barna.commons.system.OSChecker;
import org.junit.Test;

import java.io.IOException;

import static junit.framework.Assert.assertEquals;
import static junit.framework.Assert.fail;

/**
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class ByteArrayInputStreamTest {
    static final String TEST_STRING = "abc123"+ OSChecker.NEW_LINE+"321cba";

    @Test
    public void testRead(){
        java.io.ByteArrayInputStream in = new java.io.ByteArrayInputStream(TEST_STRING.getBytes());
        ByteArrayCharSequence seq = new ByteArrayCharSequence(16);
        ByteArrayInputStream stream = new ByteArrayInputStream(seq, in);

        try {
            int read = stream.read();
            assertEquals('a', read);
            read = stream.read();
            assertEquals('b', read);
            read = stream.read();
            assertEquals('c', read);
            read = stream.read();
            assertEquals('1', read);
            read = stream.read();
            assertEquals('2', read);
            read = stream.read();
            assertEquals('3', read);
            read = stream.read();
            assertEquals(OSChecker.NEW_LINE, read);
            read = stream.read();
            assertEquals('3', read);
            read = stream.read();
            assertEquals('2', read);
            read = stream.read();
            assertEquals('1', read);
            read = stream.read();
            assertEquals('c', read);
            read = stream.read();
            assertEquals('b', read);
            read = stream.read();
            assertEquals('a', read);

            read = stream.read();
            assertEquals(-1, read);
        } catch (IOException e) {
            e.printStackTrace();
            fail();
        }
    }


    @Test
    public void testReadLine(){
        java.io.ByteArrayInputStream in = new java.io.ByteArrayInputStream(TEST_STRING.getBytes());
        ByteArrayCharSequence seq = new ByteArrayCharSequence(16);
        ByteArrayInputStream stream = new ByteArrayInputStream(seq, in);

        try {
            int l = stream.readLine();
            assertEquals("abc123", seq.toString());
            assertEquals(6, l);
            l = stream.readLine();
            assertEquals("abc123321cba", seq.toString());
            assertEquals(6, l);
            l = stream.readLine();
            assertEquals(-1, l);



        } catch (IOException e) {
            e.printStackTrace();
            fail();
        }
    }


}
