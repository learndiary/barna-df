package barna.commons.io;

import barna.commons.ByteArrayCharSequence;
import org.junit.Test;

import java.io.IOException;

import static junit.framework.Assert.*;

/**
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class ThreadedIOHandlerTest {

    static final String TEST_STRING = "abc123\n321cba";

    @Test
    public void testThreadAndSimpleSameReadLine(){
        java.io.ByteArrayInputStream in_s = new java.io.ByteArrayInputStream(TEST_STRING.getBytes());
        java.io.ByteArrayInputStream in_t = new java.io.ByteArrayInputStream(TEST_STRING.getBytes());
        SimpleIOHandler simple = new SimpleIOHandler();
        ThreadedIOHandler threaded = new ThreadedIOHandler(1);

        simple.addStream(in_s);
        threaded.addStream(in_t);

        try {
            ByteArrayCharSequence read_s = simple.readLine(in_s);
            ByteArrayCharSequence read_t = threaded.readLine(in_t);
            assertEquals(read_s.toString(), read_t.toString());

            read_s = simple.readLine(in_s);
            read_t = threaded.readLine(in_t);
            assertEquals(read_s.toString(), read_t.toString());

            read_s = simple.readLine(in_s);
            read_t = threaded.readLine(in_t);

            assertNull(read_s);
            assertNull(read_t);


        } catch (IOException e) {
            e.printStackTrace();
            fail();
        }


    }
}
