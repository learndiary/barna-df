package barna.commons.thread;

import barna.commons.ByteArrayCharSequence;
import org.junit.Test;

/**
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class SyncIOHandler2Test {

    static final String TEST_STRING = "abc123\n321cba";

    @Test
    public void testThreadAndSimpleSameReadLine(){
        java.io.ByteArrayInputStream in = new java.io.ByteArrayInputStream(TEST_STRING.getBytes());
        SyncIOHandler2 h = new SyncIOHandler2(2);


        h.addStream(in);
        ByteArrayCharSequence cs = new ByteArrayCharSequence(128);
        int i = h.readLine(in, cs);
        i = h.readLine(in, cs);
        System.out.println(cs.toString());


    }

}
