package fbi.commons;

import org.junit.Test;

import static junit.framework.Assert.assertEquals;

/**
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class ByteArrayCharSequenceTest {

    @Test
    public void testTrim(){
        ByteArrayCharSequence ss = new ByteArrayCharSequence("   abc  ");
        ss.trim();
        assertEquals("abc", ss.toString());

        ss = new ByteArrayCharSequence("   abc  \t ");
        ss.trim();
        assertEquals("abc", ss.toString());
    }
}
