package fbi.commons;

import org.junit.Test;

import static junit.framework.Assert.assertEquals;
import static junit.framework.Assert.fail;

/**
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class StringUtilsTest {

    @Test
    public void testReverse(){
        assertEquals("CBA", StringUtils.reverse("ABC"));
    }

    @Test
    public void testComplement(){
        assertEquals("TGCAN-tgcanKkMmXxYy", StringUtils.complement("ACGTN-acgtnMmKkXxRr"));

        try{
            StringUtils.complement("B");
            fail();
        }catch (RuntimeException e){
            // ignore
        }
    }
    @Test
    public void testReverseComplement(){
        assertEquals("yYxXmMkKnacgt-NACGT", StringUtils.reverseComplement("ACGTN-acgtnMmKkXxRr"));

        try{
            StringUtils.complement("B");
            fail();
        }catch (RuntimeException e){
            // ignore
        }
    }
}
