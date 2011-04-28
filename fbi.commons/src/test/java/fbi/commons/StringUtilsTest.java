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

    @Test
    public void testFPrint(){
        assertEquals("10.00000", StringUtils.fprint(10, 5));
        assertEquals("10.30000", StringUtils.fprint(10.3, 5));
        assertEquals("10.00000", StringUtils.fprint(10.00000000000000000000000000003, 5));
    }
}
