package barna.commons;

import barna.commons.utils.StringUtils;
import org.junit.Test;

import java.util.Arrays;
import java.util.Collections;

import static junit.framework.Assert.assertEquals;

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
    }
    @Test
    public void testReverseComplement(){
        assertEquals("yYxXmMkKnacgt-NACGT", StringUtils.reverseComplement("ACGTN-acgtnMmKkXxRr"));
    }

    @Test
    public void testFPrint(){
        assertEquals("10.00000", StringUtils.fprint(10, 5));
        assertEquals("10.30000", StringUtils.fprint(10.3, 5));
        assertEquals("10.00000", StringUtils.fprint(10.00000000000000000000000000003, 5));
    }

    @Test
    public void testAppend() throws Exception {
        assertEquals("aAAA", StringUtils.append('A', "a", 4, false));
        assertEquals("AAAa", StringUtils.append('A', "a", 4, true));
    }

    @Test
    public void testJoin() throws Exception {
        assertEquals("aAAA", StringUtils.join("", Arrays.asList("a","A", "A",'A')));
        assertEquals("a A A A", StringUtils.join(" ", Arrays.asList("a","A", "A",'A')));
        assertEquals("", StringUtils.join(" ", Collections.EMPTY_LIST));
    }
}
