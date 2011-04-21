package fbi.commons.tools;

import org.junit.Test;

import static junit.framework.Assert.assertTrue;

/**
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class LineComparatorTest {

    @Test
    public void testNumberComparison(){
        String l1 = "01";
        String l2 = "1";
        // no numbers
        assertTrue(new LineComparator(false, "\\t", -1).compare(l1, l2) < 0);
        // with numbers
        assertTrue(new LineComparator(true, "\\t", -1).compare(l1, l2) == 0);
    }


    @Test
    public void testFieldComparison(){
        String l1 = "A\tB";
        String l2 = "B\tA";
        assertTrue(new LineComparator(false, "\\t", 0).compare(l1, l2) < 0);
        assertTrue(new LineComparator(false, "\\t", 1).compare(l1, l2) > 0);
    }

    @Test
    public void testSeparatorComparison(){
        String l1 = "AnnnB";
        String l2 = "BnnnA";
        assertTrue(new LineComparator(false, "nnn", 0).compare(l1, l2) < 0);
        assertTrue(new LineComparator(false, "nnn", 1).compare(l1, l2) > 0);
    }

    @Test
    public void testOrder(){
        String s1 = "HWUSI-EAS1678:1:1:999:1021";
        String s2 = "HWUSI-EAS1678:1:1:9999:3539";

        assertTrue(new LineComparator(false, "\t", -1).compare(s1, s2) > 0);

    }


}
