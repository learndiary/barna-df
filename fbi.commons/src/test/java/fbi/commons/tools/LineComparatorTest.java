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
        assertTrue(new LineComparator(-1, false, "\\t").compare(l1, l2) < 0);
        // with numbers
        assertTrue(new LineComparator(-1, true, "\\t").compare(l1, l2) == 0);
    }


    @Test
    public void testFieldComparison(){
        String l1 = "A\tB";
        String l2 = "B\tA";
        assertTrue(new LineComparator(0, false, "\\t").compare(l1, l2) < 0);
        assertTrue(new LineComparator(1, false, "\\t").compare(l1, l2) > 0);
    }

    @Test
    public void testSeparatorComparison(){
        String l1 = "AnnnB";
        String l2 = "BnnnA";
        assertTrue(new LineComparator(0, false, "nnn").compare(l1, l2) < 0);
        assertTrue(new LineComparator(1, false, "nnn").compare(l1, l2) > 0);
    }


}
