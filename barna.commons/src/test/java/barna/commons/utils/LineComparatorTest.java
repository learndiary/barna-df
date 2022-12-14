package barna.commons.utils;

import org.junit.Test;

import java.util.Arrays;

import static junit.framework.Assert.assertEquals;
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
    public void testFieldMergeComparison(){
        String l1 = "A\tC";
        String l2 = "A\tB";
        assertTrue(new LineComparator("\\t", 0, 1).compare(l1, l2) > 0);
        assertTrue(new LineComparator("\\t", 0).compare(l1, l2) == 0);
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

    @Test
    public void testComposites() {
        String s1 = "1\t0";
        String s2 = "1\t1";
        
        LineComparator<CharSequence> c1= new LineComparator<CharSequence>(true, "\t", 0);
        LineComparator<CharSequence> c2= new LineComparator<CharSequence>(true, "\t", 1);
        LineComparator<CharSequence> c3= new LineComparator<CharSequence>(c1)
        	.addComparator(new LineComparator<CharSequence>(c2));
        
        assertEquals(0, c1.compare(s1, s2));
        assertEquals(-1, c2.compare(s1, s2));
        assertEquals(1, c2.compare(s2, s1));
        
        assertEquals(-1, c3.compare(s1, s2));
        assertEquals(1, c3.compare(s2, s1));
    }

    @Test
    public void testCustomSplits() throws Exception {
        {
            String s = "A\tB\t\tC";
            LineComparator lc = new LineComparator("\t", 0, 1, 2, 3);
            System.out.println(Arrays.toString(lc.split(s)));
            assertTrue(Arrays.deepEquals(new String[]{"A", "B", "", "C"}, lc.split(s) ));
        }
        {
            String s = "\tA\tB\t\tC";
            LineComparator lc = new LineComparator("\t", 0, 1, 2, 3,4);
            System.out.println(Arrays.toString(lc.split(s)));
            assertTrue(Arrays.deepEquals(new String[]{"","A", "B", "", "C"}, lc.split(s) ));
        }
        {
            String s = "\tA\tB\t\tC\t";
            LineComparator lc = new LineComparator("\t", 0, 1, 2, 3,4,5);
            System.out.println(Arrays.toString(lc.split(s)));
            assertTrue(Arrays.deepEquals(new String[]{"","A", "B", "", "C", ""}, lc.split(s) ));
        }
        {
            String s = "\tA\tB\t\tC\t";
            LineComparator lc = new LineComparator("\t", 0, 1, 2);
            System.out.println(Arrays.toString(lc.split(s)));
            assertTrue(Arrays.deepEquals(new String[]{"","A", "B"}, lc.split(s) ));
        }
        {
            String s = "\tA\tB\t\tC\t";
            LineComparator lc = new LineComparator("\t", 2, 3, 4);
            System.out.println(Arrays.toString(lc.split(s)));
            assertTrue(Arrays.deepEquals(new String[]{"", "A","B","", "C"}, lc.split(s) ));
        }

    }
}
