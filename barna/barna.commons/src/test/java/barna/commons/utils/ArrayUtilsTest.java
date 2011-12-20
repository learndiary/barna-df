package barna.commons.utils;

import org.junit.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static junit.framework.Assert.*;

/**
 * Array utils tests
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class ArrayUtilsTest {

    @Test
    public void testAddSorted(){
        ArrayList<Integer> list = new ArrayList<Integer>(Arrays.asList(1, 2, 5, 6));
        assertTrue(ArrayUtils.addSorted(list, 7));
        assertEquals(5, list.size());
        assertEquals(Arrays.asList(1,2,5,6,7),list);

        assertTrue(ArrayUtils.addSorted(list, 7));
        assertEquals(6, list.size());
        assertEquals(Arrays.asList(1,2,5,6,7, 7),list);

        assertTrue(ArrayUtils.addSorted(list, 0));
        assertEquals(7, list.size());
        assertEquals(Arrays.asList(0,1,2,5,6,7,7),list);

        assertTrue(ArrayUtils.addSorted(list, 3));
        assertEquals(8, list.size());
        assertEquals(Arrays.asList(0,1,2,3,5,6,7,7),list);
    }

    @Test
    public void testAddUniqueSorted(){
        ArrayList<Integer> list = new ArrayList<Integer>(Arrays.asList(1, 2, 5, 6));
        assertTrue(ArrayUtils.addUniqueSorted(list, 7));
        assertEquals(5, list.size());
        assertEquals(Arrays.asList(1,2,5,6,7),list);

        assertFalse(ArrayUtils.addUniqueSorted(list, 7));
        assertEquals(5, list.size());
        assertEquals(Arrays.asList(1,2,5,6,7),list);

        assertTrue(ArrayUtils.addUniqueSorted(list, 0));
        assertEquals(6, list.size());
        assertEquals(Arrays.asList(0,1,2,5,6,7),list);

        assertTrue(ArrayUtils.addUniqueSorted(list, 3));
        assertEquals(7, list.size());
        assertEquals(Arrays.asList(0,1,2,3,5,6,7),list);
    }


    @Test
    public void testAddAllUniqueSortedCollection(){
        ArrayList<Integer> list = new ArrayList<Integer>(Arrays.asList(1, 2, 5, 6));
        assertTrue(ArrayUtils.addAllUniqueSorted(list, Arrays.asList(0, 3, 7), null));
        assertEquals(7, list.size());
        assertEquals(Arrays.asList(0,1,2,3,5,6,7),list);

        assertFalse(ArrayUtils.addAllUniqueSorted(list, Arrays.asList(0, 8), null));
        assertEquals(8, list.size());
        assertEquals(Arrays.asList(0,1,2,3,5,6,7,8),list);
    }

    @Test
    public void testAddAllUniqueSortedArray(){
        ArrayList<Integer> list = new ArrayList<Integer>(Arrays.asList(1, 2, 5, 6));
        assertTrue(ArrayUtils.addAllUniqueSorted(list, new Integer[]{0, 3, 7}, null));
        assertEquals(7, list.size());
        assertEquals(Arrays.asList(0,1,2,3,5,6,7),list);

        assertFalse(ArrayUtils.addAllUniqueSorted(list, new Integer[]{0, 8}, null));
        assertEquals(8, list.size());
        assertEquals(Arrays.asList(0, 1, 2, 3, 5, 6, 7, 8), list);
    }

    @Test
    public void testAddUnique(){
        ArrayList<Integer> list = new ArrayList<Integer>(Arrays.asList(1, 2, 5, 6));
        assertTrue(ArrayUtils.addUnique(list, 7));
        assertEquals(5, list.size());
        assertEquals(Arrays.asList(1,2,5,6,7),list);

        assertFalse(ArrayUtils.addUnique(list, 7));
        assertEquals(5, list.size());
        assertEquals(Arrays.asList(1,2,5,6,7),list);

        assertTrue(ArrayUtils.addUnique(list, 0));
        assertEquals(6, list.size());
        assertEquals(Arrays.asList(1,2,5,6,7,0),list);
    }

    @Test
    public void testInsert(){
        Object[] a = ArrayUtils.insert(null, 0, 0);
        assertTrue(Arrays.equals(a, new Integer[]{0}));

        Integer[] testArray = {0,1,2};
        Object[] inserted = ArrayUtils.insert(testArray, 3, 0);
        assertTrue(Arrays.equals(inserted, new Integer[]{3,0,1,2}));

        inserted = ArrayUtils.insert(testArray, 3, 1);
        assertTrue(Arrays.equals(inserted, new Integer[]{0,3,1,2}));

        inserted = ArrayUtils.insert(testArray, 3, 2);
        assertTrue(Arrays.equals(inserted, new Integer[]{0,1,3,2}));

        inserted = ArrayUtils.insert(testArray, 3, 3);
        assertTrue(Arrays.equals(inserted, new Integer[]{0,1,2,3}));

        inserted = ArrayUtils.insert(testArray, 3, 4);
        assertTrue(Arrays.equals(inserted, new Integer[]{0,1,2,3}));



        inserted = ArrayUtils.insert(testArray, 3, -1);
        assertTrue(Arrays.equals(inserted, new Integer[]{3,0,1,2}));

        inserted = ArrayUtils.insert(testArray, 3, -2);
        assertTrue(Arrays.equals(inserted, new Integer[]{0,3,1,2}));

        inserted = ArrayUtils.insert(testArray, 3, -3);
        assertTrue(Arrays.equals(inserted, new Integer[]{0,1,3,2}));

        inserted = ArrayUtils.insert(testArray, 3, -4);
        assertTrue(Arrays.equals(inserted, new Integer[]{0,1,2,3}));

        inserted = ArrayUtils.insert(testArray, 3, -5);
        assertTrue(Arrays.equals(inserted, new Integer[]{0,1,2,3}));

    }
    @Test
    public void testAdd(){
        Object[] a = ArrayUtils.add(null, new Integer(0));
        assertTrue(Arrays.equals(a, new Integer[]{0}));

        Integer[] testArray = {0,1,2};
        Object[] inserted = ArrayUtils.add(testArray, 3);
        assertTrue(Arrays.equals(inserted, new Integer[]{0,1,2, 3}));
    }


    @Test
    public void testToField(){
        Object arrayarray = ArrayUtils.toField(new int[][]{{1,2,3}, {1,2,3}});
        int[][] ia = (int[][]) arrayarray;
        assertTrue(ia.length == 2);
        assertEquals(ia[0][0],1);
        assertEquals(ia[0][1], 2);
        assertEquals(ia[0][2],3);
        assertEquals(ia[1][0],1);
        assertEquals(ia[1][1], 2);
        assertEquals(ia[1][2],3);




        Object listlist = ArrayUtils.toField(Arrays.asList(Arrays.asList(1,2,3), Arrays.asList(1,2,3)));
        Integer[][] il = (Integer[][]) listlist;
        assertTrue(il.length == 2);
        assertEquals(il[0][0],new Integer(1));
        assertEquals(il[0][1],new Integer(2));
        assertEquals(il[0][2],new Integer(3));
        assertEquals(il[1][0],new Integer(1));
        assertEquals(il[1][1],new Integer(2));
        assertEquals(il[1][2],new Integer(3));

        Object arraylist = ArrayUtils.toField(new List[]{Arrays.asList(1,2,3), Arrays.asList(1,2,3)});
        Integer[][] al = (Integer[][]) listlist;
        assertTrue(al.length == 2);
        assertEquals(al[0][0],new Integer(1));
        assertEquals(al[0][1],new Integer(2));
        assertEquals(al[0][2],new Integer(3));
        assertEquals(al[1][0],new Integer(1));
        assertEquals(al[1][1],new Integer(2));
        assertEquals(al[1][2],new Integer(3));

    }


}
