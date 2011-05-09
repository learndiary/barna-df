package fbi.commons.tools;

import java.lang.reflect.Array;
import java.util.*;


/**
 * Various array utilities.
 *
 * @author msammeth
 * @author Thasso Griebel (thasso.griebel@googlemail.com)
 */
public class ArrayUtils {

    /**
     * Searches the given element in the sorted list. If the element is not found, it is added
     * to the list. The given list MUST be sorted, otherwise results are not predictable.
     *
     * @param list    the list of elements
     * @param element the elemet
     * @param <T>     the type of elements in the list
     * @return inserted true if element was inserted successfully
     */
    public static <T> boolean addUniqueSorted(List<? extends Comparable<? super T>> list, T element) {
        return addUniqueSorted((List<T>) list, element, null);
    }

    /**
     * Searches for the given element in the sorted list. If the element is not found, it is added
     * to the list. The given list MUST be sorted, otherwise results are not predictable.
     *
     * @param list       the list
     * @param element    the element
     * @param comparator the comparator
     * @param <T>        the type of elements in the list
     * @return inserted true if element was inserted
     */
    public static <T> boolean addUniqueSorted(List<T> list, T element, Comparator<T> comparator) {
        int p = Collections.binarySearch(list, element, comparator);
        if (p < 0) {
            list.add(-(p + 1), element);
            return true;
        }
        return false;
    }

    /**
     * Searches for the given element in the sorted list. Natural order is used. The new element is added
     * to the list in the right position. If the given list is NOT empty, it MUST be sorted, otherwise
     * results are not predictable.
     *
     * @param list    the list
     * @param element the element
     * @param <T>     the type of elements in the list
     * @return inserted true if element was inserted
     */
    public static <T> boolean addSorted(List<? extends Comparable<? super T>> list, T element) {
        return addSorted((List<T>) list, element, null);
    }


    /**
     * Searches for the given element in the sorted list. The new element is added
     * to the list in the right position. If the given list is NOT empty, it MUST be sorted, otherwise
     * results are not predictable.
     *
     * @param list       the list
     * @param element    the element
     * @param comparator the comparator
     * @param <T>        the type of elements in the list
     * @return inserted true if element was inserted
     */
    public static <T> boolean addSorted(List<T> list, T element, Comparator<T> comparator) {
        int p = Collections.binarySearch(list, element, comparator);
        if (p < 0) {
            list.add(-(p + 1), element);
        } else {
            list.add(p, element);
        }
        return true;
    }


    /**
     * Adds all the elements to the given list. This assumes that the list is sorted. The results of applying
     * this method to an unsorted list are unpredictable. If the comparator is null, the element type must be comparable.
     *
     * @param list       the target list
     * @param elements   the elements
     * @param comparator the comparator (null permitted)
     * @param <T>        the type of the elements
     * @return success true if ALL elements were add successfully
     */
    public static <T> boolean addAllUniqueSorted(List<T> list, T[] elements, Comparator<T> comparator) {
        if (elements == null || elements.length == 0) {
            return true;
        }
        boolean success = true;
        for (T element : elements) {
            success &= addUniqueSorted(list, element, comparator);
        }
        return success;
    }

    /**
     * Adds all the elements to the given list. This assumes that the list is sorted. The results of applying
     * this method to an unsorted list are unpredictable. If the comparator is null, the element type must be comparable.
     *
     * @param list       the target list
     * @param elements   the elements
     * @param comparator the comparator (null permitted)
     * @param <T>        the type of the elements
     * @return success true if all elements were successfully added
     */
    public static <T> boolean addAllUniqueSorted(List<T> list, Iterable<T> elements, Comparator<T> comparator) {
        if (elements == null) {
            return true;
        }
        boolean success = true;
        for (T t : elements) {
            success &= addUniqueSorted(list, t, comparator);
        }
        return success;
    }

    /**
     * Adds the given element to the list if the list does NOT contain the element yet.
     *
     * @param list    the target list
     * @param element the element
     * @param <T>     the type of the elements in the list
     * @return success true if element was added successfully
     */
    public static <T> boolean addUnique(Collection<T> list, T element) {
        return list != null && !list.contains(element) && list.add(element);
    }


    /**
     * Inserts an object into an array at the specified position. Automatically
     * converts negative positions to insertion points as described in {@link Arrays#binarySearch(Object[], Object)} )}.
     * <p/>
     * Null is permitted for the array. In this case the array will be created.
     * <p/>
     * Note that this method always create a new array with {@code length == length + 1} and
     * copies the data.
     *
     * @param array    the target array, if null a new array will be created
     * @param element  the element to add
     * @param position the position to insert to. Negative values will translated to insertion points as described in {@link Arrays#binarySearch(Object[], Object)}
     * @return array the array
     * @throws ClassCastException in case the array and the element type do not match up
     */
    public static Object[] insert(Object[] array, Object element, int position) {
        /*
        Create array if not exist
         */
        if (array == null) {
            Object[] newA = (Object[]) Array.newInstance(element.getClass(), 1);
            newA[0] = element;
            return newA;
        }

        if (position < 0) {
            position = -(position + 1);
        }

        if (position > array.length) {
            position = array.length;
        }

        // get the array type
        Class arrayClass = array.getClass().getComponentType();
        // create a new array
        Object[] newA = (Object[]) Array.newInstance(arrayClass, array.length + 1);


        // COPY
        // copy the front part
        if (position > 0) {
            System.arraycopy(array, 0, newA, 0, position);
        }
        // set the element
        newA[position] = element;
        // copy the end part
        if (position < newA.length - 1) {
            System.arraycopy(array, position, newA, position + 1, array.length - position);
        }

        return newA;
    }

    /**
     * Add the element to the given array by extending the array
     *
     * @param array   the array
     * @param element the element
     * @return array the extended array
     */
    public static Object[] add(Object[] array, Object element) {
        if (array == null) {
            return insert(array, element, 0);
        } else {
            return insert(array, element, array.length);
        }
    }


    /**
     * Converts an eventually highdimensional Array or Vector to an Array.
     *
     * @param base the source
     * @return array resulting object
     */
    public static Object toField(Object base) {

        if (!(base instanceof Collection) && !(base instanceof Object[])) {
            return base;
        }

        Object[] o = null;
        if (base instanceof Collection) {
            o = ((Collection) base).toArray();
        } else                         // (base instanceof Object[])
        {
            o = ((Object[]) base);
        }

        if (o.length < 1)    // empty array
        {
            return null; //cannot guess the basic class for vectors anyway
        }

        Object r = null;    // find reference
        //int x= 0;
        //while (r== null&& x< o.length) {
        for (int x = 0; x < o.length; x++) {
            Object ro = toField(o[x]);
            if (ro == null) {
                continue;
            }
            if (r == null || ro.getClass().isAssignableFrom(r.getClass())) {
                r = ro;
            }
        }
        Object[] result = null;
        if (r != null) {
            result = (Object[]) Array.newInstance(r.getClass(), o.length);
            for (int i = 0; i < o.length; i++) {
                result[i] = toField(o[i]);
            }
        }

        return result;
    }
}
