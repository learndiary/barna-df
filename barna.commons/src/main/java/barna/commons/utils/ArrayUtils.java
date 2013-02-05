/*
 * Copyright (c) 2010, Micha Sammeth
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * The names of its contributors may be not used to endorse or promote
 *       products derived from this software without specific prior written
 *       permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL MICHA SAMMETH BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

package barna.commons.utils;

import java.lang.reflect.Array;
import java.lang.reflect.Method;
import java.util.*;

/**
 * Various array utilities.
 *
 * @author Micha Sammeth (gmicha@googlemail.com)
 * @author Thasso Griebel (thasso.griebel@googlemail.com)
 */
public class ArrayUtils {



    /**
     * Creates a copy of an integer array.
     * @param a The array to be cloned
     * @return The clone of the provided array
     */
    public static int[] duplicate(int[] a) {
        if (a== null)
            return null;
        int[] d= new int[a.length];
        for (int i = 0; i < d.length; i++)
            d[i]= a[i];
        return d;
    }

    /**
     * Creates a copy of a double array.
     * @param a The array to be cloned
     * @return The clone of the provided array
     */
    public static double[] duplicate(double[] a) {
        if (a== null)
            return null;
        double[] d= new double[a.length];
        for (int i = 0; i < d.length; i++)
            d[i]= a[i];
        return d;
    }

    /**
     * Creates a copy of an Object array.
     * @param a The array to be cloned
     * @return The clone of the provided array
     */
    public static Object[] duplicate(Object[] a) {
        if (a== null)
            return null;
        if (a.length== 0)
            return new Object[0];
        Object[] d= (Object[]) Array.newInstance(a[0].getClass(), a.length);
        for (int i = 0; i < d.length; i++)
            d[i]= a[i];
        return d;
    }

    /**
     * Converts arrays of the primitive types <code>int</code> or <code>double</code>
     * into arrays of the corresponding wrapper classes, <code>Integer</code> respectively
     * <code>Double</code>.
     * @param inA the array to be converted
     * @return The converted array with objects of the respective wrapper class representing
     * the provided primitive values.
     */
    public static Object[] primitiveToWrapperFieldDistinguishable(Object inA) {
        if (inA instanceof int[]) {
            int[] in= (int[]) inA;
            Integer[] out= new Integer[in.length];
            for (int i = 0; i < out.length; i++)
                out[i]= new Integer(in[i]);
            return out;
        } else if (inA instanceof double[]) {
            double[] in= (double[]) inA;
            Double[] out= new Double[in.length];
            for (int i = 0; i < out.length; i++)
                out[i]= new Double(in[i]);
            return out;
        }
        return (Object[]) inA;
    }

    /**
     * Sort two vectors according to the natural ordering of the first one.
     * @param primSort primary vector that determines the ordering
     * @param restSort secondary vector that gets permuted according to the permutations
     *                 necessary to order the elements of the primary array
     */
    public static void synchroneousSort(Object primSort, Vector restSort) {
        Object[] primO= primitiveToWrapperFieldDistinguishable(primSort);
        if (primO== null|| primO.length< 2)
            return;

        HashMap<Object,Integer> refMap= new HashMap<Object,Integer>(primO.length);
        for (int i = 0; i < primO.length; i++)
            refMap.put(primO[i], new Integer(i));

        java.util.Arrays.sort(primO);

        // sort others
        for (int j = 0; j < restSort.size(); j++) {
            if (restSort.elementAt(j) instanceof int[]) {
                int[] array= (int[]) restSort.elementAt(j);
                int[] arrayOld= duplicate(array);
                for (int i = 0; i < arrayOld.length; i++)
                    array[i]= arrayOld[refMap.get(primO[i]).intValue()];
            } else if (restSort.elementAt(j) instanceof double[]) {
                double[] array= (double[]) restSort.elementAt(j);
                double[] arrayOld= duplicate(array);
                for (int i = 0; i < arrayOld.length; i++)
                    array[i]= arrayOld[refMap.get(primO[i]).intValue()];
            } else {
                Object[] array= (Object[]) restSort.elementAt(j);
                Object[] arrayOld= duplicate(array);
                for (int i = 0; i < arrayOld.length; i++)
                    array[i]= arrayOld[refMap.get(primO[i]).intValue()];
            }
        }

        // convert prim sort
        if (primSort instanceof int[]) {
            try {
                Method m= primO[0].getClass().getMethod("intValue", null);
                int[] out= (int[]) primSort;
                for (int i = 0; i < out.length; i++)
                    out[i]= ((java.lang.Integer) m.invoke(primO[i], null)).intValue();
            } catch (Exception e) {
                e.printStackTrace();
            }
        } else if (primSort instanceof double[]) {
            try {
                Method m= primO[0].getClass().getMethod("doubleValue", null);
                double[] out= (double[]) primSort;
                for (int i = 0; i < out.length; i++)
                    out[i]= ((java.lang.Double) m.invoke(primO[i], null)).doubleValue();
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
    }

    /**
     * Wrapper to sort two integer arrays according to the natural ordering of the first
     * one.
     * @param primSort primary array of integer that determines the ordering
     * @param restSort secondary array that gets permuted according to the permutations
     *                 necessary to order the elements of the primary array
     */
    public static void synchroneousSort(int[] primSort, int[] restSort) {
        Vector v= new Vector();
        v.add(restSort);
        synchroneousSort(primSort, v);
    }


    /**
     * Converts values returned by search methods to indices where to insert
     * a corresponding value.
     * @param p
     * @return The index as is, if positive, otherwise the insertion point converted to
     * a positive index.
     */
    public static int convertInsertionPoint(int p) {
        if (p< 0)
            p= (p+1)* (-1);
        return p;
    }


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
