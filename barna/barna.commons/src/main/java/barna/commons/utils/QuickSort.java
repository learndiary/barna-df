/*
 * This file is part of the Flux Library.
 *
 * The code of the Flux Library may be freely distributed and modified under the terms of the
 * European Union Public Licence (EUPL) published on the web site <http://www.osor.eu/eupl/european-union-public-licence-eupl-v.1.1>.
 * Copyright for the code is held jointly by the individual authors, who should be listed
 * in @author doc comments. According to Article 5 and Article 11 of the EUPL, publications that
 * include results produced by the Flux Library are liable to reference the Work,
 * see the Flux Library homepage <http://flux.sammeth.net> for more information.
 */

package barna.commons.utils;

import java.util.Comparator;
import java.util.List;

/**
 * Simple quicksort implementation
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class QuickSort<T extends Comparable> {

    /**
     * Sort the array
     *
     * @param a array the array
     */
    public static void sort(Comparable[] a) {
        sort(a, 0, a.length - 1);
    }

    /**
     * Sort the list
     *
     * @param a list the list to sort
     */
    public static <T extends Comparable> void sort(List<T> a) {
        sort(a, 0, a.size() - 1);
    }

    /**
     * Sort the list
     *
     * @param a list the list to sort
     * @param comparator the comparator
     */
    public static <T> void sort(List<T> a, Comparator<T> comparator) {
        sort(a, 0, a.size() - 1, comparator);
    }


    /**
     * Quicksort the sub array from lo to hi
     *
     * @param a elements
     * @param lo lo index
     * @param hi high index
     */
    private static void sort(Comparable[] a, int lo, int hi) {
        if (hi <= lo) return;
        int j = partition(a, lo, hi);
        sort(a, lo, j-1);
        sort(a, j+1, hi);
    }
    /**
     * Quicksort the sub list from lo to hi
     *
     * @param a elements
     * @param lo lo index
     * @param hi high index
     */
    private static <T extends Comparable> void sort(List<T> a, int lo, int hi) {
        if (hi <= lo) return;
        int j = partition(a, lo, hi);
        sort(a, lo, j-1);
        sort(a, j+1, hi);
    }

    /**
     * Quicksort the sub list from lo to hi
     *
     * @param a elements
     * @param lo lo index
     * @param hi high index
     */
    private static <T> void sort(List<T> a, int lo, int hi, Comparator<T> comparator) {
        if (hi <= lo) return;
        int j = partition(a, lo, hi, comparator);
        sort(a, lo, j-1, comparator);
        sort(a, j+1, hi, comparator);
    }

    /**
     * Find and return the index j in the range from lo to hi such
     * that a[lo .. j-1] <= a[j] <= a[j+1 .. hi]
     *
     * @param a array to sort
     * @param lo lo index
     * @param hi high index
     */
    private static int partition(Comparable[] a, int lo, int hi) {
        int i = lo;
        int j = hi + 1;
        //Comparable v = a[lo];
        Comparable v = a[(hi-lo)/2];
        while (true) {
            // find item on lo to swap
            while (less(a[++i], v))
                if (i == hi) break;
            // find item on hi to swap
            while (less(v, a[--j]))
                if (j == lo) break;      // redundant since a[lo] acts as sentinel
            // check if pointers cross
            if (i >= j) break;
            exch(a, i, j);
        }
        // put v = a[j] into position
        exch(a, lo, j);
        // with a[lo .. j-1] <= a[j] <= a[j+1 .. hi]
        return j;
    }
    /**
     * Find and return the index j in the range from lo to hi such
     * that a[lo .. j-1] <= a[j] <= a[j+1 .. hi]
     *
     * @param a array to sort
     * @param lo lo index
     * @param hi high index
     */
    private static <T extends Comparable> int partition(List<T> a, int lo, int hi) {
        int i = lo;
        int j = hi + 1;
        //T v = a.get(lo);
        //use middle elemetn as pivot
        T v = a.get((hi-lo)/2);
        while (true) {
            // find item on lo to swap
            while (less(a.get(++i), v))
                if (i == hi) break;
            // find item on hi to swap
            while (less(v, a.get(--j)))
                if (j == lo) break;      // redundant since a[lo] acts as sentinel
            // check if pointers cross
            if (i >= j) break;
            exch(a, i, j);
        }
        // put v = a[j] into position
        exch(a, lo, j);
        // with a[lo .. j-1] <= a[j] <= a[j+1 .. hi]
        return j;
    }
    /**
     * Find and return the index j in the range from lo to hi such
     * that a[lo .. j-1] <= a[j] <= a[j+1 .. hi]
     *
     * @param a array to sort
     * @param lo lo index
     * @param hi high index
     */
    private static <T> int partition(List<T> a, int lo, int hi, Comparator<T> comparator) {
        int i = lo;
        int j = hi + 1;
        T v = a.get(lo);
        while (true) {
            // find item on lo to swap
            while (less(a.get(++i), v, comparator))
                if (i == hi) break;
            // find item on hi to swap
            while (less(v, a.get(--j), comparator))
                if (j == lo) break;      // redundant since a[lo] acts as sentinel
            // check if pointers cross
            if (i >= j) break;
            exch(a, i, j);
        }
        // put v = a[j] into position
        exch(a, lo, j);
        // with a[lo .. j-1] <= a[j] <= a[j+1 .. hi]
        return j;
    }



   /***********************************************************************
    *  Helper sorting functions
    ***********************************************************************/

    /**
     * Return true if is {@code v < w}
     *
     * @param v element
     * @param w element
     * @return less true if {@code v < w}
     */
    private static boolean less(Comparable v, Comparable w) {
        return (v.compareTo(w) < 0);
    }
    /**
     * Return true if is {@code v < w}
     *
     * @param v element
     * @param w element
     * @return less true if {@code v < w}
     */
    private static <T> boolean less(T v, T w, Comparator<T> comparator) {
        return (comparator.compare(v, w) < 0);
    }

    /**
     * Switch elements at i and j
     *
     * @param a the array
     * @param i index
     * @param j index
     */
    private static void exch(Object[] a, int i, int j) {
        Object swap = a[i];
        a[i] = a[j];
        a[j] = swap;
    }
    /**
     * Switch elements at i and j
     *
     * @param a the array
     * @param i index
     * @param j index
     */
    private static void exch(List a, int i, int j) {
        Object swap = a.get(i);
        a.set(i, a.get(j));
        a.set(j,swap);
    }

  // test client
    public static void main(String[] args) {

        // generate array of N random reals between 0 and 1
        int N = 100;
        Double[] a = new Double[N];
        for (int i = 0; i < N; i++) {
            a[i] = Math.random();
        }

        // sort the array
        QuickSort.sort(a);

        // display results
        for (int i = 0; i < N; i++) {
            System.out.println(a[i]);
        }
    }


}
