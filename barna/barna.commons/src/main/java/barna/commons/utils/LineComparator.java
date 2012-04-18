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

import java.util.*;

/**
 * Compares strings and supports field splitting and number parsing. The implementations default is NOT to
 * cache any values, but if you intend to use this in scenarios where the field values will be accessed
 * multiple times, you might want to provide a cache using {@link #setCache(java.util.Map)}.
 * <p>
 *     Also note that you can speed up the line splitting gif your separator is only a single character. In
 *     that case a manual split is performed and regexp usage is avoided.
 * </p>
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class LineComparator<T extends CharSequence> implements Comparator<T> {
    /**
     * Split types
     */
    private static enum SplitType {
        Regexp, Character
    }

    /**
     * The fields to use for comparison
     */
    private int[] field = new int[]{0};
    /**
     * Parse the field to a numerical value
     */
    private final boolean numerical;
    /**
     * Field separator
     */
    private String separator = "\t";

    /**
     * How to split the string
     */
    private SplitType splitType = null;

    /**
     * Possible list of sub comparators that are used if
     * this comparators result is 0
     */
    private List<Comparator<T>> subComparators;
    /**
     * Optional parent comparator
     */
    private Comparator<T> delegate;

    /**
     * Optional cache that caches field values
     */
    private Map<T, Object> cache;


    /**
     * Copy constructor
     *
     * @param copy the source
     */
    public LineComparator(LineComparator<T> copy) {
        super();
        this.field = new int[copy.field.length];
        System.arraycopy(copy.field, 0, this.field, 0, copy.field.length);
        this.numerical = copy.numerical;

        this.separator = copy.separator;
        if (copy.subComparators != null) {
            this.subComparators = new ArrayList<Comparator<T>>();
            for (Comparator<? extends CharSequence> cc : copy.subComparators) {
                if (cc instanceof LineComparator) {
                    this.subComparators.add(new LineComparator((LineComparator<T>) cc));
                } else {
                    this.subComparators.add((Comparator<T>) cc);
                }
            }
        }
        this.delegate = copy.delegate;
        if(copy.getCache() != null){
            setCache(new HashMap<T, Object>());
        }
    }

    /**
     * Create a new line comparator. If the given field is {@code < 0}, the line is not splitted. If numerical
     * is set, the field (or the complete line) is parsed to a number. The separator is used to split
     * the line. It is used as regular expression.
     *
     * @param numerical numerical comparison
     * @param separator the separator used to split fields
     * @param field     the field (split using th separator as regular expression)
     */
    public LineComparator(boolean numerical, String separator, int field) {
        this.field = new int[]{field};
        this.numerical = numerical;
        this.separator = separator;
    }

    /**
     * Create a new line comparator. The given fields are merged and treated as one single concatenated field
     * in comparison
     *
     * @param separator the separator used to split fields
     * @param fields    the fields (split using th separator as regular expression)
     */
    public LineComparator(String separator, int... fields) {
        this.field = fields;
        this.numerical = false;
        this.separator = separator;
    }

    /**
     * Create a new line comparator. The given fields are merged and treated as one single concatenated field
     * in comparison
     *
     * @param separator the separator used to split fields
     * @param fields    the fields (split using th separator as regular expression)
     */
    public LineComparator(String separator, Integer[] fields) {
        this.field = new int[fields.length];
        for (int i = 0; i < fields.length; i++) {
            this.field[i] = fields[i];
        }
        this.numerical = false;
        this.separator = separator;
    }

    /**
     * Create a line comparator stub that uses a specific separator
     * and wraps around another comparator to which comparison
     * calls are delegated.
     *
     * @param comparator parent comparator
     */
    public LineComparator(Comparator<T> comparator) {
        this.delegate = comparator;
        this.numerical = false; // dummy init to make compiler happy
    }


    /**
     * Add a sub comparator that is used if this comparator would return 0
     *
     * @param comparator the sub comparator
     * @return <code>this</code> instance with the new child comparator
     */
    public LineComparator<T> addComparator(Comparator<T> comparator) {
        if (comparator != null) {
            if (subComparators == null) {
                subComparators = new ArrayList<Comparator<T>>();
            }
            subComparators.add(comparator);
        }
        return this;
    }

    /**
     * Get the cached value for s
     *
     * @param s the key
     * @return value the cached value or null
     */
    public Object get(T s) {
        if (cache != null) {
            return cache.get(s);
        }
        return null;
    }

    /**
     * Put a value in the cache
     *
     * @param s the key
     * @param o the value
     */
    public void put(T s, Object o) {
        if (cache != null) {
            cache.put(s, o);
        }
    }

    public int compare(final T o1, final T o2) {

        // check for empty string
        if (o1.length() == 0 || o2.length() == 0) {
            return o1.toString().compareTo(o2.toString());
        }
        int result = 0;
        if (delegate == null) {
            String s1 = o1.toString();
            String s2 = o2.toString();

            Object c1 = get(o1);
            Object c2 = get(o2);

            boolean gotResult = false;
            if (numerical && c1 != null && c2 != null) {
                result = Double.compare((Double) c1, (Double) c2);
                gotResult = true;
            } else {
                if (c1 != null) s1 = c1.toString();
                if (c2 != null) s2 = c2.toString();
            }


            if (!gotResult) {
                // one field specified
                if (field.length == 1 && field[0] >= 0) {
                    // split fields
                    if (c1 == null) {
                        s1 = split(o1.toString())[field[0]];//o1.toString().split(separator)[field[0]];
                        put(o1, s1);
                    }
                    if (c2 == null) {
                        s2 = split(o2.toString())[field[0]];//o2.toString().split(separator)[field[0]];
                        put(o2, s2);
                    }
                } else if (field.length > 1) {
                    // merge multiple fields
                    // merge o1 fields
                    if (c1 == null) {
                        String[] o1_split = split(o1.toString());//o1.toString().split(separator);
                        for (int i : field) {
                            s1 += o1_split[i];
                        }
                        put(o1, s1);
                    }

                    if (c2 == null) {
                        // merge o2 fields
                        String[] o2_split = split(o2.toString());//o2.toString().split(separator);
                        for (int i : field) {
                            s2 += o2_split[i];
                        }
                        put(o2, s2);
                    }
                }


                if (numerical) {
                    try {
                        // try integer first
                        int i1 = Integer.parseInt(s1);
                        int i2 = Integer.parseInt(s2);
                        put(o1, new Double(i1));
                        put(o2, new Double(i2));
                        result = i1 - i2;
                    } catch (NumberFormatException e) {
                        // ok integer failed ... lets try double
                        double d1 = Double.parseDouble(s1);
                        double d2 = Double.parseDouble(s2);
                        put(o1, new Double(d1));
                        put(o2, new Double(d2));
                        result = Double.compare(d1, d2);
                    }
                } else {
                    result = s1.compareTo(s2);
                }
            }
        } else {
            result = delegate.compare(o1, o2);
        }
        /*
        If result is 0 check the sub comparators
         */
        if (result == 0 && subComparators != null) {
            for (Comparator<T> subComparator : subComparators) {
                int sub = subComparator.compare(o1, o2);
                if (sub != 0) {
                    return sub;
                }
            }
        }
        return result;
    }

    public void reset() {
        if (cache != null) {
            cache.clear();
        }
        if (subComparators != null) {
            for (Comparator<? extends CharSequence> subComparator : subComparators) {
                if (subComparator instanceof LineComparator) {
                    ((LineComparator) subComparator).reset();
                }
            }
        }
    }

    /**
     * Returns all subcomparators currently registered for this
     * <code>LineComparator</code>.
     *
     * @return all subcomparators currently registered
     */
    public List<Comparator<T>> getSubComparators() {
        return subComparators;
    }


    /**
     * Split the given string using a manual character comparison split
     * or a regexp, based on the separator
     *
     * @param string the string to split
     * @return splits the splitted string
     */
    protected String[] split(String string) {
        if (splitType == null) {
            // find out how to split
            if (separator.length() == 1) {
                splitType = SplitType.Character;
            } else {
                splitType = SplitType.Regexp;
            }
        }

        String[] split = null;
        if (splitType == SplitType.Character) {
            ArrayList<String> splits = new ArrayList<String>();
            // max field to get
            int max = -1;
            for (int i : field) {
                max = Math.max(i + 1, max);
            }

            char sc = separator.charAt(0);
            int s = 0;
            int field = 0;
            for (int i = 0; i < string.length(); i++) {
                if (string.charAt(i) == sc) {
                    // split
                    if (s == i) splits.add("");
                    else splits.add(string.substring(s, i));
                    s = i + 1;
                    field++;
                    if (splits.size() >= max || field >= max) break;
                }
            }
            // add final
            if (splits.size() < max) {
                splits.add(string.substring(s));
            }
            split = splits.toArray(new String[splits.size()]);
        }else{
            split = string.split(separator);
        }

        return split;
    }

    /**
     * Get the current cache or null
     * @return cache cache or null
     */
    public Map<T, Object> getCache() {
        return cache;
    }

    /**
     * Set the cache. The implementation does NOT cache by default, you
     * have to explicitly provide a cache here to enable field value caching
     *
     * @param cache the cache
     */
    public void setCache(Map<T, Object> cache) {
        this.cache = cache;
    }
}

