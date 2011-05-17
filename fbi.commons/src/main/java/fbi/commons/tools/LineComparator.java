package fbi.commons.tools;

import java.util.*;

/**
 * Compares strings and supports field splitting and number parsing.
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class LineComparator implements Comparator<String> {
    /**
     * The fields to use for comparison
     */
    private int[] field = new int[]{0};
    /**
     * Parse the field to a numerical value
     */
    private boolean numerical;
    /**
     * Field separator
     */
    private String separator = "\\t";

    /**
     * Possible list of sub comparators that are used if
     * this comparators result is 0
     */
    private List<Comparator<String>> subComparators;
    /**
     * Optional parent comparator
     */
    private Comparator<String> parent;

    /**
     * Cache separator splits
     */
    private Map<String, Object> splitCache = new HashMap<String, Object>();


    /**
     * Copy constructor
     *
     * @param copy the source
     */
    public LineComparator(LineComparator copy){
        super();
        this.field = new int[copy.field.length];
        System.arraycopy(copy.field, 0, this.field,0, copy.field.length);

        this.separator = copy.separator;
        if(copy.subComparators != null){
            this.subComparators = new ArrayList<Comparator<String>>();
            for (Comparator<String> cc : copy.subComparators) {
                if(cc instanceof LineComparator){
                    this.subComparators.add(new LineComparator((LineComparator)cc));
                }else{
                    this.subComparators.add(cc);
                }
            }
        }
        this.parent = copy.parent;
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
     * Create a line comparator that wraps around another comparator
     *
     * @param comparator parent comparator
     */
    public LineComparator(Comparator<String> comparator) {
        this.parent = comparator;
    }


    /**
     * Add a sub comparator that is used if this comparator would returns 0
     *
     * @param comparator the sub comparator
     */
    public void addComparator(Comparator<String> comparator) {
        if (comparator != null) {
            if (subComparators == null) {
                subComparators = new ArrayList<Comparator<String>>();
            }
            subComparators.add(comparator);
        }
    }


    public int compare(final String o1, final String o2) {
        // check for empty string
        if (o1.length() == 0 || o2.length() == 0) {
            return o1.compareTo(o2);
        }
        int result = 0;
        if (parent == null) {
            String s1 = o1;
            String s2 = o2;

            Object c1 = splitCache.get(o1);
            Object c2 = splitCache.get(o2);

            boolean gotResult = false;
            if(numerical && c1 != null && c2 != null){
                result = Double.compare((Double)c1, (Double)c2);
                gotResult = true;
            }else{
                if(c1 != null) s1 = c1.toString();
                if(c2 != null) s2 = c2.toString();
            }


            if(!gotResult){
                // one field specified
                if (field.length == 1 && field[0] >= 0) {
                    // split fields
                    if(c1 == null){
                        s1 = o1.split(separator)[field[0]];
                        splitCache.put(o1, s1);
                    }
                    if(c2 == null){
                        s2 = o2.split(separator)[field[0]];
                        splitCache.put(o2, s2);
                    }
                } else if (field.length > 1) {
                    // merge multiple fields
                    // merge o1 fields
                    if(c1 == null){
                        String[] o1_split = o1.split(separator);
                        for (int i : field) {
                            s1 += o1_split[i];
                        }
                        splitCache.put(o1, s1);
                    }

                    if(c2 == null){
                        // merge o2 fields
                        String[] o2_split = o2.split(separator);
                        for (int i : field) {
                            s2 += o2_split[i];
                        }
                        splitCache.put(o2, s2);
                    }
                }


                if (numerical) {
                    try {
                        // try integer first
                        int i1 = Integer.parseInt(s1);
                        int i2 = Integer.parseInt(s2);
                        splitCache.put(o1, new Double(i1));
                        splitCache.put(o2, new Double(i2));
                        result = i1 - i2;
                    } catch (NumberFormatException e) {
                        // ok integer failed ... lets try double
                        double d1 = Double.parseDouble(s1);
                        double d2 = Double.parseDouble(s2);
                        splitCache.put(o1, new Double(d1));
                        splitCache.put(o2, new Double(d2));
                        result = Double.compare(d1, d2);
                    }
                } else {
                    result = s1.compareTo(s2);
                }
            }
        } else {
            result = parent.compare(o1, o2);
        }
        /*
        If result is 0 check the sub comparators
         */
        if (result == 0 && subComparators != null) {
            for (Comparator<String> subComparator : subComparators) {
                int sub = subComparator.compare(o1, o2);
                if (sub != 0) {
                    return sub;
                }
            }
        }
        return result;
    }

    public void reset(){
        splitCache.clear();
        if(subComparators != null){
            for (Comparator<String> subComparator : subComparators) {
                if(subComparator instanceof LineComparator){
                    ((LineComparator)subComparator).reset();
                }
            }
        }
    }
}

