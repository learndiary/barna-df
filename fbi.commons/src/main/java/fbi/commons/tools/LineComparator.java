package fbi.commons.tools;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

/**
 * Compares strings and supports field splitting and number parsing.
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
class LineComparator implements Comparator<String> {
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
     * Create a new line comparator. If the given field is {@code < 0}, the line is not splitted. If numerical
     * is set, the field (or the complete line) is parsed to a number. The separator is used to split
     * the line. It is used as regular expression.
     *
     * @param numerical numerical comparison
     * @param separator the separator used to split fields
     * @param field     the field (split using th separator as regular expression)
     */
    LineComparator(boolean numerical, String separator, int field) {
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
    LineComparator(String separator, int... fields) {
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
    LineComparator(String separator, Integer[] fields) {
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
    LineComparator(Comparator<String> comparator) {
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

            // one field specified
            if (field.length == 1 && field[0] >= 0) {
                // split fields
                s1 = o1.split(separator)[field[0]];
                s2 = o2.split(separator)[field[0]];
            } else if (field.length > 1) {
                // merge multiple fields
                // merge o1 fields
                String[] o1_split = o1.split(separator);
                for (int i : field) {
                    s1 += o1_split[i];
                }
                // merge o2 fields
                String[] o2_split = o2.split(separator);
                for (int i : field) {
                    s2 += o2_split[i];
                }
            }

            if (numerical) {
                try {
                    // try integer first
                    result = Integer.parseInt(s1) - Integer.parseInt(s2);
                } catch (NumberFormatException e) {
                    // ok integer failed ... lets try double
                    result = Double.compare(Double.parseDouble(s1), Double.parseDouble(s2));
                }
            } else {
                result = s1.compareTo(s2);
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
}

