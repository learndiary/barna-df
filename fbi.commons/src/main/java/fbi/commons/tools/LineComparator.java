package fbi.commons.tools;

import java.util.Comparator;

/**
 *
 * Compares strings and supports field splitting and number parsing.
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
class LineComparator implements Comparator<String>{
    /**
     * The field to use for comparison
     */
    private int field = 0;
    /**
     * Parse the field to a numerical value
     */
    private boolean numerical;
    /**
     * Field separator
     */
    private String separator = "\\t";

    /**
     * Create a new line comparator. If the given field is {@code < 0}, the line is not splitted. If numerical
     * is set, the field (or the complete line) is parsed to a number. The separator is used to split
     * the line. It is used as regular expression.
     *
     * @param field the field (split using th separator as regular expression)
     * @param numerical numerical comparison
     * @param separator the separator used to split fields
     */
    LineComparator(int field, boolean numerical, String separator) {
        this.field = field;
        this.numerical = numerical;
        this.separator = separator;
    }

    public int compare(String o1, String o2) {
        if(field >= 0){
            // split fields
            o1 = o1.split(separator)[field];
            o2 = o2.split(separator)[field];
        }

        // trim it
        o1 = o1.trim();
        o2 = o2.trim();

        // no field sep
        if(numerical){
            try{
                // try integer first
                return Integer.parseInt(o1) - Integer.parseInt(o2);
            }catch (NumberFormatException e){
                // ok integer failed ... lets try double
                return Double.compare(Double.parseDouble(o1),Double.parseDouble(o2));
            }
        }
        return o1.compareTo(o2);
    }
}
