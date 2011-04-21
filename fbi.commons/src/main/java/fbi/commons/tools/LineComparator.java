package fbi.commons.tools;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

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
     * Possible list of sub comparators that are used if
     * this comparators result is 0
     */
    private List<Comparator<String>> subComparators;

    /**
     * Create a new line comparator. If the given field is {@code < 0}, the line is not splitted. If numerical
     * is set, the field (or the complete line) is parsed to a number. The separator is used to split
     * the line. It is used as regular expression.
     *
     * @param numerical numerical comparison
     * @param separator the separator used to split fields
     * @param field the field (split using th separator as regular expression)
     */
    LineComparator(boolean numerical, String separator, int field) {
        this.field = field;
        this.numerical = numerical;
        this.separator = separator;
    }

    /**
     * Add a sub comparator that is used if this comparator would returns 0
     *
     * @param comparator the sub comparator
     */
    public void addComparator(Comparator<String> comparator){
        if(comparator != null){
            if(subComparators == null) subComparators = new ArrayList<Comparator<String>>();
            subComparators.add(comparator);
        }
    }

    public int compare(final String o1, final String o2) {
        String s1 = o1;
        String s2 = o2;
        if(field >= 0){
            // split fields
            s1 = o1.split(separator)[field];
            s2 = o2.split(separator)[field];
        }

        // no field sep
        int result = 0;

        if(numerical){
            try{
                // try integer first
                result =  Integer.parseInt(o1) - Integer.parseInt(o2);
            }catch (NumberFormatException e){
                // ok integer failed ... lets try double
                result =  Double.compare(Double.parseDouble(o1),Double.parseDouble(o2));
            }
        }else{
            result =  o1.compareTo(o2);
        }
        /*
        If result is 0 check the sub comparators
         */
        if(result == 0 && subComparators != null){
            for (Comparator<String> subComparator : subComparators) {
                int sub = subComparator.compare(o1, o2);
                if(sub != 0) return sub;
            }
        }
        return result;
    }
}
