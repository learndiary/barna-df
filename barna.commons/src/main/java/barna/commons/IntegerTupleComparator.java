package barna.commons;

import java.util.Comparator;
import java.util.StringTokenizer;

/**
 * Created with IntelliJ IDEA.
 * User: micha
 * Date: 10/6/12
 * Time: 10:22 AM
 * To change this template use File | Settings | File Templates.
 */
public class IntegerTupleComparator implements Comparator {
    static final IntegerTupleComparator defaultIntegerTupleComparator= new IntegerTupleComparator();

    String delimiters= null;

    public String getDelimiters(String sample, boolean forceAdd) {
        if (delimiters == null|| forceAdd) {
            StringBuffer delimBuf= new StringBuffer();
            if (delimiters!= null)
                delimBuf.append(delimiters);
            String[] d= sample.split("\\d+");
            for (int i = 0; i < d.length; i++) {
                if (d[i].length()> 0&& delimBuf.indexOf(d[i])< 0)
                    delimBuf.append(d[i]);
            }
            delimiters= delimBuf.toString();
        }
        return delimiters;
    }

    public int compare(Object arg0, Object arg1) {

        StringTokenizer toki0= new StringTokenizer((String) arg0, getDelimiters((String) arg0, false));
        StringTokenizer toki1= new StringTokenizer((String) arg1, getDelimiters((String) arg1, false));

        int min= java.lang.Math.min(toki0.countTokens(), toki1.countTokens());
        while (toki0.hasMoreTokens()&& toki1.hasMoreTokens()) {
            String next0= toki0.nextToken();
            int x= 0, y= 0;
            try {
                x= java.lang.Integer.parseInt(next0);
            } catch (Exception e) {
                getDelimiters(next0, true);
                return compare(arg0,arg1);
            }
            String next1= toki1.nextToken();
            try {
                y= java.lang.Integer.parseInt(next1);
            } catch (Exception e) {
                getDelimiters(next1, true);
                return compare(arg0,arg1);
            }
            if (x< y)
                return -1;
            if (x> y)
                return 1;
        }
        // equal over min length
        if (toki0.countTokens()< toki1.countTokens())
            return -1;
        if (toki0.countTokens()> toki1.countTokens())
            return 1;
        return 0;
    }

    public static IntegerTupleComparator getDefaultIntegerTupleComparator() {
        return defaultIntegerTupleComparator;
    }

}
