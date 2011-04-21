package fbi.commons.tools;

import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.List;

/**
 * Sorter
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class Sorter {

    private InputStream in;
    private OutputStream out;
    private boolean silent;
    private String separator = "\t";
    private List<Integer> fields = new ArrayList<Integer>();
    private List<Boolean> numerics = new ArrayList<Boolean>();

    private Sorter(InputStream in, OutputStream out, boolean silent) {
        this.in = in;
        this.out = out;
        this.silent = silent;
    }

    public Sorter separator(String separator){
        this.separator = separator;
        return this;
    }

    public Sorter field(int field, boolean numeric){
        this.fields.add(field);
        this.numerics.add(numeric);
        return this;
    }

    public void sort() throws IOException {
        UnixStreamSorter s = new UnixStreamSorter(silent);
        LineComparator comparator = new LineComparator(numerics.size() > 0 ? numerics.get(0) : false, separator, fields.size() > 0 ? fields.get(0) : -1);
        if(fields.size() > 1){
            for (int i = 1; i < fields.size(); i++) {
                comparator.addComparator(new LineComparator(numerics.get(i), separator, fields.get(i)));
            }
        }
        s.setLineComparator(comparator);
        s.sort(in, out, -1, false, separator);
    }



    public static Sorter create(InputStream in, OutputStream out, boolean silent){
        return new Sorter(in, out, silent);
    }

}
