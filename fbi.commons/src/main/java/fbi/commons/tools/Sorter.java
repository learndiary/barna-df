package fbi.commons.tools;

import fbi.commons.io.DevNullOutputStream;

import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

/**
 * The Sorter follows the Builder pattern. Use the {@link #create(java.io.InputStream, java.io.OutputStream, boolean)}
 * method to get a new sorter instance. You can then call methods to configure the sorter. For example, call {@link #field(int, boolean)}
 * to specify a field that is used for sorting. You can call the {@code field} methods multiple times to add fields
 * to the consecutive sort order, i.e., if two values are equal in the first field, the next field is used.
 * <p>
 * You can start sorting by calling {@link #sort()}. To create a background task that performs sorting until
 * no more data are available or it is canceled, use {@link #sortInBackground()}. The returned feature is already submitted
 * and running. To wait until it is finished, call {@link java.util.concurrent.Future#get()}.
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class Sorter {

    /**
     * The thread executor
     */
    private static ExecutorService executor;

    static{
        /*
        Initialize executor
         */
        executor = Executors.newCachedThreadPool();
    }
    /**
     * The input stream
     */
    private InputStream in;
    /**
     * The output stream
     */
    private OutputStream out;
    /**
     * Be silent
     */
    private boolean silent;
    /**
     * The separator character
     */
    private String separator = "\t";
    /**
     * List of line comparators
     */
    private List<LineComparator> comparators = new ArrayList<LineComparator>();
    /**
     * List of interceptors
     */
    private List<Interceptable.Interceptor<String>> interceptors;


    /**
     * INTERNAL : use the {@link #create(java.io.InputStream, java.io.OutputStream, boolean)} method
     * to get an instance
     *
     * @param in the input stream
     * @param out the output stream
     * @param silent be silent
     */
    private Sorter(InputStream in, OutputStream out, boolean silent) {
        this.in = in;
        this.out = out;
        this.silent = silent;

        if(this.out == null){
            // create dev null output
            this.out = new DevNullOutputStream();
        }


    }

    /**
     * Set the field separator
     *
     * @param separator field separator
     * @return sorter this sorter
     */
    public Sorter separator(String separator){
        this.separator = separator;
        return this;
    }

    /**
     * Sort by given field
     *
     * @param field the field index
     * @param numeric is the field numeric
     * @return sorter this sorter
     */
    public Sorter field(int field, boolean numeric){
        comparators.add(new LineComparator(numeric, separator, field));
        return this;
    }

    /**
     * Sort by merged fields. All specified fields are merged and teh concatenated string is compared
     *
     * @param fields the fields
     * @return sorter this sorter
     */
    public Sorter field(int...fields){
        comparators.add(new LineComparator(separator, fields));
        return this;
    }

    /**
     * Use a custom comparator
     *
     * @param comparator the comparator
     * @return sorter this sorter
     */
    public Sorter field(Comparator<String> comparator){
        if(comparator == null) throw new NullPointerException("Null comparator is not permitted");
        comparators.add(new LineComparator(comparator));
        return this;
    }

    /**
     * Add an interceptor to the sorter.
     *
     * @param interceptor the interceptor
     * @return sorter this sorter
     */
    public Sorter addInterceptor(Interceptable.Interceptor<String> interceptor){
        if(interceptor == null) throw new NullPointerException();
        if(this.interceptors == null) this.interceptors = new ArrayList<Interceptable.Interceptor<String>>();
        this.interceptors.add(interceptor);
        return this;

    }

    /**
     * Perform the sort
     *
     * @throws IOException in case of any IO errors
     */
    public void sort() throws IOException {
        StreamSorter s = createSorter();
        s.sort(in, out);
    }

    /**
     * Submits a new background task and returns the created feature.
     * This creates a background task that performs sorting until
     * no more data are available or it is canceled, use {@link #sortInBackground()}.
     * The returned feature is already submitted and running. To wait until it is finished,
     * call {@link java.util.concurrent.Future#get()}, to cancel use {@link Future#cancel(boolean)}.
     *
     * @return feature the submitted feature
     */
    public Future sortInBackground(){
        final StreamSorter s = createSorter();
        final InputStream input = in;
        final OutputStream output = out;
        final String sep = separator;
        return executor.submit(new Callable<Object>() {
            public Object call() throws Exception {
                s.sort(input, output);
                return null;
            }
        });
    }

    /**
     * Create a new sorter
     *
     * @param in the input stream
     * @param out the output stream
     * @param silent be silent
     * @return sorter the sorter
     */
    public static Sorter create(InputStream in, OutputStream out, boolean silent){
        return new Sorter(in, out, silent);
    }

    /**
     * Create an instance of the actual sorter implementation
     *
     * @return streamSorter the stream sorter
     */
    protected StreamSorter createSorter() {
        UnixStreamSorter s = new UnixStreamSorter(silent, -1, false, separator);
        LineComparator comparator = null;
        if(comparators.size() == 0){
            comparator = new LineComparator(false, separator, -1);
        }else{
            comparator = comparators.get(0);
            for (int i = 1; i < comparators.size(); i++) {
                  comparator.addComparator(comparators.get(i));
            }
        }
        s.setLineComparator(comparator);

        if(interceptors != null){
            for (Interceptable.Interceptor<String> interceptor : interceptors) {
                s.addInterceptor(interceptor);
            }
        }
        return s;
    }


}
