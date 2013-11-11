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

package barna.io;

import barna.commons.Execute;
import barna.commons.io.DevNullOutputStream;
import barna.commons.launcher.Flux;
import barna.commons.utils.Interceptable;
import barna.commons.utils.LineComparator;

import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
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
 * <p/>
 * You can start sorting by calling {@link #sort()}. To create a background task that performs sorting until
 * no more data are available or it is canceled, use {@link #sortInBackground()}. The returned feature is already submitted
 * and running. To wait until it is finished, call {@link java.util.concurrent.Future#get()}.
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class Sorter {

	public static final String SEPARATOR_DEFAULT= "\t";
	
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
    private String separator = SEPARATOR_DEFAULT;
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
     * @param in     the input stream
     * @param out    the output stream
     * @param silent be silent
     */
    private Sorter(InputStream in, OutputStream out, boolean silent, String separator) {
        this.in = in;
        this.out = out;
        this.silent = silent;
        this.separator = separator;
        
        if (this.out == null) {
            // create dev null output
            this.out = new DevNullOutputStream();
        }


    }

    /**
     * Sort by given field
     *
     * @param field   the field index
     * @param numeric is the field numeric
     * @return sorter this sorter
     */
    public Sorter field(int field, boolean numeric) {
    	return field(this.separator, field, numeric);
    }
    
    /**
     * Sort by given field using a specific separator
     *
     * @param separator custom separator string
     * @param field   the field index
     * @param numeric is the field numeric
     * @return sorter this sorter
     */
    public Sorter field(String separator, int field, boolean numeric) {
        addComparator(new LineComparator(numeric, separator, field));
        return this;
    }

    /**
     * Sort by merged fields. All specified fields are merged and teh concatenated string is compared
     *
     * @param fields the fields
     * @return sorter this sorter
     */
    public Sorter field(int... fields) {
    	return field(this.separator, fields);
    }
    
    /**
     * Sort by merged fields. All specified fields are merged and teh concatenated string is compared
     * employing a custom separator
     *
     * @param separator custom separator
     * @param fields the fields
     * @return sorter this sorter
     */
    public Sorter field(String separator, int... fields) {
        addComparator(new LineComparator(separator, fields));
        return this;
    }

    /**
     * Use a custom comparator and the default separator
     *
     * @param comparator the comparator
     * @return sorter this sorter
     */
    public Sorter field(Comparator<? extends CharSequence> comparator) {
    	return field(this.separator, comparator);
    }
    
    /**
     * Use a custom comparator and a custom separator
     *
     * @param separator custom separator
     * @param comparator the comparator
     * @return sorter this sorter
     */
    public Sorter field(String separator, Comparator<? extends CharSequence> comparator) {
        if (comparator == null) {
            throw new NullPointerException("Null comparator is not permitted");
        }
        addComparator(new LineComparator(comparator));
        return this;
    }

    /**
     * Add an interceptor to the sorter.
     *
     * @param interceptor the interceptor
     * @return sorter this sorter
     */
    public Sorter addInterceptor(Interceptable.Interceptor<String> interceptor) {
        if (interceptor == null) {
            throw new NullPointerException();
        }
        if (this.interceptors == null) {
            this.interceptors = new ArrayList<Interceptable.Interceptor<String>>();
        }
        this.interceptors.add(interceptor);
        return this;

    }

    /**
     * Perform the sort
     *
     * @throws IOException in case of any IO errors
     */
    public void sort() throws IOException {
        StreamSorter s = createSorter(-1);
        s.sort(in, out);
    }

    /**
     * Perform the sort
     *
     * @throws IOException in case of any IO errors
     */
    public void sort(long fileSize) throws IOException {
        StreamSorter s = createSorter(fileSize);
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
    public Future sortInBackground() {
        final StreamSorter s = createSorter(-1);
        final InputStream input = in;
        final OutputStream output = out;
        final String sep = separator;

        Execute.initialize(Flux.THREADS);
        Future<Object> job = Execute.getExecutor().submit(new Callable<Object>() {
            public Object call() throws Exception {
                s.sort(input, output);
                return null;
            }
        });
        return job;
    }

    /**
     * Create a new sorter
     *
     * @param in     the input stream
     * @param out    the output stream
     * @param silent be silent
     * @param separator the separator used to tokenize the line
     * @return sorter the sorter
     */
    public static Sorter create(InputStream in, OutputStream out, boolean silent, String separator) {
    	if (separator== null)
    		separator= SEPARATOR_DEFAULT;
        return new Sorter(in, out, silent, separator);
    }

    /**
     * Create an instance of the actual sorter implementation
     *
     * @param fileSize optional current file size
     * @return streamSorter the stream sorter
     */
    protected StreamSorter createSorter(long fileSize) {
        UnixStreamSorter s = new UnixStreamSorter(silent, -1, false, separator);
        s.setFileSize(fileSize);
        LineComparator comparator = null;
        if (comparators.size() == 0) {
            comparator = new LineComparator(false, separator, -1);
        } else {
            comparator = comparators.get(0);
            for (int i = 1; i < comparators.size(); i++) {
                comparator.addComparator(comparators.get(i));
            }
        }
        s.setLineComparator(comparator);

        if (interceptors != null) {
            for (Interceptable.Interceptor<String> interceptor : interceptors) {
                s.addInterceptor(interceptor);
            }
        }
        return s;
    }

    /**
     * INTERNAL: add a line comparator and enable caching
     *
     * @param comparator the comparator
     */
    protected void addComparator(LineComparator comparator){
        comparator.setCache(new HashMap());
        this.comparators.add(comparator);
    }


}
