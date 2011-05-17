package fbi.commons.tools;


import fbi.commons.Log;
import fbi.commons.file.FileHelper;

import java.io.*;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.PriorityQueue;
import java.util.concurrent.*;

/**
 * Implements an external R-Way merge sort by creating sorted chunks, writing them to disk and
 * then merge-sort the chunks. This also supports intercepting lines and use the modified version
 * of the line.
 *
 * @author Thasso Griebel (thasso.griebel@googlemail.com)
 */
public class UnixStreamSorter implements StreamSorter, Interceptable<String> {
    /**
     * OS dependent line separator
     */
    private static final String LINE_SEP = System.getProperty("line.separator");
    /**
     * Maximum number of chunks that are sorted in one run
     */
    private int sortChunks = 4;
    /**
     * Maximum memory to use per chunk
     */
    private long memoryBound;
    /**
     * Print status
     */
    private boolean silent = true;
    /**
     * The line Comparator
     */
    private LineComparator lineComparator;
    /**
     * List of interceptors
     */
    private List<Interceptor<String>> interceptors;
    /**
     * File size of the file to sort (optional, disable with -1)
     */
    private long fileSize;

    /**
     * Create a new sorter that uses a {@code ~10%} of heapspace as chunk size
     *
     * @param field the field
     * @param numeric is the field numeric
     * @param fieldSeparator the field separator
     */
    public UnixStreamSorter(int field, boolean numeric, String fieldSeparator) {
        this((long) (Runtime.getRuntime().maxMemory() * 0.1), field, numeric, fieldSeparator);
    }

    /**
     * Create a new sorter and control if it should print status information
     *
     * @param silent be silent
     * @param field the field
     * @param numeric is the field numeric
     * @param fieldSeparator the field separator
     */
    public UnixStreamSorter(boolean silent, int field, boolean numeric, String fieldSeparator) {
        this((long) (Runtime.getRuntime().maxMemory() * 0.1), silent, field, numeric, fieldSeparator);
    }

    /**
     * Create a new sorter with given chunk size limit
     *
     * @param memoryBound the maximum chunk size in bytes
     * @param field the field
     * @param numeric is the field numeric
     * @param fieldSeparator the field separator

     */
    public UnixStreamSorter(long memoryBound, int field, boolean numeric, String fieldSeparator) {
        this(memoryBound, true, field, numeric, fieldSeparator);
    }

    /**
     * Create a new sorter and control its memory bound and the if it should print status information
     *
     * @param memoryBound the memory bound (must be {@code > 0}
     * @param silent      be silent
     * @param field the field
     * @param numeric is the field numeric
     * @param fieldSeparator the field separator

     */
    public UnixStreamSorter(long memoryBound, boolean silent, int field, boolean numeric, String fieldSeparator) {
        if (memoryBound <= 0) {
            throw new IllegalArgumentException("You have to allow memory chunk size > 0");
        }
        this.memoryBound = memoryBound;
        this.silent = silent;
        lineComparator = new LineComparator(numeric, fieldSeparator, field);
    }

    public void sort(InputStream input, OutputStream output) throws IOException {
        LineComparator comparator = getLineComparator();
        List<SorterFile> files = divide(input, comparator, memoryBound);

        /*
         * make sure we open at most SORT_CHUNK files in parallel
         */
        while (sortChunks >= 2 && files.size() > sortChunks) {
            if (Thread.interrupted()) {
                break;
            }
            // sort chunk
            ArrayList<SorterFile> chunks = new ArrayList<SorterFile>(files.subList(0, sortChunks));

            // create a temp file where we put the result of this chunk sort
            File chunk = FileHelper.createTempFile("chunk", ".srt");
            chunk.deleteOnExit();
            FileOutputStream out = new FileOutputStream(chunk);
            int lines = mergeFiles(chunks, out, comparator);

            // add chunk to list and remove the rest from the list
            files.add(new SorterFile(chunk, lines));
            files.removeAll(chunks);
        }


        // final merge
        mergeFiles(files, output, comparator);

        // make sure in and out are closed
        output.close();
        input.close();
    }


    /**
     * Divide the input stream content into sorted chunks.
     *
     * @param input       the input
     * @param comparator  the comparator
     * @param memoryBound upper bound for memory per chunk
     * @return files files create by the divider
     * @throws java.io.IOException in case of errors
     */
    private List<SorterFile> divide(InputStream input, LineComparator comparator, long memoryBound) throws IOException {
        BufferedReader reader = new BufferedReader(new InputStreamReader(input), 10 * 1024);
        String line = null;
        final List<SorterFile> files = new ArrayList<SorterFile>();

        int blocks = (int) ((fileSize / memoryBound) + 1);
        if (!silent) {
            Log.progressStart("\tdividing input to ~"+blocks + " blocks ");
        }


        try {
            // counters
            int bytes = 0;
            int separatorLength = LINE_SEP.length();
            List<String> lines = new ArrayList<String>((int) (memoryBound / 512));

            while ((line = reader.readLine()) != null) {
                if (Thread.interrupted()) {
                    return null;
                }

                // add to list
                lines.add(line);
                bytes += line.length() + separatorLength;
                // check memory
                if (bytes >= memoryBound) {
                    // write sorted chunk to temp file and reset
                    sortAndWriteTempFile(lines, comparator, files);
                    if(!silent && fileSize > 0){
                        Log.progress(files.size(), blocks);
                    }
                    bytes = 0;
                    lines.clear();
                }
            }

            // add the last file
            if (bytes > 0) {
                sortAndWriteTempFile(lines, comparator, files);
                if(!silent && fileSize > 0){
                    Log.progress(files.size(), blocks);
                }
            }

            if(!silent){
                Log.progressFinish("Done", true);
            }

        } finally {
            try {
                input.close();
                reader.close();
            } catch (IOException e) {
            }
        }
        return files;
    }

    /**
     * Write the given lines to a temp file and return the file
     *
     * @param lines      the lines
     * @param comparator the comparator
     * @param files      @return file the created file
     * @return file the file
     * @throws IOException in case of an error
     */
    private SorterFile sortAndWriteTempFile(List<String> lines, LineComparator comparator, final List<SorterFile> files) throws IOException {
        // sort the chunk
        Collections.sort(lines, comparator);
        comparator.reset();

        // write the file
        File file = FileHelper.createTempFile("sort", ".srt");
        file.deleteOnExit();
        BufferedWriter writer = new BufferedWriter(new FileWriter(file), 10 * 1024);
        for (String line : lines) {
            if (Thread.interrupted()) {
                return null;
            }
            writer.write(line);
            writer.write(LINE_SEP);
        }
        writer.flush();
        writer.close();



        // add to list of files
        SorterFile sorterFile = new SorterFile(file, lines.size());
        files.add(sorterFile);
        return sorterFile;
    }

    /**
     * Merge the content of the given files and write the sorted output to the stream.
     * The key assumption is that the file contents are sorted already, so this is essentially a
     * merge sort step. NOTE: this deletes the given files aver the merge !
     *
     * @param files      the files
     * @param output     the output stream
     * @param comparator the comparator
     * @throws IOException in case of any errors
     * @return lines number of lines written to file
     */
    private int mergeFiles(List<SorterFile> files, OutputStream output, LineComparator comparator) throws IOException {
        // create a queue and add
        PriorityQueue<CachedFileReader> queue = new PriorityQueue<CachedFileReader>();

        ExecutorService exec = Executors.newFixedThreadPool(sortChunks + 1);

        int totalLines = 0;
        // add files
        for (SorterFile file : files) {
            CachedFileReader cc = CachedFileReader.create(file.getFile(), comparator, exec);
            totalLines += file.getLines();
            if (cc != null) {
                queue.add(cc);
            }
        }
        BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(output));
        int lines= 0;

        if(!silent){
            Log.progressStart("\tmerging " + files.size() + " files");
        }
        // now iterate until everything is written
        while (queue.size() > 0) {
            if (Thread.interrupted()) {
                break;
            }
            CachedFileReader next = queue.poll();
            try {
                String line = next.pop();
                if (interceptors != null) {
                    for (Interceptor<String> interceptor : interceptors) {
                        line = interceptor.intercept(line);
                    }
                }
                writer.write(line);
                writer.write(LINE_SEP);
                lines++;
                if(!silent){
                    Log.progress(lines, totalLines);
                }


                // check if there is more in this queue
                if (next.peek() != null) {
                    queue.add(next);
                }
            } catch (Exception e) {
                Log.error("Error while sorting chunks : " + e.getMessage(), e);
                break;
            }
        }

        exec.shutdownNow();

        // make sure all reader are closed
        for (CachedFileReader reader : queue) {
            reader.close();
        }
        // delete the temp files
        for (SorterFile file : files) {
            file.getFile().delete();
        }
        writer.flush();
        writer.close();
        return lines;
    }

    /**
     * Get the lie comparator or null
     *
     * @return comparator the line comparator or null
     */
    LineComparator getLineComparator() {
        return lineComparator;
    }

    /**
     * Set the line comparator
     *
     * @param lineComparator the line comparator
     */
    void setLineComparator(LineComparator lineComparator) {
        this.lineComparator = lineComparator;
    }

    public void addInterceptor(Interceptor<String> stringInterceptor) {
        if (stringInterceptor == null) {
            return;
        }
        if (interceptors == null) {
            interceptors = new ArrayList<Interceptor<String>>();
        }
        if (!interceptors.contains(stringInterceptor)) {
            interceptors.add(stringInterceptor);
        }

    }

    /**
     * Optional method to set the file size to get mor informative progress information
     *
     * @param fileSize the filesize
     */
    public void setFileSize(final long fileSize) {
        this.fileSize = fileSize;
    }

    /**
     * Helper class that caches the first line of a file and is able to continue reading if
     * the line is popped. Used during the merge step.
     */
    static class CachedFileReader implements Comparable<CachedFileReader> {
        private BufferedReader reader;
        private LineComparator comparator;
        private ExecutorService exec;
        private int read;
        private Future<List<String>> current;
        private Future<List<String>> next;
        private List<String> lines = new ArrayList<String>();
        private int position = 0;
        private int readahead = 1000;

        /**
         * INTERNAL
         *
         * @param reader     the reader
         * @param comparator the comparator
         * @param exec the executor
         */
        private CachedFileReader(BufferedReader reader, LineComparator comparator, final ExecutorService exec) {
            this.reader = reader;
            this.comparator = comparator;
            this.exec = exec;
            comparator.reset();
            current = fetch();
        }

        /**
         * Return the current first line
         *
         * @return line the current line
         */
        String peek() {
            try {
                if(current == null || current.get() == null || current.get().size() == 0 || current.get().size() < position ) return null;
                return current.get().get(position);
            } catch (Exception e) {
                e.printStackTrace();
            }
            return null;
        }

        /**
         * Return the current first line and reads the next one.
         * If there are no more lines, this return null.
         *
         * @return line current line or null
         */
        String pop() throws ExecutionException, InterruptedException {
            if (current == null || current.get() == null || current.get().size() == 0) {
                close();
                return null;
            }
            read++;
            String last = current.get().get(position++);
            if(last == null){
                current = null;
                close();
            }else if(next == null && position >= current.get().size() / 2){
                next = fetch();
            }else if(position >= current.get().size()){
                current = next;
                position = 0;
                comparator.reset();
            }
            return last;
        }

        /**
         * Compares by the current line
         *
         * @param o the other reader
         * @return c comparison result
         */
        public int compareTo(CachedFileReader o) {
            return comparator.compare(peek(), o.peek());
        }

        /**
         * Use for emergencies and close the stream
         */
        public void close() {
            try {
                reader.close();
            } catch (IOException e) {
                // ignore
            }
        }

        private Future<List<String>> fetch(){
            return exec.submit(new Callable<List<String>>() {
                @Override
                public List<String> call() throws Exception {

                    List<String> lines = new ArrayList<String>();
                    String l = null;
                    int c = 0;
                    while(c++ < readahead && (l = reader.readLine()) != null){
                        lines.add(l);
                    }
                    return lines;
                }
            });
        }


        /**
         * Create a CachedFileReader for the given file iff the file has content
         *
         *
         * @param file       the file
         * @param comparator the comparator
         * @param exec
         * @return reader the reader or null if the file does not exist or has no content
         */
        static CachedFileReader create(File file, LineComparator comparator, final ExecutorService exec) {
            BufferedReader reader = null;
            try {
                reader = new BufferedReader(new FileReader(file));
            } catch (FileNotFoundException e) {
                e.printStackTrace();
                return null;
            }
            // read first line
            return new CachedFileReader(reader, comparator, exec);
        }
    }

    /**
     * Stores one divided file and its number of entries
     */
    private static class SorterFile{
        /**
         * The file
         */
        private File file;
        /**
         * The lines
         */
        private int lines;

        /**
         * Create a new instance
         *
         * @param file the file
         * @param lines the lines
         */
        private SorterFile(final File file, final int lines) {
            this.file = file;
            this.lines = lines;
        }

        /**
         * Get the file
         * @return file the file
         */
        public File getFile() {
            return file;
        }

        /**
         * Get the number of lines written to the file
         *
         * @return lines number of lines written to the file
         */
        public int getLines() {
            return lines;
        }
    }
}
