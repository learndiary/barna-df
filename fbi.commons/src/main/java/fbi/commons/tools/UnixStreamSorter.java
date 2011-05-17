package fbi.commons.tools;


import fbi.commons.Log;
import fbi.commons.file.FileHelper;

import java.io.*;
import java.util.*;
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
    private int sortChunks = 16;
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
        this((long) (Runtime.getRuntime().maxMemory() * 0.001), field, numeric, fieldSeparator);
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
        this((long) (Runtime.getRuntime().maxMemory() * 0.001), silent, field, numeric, fieldSeparator);
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


        // split th input
        List<SorterFile> files = divide(input, comparator, memoryBound);



        int totalLines = 0;
        for (SorterFile file : files) {
            totalLines += file.getLines();
        }

        int c = 2;
        int m = files.size();
        while(m > sortChunks){
            c++;
            m = (int) Math.ceil(m / sortChunks);
        }
        int currentMerges = 0;
        totalLines = c * totalLines;
        ExecutorService exec = Executors.newFixedThreadPool(4);


        if(!silent){
            Log.progressStart("\tmerging ~" + files.size() + " blocks");
        }
        /*
         * make sure we open at most SORT_CHUNK files in parallel
         */
        while (files.size() > 0) {
            if (Thread.interrupted()) {
                break;
            }
            try {
                if(files.size() <= sortChunks){
                    // final run
                    mergeFiles(files, exec, output, totalLines, currentMerges);
                    files.clear();
                }else{
                    // intermediate run
                    List<SorterFile> next = mergeFiles(files, exec, null, totalLines, currentMerges);
                    for (SorterFile sorterFile : next) {
                        currentMerges += sorterFile.getLines();
                    }
                    files.clear();
                    files.addAll(next);
                }
            } catch (Exception e) {
                e.printStackTrace();
            }

        }
        exec.shutdownNow();

        // make sure in and out are closed
        output.close();
        input.close();
        if(!silent){
            Log.progressFinish("Done", true);
        }
    }


    protected List<SorterFile> mergeFiles(List<SorterFile> files, ExecutorService exec, final OutputStream target, final int numMerges, int currentMerges) throws Exception {
        List<Future<SorterFile>> jobs = new ArrayList<Future<SorterFile>>();
        while(files.size() > 0){
            final ArrayList<SorterFile> chunks = new ArrayList<SorterFile>(files.subList(0, Math.min(files.size(), sortChunks)));
            final int finalCurrentMerges = currentMerges;
            jobs.add(exec.submit(new Callable<SorterFile>() {
                @Override
                public SorterFile call() throws Exception {
                    // create a temp file where we put the result of this chunk sort
                    OutputStream out = target;
                    File chunk = null;
                    if(out == null){

                        chunk = FileHelper.createTempFile("chunk", ".srt");
                        chunk.deleteOnExit();
                        out = new FileOutputStream(chunk);
                    }
                    int lines = mergeFiles(chunks, out, new LineComparator(getLineComparator()), target != null, finalCurrentMerges, numMerges);
                    out.flush();
                    out.close();
                    return new SorterFile(chunk, lines, chunks);
                }
            }));
            files.removeAll(chunks);
        }
        // wait for the jobs
        List<SorterFile> result = new ArrayList<SorterFile>();

        while(jobs.size() > 0){
            Iterator<Future<SorterFile>> iterator = jobs.iterator();
            while(iterator.hasNext()){
                Future<SorterFile> job = iterator.next();
                if(job.isDone()){
                    SorterFile sorterFile = job.get();
                    result.add(sorterFile);
                    iterator.remove();
                    currentMerges +=sorterFile.getLines();
                    if(!silent){
                        Log.progress(currentMerges, numMerges);
                    }
                }
            }
        }
        return result;
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
    private List<SorterFile> divide(InputStream input, final LineComparator comparator, long memoryBound) throws IOException {
        BufferedReader reader = new BufferedReader(new InputStreamReader(input), 10 * 1024);
        String line = null;
        final List<SorterFile> files = new ArrayList<SorterFile>();

        int blocks = (int) ((fileSize / memoryBound) + 1);
        if (!silent) {
            Log.progressStart("\tdividing input to ~"+blocks + " blocks ");
        }

        ExecutorService exec = Executors.newFixedThreadPool(4);
        try {
            // counters
            int bytes = 0;
            int separatorLength = LINE_SEP.length();
            final List<String> lines = new ArrayList<String>((int) (memoryBound / 512));

            final List<Future<SorterFile>> jobs = new ArrayList<Future<SorterFile>>();

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

                    if(!silent && fileSize > 0){
                        Log.progress(files.size(), blocks);
                    }

                    final ArrayList<String> jobLines = new ArrayList<String>(lines);
                    final LineComparator jobComparator = new LineComparator(comparator);
                    jobs.add(exec.submit(new Callable<SorterFile>() {
                        @Override
                        public SorterFile call() throws Exception {
                            SorterFile sorterFile = sortAndWriteTempFile(jobLines, jobComparator);
                            jobComparator.reset();
                            synchronized (jobs){
                                jobs.notifyAll();
                            }
                            return sorterFile;
                        }
                    })
                    );
                    bytes = 0;
                    lines.clear();
                    System.gc();


                    // wait for fobs
                    while (jobs.size() >= 4){
                        List<Future<SorterFile>> remove = new ArrayList<Future<SorterFile>>();
                        for (Future<SorterFile> job : jobs) {
                            if(job.isDone()){
                                try {
                                    files.add(job.get());
                                    remove.add(job);
                                } catch (Exception e) {
                                    e.printStackTrace();
                                }
                            }
                        }
                        if(remove.size() > 0){
                            jobs.removeAll(remove);
                        }else{
                            try {
                                synchronized (jobs){
                                    jobs.wait(10000);
                                }
                            } catch (InterruptedException e) {
                                e.printStackTrace();
                            }
                        }
                    }
                }
            }

            // wait for the rest of the jobs
            for (Future<SorterFile> job : jobs) {
                try {
                    files.add(job.get());
                } catch (Exception e) {
                    e.printStackTrace();
                }
            }

            // add the last file
            if (bytes > 0) {
                files.add(sortAndWriteTempFile(lines, comparator));
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
                exec.shutdownNow();
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
     * @return file the file
     * @throws IOException in case of an error
     */
    private SorterFile sortAndWriteTempFile(List<String> lines, LineComparator comparator) throws IOException {
        // sort the chunk
        Collections.sort(lines, comparator);
        //QuickSort.sort(lines, comparator);
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
        return new SorterFile(file, lines.size());
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
    private int mergeFiles(List<SorterFile> files, OutputStream output, LineComparator comparator, boolean status, int current, int total) throws IOException {
        // create a queue and add
        PriorityQueue<CachedFileReader> queue = new PriorityQueue<CachedFileReader>();

        // add files
        for (SorterFile file : files) {
            CachedFileReader cc = CachedFileReader.create(file.getFile(), new LineComparator(comparator));
            if (cc != null) {
                queue.add(cc);
            }
        }
        BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(output));
        int lines= 0;

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
                if(!silent && status){
                    Log.progress(++current, total);
                }

                // check if there is more in this queue
                String peek = next.peek();
                if (peek != null) {
                    queue.add(next);
                }else {
                    next.close();
                }
            } catch (Exception e) {
                Log.error("Error while sorting chunks : " + e.getMessage(), e);
                break;
            }
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
        private List<String> lines = new ArrayList<String>();
        private int position = 0;
        private int readahead = 1000;
        private boolean initialized = false;

        /**
         * INTERNAL
         *
         * @param reader     the reader
         * @param comparator the comparator
         */
        private CachedFileReader(BufferedReader reader, LineComparator comparator) {
            this.reader = reader;
            this.comparator = comparator;
            comparator.reset();
        }

        /**
         * Return the current first line
         *
         * @return line the current line
         */
        String peek() {
            if (!initialized) {
                read();
                initialized = true;
            }

            if(lines == null || lines.size() == 0 ) return null;

            // read next chunk
            if(position >= lines.size()){
                read();
            }

            // there was nothing more to read
            if(lines.size() == 0){
                lines = null;
                return null;
            }

            return lines.get(position);
        }

        /**
         * Return the current first line and reads the next one.
         * If there are no more lines, this return null.
         *
         * @return line current line or null
         */
        String pop() throws ExecutionException, InterruptedException {
            if (!initialized) {
                read();
                initialized = true;
            }

            if (lines == null || lines.size() == 0) {
                lines = null;
                close();
                return null;
            }

            // read next chunk
            if(position >= lines.size()){
                read();
            }

            // there was nothing more to read
            if(lines.size() == 0){
                lines = null;
                return null;
            }
            return lines.get(position++);
        }

        private void read(){
            String l = null;
            int c = 0;
            position = 0;
            lines.clear();
            comparator.reset();
            try {
                while(c++ < readahead && (l = reader.readLine()) != null){
                    lines.add(l);
                }
            } catch (IOException e) {
                e.printStackTrace();
                lines = null;
            }
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


        /**
         * Create a CachedFileReader for the given file iff the file has content
         *
         *
         * @param file       the file
         * @param comparator the comparator
         * @return reader the reader or null if the file does not exist or has no content
         */
        static CachedFileReader create(File file, LineComparator comparator) {
            BufferedReader reader = null;
            try {
                reader = new BufferedReader(new FileReader(file));
            } catch (FileNotFoundException e) {
                e.printStackTrace();
                return null;
            }
            // read first line
            return new CachedFileReader(reader, comparator);
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
         * source files after merge
         */
        private List<SorterFile> sourceFiles;


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
         * Create a new instance
         *
         * @param file the file
         * @param lines number of lines
         * @param sourceFiles the source files
         */
        private SorterFile(final File file, final int lines, final List<SorterFile> sourceFiles) {
            this.file = file;
            this.lines = lines;
            this.sourceFiles = sourceFiles;
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

        /**
         * Get the source files or null
         *
         * @return source the source files or null
         */
        public List<SorterFile> getSourceFiles() {
            return sourceFiles;
        }
    }



}
