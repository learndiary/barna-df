package fbi.commons.tools;


import fbi.commons.Log;
import fbi.commons.file.FileHelper;

import java.io.*;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.PriorityQueue;

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
     * Create a new sorter that uses a {@code ~25%} of heapspace as chunk size
     */
    public UnixStreamSorter(int field, boolean numeric, String fieldSeparator) {
        this((long) (Runtime.getRuntime().maxMemory() * 0.25), field, numeric, fieldSeparator);
    }

    /**
     * Create a new sorter and control if it should print status information
     *
     * @param silent be silent
     */
    public UnixStreamSorter(boolean silent, int field, boolean numeric, String fieldSeparator) {
        this((long) (Runtime.getRuntime().maxMemory() * 0.25), silent, field, numeric, fieldSeparator);
    }

    /**
     * Create a new sorter with given chunk size limit
     *
     * @param memoryBound the maximum chunk size in bytes
     */
    public UnixStreamSorter(long memoryBound, int field, boolean numeric, String fieldSeparator) {
        this(memoryBound, true, field, numeric, fieldSeparator);
    }

    /**
     * Create a new sorter and control its memory bound and the if it should print status information
     *
     * @param memoryBound the memory bound (must be {@code > 0}
     * @param silent      be silent
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
        if (!silent) {
            Log.progressStart("Sorting");
        }

        LineComparator comparator = getLineComparator();
        List<File> files = divide(input, comparator, memoryBound);
        /*
         * make sure we open at most SORT_CHUNK files in parallel
         */
        while (sortChunks >= 2 && files.size() > sortChunks) {
            if (Thread.interrupted()) {
                break;
            }
            Log.progressStart("Merging Chunk");
            // sort chunk
            ArrayList<File> chunks = new ArrayList<File>(files.subList(0, sortChunks));

            // create a temp file where we put the result of this chunk sort
            File chunk = FileHelper.createTempFile("chunk", ".srt");
            chunk.deleteOnExit();
            FileOutputStream out = new FileOutputStream(chunk);
            mergeFiles(chunks, out, comparator);

            // add chunk to list and remove the rest from the list
            files.add(chunk);
            files.removeAll(chunks);
            Log.progressFinish("OK", true);
        }


        // final merge
        mergeFiles(files, output, comparator);

        // make sure in and out are closed
        output.close();
        input.close();
        if (!silent) {
            Log.progressFinish("Done", true);
        }
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
    private List<File> divide(InputStream input, LineComparator comparator, long memoryBound) throws IOException {
        BufferedReader reader = new BufferedReader(new InputStreamReader(input), 10 * 1024);
        String line = null;
        final List<File> files = new ArrayList<File>();
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
                    bytes = 0;
                    lines.clear();
                }
            }

            // add the last file
            if (bytes > 0) {
                sortAndWriteTempFile(lines, comparator, files);
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
    private File sortAndWriteTempFile(List<String> lines, LineComparator comparator, final List<File> files) throws IOException {
        // sort the chunk
        Collections.sort(lines, comparator);

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
        files.add(file);

        return file;
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
     */
    private void mergeFiles(List<File> files, OutputStream output, LineComparator comparator) throws IOException {
        // create a queue and add
        PriorityQueue<CachedFileReader> queue = new PriorityQueue<CachedFileReader>();
        // add files
        for (File file : files) {
            CachedFileReader cc = CachedFileReader.create(file, comparator);
            if (cc != null) {
                queue.add(cc);
            }
        }
        BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(output));

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

                // check if there is more in this queue
                if (next.peek() != null) {
                    queue.add(next);
                }
            } catch (IOException e) {
                Log.error("Error while sorting chunks : " + e.getMessage(), e);
                break;
            }
        }

        // make sure all reader are closed
        for (CachedFileReader reader : queue) {
            reader.close();
        }
        // delete the temp files
        for (File file : files) {
            file.delete();
        }

        // close the writer
        if (writer != null) {
            writer.flush();
            writer.close();
        }

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
     * Helper class that caches the first line of a file and is able to continue reading if
     * the line is popped. Used during the merge step.
     */
    static class CachedFileReader implements Comparable<CachedFileReader> {
        private BufferedReader reader;
        private String currentLine = null;
        private LineComparator comparator;
        private int read;

        /**
         * INTERNAL
         *
         * @param reader     the reader
         * @param firstLine  the first line
         * @param comparator the comparator
         */
        private CachedFileReader(BufferedReader reader, String firstLine, LineComparator comparator) {
            this.reader = reader;
            this.currentLine = firstLine;
            this.comparator = comparator;
        }

        /**
         * Return the current first line
         *
         * @return line the current line
         */
        String peek() {
            return currentLine;
        }

        /**
         * Return the current first line and reads the next one.
         * If there are no more lines, this return null.
         *
         * @return line current line or null
         */
        String pop() {
            if (currentLine == null) {
                return null;
            }
            read++;
            String last = currentLine;
            try {
                currentLine = reader.readLine();

                // close the reader if no more lines are available
                if (currentLine == null) {
                    close();
                }
            } catch (IOException e) {
                close();
                Log.error("Error reading line from chunk file while sorting : " + e.getMessage(), e);
                return null;
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


        /**
         * Create a CachedFileReader for the given file iff the file has content
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
            try {
                String currentLine = reader.readLine();
                if (currentLine == null) {
                    try {
                        reader.close();
                    } catch (IOException e1) {
                    }
                    return null;
                }
                return new CachedFileReader(reader, currentLine, comparator);
            } catch (IOException e) {
                Log.error("Error while reading sort chunk : " + e.getMessage(), e);
                try {
                    reader.close();
                } catch (IOException e1) {
                }
            }
            return null;

        }
    }
}
