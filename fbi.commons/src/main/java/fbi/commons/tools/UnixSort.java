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
 * then merge-sort the chunks.
 *
 * @author Thasso Griebel (thasso.griebel@googlemail.com)
 */
public class UnixSort implements StreamSorter{
    /**
     * OS dependent line separator
     */
    private static final String LINE_SEP = System.getProperty("line.separator");
    /**
     * Maxumim number of chunks that are sorted in one run
     */
    private static final int SORT_CHUNKS = 16;
    /**
     * Maximum memory to use per chunk
     */
    private long memoryBound;

    /**
     * Create a new sorter that uses a {@code ~25%} of heapspace as chunk size
     */
    public UnixSort() {
        this((long) (Runtime.getRuntime().maxMemory() * 0.25));
    }

    /**
     * Create a new sorter with given chunk size limit
     *
     * @param memoryBound the maximum chunk size in bytes
     */
    public UnixSort(long memoryBound) {
        if(memoryBound <= 0) throw new IllegalArgumentException("You have to allow memory chunk size > 0");
        this.memoryBound = memoryBound;
    }

    public void sort(InputStream input, OutputStream output, int field, boolean numeric, String fieldSeparator) throws IOException {

        LineComparator comparator = new LineComparator(field, numeric, fieldSeparator);

        List<File> files =  divide(input, comparator, memoryBound);

        /*
         * make sure we open at most SORT_CHUNK files in parallel
         */
        while(SORT_CHUNKS >= 2 && files.size() > SORT_CHUNKS){
            // sort chunk
            ArrayList<File> chunks = new ArrayList<File>(files.subList(0, SORT_CHUNKS));

            // create a temp file where we put the result of this chunk sort
            try{
                File chunk = FileHelper.createTempFile("chunk", ".srt");
                //chunk.deleteOnExit();
                FileOutputStream out = new FileOutputStream(chunk);
                mergeFiles(chunks, out, comparator);


                // add chunk to list and remove the rest from the list
                files.add(chunk);
                files.removeAll(chunks);
            }catch (IOException tooManyFiles){
                if(tooManyFiles.getMessage().equals("Too many open files")){
                    System.err.println("Too many open files at " + files.size());
                    throw tooManyFiles;
                }
            }
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
     * @param input the input
     * @param comparator the comparator
     * @param memoryBound upper bound for memory per chunk
     * @return files files create by the divider
     * @throws java.io.IOException in case of errors
     */
    private List<File> divide(InputStream input, LineComparator comparator, long memoryBound) throws IOException {
        BufferedReader reader = new BufferedReader(new InputStreamReader(input));
        String line = null;
        List<File> files= new ArrayList<File>();
        try {

            // counters
            int bytes =0;
            List<String> lines = new ArrayList<String>();
            while( (line = reader.readLine()) != null){
                // add to list
                lines.add(line);
                bytes += line.length() + LINE_SEP.length(); // assume one byte for separator
                // check memory
                if(bytes >= memoryBound){
                    // sort the chunk
                    Collections.sort(lines, comparator);
                    // write sorted chunk to temp file and reset
                    files.add(writeTempFile(lines));
                    bytes = 0;
                    lines.clear();
                }
            }

            // add the last file
            Collections.sort(lines, comparator);
            files.add(writeTempFile(lines));

        } catch (IOException e) {
            throw e;
        }finally {
            try {reader.close();} catch (IOException e) {}
        }
        return files;
    }

    /**
     * Write the given lines to a temp file and return the file
     *
     * @param lines the lines
     * @return file the created file
     * @throws IOException in case of an error
     */
    private File writeTempFile(List<String> lines) throws IOException {
        File file = FileHelper.createTempFile("sort", ".srt");
        //file.deleteOnExit();
        BufferedWriter writer = new BufferedWriter(new FileWriter(file));
        for (String line : lines) {
            writer.write(line);
            writer.write(LINE_SEP);
        }
        writer.flush();
        writer.close();
        return file;
    }

    /**
     * Merge the content of the given files and write the sorted output to the stream.
     * The key assumption is that the file contents are sorted already, so this is essentially a
     * merge sort step. NOTE: this deletes the given files aver the merge !
     *
     * @param files the files
     * @param output the output stream
     * @param comparator the comparator
     * @throws IOException in case of any errors
     */
    private void mergeFiles(List<File> files, OutputStream output, LineComparator comparator) throws IOException {
        // create a queue and add
        PriorityQueue<CachedFileReader> queue = new PriorityQueue<CachedFileReader>();
        // add files
        for (File file : files) {
            CachedFileReader cc = CachedFileReader.create(file, comparator);
            if(cc != null){
                queue.add(cc);
            }
        }

        BufferedWriter writer = null;

        // now iterate until everything is written
        while(queue.size() > 0){
            CachedFileReader next = queue.poll();

            try {
                if(writer == null){
                    writer = new BufferedWriter(new OutputStreamWriter(output));
                }
                String line = next.pop();
                if(line == null || line.isEmpty()){
                    Log.error("Line is null or empty ?");
                }
                writer.write(line);
                writer.write(LINE_SEP);

                // check if there is more in this queue
                if(next.peek() != null) queue.add(next);

            } catch (IOException e) {
                Log.error("Error while sorting chunks : "  + e.getMessage(), e);

                for (CachedFileReader reader : queue) {
                    reader.close();
                }
                break;
            }
        }

        // make sure all reader are closed
        for (CachedFileReader reader : queue) {
            reader.close();
        }


//        for (File file : files) {
//            file.delete();
//        }

        if(writer != null){
            writer.flush();
            writer.close();
        }

    }

    /**
     * Helper class that caches the first line of a file and is able to continue reading if
     * the line is popped. Used during the merge step.
     */
    static class CachedFileReader implements Comparable<CachedFileReader>{
        private BufferedReader reader;
        private String currentLine = null;
        private LineComparator comparator;

        /**
         * INTERNAL
         *
         * @param reader the reader
         * @param firstLine the first line
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
        String peek(){
            return currentLine;
        }

        /**
         * Return the current first line and reads the next one.
         * If there are no more lines, this return null.
         *
         * @return line current line or null
         */
        String pop(){
            if(currentLine == null){
                // close
                try {reader.close();} catch (IOException e1) {}
                return null;
            }
            String last = currentLine;
            try {
                currentLine = reader.readLine();
            } catch (IOException e) {
                Log.error("Error reading line from chunk file while sorting : "+ e.getMessage(), e);
                try {reader.close();} catch (IOException e1) {}
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
         * @param file the file
         * @param comparator the comparator
         * @return reader the reader or null if the file does not exist or has no content
         */
        static CachedFileReader create(File file, LineComparator comparator){
            BufferedReader reader = null;
            try {
                reader = new BufferedReader(new FileReader(file));
            } catch (FileNotFoundException e) {
                return null;
            }

            // read first line
            try {
                String currentLine = reader.readLine();
                if(currentLine == null){
                    try {reader.close();} catch (IOException e1) {}
                    return null;
                }
                return new CachedFileReader(reader, currentLine, comparator);
            } catch (IOException e) {
                Log.error("Error while reading sort chunk : " + e.getMessage(), e);
                try {reader.close();} catch (IOException e1) {}
            }
            return null;

        }

    }
}
