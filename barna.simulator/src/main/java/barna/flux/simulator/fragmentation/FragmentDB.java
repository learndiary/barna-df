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

package barna.flux.simulator.fragmentation;

import barna.commons.ByteArrayCharSequence;
import barna.commons.io.ByteArrayInputStream;
import barna.commons.io.RandomAccessInputStream;
import barna.commons.log.Log;
import barna.commons.system.OSChecker;
import barna.commons.utils.StringUtils;
import barna.io.FileHelper;
import jdbm.PrimaryTreeMap;
import jdbm.RecordManager;
import jdbm.RecordManagerFactory;

import java.io.*;
import java.util.Iterator;

/**
 * the fragment DB allows to iterate over fragments with the same ID.
 * You have to manually call {@link #createIndex()}  first, before you can
 * use {@link #getEntries(String)} for a given id. This {@link #getEntries(String)}
 * method will never return null. You can use {@link #containsKey(String)}  to check
 * for existence.
 *
 * @author Thasso Griebel
 */
public class FragmentDB {
    /**
     * Index record name
     */
    public static final String RECORD_NAME = "fragmentpositions";
    /**
     * The original library file
     */
    private File libraryFile;
    /**
     * The index database file
     */
    private File indexFile;
    /**
     * The index store
     */
    private PrimaryTreeMap<String, Entry> index;
    /**
     * Total number of unique (base on id) entries
     */
    private long numberOfEntries;
    /**
     * Total number od scanned lines
     */
    private long numberOfLines;
    /**
     * Total number of fragments including duplicates
     */
    private long numberOfFragments;
    /**
     * The random access file
     */
    private RandomAccessFile access;
    /**
     * The input stream on the random access file
     */
    private RandomAccessInputStream inputStream;
    /**
     * The byte stream to read lines
     */
    private ByteArrayInputStream byteStream;

    /**
     * Log status messages
     */
    private boolean printStatus = false;
    /**
     * Initialize status
     */
    private boolean initializes = false;
    /**
     * Close status
     */
    private boolean closed = false;
    /**
     * The record manager
     */
    private RecordManager recman;

    /**
     * Create a new instance for a given library file
     *
     * @param libraryFile the library file
     */
    public FragmentDB(File libraryFile) {
        this.libraryFile = libraryFile;
    }

    /**
     * true if status messages are printed
     *
     * @return status true if printed
     */
    public boolean isPrintStatus() {
        return printStatus;
    }

    /**
     * Enable/disable printing of status messages
     *
     * @param printStatus print status messages
     */
    public void setPrintStatus(boolean printStatus) {
        this.printStatus = printStatus;
    }

    /**
     * Create the index from the specified library file
     *
     * @throws IOException in case of any IO errors
     */
    public void createIndex() throws IOException {

        // create index
        indexFile = FileHelper.createTempFile("fragment-", "-index.db");
        recman = RecordManagerFactory.createRecordManager(indexFile.getAbsolutePath());
        index = recman.treeMap(RECORD_NAME);

        BufferedReader libFileReader = null;
        try {
            if (isPrintStatus()) {
                Log.message("\tInitializing Fragment Index");
                Log.progressStart("Indexing");
            }

            // read the sorted file and put it in a zip form
            String newlinecharacters = null;
            try{
                newlinecharacters = FileHelper.guessFileSep(libraryFile);
            }catch(Exception e){
                newlinecharacters = OSChecker.NEW_LINE;
            }

            int lineSeparatorLength =  newlinecharacters.length();
            libFileReader = new BufferedReader(new FileReader(libraryFile));

            long totalSize = libraryFile.length();
            long currentPosition = 0;
            numberOfFragments = 0;
            numberOfLines = 0;
            numberOfEntries = 0;


            ByteArrayCharSequence cs = new ByteArrayCharSequence(300);
            long entryStart = 0;
            long entryLength = 0;

            ByteArrayCharSequence nextID = null;
            ByteArrayCharSequence currentID = null;
            String line = null;
            while ((line = libFileReader.readLine()) != null) {
                cs.set(line);
                cs.resetFind();
                nextID = cs.getToken(2);

                // initialize current
                if (currentID == null) {
                    currentID = nextID.cloneCurrentSeq();
                }

                if (!nextID.equals(currentID)) {
                    // add a new entry
                    numberOfEntries++;

                    Entry entry = new Entry(entryStart, entryLength);
                    index.put(currentID.toString(), entry);
                    if (numberOfEntries % 10000 == 0) {
                        recman.commit();
                    }
                    // reset for next entry
                    entryStart = entryStart + entryLength;
                    entryLength = 0;
                    currentID = nextID.cloneCurrentSeq();
                }


                numberOfLines++;
                currentPosition += cs.length() + lineSeparatorLength; // one for the missing newline
                entryLength += cs.length() + lineSeparatorLength;
                if (isPrintStatus()) {
                    Log.progress(currentPosition, totalSize);
                }
                int dups = cs.getTokenInt(3);
                numberOfFragments += Math.max(1, dups);

            }

            // add final entry
            numberOfEntries++;
            Entry entry = new Entry(entryStart, entryLength);
            try {
                index.put(currentID.toString(), entry);
            } catch (Exception e) {
                Log.error("Error while creating Fragment Index: " + e.getMessage());
                throw new RuntimeException("Error while creating Fragment Index: " + e.getMessage(), e);
            }
            recman.commit();


            if(isPrintStatus()){
                Log.progressFinish(StringUtils.OK, true);
                Log.message("\t" + numberOfLines + " lines indexed (" + numberOfFragments + " fragments, "+numberOfEntries+" entries)");
            }
        } finally {
            libFileReader.close();
        }
        initializes = true;
    }

    /**
     * Iterate over all fragment IDs
     *
     * @return ids iterable over all indexed ids
     */
    public Iterable<String> fragmentIDs() {
        if (!initializes) {
            throw new RuntimeException("Index is not initialzed ! Call createIndex() first !");
        }
        if (closed) {
            throw new RuntimeException("Index is closed!");
        }
        return index.keySet();
    }

    /**
     * Returns true if the index contains the given key
     *
     * @param id the id
     * @return contained true if contained
     */
    public boolean containsKey(String id) {
        if (!initializes) {
            throw new RuntimeException("Index is not initialzed ! Call createIndex() first !");
        }
        if (closed) {
            throw new RuntimeException("Index is closed!");
        }
        return index.containsKey(id);
    }

    /**
     * Get an iterator over the fragments. The iterable with lazily load the fragments from disk
     *
     * @param id the id of the fragment
     * @return fragments the iterable over the fragments
     * @throws IOException in case of any initial errors
     */
    public Iterable<ByteArrayCharSequence> getEntries(String id) throws IOException {
        if (!initializes) {
            throw new RuntimeException("Index is not initialzed ! Call createIndex() first !");
        }
        if (closed) {
            throw new RuntimeException("Index is closed!");
        }

        Entry entry = index.get(id);
        if (entry == null) return new FragmentIterable(-1);

        if (access == null) {
            access = new RandomAccessFile(this.libraryFile, "r");
            inputStream = new RandomAccessInputStream(access, false, 0);
            byteStream = new ByteArrayInputStream(inputStream);
        }

        // position the stream
        inputStream.position(entry.start);
        return new FragmentIterable(entry.length);
    }


    /**
     * Close the index. After the index is closed, you can not access it anymore
     * and access tries will result in an Exception thrown.
     */
    public void close() {
        try {
            if (inputStream != null) {
                inputStream.close();
            }
            if (byteStream != null) {
                byteStream.close();
            }
            if (access != null) {
                access.close();
            }
            recman.close();
            if (index!=null && indexFile.exists()) {
                File dir = indexFile.getParentFile();
                File[] files = dir.listFiles(new FilenameFilter() {
                    @Override
                    public boolean accept(File dir, String name) {
                        return name.contains(indexFile.getName());
                    }
                });
                for (File f : files) {
                    f.delete();
                }
            }
        } catch (IOException e) {
            Log.error("Unable to close Fragment DB!");
        } finally {
            closed = true;
        }

    }

    /**
     * Get the number of unique fragments based on the ID
     *
     * @return uniqueFragments number of fragments
     */
    public long getNumberOfEntries() {
        return numberOfEntries;
    }

    /**
     * Get the number of scanned lines.
     *
     * @return lines number of scanned lines
     */
    public long getNumberOfLines() {
        return numberOfLines;
    }

    /**
     * Number of fragments including duplication counts
     *
     * @return fragments number fragments including duplication counts
     */
    public long getNumberOfFragments() {
        return numberOfFragments;
    }

    /**
     * Simple serializable index entry
     */
    private static class Entry implements Serializable {
        long start;
        long length;

        Entry() {
        }

        Entry(long start, long length) {
            this.start = start;
            this.length = length;
        }
    }

    /**
     * Internal iterator implementation to extract the
     * fragments from a random access stream
     */
    class FragmentIterable implements Iterable<ByteArrayCharSequence>, Iterator<ByteArrayCharSequence> {
        /**
         * Internal cache
         */
        private ByteArrayCharSequence next;
        /**
         * The current instance
         */
        private ByteArrayCharSequence current;
        /**
         * Has next indicator
         */
        private boolean hasNext = true;
        /**
         * The limit indicator
         */
        private long limit;
        /**
         * Read bytes
         */
        private long bytes = 0;


        FragmentIterable(long limit) throws IOException {
            this.limit = limit;
            if (limit> 0) {
            	next = byteStream.getSequence();
            	current = new ByteArrayCharSequence(300);
            	// read the first entry
            	next.clear();
            	readNext();
            }
        }

        @Override
        public Iterator<ByteArrayCharSequence> iterator() {
            return this;
        }

        @Override
        public boolean hasNext() {
            return bytes <= limit && hasNext;
        }

        @Override
        public ByteArrayCharSequence next() {
            if (!hasNext) return null;
            current.set(next);
            try {
                next.clear();
                readNext();
            } catch (IOException e) {
                throw new RuntimeException(e);
            }
            // read the 
            return current;
        }

        @Override
        public void remove() {
            throw new UnsupportedOperationException();
        }

        /**
         * Read the next entry
         *
         * @return entry the next entry
         * @throws IOException in case of an error
         */
        private boolean readNext() throws IOException {
            int i = byteStream.readLine();
            if (i == -1 || limit < 0) {
                // end
                hasNext = false;
                return false;
            }
            bytes += (i + 1);
            return true;
        }
    }
}
