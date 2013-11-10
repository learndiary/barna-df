package barna.io.sam;

import barna.io.FileHelper;
import net.sf.samtools.*;

import java.io.*;
import java.util.concurrent.Callable;

/**
 * Wrapper class for SAM file format specific constants.
 * User: micha
 */
public class SAMConstants {

    /**
     * A single factory for file I/O.
     */
    private static SAMFileWriterFactory factory= null;

    /**
     * Lazily create a factory with the default temp dir.
     * @return a factory to create writer instances
     */
    public static SAMFileWriterFactory getFactory() {
        if (factory== null) {
            factory= new SAMFileWriterFactory();
            factory.setTempDirectory(FileHelper.tempDirectory);
        }

        return factory;
    }

    public static class SAMSorter implements Callable<Integer> {

        /**
         * Flag which is <code>true</code> if the output is sorted by position, and <code>false</code> if the output is
         * sorted by query name (=readID).
         */
        private boolean sortPosition= true;

        /**
         * Flag indicating whether <code>this</code> thread opened the input stream.
         */
        private boolean iStreamOwner= false;
        /**
         * Flag indicating whether <code>this</code> thread opened the output stream.
         */
        private boolean oStreamOwner= false;
        /**
         * Reader from which is read.
         */
        private SAMFileReader reader= null;
        /**
         * Writer to which is written.
         */
        private SAMFileWriter writer= null;

        /**
         * Constructor that initializes the streams and sets the flag according to which the output is sorted.
         * @param iStream stream from which is read
         * @param oStream stream to which is written
         * @param sortPosition <code>true</code> if sorting by position, <code>false</code> if sorting by query name
         *                     (=readID)
         */
        public SAMSorter(InputStream iStream, OutputStream oStream, boolean sortPosition) {
            this.sortPosition= sortPosition;
            reader= new SAMFileReader(iStream, false);
            SAMFileHeader header= reader.getFileHeader();
            header.setSortOrder(sortPosition? SAMFileHeader.SortOrder.coordinate: SAMFileHeader.SortOrder.queryname);
            writer= getFactory().makeSAMWriter(header, false, oStream);
        }

        /**
         * Constructor that initializes the streams and sets the flag according to which the output is sorted.
         * Optional flags allow to close both streams after processing has finished.
         * @param iStream stream from which is read
         * @param oStream stream to which is written
         * @param sortPosition <code>true</code> if sorting by position, <code>false</code> if sorting by query name
         *                     (=readID)
         */
        public SAMSorter(InputStream iStream, boolean iStreamOwner, OutputStream oStream, boolean oStreamOwner, boolean sortPosition) {
            this.sortPosition= sortPosition;
            reader= new SAMFileReader(iStream, false);
            SAMFileHeader header= reader.getFileHeader();
            header.setSortOrder(sortPosition? SAMFileHeader.SortOrder.coordinate: SAMFileHeader.SortOrder.queryname);
            writer= getFactory().makeSAMWriter(header, false, oStream);
            this.iStreamOwner= iStreamOwner;
            this.oStreamOwner= oStreamOwner;
        }

        /**
         * Constructor for reading from a file and writing to a stream.
         * @param inFile file from which is read
         * @param oStream stream to which is written
         * @param sortPosition <code>true</code> if sorting by position, <code>false</code> if sorting by query name
         *                     (=readID)
         */
        public SAMSorter(File inFile, OutputStream oStream, boolean sortPosition) {
            this.sortPosition= sortPosition;
            reader= new SAMFileReader(inFile, false);
            SAMFileHeader header= reader.getFileHeader();
            header.setSortOrder(sortPosition? SAMFileHeader.SortOrder.coordinate: SAMFileHeader.SortOrder.queryname);
            writer= getFactory().makeSAMWriter(header, false, oStream);
            iStreamOwner= true;
        }

        /**
         * Constructor for reading from a file and writing to a stream, which can be closed.
         * @param inFile file from which is read
         * @param oStream stream to which is written
         * @param sortPosition <code>true</code> if sorting by position, <code>false</code> if sorting by query name
         *                     (=readID)
         */
        public SAMSorter(File inFile, OutputStream oStream, boolean oStreamOwner, boolean sortPosition) {
            this.sortPosition= sortPosition;
            reader= new SAMFileReader(inFile, false);
            SAMFileHeader header= reader.getFileHeader();
            header.setSortOrder(sortPosition? SAMFileHeader.SortOrder.coordinate: SAMFileHeader.SortOrder.queryname);
            writer= getFactory().makeSAMWriter(header, false, oStream);
            iStreamOwner= true;
            this.oStreamOwner= oStreamOwner;
        }

        /**
         * Constructor for reading from a stream and writing to a file.
         * @param iStream stream from which is read
         * @param outFile file to which is written
         * @param sortPosition <code>true</code> if sorting by position, <code>false</code> if sorting by query name
         *                     (=readID)
         */
        public SAMSorter(InputStream iStream, File outFile, boolean sortPosition) {
            this.sortPosition= sortPosition;
            reader= new SAMFileReader(iStream, false);
            SAMFileHeader header= reader.getFileHeader();
            header.setSortOrder(sortPosition? SAMFileHeader.SortOrder.coordinate: SAMFileHeader.SortOrder.queryname);
            writer= getFactory().makeSAMWriter(header, false, outFile);
            oStreamOwner= true;
        }

        /**
         * Constructor for reading from a stream, which can optionally be closed after processing, and writing to a file.
         * @param iStream stream from which is read
         * @param outFile file to which is written
         * @param sortPosition <code>true</code> if sorting by position, <code>false</code> if sorting by query name
         *                     (=readID)
         */
        public SAMSorter(InputStream iStream, boolean iStreamOwner, File outFile, boolean sortPosition) {
            this.sortPosition= sortPosition;
            reader= new SAMFileReader(iStream, false);
            SAMFileHeader header= reader.getFileHeader();
            header.setSortOrder(sortPosition? SAMFileHeader.SortOrder.coordinate: SAMFileHeader.SortOrder.queryname);
            writer= getFactory().makeSAMWriter(header, false, outFile);
            this.iStreamOwner= iStreamOwner;
            oStreamOwner= true;
        }

        /**
         * Constructor for reading from a file and writing to a file.
         * @param inFile file from which is read
         * @param outFile file to which is written
         * @param sortPosition <code>true</code> if sorting by position, <code>false</code> if sorting by query name
         *                     (=readID)
         */
        public SAMSorter(File inFile, File outFile, boolean sortPosition) {
            this.sortPosition= sortPosition;
            reader= new SAMFileReader(inFile, false);
            SAMFileHeader header= reader.getFileHeader();
            header.setSortOrder(sortPosition? SAMFileHeader.SortOrder.coordinate: SAMFileHeader.SortOrder.queryname);
            writer= getFactory().makeSAMWriter(header, false, outFile);
            iStreamOwner= true;
            oStreamOwner= true;
        }

        /**
         * Copies bytes from the input to the output.
         * @return the number of lines copied
         * @throws Exception if an I/O Exception occurs
         */
        @Override
        public Integer call() throws Exception {

            SAMRecordIterator iter= reader.iterator();
            int n= 0;
            while (iter.hasNext()) {
                writer.addAlignment(iter.next());
                ++n;
            }

            if (iStreamOwner)
                reader.close();
            if (oStreamOwner)
                writer.close();

            return n;
        }
    }


    /**
     * Number of reported alignments that contains the query in the current record.
     * Although optional, the <code>NM</code> tag should be present according to the
     * specification.
     * Example: <code>NM:i:2</code> corresponds to 2 mismatches.
     */
    public static String SAM_OPTION_NH= "NH";

    /**
     * <code>XT:A:U</code> identifies unique mappings.
     * Note that tags starting with `X', `Y' and `Z' or tags containing lowercase letters
     * in either position are reserved for local use and will not be formally dened in any future version of
     * this specication.
     */
    public static String SAM_OPTION_XT= "XT";

    /**
     * Edit distance to the reference, including ambiguous bases but excluding clipping.
     */
    public static String SAM_OPTION_NM= "NM";


    /**
     * <code>XS:A:+/-</code> identifies the strand of the mapping.
     * Note that tags starting with `X', `Y' and `Z' or tags containing lowercase letters
     * in either position are reserved for local use and will not be formally dened in any future version of
     * this specication.
     */
    public static String SAM_OPTION_XS= "XS";

}
