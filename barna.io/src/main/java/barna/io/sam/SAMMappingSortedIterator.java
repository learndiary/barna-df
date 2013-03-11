package barna.io.sam;

import barna.commons.log.Log;
import barna.io.FileHelper;
import barna.io.MSIterator;
import barna.io.rna.UniversalReadDescriptor;
import barna.model.Mapping;
import barna.model.sam.SAMMapping;
import net.sf.samtools.*;

import java.io.BufferedOutputStream;
import java.io.IOException;
import java.io.PipedInputStream;
import java.io.PipedOutputStream;
import java.util.ArrayList;
import java.util.Iterator;

/**
 * @author Emilio Palumbo (emiliopalumbo@gmail.com)
 */

/**
 * MSIterator over an indexed BAM file. The class wraps a <code>SAMRecordIterator</code> wich is used to read over
 * the file.It returns records sorted following the SortOrder defined in the constructor.
 */
public class SAMMappingSortedIterator implements MSIterator<SAMMapping>{

    static final int DEFAULT_THRESHOLD = 5;
    static final boolean DEFAULT_READ_ALL = false;

    private SAMRecordIterator wrappedIterator;
    private SAMMapping mapping;
    private ArrayList<SAMMapping> mappings;
    private SAMFileHeader.SortOrder sortOrder;
    private final long maxRecordsInRam;
    private int currPos, markedPos;
    private boolean allReads;

    /**
     * Costruct an instance of the class.
     * @param iterator <code>SAMMappingIterator</code> used for iterating over the file
     * @param header SAM/BAM file header
     * @param maxRecordsInRam max number of records to be loaded in ram
     */
    public SAMMappingSortedIterator(SAMRecordIterator iterator, SAMFileHeader header, long maxRecordsInRam) {
        this(iterator,header, maxRecordsInRam, DEFAULT_READ_ALL);
    }

    /**
     * Costruct an instance of the class.
     * @param iterator <code>SAMMappingIterator</code> used for iterating over the file
     * @param header SAM/BAM file header
     * @param maxRecordsInRam max number of records to be loaded in ram
     */
    public SAMMappingSortedIterator(SAMRecordIterator iterator, SAMFileHeader header, long maxRecordsInRam, boolean allReads) {
        this.wrappedIterator = getSortedIterator(iterator, header);
        this.sortOrder = header.getSortOrder();
        this.maxRecordsInRam = maxRecordsInRam;
        this.currPos = this.markedPos = -1;
        this.allReads = allReads;
        readChunk();
    }

    /**
     * Returns the sorted iterator
     * @param iterator <code>SAMMappingIterator<code/> opened on the BAM file
     * @param header BAM file header
     * @return a sorted <code>SAMRecordIterator</code>
     */
    private SAMRecordIterator getSortedIterator(final SAMRecordIterator iterator, final SAMFileHeader header) {
        PipedInputStream pip = new PipedInputStream();

        try {
            PipedOutputStream pop = new PipedOutputStream(pip);
            final BufferedOutputStream out = new BufferedOutputStream(pop);
            final int mrec = (int) (maxRecordsInRam / 2);

            new Thread(
                    new Runnable(){
                        public void run(){
                            SAMFileWriter writer = new SAMFileWriterFactory().setTempDirectory(FileHelper.tempDirectory).setMaxRecordsInRam(mrec).makeSAMWriter(header, false, out);
                            SAMRecord rec = null;
                            while(iterator.hasNext()) {
                                rec = iterator.next();
                                writer.addAlignment(rec);
                            }
                            iterator.close();
                            writer.close();
                        }
                    }
            ).start();
        } catch (IOException e) {
            System.err.println("[SAMWriter THREAD]");
            e.printStackTrace();
        }

        SAMFileReader r = new SAMFileReader(pip);
        return r.iterator();
    }

    /**
     * Reads at maximum <code>maxRecordsInRam</code> records and load them in memory.
     */
    private void readChunk() {
        SAMRecord record;
//        SAMMapping mapping;

        if (mappings == null)
            mappings = new ArrayList<SAMMapping>();
        else {
            this.currPos = this.markedPos = -1;
            if (mapping != null && mappings.get(mappings.size()-1).equals(mapping)) {
                mapping=null;
            }
            mappings.clear();
            if (mapping!=null)
                mappings.add(mapping);
        }

        while(wrappedIterator.hasNext() && mappings.size() < maxRecordsInRam/2) {
            record = wrappedIterator.next();

            if (!allReads && record.getReadUnmappedFlag())
                continue;

            mapping = new SAMMapping(record, getSuffix(record));

            if (mappings.size()>1 && !mapping.getName().equals(mappings.get(mappings.size()-1).getName()) && (maxRecordsInRam-mappings.size())<= DEFAULT_THRESHOLD) {
                break;
            }

            mappings.add(mapping);
        }
        if (!wrappedIterator.hasNext()) {
            wrappedIterator.close();
        }
    }

    @Override
    public void mark() {
        markedPos = currPos;
    }

    @Override
    public void reset() {
        if (markedPos >-1) {
            currPos = markedPos;
            markedPos = -1;
        } else {
            currPos = 0;
        }
    }

    @Override
    public void setAtStart() {
        currPos = 0;
        markedPos = -1;
    }

    @Override
    public void clear() {
        wrappedIterator.close();
    }

    @Override
    public Iterator<Mapping> getMates(Mapping firstMate, UniversalReadDescriptor descriptor) {
        ArrayList<Mapping> mappings = new ArrayList<Mapping>();
        UniversalReadDescriptor.Attributes attr1 = null, attr2 = null;
        attr1 = getAttributes(firstMate,descriptor,attr1);
        if (attr1.flag == 2)
            return mappings.iterator();
        this.mark();
        while (this.hasNext()) {
            Mapping currentMapping = this.next();
            attr2 = getAttributes(currentMapping,descriptor,attr2);
            if (!attr1.id.equals(attr2.id))
                break;
            if (attr2 == null || attr2.flag == 1)
                continue;
            mappings.add(currentMapping);
        }
        this.reset();
        return mappings.iterator();
    }

    /**
     * Returns the attributes of a mapping given a read descriptor.
     * @param mapping the mapping
     * @param desc the read descriptor for the mapping
     * @param attributes an optional <code>Attributes</code> instance for object re-use
     * @return the instance of <code>Attributes</code>
     */
    private UniversalReadDescriptor.Attributes getAttributes(Mapping mapping, UniversalReadDescriptor desc, UniversalReadDescriptor.Attributes attributes) {

        CharSequence tag= mapping.getName();
        attributes= desc.getAttributes(tag, attributes);
        if (attributes == null) {
            Log.warn("Error in read ID: could not parse read identifier " + tag);
            return null;
        }
        if (desc.isPaired()&& attributes.flag<= 0) {
            Log.warn("Error in read ID: could not find mate in " + tag);
            return null;
        }
        if (desc.isStranded()&& attributes.strand< 0) {
            Log.warn("Error in read ID: could not find strand in " + tag);
            return null;
        }
        return attributes;
    }

    @Override
    public Iterator<SAMMapping> iterator() {
        return this;
    }

    @Override
    public boolean hasNext() {
        if (mappings==null)
            return false;
        if (currPos>=mappings.size()-1) {
            if (wrappedIterator.hasNext()) {
                readChunk();
            } else {
                wrappedIterator.close();
                return false;
            }
        }
        return true;
    }

    @Override
    public SAMMapping next() {
        return (mappings.get(++currPos));
    }

    @Override
    public void remove() {
        wrappedIterator.remove();
    }

    /**
     * Returns a suffix for the SAM record read name for paired end read. I returns "/1"
     * for the first mate and "/2" for the second.
     * @param record the SAM record
     * @return the suffix
     */
    private String getSuffix(SAMRecord record) {
        if (record.getReadPairedFlag()) {
            return record.getFirstOfPairFlag()?"/1":"/2";
        }
        return "";
    }
}
