package barna.io.sam;

import barna.io.FileHelper;
import barna.io.MSIterator;
import barna.model.Mapping;
import barna.model.rna.UniversalReadDescriptor;
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
    static final boolean DEFAULT_PRIMARY_ONLY = false;
    static final boolean DEFAULT_MATES_ONLY = false;
    static final boolean DEFAULT_UNIQUE_ONLY = false;

    private SAMRecordIterator wrappedIterator;
    private SAMMapping mapping;
    private ArrayList<SAMMapping> mappings;
    private SAMFileHeader.SortOrder sortOrder;
    private int maxRecordsInRam;
    private int currPos, markedPos;
    private boolean allReads;
    private int scoreFilter;
    private boolean primaryOnly;
    private boolean matesOnly;
    private boolean uniqueOnly;
    private SAMFileReader.ValidationStringency validationStringency = SAMFileReader.ValidationStringency.DEFAULT_STRINGENCY;

    /**
     * Costruct an instance of the class.
     * @param iterator <code>SAMMappingIterator</code> used for iterating over the file
     * @param header SAM/BAM file header
     * @param maxRecordsInRam max number of records to be loaded in ram
     */
    public SAMMappingSortedIterator(SAMRecordIterator iterator, SAMFileHeader header, int maxRecordsInRam) {
        this(iterator,header, maxRecordsInRam, DEFAULT_READ_ALL);
    }

    /**
     * Costruct an instance of the class.
     * @param iterator <code>SAMMappingIterator</code> used for iterating over the file
     * @param header SAM/BAM file header
     * @param maxRecordsInRam max number of records to be loaded in ram
     */
    public SAMMappingSortedIterator(SAMRecordIterator iterator, SAMFileHeader header, int maxRecordsInRam, boolean allReads) {
        this(iterator, header, maxRecordsInRam, allReads, -1, DEFAULT_PRIMARY_ONLY, DEFAULT_MATES_ONLY, DEFAULT_UNIQUE_ONLY,
                SAMFileReader.ValidationStringency.DEFAULT_STRINGENCY);
    }

    public SAMMappingSortedIterator(SAMRecordIterator iterator, SAMFileHeader header,
                                    int maxRecordsInRam, boolean allReads, int scoreFilter, boolean primaryOnly,
                                    boolean matesOnly, boolean uniqueOnly, SAMFileReader.ValidationStringency validationStringency) {
        this.sortOrder = header.getSortOrder();
        this.maxRecordsInRam = maxRecordsInRam;
        this.currPos = this.markedPos = -1;
        this.allReads = allReads;
        this.scoreFilter = scoreFilter;
        this.wrappedIterator = getSortedIterator(iterator, header);
        this.primaryOnly = primaryOnly;
        this.matesOnly = primaryOnly ? false : matesOnly;
        this.uniqueOnly = uniqueOnly;
        this.validationStringency = validationStringency;
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
            final int mrec = maxRecordsInRam / 2;

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
        r.setValidationStringency(this.validationStringency);
        return r.iterator();
    }

    /**
     * Reads at maximum <code>maxRecordsInRam</code> records and load them in memory.
     */
    private void readChunk() {
        SAMRecord record;

        if (mappings == null)
            mappings = new ArrayList<SAMMapping>();
        else {
            this.currPos = this.markedPos = -1;
            if (mapping != null && mappings.get(mappings.size()-1).equals(mapping)) {
                mapping=null;
            }
            mappings.clear();
            if (mapping!=null && (this.scoreFilter < 0 || mapping.getScore() >= this.scoreFilter) && (!this.primaryOnly || mapping.isPrimary()))
                mappings.add(mapping);
        }

        while(wrappedIterator.hasNext() && mappings.size() < maxRecordsInRam/2) {
            record = wrappedIterator.next();

            if (!allReads && record.getReadUnmappedFlag())
                continue;

            if (this.primaryOnly && record.getNotPrimaryAlignmentFlag())
                continue;

            mapping = new SAMMapping(record);

            if (uniqueOnly && !mapping.isUnique())
                continue;

            if (mappings.size()>1 && !mapping.getName(true).equals(mappings.get(mappings.size()-1).getName(true)) && ((maxRecordsInRam/2)-mappings.size())<= DEFAULT_THRESHOLD) {
                break;
            }
            if ((this.scoreFilter < 0 || mapping.getScore() >= this.scoreFilter)) {
                mappings.add(mapping);
            }
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
    public Iterator<Mapping> getMates(Mapping firstMate) {
        ArrayList<Mapping> mappings = new ArrayList<Mapping>();
        UniversalReadDescriptor.Attributes attr1 = null, attr2 = null;
        if ((firstMate.getMateFlag() == 2)
                || (!((SAMMapping) firstMate).isProperlyPaired()))
            return mappings.iterator();
        this.mark();
        while (this.hasNext()) {
            SAMMapping currentMapping = this.next();
            if (!currentMapping.isProperlyPaired())
                continue;
            if (!firstMate.getName(false).equals(currentMapping.getName(false)))
                break;
            if (currentMapping.getMateFlag() == 1)
                continue;
            if ((!this.matesOnly) || currentMapping.isMateOf((SAMMapping)firstMate)) {
                mappings.add(currentMapping);
                if (this.matesOnly) {
                    this.mappings.remove(currPos--);
                    break;
                }
            }
        }
        this.reset();
        return mappings.iterator();
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

    @Override
    public int size() {
        if (mappings== null)
            return 0;
        return mappings.size();
    }


}
