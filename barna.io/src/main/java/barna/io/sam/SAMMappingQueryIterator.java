package barna.io.sam;

import barna.io.MSIterator;
import barna.model.sam.SAMMapping;
import net.sf.samtools.*;

import java.io.File;
import java.util.Iterator;

/**
 * @author Emilio Palumbo (emiliopalumbo@gmail.com)
 */
public class SAMMappingQueryIterator implements MSIterator<SAMMapping> {

    private final File inputFile;
    private SAMRecordIterator wrappedIterator;
    private SAMRecord nextRecord;
    private SAMMapping nextMapping;
    private int start,end;
    private boolean marked;

    public SAMMappingQueryIterator(File inputFile, SAMRecordIterator wrappedIterator, int start, int end) {
        this.inputFile = inputFile;
        this.wrappedIterator = wrappedIterator;
        this.start = start;
        this.end = end;
        marked = false;
        getNext();
    }

    @Override
    public void mark() {
        marked = true;
    }

    @Override
    public void reset() {
        if (marked)
            marked = false;
    }

    @Override
    public void setAtStart() {
        wrappedIterator.close();
//        SAMFileReader reader = new SAMFileReader(inputFile);
//        wrappedIterator = reader.query("",start,end,true);
    }

    @Override
    public void clear() {
        wrappedIterator.close();
    }

    @Override
    public Iterator<SAMMapping> iterator() {
        return this;
    }

    @Override
    public boolean hasNext() {
        return nextMapping!=null;
    }

    @Override
    public SAMMapping next() {
        SAMMapping ret = nextMapping;
        getNext();
        return ret;
    }

    @Override
    public void remove() {
        wrappedIterator.remove();
    }

    private void getNext() {
        nextMapping = null;
        if (marked) {
            if (nextRecord!=null) {// && nextRecord.getFirstOfPairFlag()) {
                SAMFileReader reader=null;
                SAMRecord mate=null;
                if (reader==null)
                    reader = new SAMFileReader(inputFile);
                mate = reader.queryMate(nextRecord);
                if (mate != null) {
                    if (!(nextRecord.getSecondOfPairFlag() && mate.getAlignmentStart()>=start&&mate.getAlignmentStart()<end)) {
                        nextMapping = new SAMMapping(mate,getSuffix(mate));
                    }
                }
                reader.close();
            }
        }
        while (nextMapping==null&&wrappedIterator.hasNext()) {
            SAMRecord rec = wrappedIterator.next();
            if (rec.getReadUnmappedFlag()) {
                continue;
            }
            nextRecord = rec;
            nextMapping = new SAMMapping(rec,getSuffix(rec));
        }
    }

    private String getSuffix(SAMRecord record) {
        return record.getFirstOfPairFlag()?"/1":"/2";
    }
}
