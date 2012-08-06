package barna.io.sam;

import barna.io.MSIterator;
import barna.io.rna.UniversalReadDescriptor;
import barna.model.Mapping;
import barna.model.sam.SAMMapping;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

import java.io.File;
import java.util.Iterator;

/**
 * @author Emilio Palumbo (emiliopalumbo@gmail.com)
 */
public class SAMMappingQueryIterator implements MSIterator<SAMMapping> {

    private SAMRecordIterator wrappedIterator;
    private File mappingFile;
    private SAMRecord nextRecord;
    private SAMMapping nextMapping;
    private SAMRecord mate;
    private int start,end;
    private boolean marked;
    private boolean paired;

    public SAMMappingQueryIterator(File inputFile, SAMRecordIterator wrappedIterator, int start, int end, boolean isPaired) {
        this.mappingFile = inputFile;
        this.wrappedIterator = wrappedIterator;
        this.start = start;
        this.end = end;
        marked = false;
        paired = isPaired;
        getNext();
    }

    @Override
    public void mark() {
        marked = true;
    }

    @Override
    public void reset() {
        if (marked) {
            marked = false;
        }
    }

    @Override
    public void setAtStart() {
        wrappedIterator.close();
    }

    @Override
    public void clear() {
        wrappedIterator.close();
    }

    @Override
    public Iterator<Mapping> getMates(Mapping firstMapping, UniversalReadDescriptor descriptor) {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
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
        if (marked && mate !=null)
            nextMapping = new SAMMapping(mate, getSuffix(mate));
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
        while (nextMapping==null&&wrappedIterator.hasNext()) {
            SAMRecord rec;
            if (paired && nextRecord!=null && nextRecord.getFirstOfPairFlag()) {
                getMate();
            }
            rec = wrappedIterator.next();
            if (rec.getReadUnmappedFlag()) {
                continue;
            }
            nextRecord = rec;
            nextMapping = new SAMMapping(rec,getSuffix(rec));
        }
    }

    private void getMate() {
        SAMFileReader reader = new SAMFileReader(mappingFile);
        mate = reader.queryMate(nextRecord);
        reader.close();
    }

    private String getSuffix(SAMRecord record) {
        return record.getFirstOfPairFlag()?"/1":"/2";
    }
}
