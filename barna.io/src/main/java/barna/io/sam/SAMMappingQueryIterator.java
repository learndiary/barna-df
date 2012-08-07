package barna.io.sam;

import barna.io.MSIterator;
import barna.io.rna.UniversalReadDescriptor;
import barna.model.Mapping;
import barna.model.sam.SAMMapping;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

import java.io.File;
import java.util.ArrayList;
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

    @Override
    public Iterator<Mapping> getMates(Mapping firstMate, UniversalReadDescriptor descriptor) {
        ArrayList<Mapping> mates = new ArrayList<Mapping>();
        SAMFileReader reader = new SAMFileReader(mappingFile);
        SAMRecordIterator iter = reader.query(firstMate.getChromosome().toString(),firstMate.getStart()+1,firstMate.getEnd(),true);
        SAMRecord rec = null;
        while (iter.hasNext()) {
            rec = iter.next();
            if (firstMate.getName().toString().contains(rec.getReadName()))
                break;
        }
        iter.close();
        mate = reader.queryMate(rec);
        if (mate != null) {
            if (rec.getSecondOfPairFlag()&&mate.getMateAlignmentStart()==rec.getAlignmentStart()) {
                reader.close();
                return mates.iterator();
            }
            mates.add(new SAMMapping(mate, getSuffix(mate)));
        }
        reader.close();
        return mates.iterator();
    }

    private void getNext() {
        nextMapping = null;
        while (nextMapping==null&&wrappedIterator.hasNext()) {
            SAMRecord rec;
            rec = wrappedIterator.next();
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
