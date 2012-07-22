package barna.io.sam;

import barna.io.MSIterator;
import barna.model.sam.SAMMapping;
import net.sf.samtools.*;

import java.io.File;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/**
 * @author Emilio Palumbo (emiliopalumbo@gmail.com)
 */
public class SAMMappingQueryIterator implements MSIterator<SAMMapping> {

    private final File inputFile;
    private final SAMRecordIterator wrappedIterator;
    private SAMRecord nextRecord, nextMate;
    private SAMMapping nextMapping, nextMappingMate;
    int currPos, markedPos,start,end;
    private boolean returnedMapping, returnedMate;
    private int altMappingIndex, altMateIndex;

    public SAMMappingQueryIterator(File inputFile, SAMRecordIterator wrappedIterator, int start, int end) {
        this.inputFile = inputFile;
        this.wrappedIterator = wrappedIterator;
        this.start = start;
        this.end = end;
        this.currPos = this.markedPos = -1;
        this.altMappingIndex = this.altMateIndex = -1;
        getNext();
    }

    private List<SAMMapping> addRecord(SAMRecord record, List<SAMMapping> tmp) {
        SAMMapping mapping;
        String name = record.getReadName();
        if (record.getFirstOfPairFlag())
            record.setReadName(name+"/1");
        else {
            if (record.getSecondOfPairFlag())
                record.setReadName(name+"/2");
        }
        mapping = new SAMMapping(record);

        if (tmp==null)
            tmp = new ArrayList<SAMMapping>();
        tmp.add(mapping);

        if (mapping.hasAlternates()) {
            tmp.addAll(mapping.getAlternates());
        }

        return tmp;
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
    public Iterator<SAMMapping> iterator() {
        return this;
    }

    @Override
    public boolean hasNext() {
        return nextRecord!=null;
    }

    @Override
    public SAMMapping next() {
        if (!returnedMapping) {
            returnedMapping = true;
            return nextMapping;
        }
        if (nextMappingMate!=null) {
            if (!returnedMate) {
                returnedMate = true;
                return nextMappingMate;
            }
            if (nextMappingMate.hasAlternates()) {
                if (altMateIndex<nextMappingMate.getAlternates().size()-1)
                    return nextMappingMate.getAlternates().get(++altMateIndex);
            }
        }
        if (nextMapping.hasAlternates()) {
            altMateIndex = -1;
            if (altMappingIndex<nextMapping.getAlternates().size()-1)
                return nextMapping.getAlternates().get(++altMappingIndex);
        }
        getNext();
        returnedMapping = true;
        return nextMapping;
    }

    @Override
    public void remove() {
        wrappedIterator.remove();
    }

    private void getNext() {
        nextRecord = nextMate = null;
        //nextMapping = nextMappingMate = null;
        returnedMapping = false;
        returnedMate = false;
        altMappingIndex = -1;
        altMateIndex = -1;
        while (nextRecord==null&&wrappedIterator.hasNext()) {
            SAMRecord rec = wrappedIterator.next();
//            if (rec.getFirstOfPairFlag()) {
                nextRecord = rec;
                nextMapping = new SAMMapping(rec,"/1");
                break;
//            }
        }
        if (nextRecord!=null) {// && nextRecord.getFirstOfPairFlag()) {
            SAMFileReader reader=null;
            SAMRecord mate=null;
            if (reader==null)
                reader = new SAMFileReader(inputFile);
            mate = reader.queryMate(nextRecord);
            if (mate != null) {
                if (!(nextRecord.getSecondOfPairFlag() && mate.getAlignmentStart()>=start&&mate.getAlignmentStart()<end)) {
                    nextMate = mate;
                    nextMappingMate = new SAMMapping(mate,"/2");
                }
            }
            reader.close();
        }
    }
}
