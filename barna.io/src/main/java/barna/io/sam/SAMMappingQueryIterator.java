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
    private SAMMapping nextMapping, nextMappingMate, current;
    int currPos, markedPos,start,end;
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
        SAMMapping ret = current;
        getNext();
        return ret;
    }

    @Override
    public void remove() {
        wrappedIterator.remove();
    }

    private void getNext() {
        current = null;
        if (nextRecord!=null && nextMate==null) {// && nextRecord.getFirstOfPairFlag()) {
            SAMFileReader reader=null;
            SAMRecord mate=null;
            if (reader==null)
                reader = new SAMFileReader(inputFile);
            mate = reader.queryMate(nextRecord);
            if (mate != null) {
                if (!(nextRecord.getSecondOfPairFlag() && mate.getAlignmentStart()>=start&&mate.getAlignmentStart()<end)) {
                    nextMate = mate;
                    nextMappingMate = new SAMMapping(mate,getSuffix(mate));
                }
            }
            reader.close();
        }
        while (nextRecord==null&&wrappedIterator.hasNext()) {
            SAMRecord rec = wrappedIterator.next();
//            if (rec.getFirstOfPairFlag()) {
            nextRecord = rec;
            nextMapping = new SAMMapping(rec,getSuffix(rec));
//            }
        }
        if (nextMate != null) {
            if (altMateIndex == -1) {
                current = nextMappingMate;
                altMateIndex++;
                if (!nextMappingMate.hasAlternates()) {
                    altMateIndex = -1;
                }
            } else {
                current = nextMappingMate.getAlternates().get(altMateIndex++);
                if (altMateIndex==nextMappingMate.getAlternates().size()) {
                    altMateIndex = -1;
                    if (altMappingIndex==nextMapping.getAlternates().size()) {
                        altMappingIndex = altMateIndex = -1;
                        nextRecord = nextMate = null;
                    }
                }
            }
            if (altMateIndex == -1) {
                if (!nextMapping.hasAlternates()) {
                    altMappingIndex = -1;
                    nextRecord = nextMate = null;
                } else {
                    current = nextMapping.getAlternates().get(altMappingIndex++);
                }
            }
        }
        if (nextRecord!=null && altMappingIndex==-1) {
            current = nextMapping;
            altMappingIndex++;
        }
    }

    private String getSuffix(SAMRecord record) {
        return record.getFirstOfPairFlag()?"/1":"/2";
    }
}
