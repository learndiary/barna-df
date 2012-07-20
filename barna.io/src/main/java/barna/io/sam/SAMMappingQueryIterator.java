package barna.io.sam;

import barna.io.MSIterator;
import barna.model.sam.SAMMapping;
import net.sf.samtools.*;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

/**
 * @author Emilio Palumbo (emiliopalumbo@gmail.com)
 */
public class SAMMappingQueryIterator implements MSIterator<SAMMapping> {

    private final File inputFile;
    private final SAMRecordIterator wrappedIterator;
    ArrayList<SAMMapping> mappings;
    int currPos, markedPos,start,end;

    public SAMMappingQueryIterator(File inputFile, SAMRecordIterator wrappedIterator, int start, int end) {
        this.inputFile = inputFile;
        this.wrappedIterator = wrappedIterator;
        this.start = start;
        this.end = end;
        this.currPos = this.markedPos = -1;
        init();
    }

    private void init() {
        SAMRecord record;
        SAMMapping mapping;
        List<SAMMapping> tmp=null;
        SAMFileReader reader=null;
        SAMRecord mate =null;
        while(wrappedIterator.hasNext()) {
            record = wrappedIterator.next();

            if (record.getReadUnmappedFlag())
                continue;

            int mateStart = record.getMateAlignmentStart();

            if (record.getMateReferenceName().equals(record.getReferenceName()) && (mateStart>=end || mateStart<start)) {
                if (reader==null)
                    reader = new SAMFileReader(inputFile);
                mate=reader.queryMate(record);
                if (mate != null)
                    tmp = addRecord(mate, tmp);
            }

            tmp = addRecord(record, tmp);

            for (SAMMapping m : tmp) {
                if (mappings==null)
                   mappings=new ArrayList<SAMMapping>();
                if (m.getStart()>=start && m.getEnd()<=end)
                    mappings.add(m);
            }
            tmp.clear();
        }
        if (mappings!=null) {
            Collections.sort(mappings, new SAMMapping.SAMIdComparator());
        }

        wrappedIterator.close();
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
        if (mappings==null)
            return false;
        return currPos<mappings.size()-1;
    }

    @Override
    public SAMMapping next() {
        return (mappings.get(++currPos));
    }

    @Override
    public void remove() {
        wrappedIterator.remove();
    }
}
