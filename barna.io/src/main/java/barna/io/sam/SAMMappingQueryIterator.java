package barna.io.sam;

import barna.io.MSIterator;
import barna.model.sam.SAMMapping;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.SAMRecordQueryNameComparator;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

/**
 * @author Emilio Palumbo (emiliopalumbo@gmail.com)
 */
public class SAMMappingQueryIterator implements MSIterator<SAMMapping> {

    private final SAMRecordIterator wrappedIterator;
    ArrayList<SAMMapping> mappings;
    int currPos, markedPos,start,end;

    public SAMMappingQueryIterator(SAMRecordIterator wrappedIterator, int start, int end) {
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
        while(wrappedIterator.hasNext()) {
            record = wrappedIterator.next();

            if (record.getReadUnmappedFlag())
                continue;

            String name = record.getReadName();
            if (record.getMateNegativeStrandFlag())
                record.setReadName(name+"/2");
            else
                record.setReadName(name+"/1");
            mapping = new SAMMapping(record);

            if (mapping.hasAlternates()) {
                tmp = mapping.getAlternates();
            }
            if (tmp==null)
                tmp = new ArrayList<SAMMapping>();
            tmp.add(mapping);

            for (SAMMapping m : tmp) {
                if (mappings==null)
                   mappings=new ArrayList<SAMMapping>();
                mappings.add(m);
            }
            tmp.clear();
        }
        if (mappings!=null) {
            Collections.sort(mappings, new SAMMapping.SAMIdComparator());
        }

        wrappedIterator.close();
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
