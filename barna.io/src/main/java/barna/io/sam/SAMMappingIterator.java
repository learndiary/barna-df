package barna.io.sam;

import barna.io.MSIterator;
import barna.model.sam.SAMMapping;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordComparator;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.SAMRecordQueryNameComparator;

import java.util.*;

/**
 * @author Emilio Palumbo (emiliopalumbo@gmail.com)
 */
public class SAMMappingIterator implements MSIterator<SAMMapping>{

    SAMRecordIterator wrappedIterator;
    ArrayList<SAMMapping> mappings;
    int currPos, markedPos;

    String chromosome;
    int start, end;

    public SAMMappingIterator(String chromosome, int start, int end, SAMRecordIterator iterator) {
        this.chromosome = chromosome;
        this.end = end;
        this.start = start;
        this.wrappedIterator = iterator;
        this.currPos = markedPos = -1;
        init();
    }

    private void init() {
        SAMRecord record;
        SAMMapping mapping;
        List<SAMMapping> tmp = null;
//        ArrayList<SAMRecord> tmp=new ArrayList<SAMRecord>();
        while(wrappedIterator.hasNext()) {
            //mapping = new SAMMapping(wrappedIterator.next());
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
                tmp=new ArrayList<SAMMapping>();
            tmp.add(mapping);
            int i=0;
            for (SAMMapping m : tmp) {
                if (i>30)
                    break;
                if (m.getChromosome().equals(chromosome)) {
                    if (m.getStart()>=start && m.getEnd()<=end) {
                        if (mappings == null)
                            mappings = new ArrayList<SAMMapping>();
                        mappings.add(m);
                    }
                }
                i++;
            }
            tmp.clear();
        }
        //Collections.sort(tmp, new SAMRecordQueryNameComparator());
        Collections.sort(mappings, new SAMMapping.SAMIdComparator());

//        for (SAMRecord r : tmp) {
//            if (mappings==null)
//                mappings = new ArrayList<SAMMapping>();
//            SAMMapping sam = new SAMMapping(r);
//            mappings.add(sam);
//            if (sam.hasAlternates()) {
//                for (SAMMapping m : sam.getAlternates()) {
//                    if (m.getChromosome().equals(chromosome)) {
//                        if (m.getStart()>=start&&m.getEnd()<=end) {
//                            mappings.add(m);
//                        }
//                    }
//                }
//            }
//        }
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
    }

    @Override
    public Iterator<SAMMapping> iterator() {
        return this;
    }

    @Override
    public boolean hasNext() {
        return currPos<mappings.size()-1;
    }

    @Override
    public SAMMapping next() {
        return (mappings.get(++currPos));
    }

    @Override
    public void remove() {
    }
}
