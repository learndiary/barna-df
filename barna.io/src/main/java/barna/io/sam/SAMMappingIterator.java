package barna.io.sam;

import barna.io.MSIterator;
import barna.io.rna.UniversalReadDescriptor;
import barna.model.Mapping;
import barna.model.sam.SAMMapping;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

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
        this.currPos = this.markedPos = -1;
        init();
    }

    private void init() {
        SAMRecord record;
        SAMMapping mapping;
        List<SAMMapping> tmp = null;

        while(wrappedIterator.hasNext()) {
            record = wrappedIterator.next();

            if (record.getReadUnmappedFlag())
                continue;

            mapping = new SAMMapping(record, getSuffix(record));

           if (mapping.getChromosome().equals(chromosome)) {
                if (mapping.getStart()>=start && mapping.getEnd()<=end) {
                    if (mappings == null)
                        mappings = new ArrayList<SAMMapping>();
                    mappings.add(mapping);
                }
            }
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
    public Iterator<Mapping> getMates(Mapping firstMapping, UniversalReadDescriptor descriptor) {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
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

    private String getSuffix(SAMRecord record) {
        return record.getFirstOfPairFlag()?"/1":"/2";
    }}
