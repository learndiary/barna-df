package barna.io.sam;

import barna.io.MSIterator;
import barna.model.sam.SAMMapping;
import net.sf.samtools.SAMRecordIterator;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;

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
        SAMMapping mapping;
        while(wrappedIterator.hasNext()) {
            mapping = new SAMMapping(wrappedIterator.next());
            if (mapping.getChromosome().equals(chromosome)) {
                if (mapping.getStart()>=start && mapping.getEnd()<=end) {
                    if (mappings==null)
                        mappings = new ArrayList<SAMMapping>();
                    mappings.add(mapping);
                }
            }
        }
        //Arrays.sort(mappings.toArray());
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
