package barna.io.sam;

import barna.io.MSIterator;
import barna.model.sam.SAMMapping;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMRecordIterator;

import java.util.Iterator;

/**
 * @author Emilio Palumbo (emiliopalumbo@gmail.com)
 */
public class SAMMappingQueryIterator implements MSIterator<SAMMapping> {

    private final SAMRecordIterator wrappedIterator;

    public SAMMappingQueryIterator(SAMRecordIterator wrappedIterator) {
        this.wrappedIterator = wrappedIterator;
    }

    @Override
    public void mark() {
    }

    @Override
    public void reset() {
    }

    @Override
    public void setAtStart() {
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
        return wrappedIterator.hasNext();
    }

    @Override
    public SAMMapping next() {
        if (hasNext()) {
            return new SAMMapping(wrappedIterator.next());
        }
        return null;
    }

    @Override
    public void remove() {
        wrappedIterator.remove();
    }
}
