package barna.io.sam;

import barna.commons.log.Log;
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
    UniversalReadDescriptor descriptor;
    int currPos, markedPos;

    public SAMMappingIterator(SAMRecordIterator iterator, UniversalReadDescriptor descriptor) {
        this.wrappedIterator = iterator;
        this.descriptor = descriptor;
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

            if (mappings == null)
                mappings = new ArrayList<SAMMapping>();
            mappings.add(mapping);
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
    public Iterator<Mapping> getMates(Mapping firstMate, UniversalReadDescriptor descriptor) {
        ArrayList<Mapping> mappings = new ArrayList<Mapping>();
        UniversalReadDescriptor.Attributes attr1 = null, attr2 = null;
        attr1 = getAttributes(firstMate,descriptor,attr1);
        if (attr1.flag == 2)
            return mappings.iterator();
        this.mark();
        while (this.hasNext()) {
            Mapping currentMapping = this.next();
            attr2 = getAttributes(currentMapping,descriptor,attr2);
            if (!attr1.id.equals(attr2.id))
                break;
            if (attr2 == null || attr2.flag == 1)
                continue;
            mappings.add(currentMapping);
        }
        this.reset();
        return mappings.iterator();
    }

    private UniversalReadDescriptor.Attributes getAttributes(Mapping mapping, UniversalReadDescriptor desc, UniversalReadDescriptor.Attributes attributes) {

        CharSequence tag= mapping.getName();
        attributes= desc.getAttributes(tag, attributes);
        if (attributes == null) {
            Log.warn("Error in read ID: could not parse read identifier " + tag);
            return null;
        }
        if (desc.isPaired()&& attributes.flag<= 0) {
            Log.warn("Error in read ID: could not find mate in " + tag);
            return null;
        }
        if (desc.isStranded()&& attributes.strand< 0) {
            Log.warn("Error in read ID: could not find strand in " + tag);
            return null;
        }
        return attributes;
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
        if (descriptor.isPaired()) {
            char sep = descriptor.toString().charAt(descriptor.toString().indexOf("{MATE}")-1);
            return record.getFirstOfPairFlag()?sep+"1":sep+"2";
        } else {
            //to get it working also with paired-end data mapped as single end
            return record.getFirstOfPairFlag()?"/1":"/2";
        }
    }
}
