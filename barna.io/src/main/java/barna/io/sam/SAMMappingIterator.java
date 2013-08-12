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

    static final boolean DEFAULT_ALL_READS = false;
    static final boolean DEFAULT_PRIMARY_ONLY = false;
    static final boolean DEFAULT_MATES_ONLY = false;
    static final boolean DEFAULT_UNIQUE_ONLY = false;

    private SAMRecordIterator wrappedIterator;
    private ArrayList<SAMMapping> mappings;
    private int currPos, markedPos;
    private boolean allReads;
    /**
     * Filter reads by their score
     */
    private int scoreFilter;

    /**
     * Only keep primary alignments
     */
    private boolean primaryOnly;

    /**
     * Only use pairing information from BAM file
     */
    private boolean matesOnly;

    /**
     * Only keep unique alignments
     */
    private boolean uniqueOnly;

    public SAMMappingIterator(SAMRecordIterator iterator) {
        this(iterator, DEFAULT_ALL_READS);
    }

    public SAMMappingIterator(SAMRecordIterator iterator, boolean allReads) {
        this(iterator, allReads, -1, DEFAULT_PRIMARY_ONLY, DEFAULT_MATES_ONLY, DEFAULT_UNIQUE_ONLY);
    }

    public SAMMappingIterator(SAMRecordIterator iterator, boolean allReads, int scoreFilter) {
        this(iterator, allReads, scoreFilter, DEFAULT_PRIMARY_ONLY, DEFAULT_MATES_ONLY, DEFAULT_UNIQUE_ONLY);
    }

    public SAMMappingIterator(SAMRecordIterator iterator, boolean allReads, int scoreFilter, boolean primaryOnly, boolean matesOnly, boolean uniqueOnly) {
        this.wrappedIterator = iterator;
        this.currPos = this.markedPos = -1;
        this.allReads = allReads;
        this.scoreFilter = scoreFilter;
        this.primaryOnly = primaryOnly;
        this.matesOnly = primaryOnly ? false : matesOnly;
        this.uniqueOnly = uniqueOnly;
        init();
    }

    private void init() {
        SAMRecord record;
        SAMMapping mapping;
        List<SAMMapping> tmp = null;

        while(wrappedIterator.hasNext()) {
            record = wrappedIterator.next();

            if (!allReads && record.getReadUnmappedFlag())
                continue;
            if (primaryOnly && record.getNotPrimaryAlignmentFlag())
                continue;

            mapping = new SAMMapping(record, getSuffix(record));

            if (uniqueOnly && !mapping.isUnique())
                continue;

            if (mappings == null)
                mappings = new ArrayList<SAMMapping>();
            if((scoreFilter < 0 || mapping.getScore() >= this.scoreFilter)){
                mappings.add(mapping);
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
    public Iterator<Mapping> getMates(Mapping firstMate, UniversalReadDescriptor descriptor) {
        ArrayList<Mapping> mappings = new ArrayList<Mapping>();
        UniversalReadDescriptor.Attributes attr1 = null, attr2 = null;
        attr1 = getAttributes(firstMate,descriptor,attr1);
        if (attr1.flag == 2)
            return mappings.iterator();
        if (!((SAMMapping)firstMate).isProperlyPaired())
            return mappings.iterator();
        this.mark();
        while (this.hasNext()) {
            SAMMapping currentMapping = this.next();
            if (!currentMapping.isProperlyPaired())
                continue;
            attr2 = getAttributes(currentMapping,descriptor,attr2);
            if (!attr1.id.equals(attr2.id))
                break;
            if (attr2 == null || attr2.flag == 1)
                continue;
            if (!this.matesOnly || currentMapping.isMateOf((SAMMapping)firstMate))
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
        if (record.getReadPairedFlag()) {
            return record.getFirstOfPairFlag()?"/1":"/2";
        }
        return "";
    }
}
