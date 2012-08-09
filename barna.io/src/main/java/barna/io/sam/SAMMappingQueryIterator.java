package barna.io.sam;

import barna.io.rna.UniversalReadDescriptor;
import barna.model.Mapping;
import barna.model.sam.SAMMapping;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

import java.io.File;
import java.util.ArrayList;
import java.util.Iterator;

/**
 * @author Emilio Palumbo (emiliopalumbo@gmail.com)
 */
public class SAMMappingQueryIterator implements Iterator<SAMMapping> {

    private SAMRecordIterator wrappedIterator;
    private File mappingFile;
    private final UniversalReadDescriptor descriptor;
    private SAMRecord lastRecord;

    public SAMMappingQueryIterator(File inputFile, SAMRecordIterator wrappedIterator, UniversalReadDescriptor descriptor) {
        this.mappingFile = inputFile;
        this.wrappedIterator = wrappedIterator;
        this.descriptor = descriptor;
        getNext();
    }

    @Override
    public boolean hasNext() {
        return wrappedIterator.hasNext();
    }

    @Override
    public SAMMapping next() {
        SAMRecord rec = getNext();
        lastRecord = rec;
        return new SAMMapping(rec,getSuffix(rec));
    }

    @Override
    public void remove() {
        wrappedIterator.remove();
    }

    //very slow because always access to disk
    public Iterator<Mapping> getMates(Mapping firstMate, UniversalReadDescriptor descriptor) {
        ArrayList<Mapping> mates = new ArrayList<Mapping>();
        SAMFileReader reader = new SAMFileReader(mappingFile);
        SAMRecord mate = reader.queryMate(lastRecord);
        if (mate != null) {
            if (lastRecord.getSecondOfPairFlag()&&mate.getMateAlignmentStart()==lastRecord.getAlignmentStart()) {
                reader.close();
                return mates.iterator();
            }
            mates.add(new SAMMapping(mate, getSuffix(mate)));
        }
        reader.close();
        return mates.iterator();
    }

    private SAMRecord getNext() {
        SAMRecord rec = null;
        while (rec==null&&wrappedIterator.hasNext()) {
            rec = wrappedIterator.next();
            if (rec.getReadUnmappedFlag()) {
                continue;
            }
        }
        return rec;
    }

    private String getSuffix(SAMRecord record) {
        if (descriptor.isPaired()) {
            char sep = descriptor.toString().charAt(descriptor.toString().indexOf("{MATE}")-1);
            return record.getFirstOfPairFlag()?sep+"1":sep+"2";
        }
        return "";
    }
}
