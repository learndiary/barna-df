/**
 * 
 */
package barna.io.sam;

import java.io.File;
import java.io.OutputStream;
import java.lang.instrument.Instrumentation;
import java.util.Iterator;

import barna.commons.log.Log;
import barna.io.AbstractFileIOWrapper;
import barna.io.MSIterator;
import barna.io.MappingReader;
import barna.io.rna.UniversalReadDescriptor;
import barna.model.Mapping;
import barna.model.constants.Constants;
import net.sf.samtools.*;

/**
 * @author Emilio Palumbo (emiliopalumbo@gmail.com)
 */
public class SAMReader extends AbstractFileIOWrapper implements
        MappingReader {

    public static final boolean CONTAINED_DEFAULT = false;

	private Mapping[] mappings;
    private SAMFileReader reader;
    private boolean contained;
    private MSIterator iter;
    private boolean paired;

    int countAll;
    int countEntire;
    int countSplit;
    int countReads;
    int countSkippedLines;

	/**
     * Creates an instance of the reader
	 * @param inputFile the file to read
	 */
	public SAMReader(File inputFile, boolean isPaired) {
		this(inputFile, CONTAINED_DEFAULT, isPaired);
	}

    /**
     * Creates an instance of the reader
     * @param inputFile the file to read
     * @param contained flag to decide whether return mappings which are
     *                  contained/overlapping the query region
     */
    public SAMReader(File inputFile, boolean contained, boolean isPaired) {
        super(inputFile);
        reader = new SAMFileReader(this.inputFile);
        this.contained = contained;
        paired = isPaired;
    }
	
	/**
	 * Creates an instance using a specific path to a file.	 .
	 * @param absolutePath path to the file the wrapper is based on
	 */
	public SAMReader(String absolutePath, boolean isPaired) {
		this(new File(absolutePath), isPaired);
	}

    @Override
    public void read() {
    }

	@Override
	public void write() {
		// TODO Auto-generated method stub
	}

	@Override
	public boolean isApplicable() {
		if (reader==null)
            return false;
        if (!reader.isBinary())
            return false;
        if (!reader.hasIndex())
            return false;
        return true;
	}

	/**
	 * @see barna.io.IOWrapper#sort(java.io.OutputStream)
	 */
	@Override
	public void sort(OutputStream outputStream) {
        //do nothing
	}

    @Override
    public MSIterator read(String chromosome, int start, int end) {
        if (reader==null)
            reader = new SAMFileReader(this.inputFile);
        if (isApplicable())
            iter = new SAMMappingQueryIterator(inputFile, reader.query(chromosome, start, end, contained), start, end, paired);
        else
            iter = new SAMMappingIterator(chromosome, start, end, reader.iterator());
        return iter;
    }

    /**
     * @see barna.io.MappingReader#getCountReads()
     */
	@Override
	public int getCountReads() {
		return countReads;
	}

	/**
	 * @see barna.io.MappingReader#getCountMappings()
	 */
	@Override
	public int getCountMappings() {
		return countAll;
	}

	/**
	 * @see barna.io.MappingReader#getCountContinuousMappings()
	 */
	@Override
	public int getCountContinuousMappings() {
		return countEntire;
	}

    /**
     * @see barna.io.MappingReader#getCountSplitMappings()
     */
	@Override
	public int getCountSplitMappings() {
		return countSplit;
	}

	@Override
	public boolean isApplicable(UniversalReadDescriptor descriptor) {
		return true;
	}

    @Override
    public boolean close() {
        try {
            if (reader!=null) {
                reader.close();
            }
            return true;
        } catch (Exception e) {
            return false;
        }
    }

    @Override
    public void reset() {
        reader = null;
        iter = null;
    }

    @Override
    public boolean reset(String chr) {
        reader = null;
        iter=null;
        return true;
    }

	@Override
	public void scanFile() {
        if (reader == null)
            reader = new SAMFileReader(this.inputFile);

        countAll = 0; countEntire = 0; countSplit = 0; countReads = 0; countSkippedLines = 0;

        for(SAMRecord rec : reader) {
            countReads++;
            if (rec.getReadUnmappedFlag())
                ++countSkippedLines;
            else {
                ++countAll;
                if (rec.getAlignmentBlocks().size()>1) {
                    if (rec.getCigarString().contains("N"))
                        ++countSplit;
                    else
                        ++countEntire;
                }
                else
                    ++countEntire;
            }
        }

        Log.progressFinish(Constants.OK, true);

        this.close();
        this.reset();
	}

	@Override
	public int getNrInvalidLines() {
		return countSkippedLines;
	}

    @Override
    public MSIterator<Mapping> iterator() {
        return iter;
    }
}
