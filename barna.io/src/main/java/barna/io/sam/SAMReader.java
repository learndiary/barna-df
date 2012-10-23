/**
 * 
 */
package barna.io.sam;

import barna.commons.log.Log;
import barna.io.AbstractFileIOWrapper;
import barna.io.MSIterator;
import barna.io.MappingReader;
import barna.io.rna.UniversalReadDescriptor;
import barna.model.Mapping;
import barna.model.constants.Constants;
import net.sf.samtools.*;

import java.io.*;

/**
 * @author Emilio Palumbo (emiliopalumbo@gmail.com)
 */
public class SAMReader extends AbstractFileIOWrapper implements
        MappingReader {

    public static final boolean CONTAINED_DEFAULT = false;

    private SAMFileReader reader;
    private final UniversalReadDescriptor descriptor;
    private boolean contained;
    private MSIterator iter;

    int countAll;
    int countEntire;
    int countSplit;
    int countReads;
    int countSkippedLines;

	/**
     * Creates an instance of the reader
	 * @param inputFile the file to read
     * @param descriptor the descriptor to be used
	 */
	public SAMReader(File inputFile, UniversalReadDescriptor descriptor) {
		this(inputFile, CONTAINED_DEFAULT, descriptor);
	}

    /**
     * Creates an instance using a specific path to a file.	 .
     * @param absolutePath path to the file the wrapper is based on
     * @param descriptor the descriptor to be used
     */
    public SAMReader(String absolutePath, UniversalReadDescriptor descriptor) {
        this(new File(absolutePath), CONTAINED_DEFAULT, descriptor);
    }

    /**
     * Creates an instance using a specific path to a file.	 .
     * @param absolutePath path to the file the wrapper is based on
     * @param contained flag to decide whether return mappings which are
     *                  contained/overlapping the query region
     * @param descriptor the descriptor to be used
     */
    public SAMReader(String absolutePath, boolean contained, UniversalReadDescriptor descriptor) {
        this(new File(absolutePath), contained, descriptor);
    }

    /**
     * Creates an instance of the reader
     * @param inputFile the file to read
     * @param contained flag to decide whether return mappings which are
     *                  contained/overlapping the query region
     * @param descriptor the descriptor to be used
     */
    public SAMReader(File inputFile, boolean contained, UniversalReadDescriptor descriptor) {
        super(inputFile);
        reader = new SAMFileReader(this.inputFile);
        this.contained = contained;
        this.descriptor = descriptor;
    }

    @Override
    public void read() {
    }

	@Override
	public void write() {
	}

	@Override
	public boolean isApplicable() {
		if (reader==null)
            reader = new SAMFileReader(this.inputFile);
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
        //do nothing for the moment
	}

    @Override
    public MSIterator read(String chromosome, int start, int end) {
        if (reader==null) {
            reader = new SAMFileReader(this.inputFile);
        }
        if (isApplicable())
//            iter = new SAMMappingQueryIterator(inputFile, reader.query(chromosome, start, end, contained), start, end, paired);
            iter = new SAMMappingIterator(reader.query(chromosome, start, end, contained), descriptor);
        else
            throw new UnsupportedOperationException("Only indexed BAM files are currently supported!");
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
		return isApplicable();
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
        reader.getFileHeader().setSortOrder(SAMFileHeader.SortOrder.queryname);
        countAll = 0; countEntire = 0; countSplit = 0; countReads = 0; countSkippedLines = 0;
        String lastReadId = null;

        PipedInputStream pip = new PipedInputStream();

        try {
            final PipedOutputStream pop = new PipedOutputStream(pip);

            new Thread(
                new Runnable(){
                    public void run(){
                        SAMFileWriter writer = new SAMFileWriterFactory().makeSAMWriter(reader.getFileHeader(),false, pop);
                        for(final SAMRecord rec : reader) {
                            writer.addAlignment(rec);
                            try {
                                pop.flush();
                            } catch (IOException e) {
                            }
                        }
                        writer.close();
                    }
                }
            ).start();
        } catch (IOException e) {
            System.err.println("[THREAD]");
            e.printStackTrace();
        }

        SAMFileReader r = new SAMFileReader(pip);

        for(final SAMRecord rec : r) {
            if (rec.getReadUnmappedFlag()) {
                ++countSkippedLines;
            } else {
                String readId = rec.getReadName();
                if (rec.getReadPairedFlag()) {
                    readId += "/"+(rec.getFirstOfPairFlag()?1:2);
                }
                if (!readId.equals(lastReadId)) {
                    ++countReads;
                    lastReadId=readId;
                }
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
