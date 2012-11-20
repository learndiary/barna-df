/**
 * 
 */
package barna.io.sam;

import barna.commons.io.DevNullOutputStream;
import barna.commons.log.Log;
import barna.commons.system.OSChecker;
import barna.commons.utils.Interceptable;
import barna.io.AbstractFileIOWrapper;
import barna.io.MSIterator;
import barna.io.MappingReader;
import barna.io.Sorter;
import barna.io.rna.UniversalReadDescriptor;
import barna.model.Mapping;
import barna.model.constants.Constants;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;

import java.io.*;
import java.util.concurrent.Future;

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
    private boolean sortInRam;
    private int maxRecords = 500000;

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
		this(inputFile, CONTAINED_DEFAULT, descriptor, false);
	}

    /**
     * Creates an instance using a specific path to a file.	 .
     * @param absolutePath path to the file the wrapper is based on
     * @param descriptor the descriptor to be used
     */
    public SAMReader(String absolutePath, UniversalReadDescriptor descriptor) {
        this(new File(absolutePath), CONTAINED_DEFAULT, descriptor, false);
    }

    /**
     * Creates an instance using a specific path to a file.	 .
     * @param absolutePath path to the file the wrapper is based on
     * @param contained flag to decide whether return mappings which are
     *                  contained/overlapping the query region
     * @param descriptor the descriptor to be used
     */
    public SAMReader(String absolutePath, boolean contained, UniversalReadDescriptor descriptor) {
        this(new File(absolutePath), contained, descriptor, false);
    }

    /**
     * Creates an instance of the reader
     * @param inputFile the file to read
     * @param contained flag to decide whether return mappings which are
     *                  contained/overlapping the query region
     * @param descriptor the descriptor to be used
     */
    public SAMReader(File inputFile, boolean contained, UniversalReadDescriptor descriptor, boolean sortInRam) {
        super(inputFile);
        reader = new SAMFileReader(this.inputFile);
        this.contained = contained;
        this.descriptor = descriptor;
        this.sortInRam = sortInRam;
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
            throw new RuntimeException("The input must be a BAM file.");
        if (!reader.hasIndex())
            throw new RuntimeException("The input BAM file must be sorted and indexed.");
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
            reader.enableIndexCaching(true);
            reader.enableIndexMemoryMapping(false);
        }
        if (isApplicable()) {
//            iter = new SAMMappingQueryIterator(inputFile, reader.query(chromosome, start, end, contained), start, end, paired);
            if (sortInRam) {
                try {
                    iter = new SAMMappingIterator(reader.query(chromosome, start, end, contained), descriptor);
                }
                catch (OutOfMemoryError error) {
                    SAMFileHeader header =  reader.getFileHeader();
                    header.setSortOrder(SAMFileHeader.SortOrder.queryname);
                    iter = new SAMMappingSortedIterator(reader.query(chromosome, start, end, contained), descriptor, header, maxRecords);
                }
            } else {
                SAMFileHeader header =  reader.getFileHeader();
                header.setSortOrder(SAMFileHeader.SortOrder.queryname);
                iter = new SAMMappingSortedIterator(reader.query(chromosome, start, end, contained), descriptor, header, maxRecords);
            }
        }
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
		if (descriptor.equals(UniversalReadDescriptor.getDefaultDescriptor()))
            return true;
        else
            throw new RuntimeException("The read descriptor cannot be specified when using BAM input files.");
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
        if (iter != null) {
            iter.clear();
            iter = null;
        }
        if (reader != null) {
            reader.close();
            reader = null;
        }
    }

    @Override
    public boolean reset(String chr) {
        reader = null;
        iter=null;
        return true;
    }

	@Override
	public void scanFile() {
        if (reader == null) {
            reader = new SAMFileReader(this.inputFile);
        }
        //reader.getFileHeader().setSortOrder(SAMFileHeader.SortOrder.queryname);        countAll = 0; countEntire = 0; countSplit = 0; countReads = 0; countSkippedLines = 0;
        String lastReadId = null;

        BufferedReader buffy = null;
        BufferedWriter tmpWriter = null;
        PipedInputStream in = null;
        PipedOutputStream out = null;
        Future sorterFuture = null;

        int primaryMappings = 0;

        try {
            out = new PipedOutputStream();
            in = new PipedInputStream(out);
            tmpWriter= new BufferedWriter(new OutputStreamWriter(out));

            sorterFuture = Sorter.create(in, new DevNullOutputStream(), true, "\t")
                    .field(0, false)
                    .addInterceptor(new Interceptable.Interceptor<String>() {
                        String lastLine = null;

                        public String intercept(String line) {
                            if (lastLine == null || !line.equals(lastLine)) {
                                ++countReads;
                            }
                            lastLine = line;
                            return line;
                        }
                    })
                    .sortInBackground();

            for(final SAMRecord rec : reader) {
                if (rec.getReadUnmappedFlag()) {
                    ++countSkippedLines;
                } else {
                    String readId = rec.getReadName();
                    if (rec.getReadPairedFlag()) {
                        readId += "/"+(rec.getFirstOfPairFlag()?1:2);
                    }
                    if (!rec.getNotPrimaryAlignmentFlag()) {
                        ++primaryMappings;
                        tmpWriter.write(readId);
                        tmpWriter.write(OSChecker.NEW_LINE);
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

            tmpWriter.flush();
            tmpWriter.close();
            sorterFuture.get();

            Log.progressFinish(Constants.OK, true);

            this.close();
            this.reset();

        } catch (Exception e) {
            System.err.println("[Sorter THREAD]");
            e.printStackTrace();
        } finally {
            if(buffy != null)try {buffy.close();} catch (IOException e) {}
            if(tmpWriter != null)try {tmpWriter.flush();tmpWriter.close();} catch (IOException e) {}
            if(in != null)try {in.close();} catch (IOException e) {}
            if(out != null)try {out.close();} catch (IOException e) {}
            if(sorterFuture != null)sorterFuture.cancel(true);
        }

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
