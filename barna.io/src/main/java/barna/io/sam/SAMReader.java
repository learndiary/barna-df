/**
 * 
 */
package barna.io.sam;

import java.io.File;
import java.io.FileOutputStream;
import java.io.OutputStream;
import java.util.Iterator;

import barna.io.AbstractFileIOWrapper;
import barna.io.MSIterator;
import barna.io.MappingReader;
import barna.io.rna.UniversalReadDescriptor;
import barna.model.Mapping;
import net.sf.samtools.*;

/**
 * @author Emilio Palumbo (emiliopalumbo@gmail.com)
 */
public class SAMReader extends AbstractFileIOWrapper implements
        MappingReader {

    public static final boolean CONTAINED_DEFAULT = false;

	private Mapping[] mappings= null;
    private SAMFileReader reader = null;
    private boolean contained = CONTAINED_DEFAULT;

    int countAll;
    int countEntire;
    int countSplit;
    int countReads;

	/**
	 * @param inputFile
	 */
	public SAMReader(File inputFile) {
		super(inputFile);
	}

    public SAMReader(File inputFile, boolean contained) {
        super(inputFile);
        this.contained=contained;
    }
	
	/**
	 * Creates an instance using a specific path to a file 
	 * and the default line comparator.
	 * @param absolutePath path to the file the wrapper is based on
	 */
	public SAMReader(String absolutePath) {
		this(new File(absolutePath));
	}

    @Override
    public void read() {
    }

    /* (non-Javadoc)
      * @see barna.io.IOWrapper#write()
      */
	@Override
	public void write() {
		// TODO Auto-generated method stub
	}

	/* (non-Javadoc)
	 * @see barna.io.IOWrapper#isApplicable()
	 */
	@Override
	public boolean isApplicable() {
		return true;
	}

	/* (non-Javadoc)
	 * @see barna.io.IOWrapper#sort(java.io.OutputStream)
	 */
	@Override
	public void sort(OutputStream outputStream) {
//        reader.getFileHeader().setSortOrder(SAMFileHeader.SortOrder.coordinate);
//        final SAMFileWriter writer = new SAMFileWriterFactory().makeSAMOrBAMWriter(reader.getFileHeader(), false, null);
//
//        final Iterator<SAMRecord> iterator = reader.iterator();
//        while (iterator.hasNext()) {
//            writer.addAlignment(iterator.next());
//        }
//
//        reader.close();
//        writer.close();
	}

    @Override
    public MSIterator read(String chromosome, int start, int end) {
        if (reader == null)
            reader = new SAMFileReader(this.inputFile);
        if (reader.hasIndex())
            return new SAMMappingQueryIterator(inputFile, reader.query(chromosome, start, end, contained), start, end);
        else
            return new SAMMappingIterator(chromosome, start, end, reader.iterator());
    }

    /* (non-Javadoc)
    * @see barna.io.MappingReader#getCountReads()
    */
	@Override
	public int getCountReads() {
		return countReads;
	}

	/* (non-Javadoc)
	 * @see barna.io.MappingReader#getCountMappings()
	 */
	@Override
	public int getCountMappings() {
		return countAll;
	}

	/* (non-Javadoc)
	 * @see barna.io.MappingReader#getCountContinuousMappings()
	 */
	@Override
	public int getCountContinuousMappings() {
		return countEntire;
	}

	/* (non-Javadoc)
	 * @see barna.io.MappingReader#getCountSplitMappings()
	 */
	@Override
	public int getCountSplitMappings() {
		return countSplit;
	}

	/* (non-Javadoc)
	 * @see barna.io.MappingReader#isApplicable(barna.io.rna.UniversalReadDescriptor)
	 */
	@Override
	public boolean isApplicable(UniversalReadDescriptor descriptor) {
		return true;
	}

    @Override
    public boolean close() {
        return false;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public void reset() {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public boolean reset(String chr) {
        return true;  //To change body of implemented methods use File | Settings | File Templates.
    }

    /* (non-Javadoc)
      * @see barna.io.AbstractFileIOWrapper#scanFile()
      */
	@Override
	public void scanFile() {
        if (reader == null)
            reader = new SAMFileReader(this.inputFile);

        countAll= 0; countEntire= 0; countSplit= 0; countReads= 0;

        for(SAMRecord rec : reader) {
            countReads++;
            if (!rec.getReadUnmappedFlag()) {
                countAll++;
                if (rec.getAlignmentBlocks().size()>1)
                    countSplit++;
                else
                    countEntire++;
            }
        }
        reader.close();
        reader=null;
	}

	/* (non-Javadoc)
	 * @see barna.io.AbstractFileIOWrapper#getNrInvalidLines()
	 */
	@Override
	public int getNrInvalidLines() {
		// TODO Auto-generated method stub
		return 0;
	}

    @Override
    public Iterator<Mapping> iterator() {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }
}
