/**
 * 
 */
package barna.io.sam;

import java.io.File;
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

	private Mapping[] mappings= null;
    private SAMFileReader reader = null;

	/**
	 * @param inputFile
	 */
	public SAMReader(File inputFile) {
		super(inputFile);
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
		// TODO Auto-generated method stub
		return false;
	}

	/* (non-Javadoc)
	 * @see barna.io.IOWrapper#sort(java.io.OutputStream)
	 */
	@Override
	public void sort(OutputStream outputStream) {
		// TODO Auto-generated method stub

	}

    @Override
    public MSIterator read(String chromosome, int start, int end) {
        if (reader == null)
            reader = new SAMFileReader(this.inputFile);
        if (reader.hasIndex())
            return new SAMMappingQueryIterator(reader.query(chromosome,start,end, true));
        else
            return new SAMMappingIterator(chromosome, start, end, reader.iterator());
    }

    /* (non-Javadoc)
    * @see barna.io.MappingReader#getCountReads()
    */
	@Override
	public int getCountReads() {
		// TODO Auto-generated method stub
		return 0;
	}

	/* (non-Javadoc)
	 * @see barna.io.MappingReader#getCountMappings()
	 */
	@Override
	public int getCountMappings() {
		// TODO Auto-generated method stub
		return 0;
	}

	/* (non-Javadoc)
	 * @see barna.io.MappingReader#getCountContinuousMappings()
	 */
	@Override
	public int getCountContinuousMappings() {
		// TODO Auto-generated method stub
		return 0;
	}

	/* (non-Javadoc)
	 * @see barna.io.MappingReader#getCountSplitMappings()
	 */
	@Override
	public int getCountSplitMappings() {
		// TODO Auto-generated method stub
		return 0;
	}

	/* (non-Javadoc)
	 * @see barna.io.MappingReader#isApplicable(barna.io.rna.UniversalReadDescriptor)
	 */
	@Override
	public boolean isApplicable(UniversalReadDescriptor descriptor) {
		// TODO Auto-generated method stub
		return false;
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
        return false;  //To change body of implemented methods use File | Settings | File Templates.
    }

    /* (non-Javadoc)
      * @see barna.io.AbstractFileIOWrapper#scanFile()
      */
	@Override
	public void scanFile() {
		// TODO Auto-generated method stub

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
