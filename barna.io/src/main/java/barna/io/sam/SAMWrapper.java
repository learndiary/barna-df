/**
 * 
 */
package barna.io.sam;

import java.io.File;
import java.io.OutputStream;
import java.util.Iterator;
import java.util.Vector;

import barna.commons.utils.ArrayUtils;
import barna.io.AbstractFileIOWrapper;
import barna.io.MSIterator;
import barna.io.MappingReader;
import barna.io.rna.UniversalReadDescriptor;
import barna.model.Gene;
import barna.model.Mapping;
import barna.model.bed.BEDobject;
import net.sf.samtools.*;
/**
 * @author emilio
 *
 */
public class SAMWrapper extends AbstractFileIOWrapper implements
        MappingReader {

	public BEDobject[] beds= null;
	/**
	 * @param inputFile
	 */
	public SAMWrapper(File inputFile) {
		super(inputFile);
		// TODO Auto-generated constructor stub
	}
	
	/**
	 * Creates an instance using a specific path to a file 
	 * and the default line comparator.
	 * @param absolutePath path to the file the wrapper is based on
	 */
	public SAMWrapper(String absolutePath) {
		this(new File(absolutePath));
	}

	/* (non-Javadoc)
	 * @see barna.io.IOWrapper#read()
	 */
	@Override
	public void read() {		
		Vector objV= new Vector();
		
		try {	
			final SAMFileReader inputSam = new SAMFileReader(inputFile);
			final SAMFileHeader inputSamHeader =  inputSam.getFileHeader();
			
			for (final SAMRecord rec : inputSam) 
			{
				if (!rec.getReadUnmappedFlag()) {
					BEDobject bed= BEDobject.getRecycleObj(); // new BEDobject();
					
					bed.setChrom(rec.getReferenceName());				
					bed.setStart(rec.getAlignmentStart()-1);
					bed.setEnd(rec.getAlignmentEnd()-1);
					bed.setName(rec.getReadName());
					bed.setScore(rec.getMappingQuality());			
					
					if (rec.getReadNegativeStrandFlag())				
						bed.setStrand("-");
					else
						bed.setStrand("+");							
					
					/*
						bed.setThickStart();					
						bed.setThickEnd();											
						bed.setCol();								
						bed.setBlockCount());					
						bed.setBlockSizes();					 
						bed.setBlockStarts();
					*/
			
					//objV.add(bed);	
				}
	        } 
			
			
	} catch (Exception e) {
		throw new RuntimeException(e);
	}


	beds= (BEDobject[]) ArrayUtils.toField(objV);
        
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
    public MSIterator<Mapping> read(String chromosome, int start, int end) {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
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
