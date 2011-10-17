

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.Vector;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Future;

import fbi.commons.Execute;
import fbi.genome.io.FileHelper;
import fbi.genome.io.bed.BEDwrapper;
import fbi.genome.io.rna.UniversalReadDescriptor;
import fbi.genome.io.rna.UniversalReadDescriptor.Attributes;
import fbi.genome.io.state.MappingWrapperState;
import fbi.genome.model.bed.BEDobject2;

/**
 * Takes a set of reads mapped to the genome and 
 * finds nucleosomes, outputting some files for
 * visualization.
 * 
 * @author Micha Sammeth (gmicha@gmail.com)
 *
 */
public class NucleosomeFinder implements Callable<Object>{
	
	public static final byte MODE_MONONUCLEOSOME= 1;
	public static final byte MODE_DINUCLEOSOME= 2;
	public static final byte MODE_TRINUCLEOSOME= 3;

	/**
	 * The file containing the mappings.
	 */
	File fileMappings= null;
	
	/**
	 * Coverage of mapped reads separated
	 * by mapping directionality.
	 */
	File fileCoverageReads= null;
	
	/**
	 * Coverage of predicted nucleosome 
	 * centers.
	 */
	File fileCoverageDYLD= null;
	
	/**
	 * File with positioned nucleosomes.
	 */
	File fileNucleosomePos= null;
	
	/**
	 * Size the core nucleosome binds to, 
	 * estimated from distogram.
	 */
	int coreLength= 147;
	
	/**
	 * Size of the outer window ~2* expected 
	 * mono-/di-/tri-nucleosome core length.
	 */
	int windowOuterFlank= 150;	// 300nt window size
	
	/**
	 * Size of the window for kernel smoothing
	 */
	int windowKernelSize= 30;
	
	/**
	 * Half size of the inner window.
	 */
	int signalWindowFlank= 11;	// 5nt or 23nt
	

	/**
	 * Offset from annotated position to <b>read</b>
	 * start, i.e., first genomic position of sense
	 * reads
	 */
	int offsetSense= 0;
	
	/**
	 * Offset from annotated position to <b>read</b>
	 * end, i.e., last genomic position of anti-sense
	 * reads
	 */
	int offsetAntisense= 0;
	
	/**
	 * BED objects sorted by position.
	 */
	Vector<BEDobject2> bedsSortedPosition= null;
	
	/**
	 * BED objects sorted by ID.
	 */
	HashMap<BEDobject2, Integer> bedsSortedID= null;
	
	/**
	 * Positions with DYLD.
	 */
	int[] dyld= null;

	byte mode= MODE_MONONUCLEOSOME;
	UniversalReadDescriptor descriptor= null;
	boolean pairedEnd= false;
	int insertSizeMin= 0;
	int insertSizeMax= 0;
	boolean outPlusReads= false;
	boolean outMinusReads= false;
	boolean outDyld= false;
	boolean outPrediction= false;
	
	private BufferedWriter writerPlus= null;
	private BufferedWriter writerMinus= null;
	private BufferedWriter writerDyld= null;

	
	File mappingFile= null;
	
	public static void main(String[] args) {
		Execute.initialize(2);
		NucleosomeFinder myFinder= new NucleosomeFinder(new File("/Users/micha/projects/demassy/download/IP5300109chrall_sorted.bed"));
		myFinder.insertSizeMax= 200;
		myFinder.insertSizeMax= 400;
		myFinder.mode= MODE_DINUCLEOSOME;
		myFinder.coreLength= 150;
		myFinder.outDyld= true;
		Future<Object> captain= Execute.getExecutor().submit(myFinder);
		try {
			captain.get();
		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (ExecutionException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		Execute.shutdown();
	}
	
	public NucleosomeFinder(File file) {
		this.mappingFile= file;
	}
	
	public Object call() {
		
		BEDwrapper bedReader= new BEDwrapper(mappingFile);
		if (!bedReader.isApplicable()) {
			File f= FileHelper.getSortedFile(mappingFile);
			bedReader.sort(f);
			mappingFile= f;
			bedReader= new BEDwrapper(mappingFile);
			System.err.println("Sorted file at "+ mappingFile.getAbsolutePath());
		}
		
		String chr= "";
		MappingWrapperState state= null;
		
		
		for(;;chr= state.nextChr) {
			state= bedReader.read(chr, 1, Integer.MAX_VALUE);
			if (state.result== null)
				continue;
			
			double[] dyld= nextDYLD((BEDobject2[]) state.result);
			if (outDyld)
				try {
					write(chr, dyld);
				} catch (Exception e) {
					throw new RuntimeException(e);
				}
			
			// stringency
			calcStringency(dyld);
			
			if (state.state== MappingWrapperState.STATE_END_OF_FILE)
				break;
		}
		
		if (outDyld) {
			try {
				getWriterDyld().flush();
				getWriterDyld().close();
			} catch (IOException e) {
				throw new RuntimeException(e);
			}
		}
		
		return null;
	}
	
	private void calcStringency(double[] dyld) {
		
		// kernel
		smoothen(dyld, windowKernelSize);
		
		// stringency
		//stringency(dyld, windowOuterFlank);
		
	}

	private void smoothen(double[] dyld, int w) {
		
		for (int i = w; i < (dyld.length- w); i++) {
			double sum= 0;
			for (int j = -w; j <= w; j++) {
				sum+= triweight(dyld, i-j, w);
			}
		}
		
	}
	
	private double triweight(double[] dyld, int u, int w) {
		return Math.pow(1d- Math.pow(u/(double) w, 2d), 3d);
	}

	void write(String chr, double[] dyld) throws Exception {
		for (int i = 0; i < dyld.length; i++) {
			if (dyld[i]== 0)
				continue;
			getWriterDyld().write(chr+" "+ i+" "+ i+ " "+ Float.toString((float) dyld[i])+ "\n");
		}
		getWriterDyld().flush();
	}

	private BufferedWriter getWriterDyld() {
		if (writerDyld == null) {
			
			try {
				String outF= FileHelper.append(
						mappingFile.getAbsolutePath(), 
						"_dyld", 
						true, 
						"bedGraph");
				System.err.println("writing DYLD graph to "+ outF);
				writerDyld = new BufferedWriter(new FileWriter(outF));
				writerDyld.write("track type=bedGraph name=\"DYLD\" description=\"DYLD positions\" visibility=full color=200,100,0 altColor=0,100,200 priority=20\n");
			} catch (IOException e) {
				throw new RuntimeException(e);
			}
		}

		return writerDyld;
	}
	
	double[] nextDYLD(BEDobject2[] beds) {

		double[] dyld= new double[beds[beds.length- 1].getEnd()];
		
		for (int i = 0; i < beds.length; i++) {
			
			// for pairing
			if (pairedEnd) {
				if (beds[i].getStrand()< 0)
					continue;
				for (int j = i+1; j < beds.length; j++) {
					if (beds[j].getStrand()> 0
							|| beds[j].getEnd()- beds[i].getStart()< insertSizeMin
							|| (pairedEnd&& !isPaired(beds[i], beds[j])))
						continue;
					if (beds[j].getEnd()- beds[i].getStart()> insertSizeMax)
						break;
					markDYLD(dyld, beds[i], beds[j]);
				}
			} else {
				markDYLD(dyld, beds[i], null);
			}
		}
		
		return dyld;
	}

	
	
	private void markDYLD(double[] dyld, BEDobject2 bed1, BEDobject2 bed2) {
		
		// init coordinates
		int start= 0, end= 0;
		if (pairedEnd) {
			
		} else {
			if (bed2== null) {
				int insertSize= (int) Math.round((insertSizeMin+ insertSizeMax)/ 2d);
				if (bed1.getStrand()> 0) {
					start= bed1.getStart()+1 ;
					end= start+ insertSize- 1;
				} else {
					end= bed1.getEnd();
					start= end- insertSize+ 1;
				}
			} else {
				System.err.println("Not implemented 2x unpaired reads.");
				System.exit(-1);
			}
		}
		
		// mark
		switch (mode) {
		case MODE_TRINUCLEOSOME:
			// take middle of dyld
			int p= (int) Math.round((start+ end)/ 2d);
			++dyld[p];
			
		case MODE_DINUCLEOSOME:
			// take 1/2 dyld from downstream end of fragment
			p= (int) Math.round(end- (coreLength/ 2d));
			if (p>= 0&& p< dyld.length)
				++dyld[p];
			
		case MODE_MONONUCLEOSOME:
			// take half from upstream start of fragment
			p= (int) Math.round(start+ (coreLength/ 2d));
			if (p< dyld.length)
				++dyld[p];

		}
		
		
	}

	private boolean isPaired(BEDobject2 bed1, BEDobject2 bed2) {
		
		Attributes at1= descriptor.getAttributes(bed1, null);
		Attributes at2= descriptor.getAttributes(bed2, null);
		
		if (at1.id.equals(at2.id)&& at1.flag!= at2.flag)
			return true;
		return false;
	}

	void oldStuff() {
		
		BEDwrapper bedReader= null;
		
		int pos= 1;
		String chr= null;	// or ""
		MappingWrapperState state= null;
		HashMap<String, long[]> chrMap= bedReader.getMapChr();
		int mapSize= 0;
		
		if (state.state== MappingWrapperState.STATE_END_OF_CHROMOSOME) {
			pos= 2* windowOuterFlank+ 1;
			state= bedReader.read(state.nextChr, 1, pos);
			// remove all
			bedsSortedPosition.removeAllElements();
			bedsSortedID.clear();
			
		} else {
			// remove out
			int posOut= pos- (2* windowOuterFlank+ 1);
			for (int i = 0; i < bedsSortedPosition.size(); i++) {
				BEDobject2 bed= bedsSortedPosition.elementAt(i);
				if (bed.getStart()<= posOut) {
					bedsSortedPosition.remove(i--);
					bedsSortedID.remove(bed);
				}
			}
		}
		
		// add new
		BEDobject2[] beds= (BEDobject2[]) state.result;
		for (int i = 0; i < beds.length; i++) {
			bedsSortedPosition.add(beds[i]);
			Integer edgePos= (beds[i].getStrand()== 1? beds[i].getStart(): beds[i].getEnd());
			bedsSortedID.put(beds[i], edgePos);
		}

	}
	
	
}
