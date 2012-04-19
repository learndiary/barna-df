package barna.genome.tools.chipseq;


import barna.commons.Execute;
import barna.commons.cli.jsap.JSAPParameters;
import barna.commons.launcher.FluxTool;
import barna.flux.capacitor.reconstruction.Kernel;
import barna.io.FileHelper;
import barna.io.bed.BEDwrapper;
import barna.io.rna.UniversalReadDescriptor;
import barna.io.rna.UniversalReadDescriptor.Attributes;
import barna.io.state.MappingWrapperState;
import barna.model.bed.BEDobject2;
import com.martiansoftware.jsap.JSAPResult;
import com.martiansoftware.jsap.Parameter;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Vector;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Future;

/**
 * Takes a set of reads mapped to the genome and 
 * finds nucleosomes, outputting some files for
 * visualization.
 * 
 * @author Micha Sammeth (gmicha@gmail.com)
 *
 */
public class NucleosomeFinder implements FluxTool<Void> {
	
	public static final byte MODE_MONONUCLEOSOME= 1;
	public static final byte MODE_DINUCLEOSOME= 2;
	public static final byte MODE_TRINUCLEOSOME= 3;

	/**
	 * The file containing the mappings.
	 */
	File fileMappings= null;
	
	/**
	 * The file containing the parameters.
	 */
	File fileParameters= null;
	
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
	int coreLength= 150;
	
	/**
	 * Size of the outer window ~2* expected 
	 * mono-/di-/tri-nucleosome core length.
	 */
	int windowOuterFlank= 150;	// 300nt window size
	
	/**
	 * Size of the window for kernel smoothing
	 */
	int windowKernelFlank= 30;
	
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
	
	byte mode= MODE_MONONUCLEOSOME;
	UniversalReadDescriptor descriptor= null;
	boolean pairSingleReads= false;
	int insertSizeMin= 0;
	int insertSizeMax= 0;

	boolean outOccupancy= true;
	boolean outStringency= true;
	boolean outPeaks= true;
	
	private BufferedWriter writerPeaks= null;
	private BufferedWriter writerStringency= null;
	private BufferedWriter writerOccupancy= null;

	private int ctrReads= 0;
	private int ctrDylds= 0;
	
	
	File mappingFile= null;
	
	public static void main(String[] args) {
		
		Execute.initialize(2);
		
//		NucleosomeFinder myFinder= new NucleosomeFinder(new File(
//			"/Users/micha/projects/demassy/download/IP5300109chrall_sorted_chrY.bed"));
//		myFinder.insertSizeMin= 250;
//		myFinder.insertSizeMax= 400;
//		myFinder.mode= MODE_DINUCLEOSOME;
//		myFinder.coreLength= 150;
		
		NucleosomeFinder myFinder= new NucleosomeFinder(new File(
				//"/Users/micha/projects/demassy/download_new/B6+K4me3+200511_sorted_chrY_testchrx.bed"));
				"/Users/micha/projects/demassy/download_new/transfer/B6+K4me3+200511.sorted.bed_3K"));
		myFinder.descriptor= new UniversalReadDescriptor();
		myFinder.descriptor.init(UniversalReadDescriptor.getDescriptor(UniversalReadDescriptor.DESCRIPTORID_PAIRED));

		Future<Void> captain= Execute.getExecutor().submit(myFinder);
		try {
			captain.get();
		} catch (InterruptedException e) {
			e.printStackTrace();
		} catch (ExecutionException e) {
			e.printStackTrace();
		}
		Execute.shutdown();
	}
	
	/**
	 * Empty constructor to comply with FluxTool implementation.
	 */
	public NucleosomeFinder() {
	}
	
	public NucleosomeFinder(File file) {
		this.fileMappings= file;
	}
	
	public Void call() {
		
		// sort
		File sortedInput= fileMappings;
		sortedInput= BEDwrapper.getSortedFile(fileMappings, null, 
				((descriptor!= null&& descriptor.isPaired())?BEDwrapper.COMPARATOR_PAIRED_END:null));

		// distribution parameters
		ChipSeqMappingAnalyzer distanceDistr= new ChipSeqMappingAnalyzer(
				sortedInput,
				ChipSeqMappingAnalyzer.getFileOutput(fileMappings),
				descriptor);
		Future<int[]> captain= Execute.getExecutor().submit(distanceDistr);
		int[] distr= null;
		try {
			distr= captain.get();
		} catch (Exception e) {
			throw new RuntimeException(e);
		}
		insertSizeMin= distr[0];
		insertSizeMax= distr[2];
		setMode(distr[1]);
		// TODO core length via phaso-/distogram ?!
		System.err.println("peak: "+ distr[1]);
		System.err.println("insert min: "+ insertSizeMin);
		System.err.println("insert max: "+ insertSizeMax);
		System.err.println("mode: "+ mode);
		
		// let's go
		String chr= "";
		double max= -1;
		MappingWrapperState state= null;
		BEDwrapper bedReader= new BEDwrapper(sortedInput);
		bedReader.setMaxBEDObjects(1000);
		
		Vector<BEDobject2> v= new Vector<BEDobject2>();
		Attributes a= descriptor.createAttributes();
		float[] dyld= null;
		ctrReads= 0;
		ctrDylds= 0;

		for(;state== null|| state.state!= MappingWrapperState.STATE_END_OF_FILE; chr= state.nextChr) {
			
			// input
			state= bedReader.read(chr, 1, Integer.MAX_VALUE);
			if (state.result== null)
				continue;

			// concatenate with left-over
			BEDobject2[] beds= (BEDobject2[]) state.result;
			state.result= null;
			if (v.size()> 0) {
				if (descriptor.getAttributes(beds[0].getName(), a).id.equals(
						descriptor.getAttributes(v.elementAt(0).getName(), a))) {
					BEDobject2[] b= new BEDobject2[v.size()+ beds.length];
					for (int i = 0; i < v.size(); i++) 
						b[i]= v.elementAt(i);
					System.arraycopy(beds, 0, b, v.size(), beds.length);
					beds= b;
				} else
					v.removeAllElements();
			}
			
			// get dyld positions
			dyld= nextDYLD(dyld, beds, v.size());
			
			
			// save bucket with last read ID			
			if (chr.equals(state.nextChr)) {
				CharSequence id= descriptor.getAttributes(beds[beds.length- 1].getName(), a).id;
				int x= beds.length- 2;
				while (x> 0&& descriptor.getAttributes(beds[x--].getName(), a).id.equals(id));
				if (x> 0)
					x+= 2;
				else {
					if (v.size()> 0&& 
							!descriptor.getAttributes(v.elementAt(v.size()- 1).getName(), a).id.equals(id))
						v.removeAllElements();
				}
				for (int i = x; i < beds.length; i++) 
					v.add(beds[x]);
				beds= null;
				System.gc();
				
				continue;	// chr not finished yet
			}
			
			// else..
			v.removeAllElements();
			beds= null;
			System.gc();
			System.err.println("\tfound: "+ ctrReads+ " reads, "+ ctrDylds+ " dylds");
			ctrReads= 0;
			ctrDylds= 0;

			if (dyld== null)
				continue;
			
			// calculate stringency
			max= process(chr, dyld, windowKernelFlank, windowOuterFlank, max);
			dyld= null;
		}
		
		try {
			getWriterPeaks().flush();
			getWriterPeaks().close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		if (outStringency) {
			try {
				getWriterStringency().flush();
				getWriterStringency().close();
			} catch (IOException e) {
				throw new RuntimeException(e);
			}
		}
		if (outOccupancy) {
			try {
				getWriterOccupancy().flush();
				getWriterOccupancy().close();
			} catch (IOException e) {
				throw new RuntimeException(e);
			}
		}

		return null;
		
	}
	
	public void setMode(int insertMax) {
		
		byte nrNuc= (byte) Math.round(insertMax/ (float) 150);
		mode= nrNuc;

	}

	/**
	 * Smoothens dyld-positions and calculates stringency measurement per position.
	 * @param chr
	 * @param dyld
	 * @param windowKernelFlank
	 * @param windowOuterFlank
	 */
	protected double process(String chr, float[] dyld, int windowKernelFlank, int windowOuterFlank, double max) {

		// perform kernel-smoothing
		double F= 1.09d;
		double[] triweight= getTriweightKernel(windowKernelFlank, 1d/ F);
		smoothen(dyld, windowKernelFlank, triweight);
		
		// calculate stringency
		double max2= stringency(dyld, windowKernelFlank, windowOuterFlank, F, chr);
		if (max< 0)
			max= max2;
		else if (max2> max)
			System.err.println("max increased: "+ max+ " -> "+ max2);
		
		// predict
		predictDYLD(dyld, chr, max, coreLength);
		
		return max;
	}
	
	protected void predictDYLD(float[] dyld, String chr, double maxScore, int coreLength) {
		
		int w= windowOuterFlank; // coreFlank?
		
		// init
		float max= 0;
		int argmax= -1;
		float[] buf= new float[2* w+ 1];	
		System.arraycopy(dyld, 0, buf, 0, buf.length);
		for (int i = 0; i < buf.length; i++) {
			boolean peakOK= true;
			// check left
			for (int j = i-1; peakOK&& j >= 0; --j) 
				if (buf[j]> buf[i])
					peakOK= false;
				else if (buf[j]< buf[i])
					break;
			// check right
			for (int j = i+1; peakOK&& j< buf.length; ++j)
				if (buf[j]> buf[i])
					peakOK= false;
				else if (buf[j]< buf[i])
					break;
			
			if (peakOK) {
				if (buf[i]> max) {
					max= buf[i];
					argmax= i;
				}
			} else
				buf[i]= 0;
		}
		
		// scan
		int counter= 1;
		if (argmax== w) 
			writeBED(getWriterPeaks(), chr, argmax- w, argmax+ w, 
					"Nuc"+ fill(counter++, 8), (max> maxScore? 1000: 1000* max/maxScore));
		for (int i = w+ 1; i < dyld.length- w; i++) {
			if (dyld[i+w]> dyld[i+w-1])
				buf[buf.length- 1]= 0;	// eliminate former border peak, if one
			System.arraycopy(buf, 1, buf, 0, buf.length- 1);	// decrement
			
			boolean peakOK= true;
			for (int j = 1; peakOK&& j <= w; j++) 
				if (dyld[i+w-j]> dyld[i+w])
					peakOK= false;
				else if (dyld[i+w-j]< dyld[i+w])
					break;
			
			float val= peakOK? dyld[i+w]: 0;
			buf[buf.length- 1]= val;
			if (argmax== 0) {	// was at [0]
				max= 0;
				argmax= -1;
				for (int j = 0; j < buf.length; j++) {
					if (buf[j]> max) {
						max= buf[j];
						argmax= j;
					}
				}
			} else {
				if (val> max) {
					max= val;
					argmax= buf.length- 1;
				} else
					--argmax;
			}
			
			if (argmax== w) 
				writeBED(getWriterPeaks(), chr, i- coreLength/ 2, i+ coreLength/ 2, 
						"Nuc"+ fill(counter++, 8), (max> maxScore? 1000: 1000* max/maxScore));
		}

	}

	void writeBED(BufferedWriter writer, String chr, int start,
			int end, String id, double score) {
		
		try {
			writer.write(chr+"\t"+ start+"\t"+ end+ "\t"+ id+ "\t"+ Float.toString((float) score)+ "\n");
		} catch (IOException e) {
			e.printStackTrace();
		}

		
	}

	private String fill(int counter, int n) {
		StringBuilder s= new StringBuilder(Integer.toString(counter));
		for (int i = s.length(); i <= n; i++) 
			s.insert(0, "0");
		return s.toString();
	}

	protected double stringency(float[] dyld, int windowKernelFlank,
			int windowOuterFlank, double F, String chr) {
		
		double max= 0;
		float[] buf= new float[windowOuterFlank];
		System.arraycopy(dyld, 0, buf, 0, windowOuterFlank);
		for (int i = 0; i < dyld.length; i++) {
			double sum= 0;
			for (int j = Math.max(-windowOuterFlank, -i); j < 0; j++) 
				sum+= buf[Math.min(i+ j, windowOuterFlank+ j)];
			for (int j = i; j <= Math.min(i+ windowOuterFlank, dyld.length- 1); j++) 
				sum+= dyld[j];
			
			// correct window size of positions really counted
			int corr= 0;	
			if (i- windowKernelFlank< 0)
				corr+= windowKernelFlank- i;
			if (i+ windowKernelFlank>= dyld.length)
				corr+= i+ windowKernelFlank- dyld.length+ 1;
			
			double s= (dyld[i]* (windowKernelFlank- corr))/ (sum* F);
			if (s< 0|| s> 1.1)
				System.err.println("stringency out of range: "+ s);
			if (outOccupancy)
				writeBedGraph(getWriterOccupancy(), chr, i, dyld[i]);

			if (outStringency)
				writeBedGraph(getWriterStringency(), chr, i, s);
			
			// next
			if (i> windowOuterFlank) {
				System.arraycopy(buf, 1, buf, 0, buf.length- 1);
				buf[buf.length- 1]= dyld[i];
			}
			dyld[i]*= s;
			max= (dyld[i]> max? dyld[i]: max);
		}
		if (outOccupancy)
			try {
				getWriterOccupancy().flush();
			} catch (IOException e) {
				e.printStackTrace();
			}
		if (outStringency)
			try {
				getWriterStringency().flush();
			} catch (IOException e) {
				e.printStackTrace();
			}
			
		return max;
	}

	protected void smoothen(float[] dyld, int w, double[] triweight) {
		
		float[] buf= new float[w];
		System.arraycopy(dyld, 0, buf, 0, w);
		// TODO assume 0-values beyond chr edges to smoothen
		for (int i = 0; i < dyld.length; i++) {
			double sum= 0;
			for (int j = -w; j < 0; j++)
				sum+= (i< w? (i+ j< 0? 0: triweight[j+ w]* buf[i+ j]): triweight[j+ w]* buf[w+ j]);
			for (int j = 0; j <= w; j++) 
				sum+= (i+ j>= dyld.length? 0: triweight[j+ w]* dyld[i+ j]);
			
			// next
			if (i>= w) {
				System.arraycopy(buf, 1, buf, 0, buf.length- 1);	// omit leftmost
				buf[buf.length- 1]= dyld[i];
			}
			dyld[i]= (float) sum;
		}
		
	}
	
	private double[] triweightKernel= null; 
	private double[] getTriweightKernel(int w, double F) {
		if (triweightKernel == null) {
			triweightKernel = Kernel.getTriweightKernel(w, F);
		}

		return triweightKernel;
	}

	void writeBedGraph(BufferedWriter writer, String chr, double[] dyld) throws Exception {
		try {
			for (int i = 0; i < dyld.length; i++) {
				writeBedGraph(writer, chr, i, dyld[i]);
			}
			writer.flush();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	void writeBedGraph(BufferedWriter writer, String chr, int pos, double val) {
		if (val!= 0)
			try {
				writer.write(chr+" "+ pos+" "+ (pos+1)+ " "+ Float.toString((float) val)+ "\n");
			} catch (Exception e) {
				e.printStackTrace();
			}

	}

	private BufferedWriter getWriterStringency() {
		if (writerStringency == null) {
			
			try {
				String outF= FileHelper.append(
						fileMappings.getAbsolutePath(), 
						"_stringency", 
						true, 
						".bedGraph");
				System.err.println("writing stringency graph to "+ outF);
				writerStringency = new BufferedWriter(new FileWriter(outF));
				writerStringency.write("track type=bedGraph name=\""
						+ FileHelper.stripExtension(fileMappings.getName())
						+ " stringency\" description=\"nucleosome stringency\" visibility=full color=200,100,0 altColor=0,100,200 priority=20\n");
			} catch (IOException e) {
				throw new RuntimeException(e);
			}
		}

		return writerStringency;
	}
	
	private BufferedWriter getWriterOccupancy() {
		if (writerOccupancy == null) {
			try {
				String outF= FileHelper.append(
						fileMappings.getAbsolutePath(), 
						"_occupancy", 
						true, 
						".bedGraph");
				System.err.println("writing occupancy graph to "+ outF);
				writerOccupancy = new BufferedWriter(new FileWriter(outF));
				writerOccupancy.write("track type=bedGraph name=\""
						+ FileHelper.stripExtension(fileMappings.getName())
						+ " occupancy\" description=\"nucleosome occupancy\" visibility=full color=200,100,0 altColor=0,100,200 priority=20\n");
			} catch (IOException e) {
				throw new RuntimeException(e);
			}
			
		}

		return writerOccupancy;
	}
	
	private BufferedWriter getWriterPeaks() {
		if (writerPeaks == null) {
			try {
				String outF= FileHelper.append(
						fileMappings.getAbsolutePath(), 
						"_peaks", 
						true, 
						".bed");
				System.err.println("writing peak predictions to "+ outF);
				writerPeaks = new BufferedWriter(new FileWriter(outF));
				writerPeaks.write("track name=\""
						+ FileHelper.stripExtension(fileMappings.getName())
						+ " peaks\" description=\"peak predictions\" useScore=1\n");
			} catch (IOException e) {
				throw new RuntimeException(e);
			}
			
		}

		return writerPeaks;
	}
	
	float[] nextDYLD(float[]  dyld, BEDobject2[] beds, int x) {

		// get max
		int max= -1;
		for (int i = 0; i < beds.length; i++) {
			int p= beds[i].getEnd();
			if (p> max)
				max= p;
		}
		if (max< 0)
			return null;
		
		// init
		if (dyld== null)
			dyld= new float[max];
		else if (max> dyld.length) {
			float[] d= new float[max];
			System.arraycopy(dyld, 0, d, 0, dyld.length);
			dyld= d;
			System.gc();
		}
		
		// mark dyld centers
		boolean pairedEnd= descriptor!= null&& descriptor.isPaired();
		BEDobject2 bed1, bed2;
		Attributes at1= null, at2= null;
		for (int i = 0; i < beds.length; i++) {
			
			++ctrReads;
			
			// for pairing
			if (pairedEnd) {				
				at1= descriptor.getAttributes(beds[i].getName(), at1);
				if (at1.flag!= 1)
					continue;
				for (int j = i+1; j < beds.length; j++) {
					if (i< x&& j< x)
						continue;	// skip left-over
					at2= descriptor.getAttributes(beds[j].getName(), at2);
					if (!at1.id.equals(at2.id))
						break;
					if (at2.flag!= 2)
						continue;
					if (beds[i].getStart()< beds[j].getStart()) {
						bed1= beds[i];
						bed2= beds[j];
					} else {
						bed2= beds[i];
						bed1= beds[j];
					}
					if (!(bed1.getStrand()> 0&& bed2.getStrand()< 0))
						continue;	// orientation check
					int dist= bed2.getEnd()- bed1.getStart();
					if (dist> insertSizeMin&& dist< insertSizeMax)
						markDYLD(dyld, bed1, bed2);
				}
			} else {
				if (pairSingleReads) {
					if (beds[i].getStrand()< 0)
						continue;	// start pairs only in correct orientation
					bed1= beds[i];					
					for (int j = i+1; j < beds.length; j++) {
						if (i< x&& j< x)
							continue;	// skip left-over
						if (beds[j].getStrand()> 0)
							continue;
						bed2= beds[j];
						int dist= bed2.getEnd()- bed1.getStart();
						if (dist> insertSizeMax)
							break;
						markDYLD(dyld, bed1, bed2);
					}

				} else {
					if (i< x)
						continue;	// skip left-over
					markDYLD(dyld, beds[i], null);
				}
			}
		}
		
		return dyld;
	}

	
	
	void markDYLD(float[] dyld, BEDobject2 bed1, BEDobject2 bed2) {
		
		++ctrDylds;
		
		// init coordinates
		int start= 0, end= 0;
		if (bed2== null) {
			int insertSize= (int) Math.round((insertSizeMin+ insertSizeMax)/ 2d);
			if (bed1.getStrand()> 0) {
				start= bed1.getStart()+ 1;
				end= start+ insertSize- 1;
			} else {
				end= bed1.getEnd();
				start= end- insertSize+ 1;
			}
		} else {	// paired reads
			start= bed1.getStart();
			end= bed2.getEnd();
		}
		
		// mark
		int val= (bed2== null? 1: 2);
		switch (mode) {
			case MODE_TRINUCLEOSOME:
				// do both
				
			case MODE_DINUCLEOSOME:
				// take 1/2 dyld from downstream end of fragment
				int p= (int) Math.round(end- (coreLength/ 2d));
				if (p>= 0&& p< dyld.length)
					dyld[p]+= val;
				p= (int) Math.round(start+ (coreLength/ 2d));
				if (p>= 0&& p< dyld.length)
					dyld[p]+= val;
				if (mode== MODE_DINUCLEOSOME)
					break;
				
			case MODE_MONONUCLEOSOME:
				// take half from upstream start of fragment
				p= (int) Math.round(start+ ((end- start)/ 2d));	// mid
				if (p>= 0&& p< dyld.length)
					dyld[p]+= val;
		}
		
		
	}

	boolean isPaired(BEDobject2 bed1, BEDobject2 bed2) {
		
		Attributes at1= descriptor.getAttributes(bed1.getName(), null);
		Attributes at2= descriptor.getAttributes(bed2.getName(), null);
		
		if (at1.id.equals(at2.id)&& at1.flag!= at2.flag)
			return true;
		return false;
	}

    //@Cli(name = "nfind", description = "Nucleosome finder")
    @Override
    public String getName() {
        return "nfind";
    }

    @Override
    public String getDescription() {
        return "Nucleosome finder";
    }

    @Override
    public List<Parameter> getParameter() {
        ArrayList<Parameter> parameters = new ArrayList<Parameter>();
        parameters.add(JSAPParameters.flaggedParameter("input", 'i').type(File.class).help("set input mappings file (BED)").valueName("bed").required().get());
        parameters.add(JSAPParameters.flaggedParameter("descriptor", 'd').help("read descriptor").valueName("bed").required().get());
        return parameters;

    }

    @Override
    public boolean validateParameter(JSAPResult args) {
        setDescriptor(args.getString("descriptor"));
        setParameters(args.getFile("input"));
        if (fileMappings== null|| descriptor== null)
			return false;
		return true;
	}
	

    public void setParameters(File file) {
        this.fileMappings= file;
    }


    public void setDescriptor(String descriptor) {
    	this.descriptor= new UniversalReadDescriptor();
    	try {
    		this.descriptor.init(UniversalReadDescriptor.getDescriptor(descriptor));
    	} catch (Exception e) {
    		e.printStackTrace();
    		this.descriptor= null;
    	}
    }
	
	
}
