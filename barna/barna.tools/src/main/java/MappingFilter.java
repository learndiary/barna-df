import barna.commons.ByteArrayCharSequence;
import barna.commons.Execute;
import barna.commons.launcher.FluxTool;
import barna.commons.launcher.HelpPrinter;
import barna.commons.log.Log;
import barna.genome.io.BufferedIteratorDisk;
import barna.genome.io.FileHelper;
import barna.genome.io.bed.BEDwrapper;
import barna.genome.io.rna.UniversalReadDescriptor;
import barna.genome.io.rna.UniversalReadDescriptor.Attributes;
import barna.model.bed.BEDobject2;
import org.cyclopsgroup.jcli.ArgumentProcessor;
import org.cyclopsgroup.jcli.annotation.Cli;
import org.cyclopsgroup.jcli.annotation.Option;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.HashSet;
import java.util.Iterator;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Future;

/**
 * Class to filter mappings according to custom criteria.
 * @author Micha Sammeth (gmicha@gmail.com)
 *
 */
@Cli(name = "mfilter", description = "Mapping filter")
public class MappingFilter implements FluxTool<Void> {

	public static void main(String[] args) {
		
		Execute.initialize(2);

		MappingFilter myFilter= new MappingFilter();
		myFilter.setInputFile(new File(
				"/Users/micha/projects/demassy/download_new/B6+K4me3+200511_sorted_chrY_testchrx.bed"
				//"/Users/micha/projects/demassy/download_new/transfer/B6+K4me3+200511.sorted.bed"
		));
		myFilter.setDescriptor(UniversalReadDescriptor.DESCRIPTORID_PAIRED);
		myFilter.setInsertMin(60);
		myFilter.setInsertMax(250);
		
		Future<Void> captain= Execute.getExecutor().submit(myFilter);
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
	 * Default suffix for output files
	 */
	public static String SFX_FILTERED= "_filtered";
	
	/**
	 * Creates the default output file.
	 * @param in the input file
	 * @return output file
	 */
	public static File getDefaultOutputFile(File in) {
		String s= FileHelper.append(in.getAbsolutePath(), SFX_FILTERED, 
				false, null);
		return new File(s);
	}
	
	File fileIn;
	File fileOut;
	UniversalReadDescriptor descriptor;
	int insertMin= 170;
	int insertMax= 230;
	
	public Void call() throws Exception {

		// output
		if (fileOut== null)
			fileOut= getDefaultOutputFile(fileIn);
		Log.info("Output in "+ fileOut.getAbsolutePath());
		
		// sort
		File sortedInput= fileIn;
		sortedInput= BEDwrapper.getSortedFile(fileIn, null, 
				((descriptor!= null&& descriptor.isPaired())?BEDwrapper.COMPARATOR_PAIRED_END:null));

		// doit
		filter(sortedInput, fileOut);
		
		return null;
	}
	
	

	protected void filter(File sortedInput, File fileOut) {
		
		try {
			BufferedIteratorDisk iter= new BufferedIteratorDisk(sortedInput);
			BufferedWriter writer= new BufferedWriter(new FileWriter(fileOut));

			BEDobject2 bed1, bed2, bedLeft, bedRight; 
			Attributes at1= null, at2= null;
			CharSequence bucketID= null;
			HashSet<ByteArrayCharSequence> map= new HashSet<ByteArrayCharSequence>(10);
			
			Log.progressStart("Filtering");
			long fsize= sortedInput.length();
			long ctr= 0l;
			
			while(iter.hasNext()) {
				
				bed1= new BEDobject2(iter.next());
				
				ctr+= bed1.length()+ 1;
				Log.progress(ctr, fsize);

				at1= descriptor.getAttributes(bed1.getName(), at1);
				if (at1.flag!= 1)
					continue;
				
				boolean written1= false;
				iter.mark();
				
				// output bucket
				if (bucketID== null)
					bucketID= at1.id;
				else if (!at1.id.equals(bucketID)) {
					Iterator<ByteArrayCharSequence> iter2= map.iterator();
					while(iter2.hasNext()) {
						writer.write(iter2.next().toCharArray());
						writer.write("\n");
					}
					writer.flush();
					map.clear();
					bucketID= at1.id;
				}
				
				while (iter.hasNext()) {
					bed2= new BEDobject2(iter.next());
					at2= descriptor.getAttributes(bed2.getName(), at2);
					if (!at1.id.equals(at2.id))
						break;
					
					if (bed1.getStart()< bed2.getStart()) {
						bedLeft= bed1;
						bedRight= bed2;
					} else {
						bedLeft= bed2;
						bedRight= bed1;
					}
					if (!(bedLeft.getStrand()> 0&& bedRight.getStrand()< 0))
						continue;	// orientation check
					
					int dist= bedRight.getEnd()- bedLeft.getStart();					
					if (dist> insertMin&& dist< insertMax) {
						if (!written1) {
							writer.write(bed1.toCharArray());
							writer.write("\n");
							written1= true;
						}
						// don't break for + --> - - constellation
						map.add(bed2);
					}
					
				} // end inner hasNext()

				iter.reset();
				
			} // end outer iter.next()

			// last bucket
			Iterator<ByteArrayCharSequence> iter2= map.iterator();
			while(iter2.hasNext()) {
				writer.write(iter2.next().toCharArray());
				writer.write("\n");
			}
			writer.flush();
			map.clear();

			writer.close();
			
			Log.progressFinish("Done.", true);
			
		} catch (Exception e) {
			throw new RuntimeException(e);
		}
		
	}



	public boolean validateParameters(HelpPrinter printer,
			ArgumentProcessor toolArguments) {
		if (fileIn== null|| (!fileIn.exists())|| (!fileIn.canRead())
				|| descriptor== null)
			return false;
		return true;
	}

	
    @Option(name = "i", longName = "input", description = "set input mappings file (BED)", displayName = "input", required = true)
    public void setInputFile(File file) {
        this.fileIn= file;
    }

    @Option(name = "imin", longName = "insertmin", description = "set insert size minimum", displayName = "input size minimum", required = false)
    public void setInsertMin(int imin) {
        this.insertMin= imin;
    }

    @Option(name = "imax", longName = "insertmax", description = "set insert size maximum", displayName = "input size maximum", required = false)
    public void setInsertMax(int imax) {
        this.insertMax= imax;
    }

    @Option(name = "d", longName = "descriptor", description = "set read descriptor", displayName = "descriptor", required = true)
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
