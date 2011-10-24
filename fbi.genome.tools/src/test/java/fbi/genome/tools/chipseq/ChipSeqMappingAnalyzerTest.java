package fbi.genome.tools.chipseq;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Future;

import org.junit.Test;

import fbi.commons.Execute;
import fbi.commons.flux.Flux;
import fbi.genome.io.FileHelper;
import fbi.genome.io.rna.UniversalReadDescriptor;

public class ChipSeqMappingAnalyzerTest {

    public static void setUp(){
		Execute.initialize(2);
    }
	
    public static void shutdown(){
		Execute.shutdown();
    }
	
	@Test
	public void testPaired() {

		setUp();
		
		// copy input file
		File f= new File(ChipSeqMappingAnalyzerTest.class.getResource("/Paired_sorted_chrY.bed").getFile());
		
		// instantiate and run
		ChipSeqMappingAnalyzer myRun= new ChipSeqMappingAnalyzer(f);
		myRun.descriptor= new UniversalReadDescriptor();
		myRun.descriptor.init(UniversalReadDescriptor.getDescriptor(UniversalReadDescriptor.DESCRIPTORID_PAIRED));
				
		Future<int[]> captain= Execute.getExecutor().submit(myRun);
		int[] distr= null;
		try {
			distr= captain.get();
		} catch (InterruptedException e) {
			e.printStackTrace();
		} catch (ExecutionException e) {
			e.printStackTrace();
		}

		if (distr!= null)
			System.err.println("peak "+ distr[1]+ ", bounds=["+ distr[0]+ ","+ distr[2]+ "]");

		shutdown();
	}

	@Test
		public void testPairedFlux() {
			
			// copy input file
			File f= new File(ChipSeqMappingAnalyzerTest.class.getResource("/Paired_sorted_chrY.bed").getFile());
			File tmpF= null;
			try {
				tmpF= FileHelper.createTempFile(FileHelper.stripExtension(f.getName()), FileHelper.getExtension(f));
			} catch (IOException e) {
				e.printStackTrace();
			}
			FileHelper.copy(f, tmpF);
			
			
			// write par
			File parF= null;
			BufferedWriter writer= null;
			try {
				parF= FileHelper.createTempFile(ChipSeqMappingAnalyzerTest.class.getSimpleName(), ".par");
				writer= new BufferedWriter(new FileWriter(parF));
				writer.write(ChipSeqSettings.FILE_INPUT.getName()+ " "+ tmpF.getAbsolutePath()+ "\n");
				writer.write(ChipSeqSettings.READ_DESCRIPTOR.getName()+ " "+ UniversalReadDescriptor.DESCRIPTORID_PAIRED+ "\n");
			} catch (Exception e) {
				e.printStackTrace();
			} finally {
				if (writer!= null)
					try {
						writer.close();
					} catch (IOException e) {
						e.printStackTrace();
					}
			}
			
			
			// instantiate and run
			String[] args= {"-t", "mapdist", "-p", parF.getAbsolutePath()};
			Flux.main(args);
	
		}

}
