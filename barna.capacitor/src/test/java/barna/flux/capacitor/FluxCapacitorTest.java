package barna.flux.capacitor;

import barna.commons.Execute;
import barna.flux.capacitor.reconstruction.FluxCapacitor;
import barna.flux.capacitor.reconstruction.FluxCapacitorSettings;
import barna.flux.capacitor.reconstruction.FluxCapacitorSettings.AnnotationMapping;
import barna.io.FileHelper;
import barna.io.Sorter;
import barna.io.rna.UniversalReadDescriptor;
import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.*;
import java.util.concurrent.Future;
import java.util.zip.GZIPOutputStream;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

import static junit.framework.Assert.*;

public class FluxCapacitorTest {

	static final int SORTED= -1, UNSORT_GTF= 8, UNSORT_BED= 10; 
	final File GTF_SORTED= new File(getClass().getResource("/mm9_chr1_chrX.gtf").getFile());
	final File BED_SORTED= new File(getClass().getResource("/chr1_chrX.bed").getFile());
	final String subdirMappings= "mappings";
	final String subdirAnnotation= "annotation";
	final String suffixOutput= "gtf";
	final String suffixParameter= "par";
	
	protected File tmpDir= null;
	protected File anoDir= null;
	protected File mapDir= null;
	protected File outFile= null;
	protected File parFile= null;
	protected File gtfFile= null;
	protected File bedFile= null;
	
	
	
	private void initFileNames(byte compressionGTF, byte compressionBED) throws Exception {
		// set up file structure
		tmpDir= new File(System.getProperty("java.io.tmpdir"));
		anoDir= FileHelper.createTempDir(getClass().getSimpleName(), subdirAnnotation, tmpDir);
		anoDir.deleteOnExit();
		mapDir= FileHelper.createTempDir(getClass().getSimpleName(), subdirMappings, tmpDir);
		mapDir.deleteOnExit();
		outFile= File.createTempFile(getClass().getSimpleName(), suffixOutput, anoDir);
		outFile.delete();
		parFile= File.createTempFile(getClass().getSimpleName(), suffixParameter, anoDir);
		gtfFile= new File(FileHelper.append(anoDir.getAbsolutePath()+ File.separator+ 
				GTF_SORTED.getName(), null, false, FileHelper.getCompressionExtension(compressionGTF)));
		gtfFile.delete();		
		bedFile= new File(FileHelper.append(mapDir.getAbsolutePath()+ 
				File.separator+ 
				BED_SORTED.getName(), null, false, FileHelper.getCompressionExtension(compressionBED)));
		bedFile.delete();
	}
	
	protected void copy(File source, File target, int sortField, byte compression) throws Exception {
		if (sortField< 0) {
			if (compression== FileHelper.COMPRESSION_NONE)
				FileHelper.copy(source, target);
			else
				FileHelper.deflate(source, target, compression);
		} else {
			
			// init input stream
			PipedInputStream pin= new PipedInputStream();
			PipedOutputStream pout= new PipedOutputStream(pin);
			
			// init output stream
			OutputStream ostream= new FileOutputStream(target);
			if (compression== FileHelper.COMPRESSION_GZIP)
				ostream= new GZIPOutputStream(ostream);
			else if (compression== FileHelper.COMPRESSION_ZIP) {
				ZipOutputStream zstream= new ZipOutputStream(ostream);
				zstream.putNextEntry(new ZipEntry(source.getName()));
				ostream= zstream;
			}
			
			// init sorter
			Sorter s= Sorter.create(pin, ostream, true, "\\s")
				.field(sortField, false);
			Future captain= s.sortInBackground();
			
			// feed
			BufferedReader buffy= new BufferedReader(new FileReader(source));
			BufferedWriter owriter= new BufferedWriter(new OutputStreamWriter(pout));
			char[] buf= new char[1024];
			int x= -1;
			while((x= buffy.read(buf))!= -1)
				owriter.write(buf, 0, x);
			owriter.close();
			
			// close handles
			captain.get();
			if (compression== FileHelper.COMPRESSION_ZIP) 
				((ZipOutputStream) ostream).closeEntry();
			ostream.close();
		}
		
		target.deleteOnExit();
	}
	
	protected void writeParFile(boolean keepSorted) throws Exception {
		UniversalReadDescriptor descriptor= new UniversalReadDescriptor();
		descriptor.init(UniversalReadDescriptor.getDescriptor("SIMULATOR"));
		FluxCapacitorSettings settings= new FluxCapacitorSettings();
		settings.set(FluxCapacitorSettings.ANNOTATION_FILE, 
				new File(gtfFile.getAbsolutePath()));
		settings.set(FluxCapacitorSettings.MAPPING_FILE, 
				new File(bedFile.getAbsolutePath()));
		settings.set(FluxCapacitorSettings.READ_DESCRIPTOR, 
				descriptor);
		settings.set(FluxCapacitorSettings.SORT_IN_RAM, 
				false);
		settings.set(FluxCapacitorSettings.KEEP_SORTED_FILES, 
				keepSorted);
		settings.set(FluxCapacitorSettings.ANNOTATION_MAPPING, 
				AnnotationMapping.PAIRED);
		settings.set(FluxCapacitorSettings.STDOUT_FILE, 
				outFile);
//		settings.set(FluxCapacitorSettings.COVERAGE_STATS, 
//				true);
//		settings.set(FluxCapacitorSettings.COVERAGE_FILE, 
//				new File(anoDir.getAbsolutePath()+ File.separator+ getClass().getSimpleName()+ "_coverage.txt"));
		
		BufferedWriter buffy= new BufferedWriter(new FileWriter(parFile));
		buffy.write(settings.toString());
		buffy.close();
		
		parFile.deleteOnExit();
	}

	protected void runCapacitor() throws Exception{
		FluxCapacitor capacitor= new FluxCapacitor();
		capacitor.setFile(parFile);
		Future<Void> captain= Execute.getExecutor().submit(capacitor);
		captain.get();
		outFile.deleteOnExit();
	}

	protected void initFiles(byte compressionGTF, int sortGTF, boolean writeProtectGTF,   
			byte compressionBED, int sortBED, boolean writeProtectBED, boolean keepSorted) {
		
		try {
			
			// file names
			initFileNames(compressionGTF, compressionBED);
	
			// put files in location
			copy(GTF_SORTED, gtfFile, sortGTF, compressionGTF);
			copy(BED_SORTED, bedFile, sortBED, compressionBED);
			writeParFile(keepSorted);
			
			// folder rights
			if (writeProtectGTF)
				anoDir.setReadOnly();
			if (writeProtectBED)
				mapDir.setReadOnly();
			
		} catch (Exception e) {
			throw new RuntimeException(e);
		}
	}
	
	@BeforeClass
	public static void initExecuter() {
		Execute.initialize(2);
	}
	
	@AfterClass
	public static void shutdownExecuter() {
		Execute.shutdown();
	}
	
	@Test
	public void testIOflatSortedWritableGTFflatSortedWritableBEDnoKeep() {

		try {					
			initFiles(
					// GTF: compressed, sorted, readOnly
					FileHelper.COMPRESSION_NONE, 
					SORTED,
					false,
					// BED: compressed, sorted, readOnly
					FileHelper.COMPRESSION_NONE,
					SORTED, 
					false,
					// keep sorted
					false);

			runCapacitor();
			
			// check
			assertTrue(gtfFile.exists());
			assertTrue(bedFile.exists());
			assertTrue(outFile.exists());
			String[] files= anoDir.list();
			assertTrue(files.length== 3);	// annotation+ parameter+ output
			files= mapDir.list();
			assertTrue(files.length== 1);	// mapping file only
			
		} catch (Exception e) {
			throw new RuntimeException(e);
		} finally {
			FileHelper.rmDir(mapDir);
			FileHelper.rmDir(anoDir);
		}

	}

	@Test
	public void testIOgzippedSortedWritableGTFflatSortedWritableBEDnoKeep() {

		try {					
			initFiles(
					// GTF: compressed, sorted, readOnly
					FileHelper.COMPRESSION_GZIP, 
					SORTED,
					false,
					// BED: compressed, sorted, readOnly
					FileHelper.COMPRESSION_NONE,
					SORTED, 
					false,
					// keep sorted
					false);

			runCapacitor();
			
			// check
			assertTrue(gtfFile.exists());
			assertTrue(bedFile.exists());
			assertTrue(outFile.exists());
			String[] files= anoDir.list();
			assertTrue(files.length== 3);	// annotation+ parameter+ output
			files= mapDir.list();
			assertTrue(files.length== 1);	// mapping file only
			
		} catch (Exception e) {
			throw new RuntimeException(e);
		} finally {
			FileHelper.rmDir(mapDir);
			FileHelper.rmDir(anoDir);
		}

	}

	@Test
	public void testRPKMwithNrReadsMapped() {
	
		try {					
			initFiles(
					// GTF: compressed, sorted, readOnly
					FileHelper.COMPRESSION_NONE, 
					SORTED,
					false,
					// BED: compressed, sorted, readOnly
					FileHelper.COMPRESSION_NONE,
					SORTED, 
					false,
					// keep sorted
					false);
			runCapacitor();
			File out1= outFile;
			
			
			initFiles(
					// GTF: compressed, sorted, readOnly
					FileHelper.COMPRESSION_NONE, 
					SORTED,
					false,
					// BED: compressed, sorted, readOnly
					FileHelper.COMPRESSION_NONE,
					SORTED, 
					false,
					// keep sorted
					false);
			BufferedWriter buffy= new BufferedWriter(new FileWriter(parFile, true));
			try {
				buffy.write(FluxCapacitorSettings.NR_READS_MAPPED.getName()+" "+
						Integer.toString(7897));
				buffy.close();
			} catch (Exception e) {
				throw new RuntimeException(e);
			}
			runCapacitor();

			// check
			try{
				BufferedReader b1= new BufferedReader(new FileReader(out1)), 
						b2= new BufferedReader(new FileReader(outFile));
				String s1, s2;
				while ((s1= b1.readLine())!= null&& (s2= b2.readLine())!= null) {
					System.err.println(s1);
					assertEquals(s1, s2);
				}
				assertFalse(b1.ready());
				assertFalse(b2.ready());
				
			} catch (Exception e) {
				throw new RuntimeException(e);
			}
			
		} catch (Exception e) {
			throw new RuntimeException(e);
		} finally {
			FileHelper.rmDir(mapDir);
			FileHelper.rmDir(anoDir);
		}
	
	}

}
