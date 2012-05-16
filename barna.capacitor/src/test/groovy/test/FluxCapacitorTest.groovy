	
package test

import barna.commons.Execute
import barna.flux.capacitor.reconstruction.FluxCapacitorSettings
import barna.flux.capacitor.reconstruction.FluxCapacitorSettings.AnnotationMapping
import barna.io.FileHelper
import barna.io.Sorter
import barna.io.rna.UniversalReadDescriptor
import java.util.concurrent.Future
import java.util.zip.GZIPOutputStream
import java.util.zip.ZipEntry
import java.util.zip.ZipOutputStream
import org.junit.AfterClass
import org.junit.BeforeClass
import org.junit.Test
import static junit.framework.Assert.assertTrue

/**
 * 
 * @author Thasso Griebel (thasso.griebel@gmail.com)
 */

class FluxCapacitorTest{

	static final int SORTED= -1, UNSORT_GTF= 2, UNSORT_BED= 9;
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
		anoDir= FileHelper.createTempDir(getClass().getSimpleName(), subdirAnnotation, null);
		anoDir.deleteOnExit();
		mapDir= FileHelper.createTempDir(getClass().getSimpleName(), subdirMappings, null);
		mapDir.deleteOnExit();
		outFile= File.createTempFile(getClass().getSimpleName(), suffixOutput, null);	// for checking later on write protection on anoDir
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
				try {
					((ZipOutputStream) ostream).closeEntry();
				} catch (Exception e) {
					; // :)
				}
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
		if (tmpDir!= null)
			settings.set(FluxCapacitorSettings.TMP_DIR, tmpDir);
		
		BufferedWriter buffy= new BufferedWriter(new FileWriter(parFile));
		buffy.write(settings.toString());
		buffy.close();
		
		parFile.deleteOnExit();
	}

	protected String runCapacitor() throws Exception{
		
		// start process
        /**
         * Get rid of any IDEA jars files in the classpath
         */
        def classpath = System.getProperty("java.class.path")
        def cp  = []
        // remove libs containing spaces in the name, eg "blabla IDEA bla.jar"
        classpath.split(":").each {e->
            if(!e.contains(" ")) cp << e
        }
        def cpp= cp.join(":")
        String cmd= "java -cp "+ cpp

		if (tmpDir!= null)
			cmd+= " -Dflux.io.deny.tmpdir=yes"
		cmd+= " -Xmx1G barna.commons.launcher.Flux -t capacitor -p "+parFile.getAbsolutePath()

        System.out.println("Try executing : " + cmd);
		Process process= cmd.execute()
		process.waitFor()

		// internal start
		//		FluxCapacitor capacitor= new FluxCapacitor();
		//		capacitor.setFile(parFile);
		//		Future<Void> captain= Execute.getExecutor().submit(capacitor);
		//		captain.get();
		//		outFile.deleteOnExit();
		
		String stderr= "STDERR: ${process.err.text}"
		return stderr;
		
		
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
			
	static final String[] STDERR_MAPPED= ["8009","8044"]
	static final String[] STDERR_ACCESS_DENIED= ["access denied"]
	
	void assertFiles(int nrFilesInGTF, int nrFilesInBED, String stderr, String[] occurrences) {

        System.out.println("Error output : " + stderr);
		
		// current stderr output
		//		[INFO] 	7897 reads, 8009 mappings: R-factor 1.0141826
		//		[INFO] 	6145 entire, 1864 split mappings (2.3273816%)
		//
		//		1 single transcript loci
		//		8009 mappings in file
		//		566 mappings fall in single transcript loci
		//		283 mappings map to annotation
		//		1172 mappings in annotation-mapped pairs
		//		40,75 min/max read length
		//
		//		8009 mappings read from file
		//		8044 mapping pairs map to annotation
		//		0 pairs without tx evidence
		//		208 pairs in wrong orientation
		//		0 single mappings forced
		//
		//		4 transcripts, 4 detected
				
		//		println "EXIT: ${process.exitValue()}"
		//		println "STDOUT: ${process.in.text}"
		//		// alternative to combine wait and stdout
		//		println "ls ${parameter}".execute().text
		for (int i = 0; i < occurrences.length; i++) {
			assertTrue(stderr.contains(occurrences[i]))
		}
		
		
		
		
		assertTrue(gtfFile.exists());
		assertTrue(bedFile.exists());
		if (occurrences!= STDERR_ACCESS_DENIED)
			assertTrue(outFile.exists());
		String[] files= anoDir.list();
		assertTrue(files.length== nrFilesInGTF);	// annotation+ parameter+ output
		files= mapDir.list();
		assertTrue(files.length== nrFilesInBED);	// mapping file only
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

			String stderr= runCapacitor();
			
			assertFiles(2, 1, stderr, STDERR_MAPPED);
			
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

			String stderr= runCapacitor();
			
			assertFiles(2, 1, stderr, STDERR_MAPPED);
			
		} catch (Exception e) {
			throw new RuntimeException(e);
		} finally {
			FileHelper.rmDir(mapDir);
			FileHelper.rmDir(anoDir);
		}

	}
	
    @Test
	public void testIOzippedSortedWritableGTFflatSortedWritableBEDnoKeep() {
	
		try {
			initFiles(
					// GTF: compressed, sorted, readOnly
					FileHelper.COMPRESSION_ZIP,
					SORTED,
					false,
					// BED: compressed, sorted, readOnly
					FileHelper.COMPRESSION_NONE,
					SORTED,
					false,
					// keep sorted
					false);
	
			String stderr= runCapacitor();
			
			assertFiles(2, 1, stderr, STDERR_MAPPED);
			
		} catch (Exception e) {
			throw new RuntimeException(e);
		} finally {
			FileHelper.rmDir(mapDir);
			FileHelper.rmDir(anoDir);
		}
	
	}
	
	@Test
	public void testIOflatSortedWritableGTFgzippedSortedWritableBEDnoKeep() {

		try {
			initFiles(
					// GTF: compressed, sorted, readOnly
					FileHelper.COMPRESSION_NONE,
					SORTED,
					false,
					// BED: compressed, sorted, readOnly
					FileHelper.COMPRESSION_GZIP,
					SORTED,
					false,
					// keep sorted
					false);

			String stderr= runCapacitor();
			
			assertFiles(2, 1, stderr, STDERR_MAPPED);
			
		} catch (Exception e) {
			throw new RuntimeException(e);
		} finally {
			FileHelper.rmDir(mapDir);
			FileHelper.rmDir(anoDir);
		}

	}

	@Test
	public void testIOflatSortedWritableGTFzippedSortedWritableBEDnoKeep() {

		try {
			initFiles(
					// GTF: compressed, sorted, readOnly
					FileHelper.COMPRESSION_NONE,
					SORTED,
					false,
					// BED: compressed, sorted, readOnly
					FileHelper.COMPRESSION_ZIP,
					SORTED,
					false,
					// keep sorted
					false);

			String stderr= runCapacitor();
			
			assertFiles(2, 1, stderr, STDERR_MAPPED);
			
		} catch (Exception e) {
			throw new RuntimeException(e);
		} finally {
			FileHelper.rmDir(mapDir);
			FileHelper.rmDir(anoDir);
		}

	}

	@Test
	public void testIOgzippedSortedWritableGTFzippedSortedWritableBEDnoKeep() {

		try {
			initFiles(
					// GTF: compressed, sorted, readOnly
					FileHelper.COMPRESSION_GZIP,
					SORTED,
					false,
					// BED: compressed, sorted, readOnly
					FileHelper.COMPRESSION_ZIP,
					SORTED,
					false,
					// keep sorted
					false);

			String stderr= runCapacitor();
			
			assertFiles(2, 1, stderr, STDERR_MAPPED);
			
		} catch (Exception e) {
			throw new RuntimeException(e);
		} finally {
			FileHelper.rmDir(mapDir);
			FileHelper.rmDir(anoDir);
		}

	}

	@Test
	public void testIOzippedSortedWritableGTFgzippedSortedWritableBEDnoKeep() {

		try {
			initFiles(
					// GTF: compressed, sorted, readOnly
					FileHelper.COMPRESSION_ZIP,
					SORTED,
					false,
					// BED: compressed, sorted, readOnly
					FileHelper.COMPRESSION_GZIP,
					SORTED,
					false,
					// keep sorted
					false);

			String stderr= runCapacitor();
			
			assertFiles(2, 1, stderr, STDERR_MAPPED);
			
		} catch (Exception e) {
			throw new RuntimeException(e);
		} finally {
			FileHelper.rmDir(mapDir);
			FileHelper.rmDir(anoDir);
		}

	}

	@Test
	public void testIOflatUnsortedWritableGTFflatSortedWritableBEDnoKeep() {

		try {
			initFiles(
					// GTF: compressed, sorted, readOnly
					FileHelper.COMPRESSION_NONE,
					UNSORT_GTF,
					false,
					// BED: compressed, sorted, readOnly
					FileHelper.COMPRESSION_NONE,
					SORTED,
					false,
					// keep sorted
					false);

			String stderr= runCapacitor();
			
			assertFiles(2, 1, stderr, STDERR_MAPPED);
			
		} catch (Exception e) {
			throw new RuntimeException(e);
		} finally {
			FileHelper.rmDir(mapDir);
			FileHelper.rmDir(anoDir);
		}

	}

	@Test
	public void testIOflatSortedWritableGTFflatUnsortedWritableBEDnoKeep() {

		try {
			initFiles(
					// GTF: compressed, sorted, readOnly
					FileHelper.COMPRESSION_NONE,
					SORTED,
					false,
					// BED: compressed, sorted, readOnly
					FileHelper.COMPRESSION_NONE,
					UNSORT_BED,
					false,
					// keep sorted
					false);

			String stderr= runCapacitor();
			
			assertFiles(2, 1, stderr, STDERR_MAPPED);
			
		} catch (Exception e) {
			throw new RuntimeException(e);
		} finally {
			FileHelper.rmDir(mapDir);
			FileHelper.rmDir(anoDir);
		}

	}

	@Test
	public void testIOflatUnSortedWritableGTFflatUnsortedWritableBEDkeep() {

		try {
			initFiles(
					// GTF: compressed, sorted, readOnly
					FileHelper.COMPRESSION_NONE,
					UNSORT_GTF,
					false,
					// BED: compressed, sorted, readOnly
					FileHelper.COMPRESSION_NONE,
					UNSORT_BED,
					false,
					// keep sorted
					true);

			String stderr= runCapacitor();
			
			assertFiles(3, 2, stderr, STDERR_MAPPED);
			
		} catch (Exception e) {
			throw new RuntimeException(e);
		} finally {
			FileHelper.rmDir(mapDir);
			FileHelper.rmDir(anoDir);
		}

	}

	@Test
	public void testIOflatUnSortedReadOnlyGTFflatUnsortedReadOnlyBEDkeep() {

		try {
			initFiles(
					// GTF: compressed, sorted, readOnly
					FileHelper.COMPRESSION_NONE,
					UNSORT_GTF,
					true,
					// BED: compressed, sorted, readOnly
					FileHelper.COMPRESSION_NONE,
					UNSORT_BED,
					true,
					// keep sorted
					true);

			String stderr= runCapacitor();
			
			assertFiles(2, 1, stderr, STDERR_MAPPED);
			
		} catch (Exception e) {
			throw new RuntimeException(e);
		} finally {
			FileHelper.rmDir(mapDir);
			FileHelper.rmDir(anoDir);
		}

	}
	
	@Test
	public void testIOzippedUnSortedReadOnlyGTFgzippedUnsortedReadOnlyBEDnoKeep() {

		try {
			initFiles(
					// GTF: compressed, sorted, readOnly
					FileHelper.COMPRESSION_ZIP,
					UNSORT_GTF,
					true,
					// BED: compressed, sorted, readOnly
					FileHelper.COMPRESSION_GZIP,
					UNSORT_BED,
					true,
					// keep sorted
					false);

			String stderr= runCapacitor();
			
			assertFiles(2, 1, stderr, STDERR_MAPPED);
			
		} catch (Exception e) {
			throw new RuntimeException(e);
		} finally {
			FileHelper.rmDir(mapDir);
			FileHelper.rmDir(anoDir);
		}

	}

	@Test
	public void testTmpDir() {
		try {
			
			// create tmp dir
			tmpDir= new File(System.getProperty("java.io.tmpdir"));// FileHelper.createTempDir(getClass().getSimpleName(), "aTempDir", null);
			
			initFiles(
					// GTF: compressed, sorted, readOnly
					FileHelper.COMPRESSION_NONE,
					SORTED,
					true,
					// BED: compressed, sorted, readOnly
					FileHelper.COMPRESSION_NONE,
					SORTED,
					true,
					// keep sorted
					false);

			String stderr= runCapacitor();
			
			assertFiles(2, 1, stderr, STDERR_ACCESS_DENIED);
			
		} catch (Exception e) {
			throw new RuntimeException(e);
		} finally {
			FileHelper.rmDir(mapDir);
			FileHelper.rmDir(anoDir);
		}
	}

}