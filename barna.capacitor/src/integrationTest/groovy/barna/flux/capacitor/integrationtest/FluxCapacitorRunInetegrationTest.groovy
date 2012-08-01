	
package barna.flux.capacitor.integrationtest

import barna.commons.Execute
import barna.commons.system.OSChecker
import barna.flux.capacitor.reconstruction.FluxCapacitorSettings
import barna.flux.capacitor.reconstruction.FluxCapacitorSettings.AnnotationMapping
import barna.io.FileHelper
import barna.io.Sorter
import barna.io.rna.UniversalReadDescriptor
import org.junit.AfterClass
import org.junit.BeforeClass
import org.junit.Test

import java.util.concurrent.Future
import java.util.zip.GZIPOutputStream
import java.util.zip.ZipEntry
import java.util.zip.ZipOutputStream

import static junit.framework.Assert.assertTrue
import static org.junit.Assert.fail
import org.junit.Before
import org.junit.After

/**
 * 
 * @author Thasso Griebel (thasso.griebel@gmail.com)
 */

class FluxCapacitorRunInetegrationTest {

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

    static String executable

    @BeforeClass
    public static void setUp(){
        executable = System.getProperty("dist.exe")
        if(executable == null){
            fail("No capacitor executable specified")
        }
        Execute.initialize(2);

    }

    @AfterClass
    public static void shutdownExecuter() {
        Execute.shutdown();
    }




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
	


	protected String runCapacitor() throws Exception{

        def pb = new ProcessBuilder()
        pb.environment().put("FLUX_MEM", "1G")
        if (tmpDir != null){
            pb.environment().put("JAVA_OPTS", "-Dflux.io.deny.tmpdir=yes")
        }
        def cmd = [executable, "-p", parFile.getAbsolutePath()]
        if (OSChecker.isWindows()) {
            cmd = ["cmd", "/c", executable, "-p", parFile.getAbsolutePath()]
        }
        def process = pb.directory(tmpDir != null ? tmpDir : parFile.getParentFile())
                .redirectErrorStream(true)
                .command(cmd)
                .start()
        String output = process.inputStream.text
        process.waitFor()
		return output;
		
		
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
			
	static final String[] STDERR_MAPPED= ["8009","8192"]
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
		assertTrue("""Number of gtf files does not match!
Expected : ${files.length}
Found    : ${nrFilesInGTF}
List of Files : ${files.join(", ")}
""", files.length== nrFilesInGTF);	// annotation+ parameter+ output
		files= mapDir.list();
		assertTrue(files.length== nrFilesInBED);	// mapping file only
	}

    File currentTestDirectory = null
    @Before
    public void setUpTest(){
        currentTestDirectory = FileHelper.createTempDir("cap-integration", "", null)
    }

    @After
    public void cleanup(){
        if(currentTestDirectory != null){
            FileHelper.rmDir(currentTestDirectory)
        }
    }

    @Test
    public void testIOflatSortedWritableGTFflatSortedWritableBEDnoKeep_new() {
        FluxCapacitorRunner.createTestDir(currentTestDirectory, [
               "ANNOTATION_FILE" : GTF_SORTED,
               "MAPPING_FILE" : BED_SORTED,
        ])
        String stderr= runCapacitor();
        assertFiles(2, 1, stderr, STDERR_MAPPED);
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

        println "Hello"
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