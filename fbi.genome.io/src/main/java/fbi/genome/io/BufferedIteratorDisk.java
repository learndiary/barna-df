package fbi.genome.io;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.PipedInputStream;
import java.io.PipedOutputStream;
import java.util.Comparator;
import java.util.Iterator;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Future;

import fbi.commons.ByteArrayCharSequence;
import fbi.commons.CharsequenceComparator;
import fbi.commons.Execute;
import fbi.commons.Log;
import fbi.commons.tools.Sorter;
import fbi.genome.model.bed.BEDobject2;

/**
 * A class implementing the <code>BufferedBEDiterator</code> interface
 * by reading data from a file on disk.
 * 
 * @author Micha Sammeth (gmicha@gmail.com)
 * @see BufferedIteratorMemory
 */
public class BufferedIteratorDisk implements Bufferediterator {

	/**
	 * Class for capsulating file copy processes: stream-to-file, 
	 * file-to-stream, and stream-to-stream. The actual process of 
	 * copying data from the source to the target can be executed
	 * in parallel by implementing the <code>Callable</code>
	 * interface.
	 * 
	 * @author Micha Sammeth (gmicha@gmail.com)
	 */
	static class Copier implements Callable<Void> {
		
		/**
		 * Default size of the buffer used for copying data from 
		 * the source to the target.
		 */
		public static int DEFAULT_BUFFER_SIZE= 100;
		
		/**
		 * A source stream.
		 */
		InputStream inputStream;
		/**
		 * A target stream.
		 */
		OutputStream outputStream;
		/**
		 * A source file.
		 */
		File inputFile;
		/**
		 * A target file.
		 */
		File outputFile;
		/**
		 * The actual size of the buffer that is used to copy from
		 * the source to the target.
		 */
		int bufferSize;
		
		/**
		 * Constructor to copy from a source stream to a target stream,
		 * employing a default buffer size.
		 * @param inputStream the source stream
		 * @param outputStream the target stream
		 */
		public Copier(InputStream inputStream, OutputStream outputStream) {
			this(inputStream, outputStream, DEFAULT_BUFFER_SIZE);
		}
		/**
		 * Constructor to copy from a source file to a target stream,
		 * employing a default buffer size.
		 * @param inputFile the source file
		 * @param outputStream the target stream
		 */
		public Copier(File inputFile, OutputStream outputStream) {
			this(inputFile, outputStream, DEFAULT_BUFFER_SIZE);
		}
		/**
		 * Constructor to copy from a source stream to a target file,
		 * employing a default buffer size.
		 * @param inputStream the source stream
		 * @param outputFile the target file
		 */
		public Copier(InputStream inputStream, File outputFile) {
			this(inputStream, outputFile, DEFAULT_BUFFER_SIZE);
		}
		/**
		 * Constructor to copy from a source stream to a target stream,
		 * employing a custom buffer size.
		 * @param inputStream the source stream
		 * @param outputStream the target stream
		 * @param bufferSize the size of the buffer used for copying
		 */
		public Copier(InputStream inputStream, OutputStream outputStream, int bufferSize) {
			this.inputStream= inputStream;
			this.outputStream= outputStream;
			this.bufferSize= bufferSize;
		}
		
		/**
		 * Constructor to copy from a source file to a target stream,
		 * employing a custom buffer size.
		 * @param inputFile the source file
		 * @param outputStream the target stream
		 * @param bufferSize the size of the buffer used for copying
		 */
		public Copier(File inputFile, OutputStream outputStream, int bufferSize) {
			this.inputFile= inputFile;
			this.outputStream= outputStream;
			this.bufferSize= bufferSize;
		}
		
		/**
		 * Constructor to copy from a source stream to a target file,
		 * employing a custom buffer size.
		 * @param inputStream the source stream
		 * @param outputFile the target file
		 * @param bufferSize the size of the buffer used for copying
		 */
		public Copier(InputStream inputStream, File outputFile, int bufferSize) {
			this.inputStream= inputStream;
			this.outputFile= outputFile;
			this.bufferSize= bufferSize;
		}
		
		/**
		 * Determines the nature of source and target (i.e., whether 
		 * streams or files have been provided), and implements the
		 * copy process as reading/writing blocks of 
		 * <code>bufferSize</code> from the source to the target, 
		 * respectively.
		 * @throws Exception
		 */
		@Override
		public Void call() throws Exception {
			
			if (inputFile!= null) 
				inputStream= new FileInputStream(inputFile);
			if (outputFile!= null) 
				outputStream= new FileOutputStream(outputFile);
			
			byte[] b= new byte[bufferSize];
			while(inputStream.read(b)>= 0) {
				int k= inputStream.read(b);
				outputStream.write(b, 0, k);
			}
			
			outputStream.flush();
			if (inputFile!= null) 
				inputStream.close();
			if (outputFile!= null) 
				outputStream.close();
			
			return null;
		}
	}
	
	/**
	 * Inner comparator for string comparisons emplyoing the 
	 * <code>CharSequence</code> interface. 
	 */
	Comparator<CharSequence> comparator;
	/**
	 * The temporary (sorted!) file that is iterated by 
	 * <code>this</code> instance. If input data is provided
	 * unsorted and/or in form of a stream, the file is created.
	 * In this case, a new file is create in the current 
	 * temporary directory if no other is specified by 
	 * <code>directory</code>, involving an optional prefix. 
	 * @see #prefix
	 * @see #directory
	 */
	File tmpFile;
	/**
	 * <code>Future</code> object which knows when the sorted file
	 * is available. 
	 * @see #init()
	 */
	Future captain; // allow flexible types of different Callables
	/**
	 * The directory where the sorted file is created.
	 * @see tmpFile
	 * @see prefix 
	 */
	File directory;
	/**
	 * The prefix of the sorted file.
	 * @see directory
	 * @see tmpFile
	 */
	String prefix;
	/**
	 * The suffix of the sorted file.
	 * @see tmpFile
	 * @see prefix
	 * @see directory
	 */
	String suffix=".bed";
	/**
	 * A file with BED lines provided as input. 
	 */
	File inputFile;
	/**
	 * A stream with BED lines provided as input.
	 */
	InputStream inputStream;
	/**
	 * Flag indicated whether the BED lines from 
	 * the input are already sorted.
	 */
	boolean sorted= false;
	/**
	 * Internal reader instance.
	 */
	BufferedBACSReader reader;
	
	/**
	 * Number of bytes used for the reader buffer.
	 * @see #reader
	 */
	int capacity;
	
	/**
	 * An instance to re-use the <code>byte[]</code>
	 * representing the BED lines iterated.
	 */
	ByteArrayCharSequence cs= null;
	
	/**
	 * Flag indicating whether the <code>init()</code>
	 * method has already been invoked.
	 */
	boolean inited= false;
	
	/**
	 * Creates an instance with BED lines read from a stream. As no 
	 * comparator is provided, <code>CharsequenceComparator</code>
	 * is used as default.
	 * @param istream the input stream
	 * @param sorted flag to indicate whether BED lines are sorted
	 */
	public BufferedIteratorDisk(InputStream istream, boolean sorted) {
		this(istream, sorted, CharsequenceComparator.DEFAULT_CHARSEQUENCE_COMPARATOR, -1, null, null);
	}
	
	/**
	 * Creates an instance with BED lines read from a stream and 
	 * sorted according to a custom <code>Comparator</code>.
	 * @param istream the input stream
	 * @param sorted flag to indicate whether BED lines are sorted
	 * @param comparator the comparator that is used for sorting
	 */
	public BufferedIteratorDisk(InputStream istream, boolean sorted, Comparator<CharSequence> comparator) {
		this(istream, sorted, comparator, -1, null, null);
	}
	
	/**
	 * Creates an instance with BED lines read from a stream and 
	 * sorted according to a custom <code>Comparator</code>. 
	 * Additionally the capacity of the buffer used for reading can 
	 * be declared (>0).
	 * @param istream the input stream
	 * @param sorted flag to indicate whether BED lines are sorted
	 * @param comparator the comparator that is used for sorting
	 * @param capacity value specifying the capacity of the reader
	 * @see #reader
	 */
	public BufferedIteratorDisk(InputStream istream, boolean sorted, Comparator<CharSequence> comparator, int capacity) {
		this(istream, sorted, comparator, capacity, null, null);
	}
	
	/**
	 * Creates an instance with BED lines read from a stream and 
	 * sorted according to a custom <code>Comparator</code>. 
	 * Additionally the capacity of the buffer used for reading (>0),
	 * and the prefix of the sorted file can be declared. 
	 * @param istream the input stream
	 * @param sorted flag to indicate whether BED lines are sorted
	 * @param comparator the comparator that is used for sorting
	 * @param capacity value specifying the capacity of the reader
	 * @param prefix the prefix of the sorted file's name
	 * @see #reader
	 * @see #tmpFile
	 */
	public BufferedIteratorDisk(InputStream istream, boolean sorted, Comparator<CharSequence> comparator, String prefix) {
		this(istream, sorted, comparator, -1, prefix, null);
	}
	
	/**
	 * Creates an instance with BED lines read from a stream and 
	 * sorted according to a custom <code>Comparator</code>. 
	 * Additionally the capacity of the buffer used for reading (>0),
	 * and directory and prefix of the sorted file can be declared. 
	 * @param istream the input stream
	 * @param sorted flag to indicate whether BED lines are sorted
	 * @param comparator the comparator that is used for sorting
	 * @param capacity value specifying the capacity of the reader
	 * @param prefix the prefix of the sorted file's name
	 * @param directory folder where the sorted file is created
	 * @see #reader
	 * @see #tmpFile
	 * @see #directory
	 */
	public BufferedIteratorDisk(InputStream istream, boolean sorted, Comparator<CharSequence> comparator, int capacity, String prefix, File directory) {		
		this.inputStream= istream;
		this.sorted= sorted;
		this.comparator= comparator;
		this.prefix= prefix;
		this.directory= directory;
		this.capacity= capacity;
	}
	
	/**
	 * Creates an instance with BED lines read from a file.
	 * @param inputFile the input file
	 * @param sorted flag to indicate whether BED lines are sorted
	 */
	public BufferedIteratorDisk(File inputFile, boolean sorted) {
		this(inputFile, sorted, -1, null, null);
	}
	
	/**
	 * Creates an instance with BED lines read from a file. 
	 * Additionally the number of bytes used for the read-buffer
	 * can be declared.
	 * @param inputFile the input file
	 * @param sorted flag to indicate whether BED lines are sorted
	 * @param capacity number of bytes used for the reading buffer
	 */
	public BufferedIteratorDisk(File inputFile, boolean sorted, int capacity) {
		this(inputFile, sorted, capacity, null, null);
	}
	
	/**
	 * Creates an instance with BED lines read from a file. 
	 * Additionally the number of bytes used for the read-buffer,
	 * and the sorted file's name and directory can be declared.
	 * @param inputFile the input file
	 * @param sorted flag to indicate whether BED lines are sorted
	 * @param capacity number of bytes used for the reading buffer
	 * @param prefix prefix of the sorted file's name
	 * @param directory folder where the sorted file is created 
	 */
	public BufferedIteratorDisk(File inputFile, boolean sorted, int capacity, String prefix, File directory) {
		this.inputFile= inputFile;
		this.sorted= sorted;
		this.prefix= prefix;
		this.directory= directory;
		this.capacity= capacity;
	}
	
	/**
	 * Creates a <code>File</code> instance according to the provided 
	 * parameters <code>prefix</code> and <code>directory</code>. The
	 * suffix is fixed to &quot;.bed&quot;
	 * @return a file handle of the temporary file that is iterated 
	 */
	protected File createTmpFile() {
		String pfx= (prefix== null?this.getClass().getName() : prefix);
		try {
			if (directory== null)
				return File.createTempFile(pfx, suffix);
			else
				return File.createTempFile(pfx, suffix, directory);
		} catch (IOException e) {
			Log.error("Couldn't create temporary file "+ pfx+ "*"+ suffix+
					" in "+ (directory== null? System.getProperty("java.io.tmpdir"): directory.getAbsolutePath()));
			e.printStackTrace(Log.logStream);
			return null;
		}
	}
	
	/**
	 * Obtains the sorted file the <code>BEDiteratorDisk</code> instance
	 * is based on. <b>Note</b>: intrinsically invokes <code>init()</code>
	 * if no explicit call has incurred yet.
	 * @return the sorted file on which <code>this</code> iterator is
	 * based on
	 * @throws ExecutionException
	 * @throws IOException
	 * @see #init()
	 */
	protected File getTmpFile() throws ExecutionException, IOException {
		
		if (!inited) 
			init();
		
		while (captain != null) {
			try {
				captain.get();
				captain= null;
			} catch (InterruptedException e) {
				; // :)
			}
		}

		return tmpFile;
	}
	
	/**
	 * Produces a sorted file for subsequent iteration from the 
	 * correspondingly provided input (i.e., streams or files, 
	 * sorted or unsorted). Possible copying and sorting processes
	 * are executed in parallel, the result is availble as soon 
	 * as the corresponding <code>Future</object> returns.<br>  
	 * <b>Note</b>: the method has to be called before iteration 
	 * can start. It is recommended to call <code>init()</code>
	 * directly after the constructor. Otherwise, an intrinsic 
	 * call will ensure correct initialization by the first time
	 * <code>hasNext()</code> or <code>next()</code> are invoked.
	 * @throws FileNotFoundException
	 * @throws IOException
	 * @see #getTmpFile()
	 * @see #captain
	 */
	public void init() throws FileNotFoundException, IOException {
		if (inputStream== null&& inputFile== null) 
			Log.error("No input data");
		
		if (sorted) {
			// save sorted stream to file
			if (inputFile== null) {				
				tmpFile= createTmpFile();
				Copier copy= new Copier(inputStream, tmpFile);
				this.captain= Execute.getExecutor().submit(copy);
				
			} else {
				tmpFile= inputFile;
				Callable<Void> callme= new Callable<Void>() {
					@Override
					public Void call() throws Exception {
						return null;
					}
				};
				this.captain= Execute.getExecutor().submit(callme);
			}
			
		} else {	// unsorted
			
			if (inputFile== null) {
				tmpFile= createTmpFile();
				FileOutputStream fos= new FileOutputStream(tmpFile);
				Sorter s= Sorter.create(inputStream, fos, true)
					.field(comparator);
				this.captain= s.sortInBackground();

			} else {
				
				tmpFile= createTmpFile();
				Callable<Void> callme= new Callable<Void>() {
					@Override
					public Void call() throws Exception {
						PipedOutputStream pout= new PipedOutputStream();
						PipedInputStream pin= new PipedInputStream(pout);
						Copier copy= new Copier(inputFile, pout);
						Future future1= Execute.getExecutor().submit(copy);
						FileOutputStream fos= new FileOutputStream(tmpFile);
						Sorter s= Sorter.create(pin, fos, true)
							.field(comparator);
						Future future2= s.sortInBackground();
						future1.get();
						future2.get();
						return null;
					}
				};
				this.captain= Execute.getExecutor().submit(callme); 
			}
		}
		inited= true;
	}
	
	/**
	 * Returns the current reader or instantiates 
	 * @param pos position in bytes where the obtained 
	 * reader will start to read <i>iff</i> pos>0.
	 * @return a <code>BufferedBACSReader</code> instance
	 * @throws FileNotFoundException 
	 * @throws IOException
	 */
	private BufferedBACSReader getReader(long pos) throws FileNotFoundException, IOException {
		if (reader == null) {
				
			FileInputStream fis= null;
			
			try {
				fis= new FileInputStream(getTmpFile());
			} catch (ExecutionException e) {
				throw new RuntimeException(e);
			}
			if (pos> 0)
				fis.skip(pos);
			if (capacity> 0)
				reader = new BufferedBACSReader(fis, capacity);
			else
				reader = new BufferedBACSReader(fis);
		}

		return reader;
	}
	
	/**
	 * Returns the current instance implementing the 
	 * <code>Iterator</code> interface.
	 * @return <code>this</code> instance
	 */
	@Override
	public Iterator<BEDobject2> iterator() {		
		return this;
	}

	/**
	 * Determines whether there are more elements by checking 
	 * the number of bytes available for reading in the buffer.
	 * @return <code>true</code> if there are more characters
	 * available, <code>false</code> otherwise
	 * @see fbi.genome.io.BufferedBACSReader#available()
	 */
	@Override
	public boolean hasNext() {
		
		try {
			reader= getReader(-1);
			return (reader.available()> 0);
		} catch (Exception e) {
			throw new RuntimeException(e);
		}
		
	}

	/**
	 * Obtains the next element of <code>BEDobject2</code> instances
	 * from the iterator.
	 * @return the next BED line of this iterator
	 */
	@Override
	public BEDobject2 next() {
		
		try {
			reader= getReader(-1);
			cs= reader.readLine(cs);
			// bedObject2 clones byte[]
			return new BEDobject2(cs);
		} catch (Exception e) {
			throw new RuntimeException(e);
		}
		
	}

	/**
	 * Removing elements is not supported as data is kept in
	 * an underlying file. The method does absolutely nothing
	 * and is just kept to implement the interface.
	 */
	@Override
	public void remove() {
		// TODO Auto-generated method stub
	}

	/**
	 * Marks the current position as target for re-positioning
	 * later on. 
	 * @see #reset()
	 */
	@Override
	public void mark() {
		try {
			reader= getReader(-1);
		} catch (Exception e) {
			throw new RuntimeException(e);
		}
		
		reader.mark();
	}

	/**
	 * Resets the reader to the last marked position.
	 * If marked position is still in the reader buffer,
	 * it directly gets repositioned accordingly. Otherwise,
	 * the stream is closed and re-opened, sweeping directly
	 * to the marked position.
	 * @see #mark() 
	 */
	@Override
	public void reset() {
		
		try {
			reader= getReader(-1);
		} catch (Exception e) {
			throw new RuntimeException(e);
		}
		long pos= reader.reset();
		if (pos>= 0) {
			try {
				reader.close();
				reader= null;
				getReader(pos);
			} catch (FileNotFoundException e) {
				Log.error(e+ " reading from file "+ tmpFile.getName());
			} catch (IOException e) {
				Log.error(e+ "Error reading from file "+ tmpFile.getName());
			}
		}
			
	}

	/**
	 * Determines the number of elements left in the 
	 * iterator, starting to count at the current
	 * position.<br>
	 * <b>Warning</b>: time complexity of the method
	 * scales linearly with the elements that are 
	 * left and have to be read from the file. 
	 * @return the number of elements from the current
	 * reading position until the end
	 */
	public int countRemainingElements() {
	
		mark();
		int n= 0;
		while (hasNext()) {
			next();
			++n;
		}
		reset();
		
		return n;
	}
	
}