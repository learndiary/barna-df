/*
 * Copyright (c) 2010, Micha Sammeth
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * The names of its contributors may be not used to endorse or promote
 *       products derived from this software without specific prior written
 *       permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL MICHA SAMMETH BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

package barna.io;

import barna.commons.ByteArrayCharSequence;
import barna.commons.Execute;
import barna.commons.log.Log;

import java.io.*;
import java.util.Comparator;
import java.util.Iterator;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Future;

/**
 * A class implementing the <code>BufferedBEDiterator</code> interface
 * by reading data from a file on disk.
 * 
 * @author Micha Sammeth (gmicha@gmail.com)
 * @see BufferedIteratorRAM
 */
public class BufferedIteratorDisk implements BufferedIterator {

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
	 * A file with BED lines provided as input. 
	 */
	File inputFile;
	/**
	 * A stream with BED lines provided as input.
	 */
	InputStream inputStream;
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
	 * Creates an instance with <i>sorted</i> BED lines read from a stream
	 * and written to the intermediate file.
	 * 
	 * @param istream input stream with <b>sorted</b> lines
	 * @param tmpFile temporary file storing the content of the input stream
	 */
	public BufferedIteratorDisk(InputStream istream, File tmpFile) {
		this(istream, tmpFile, null, -1);
	}
	
	/**
	 * Creates an instance with <i>unsorted</i> BED lines read from a stream
	 * and sorted to the intermediate file employing the comparator provided.
	 * 
	 * @param istream input stream with <b>unsorted</b> lines
	 * @param tmpFile temporary file storing the content of the input stream
	 * @param comparator rules of comparison
	 */
	public BufferedIteratorDisk(InputStream istream, File tmpFile, Comparator<CharSequence> comparator) {
		this(istream, tmpFile, comparator, -1);
	}
	
	/**
	 * Creates an instance with <i>unsorted</i> BED lines read from a stream
	 * and sorted to the intermediate file employing the comparator provided
	 * and a certain capacity in bytes for the reading buffer.
	 * 
	 * @param istream input stream with <b>unsorted</b> lines
	 * @param tmpFile temporary file storing the content of the input stream
	 * @param comparator rules of comparison
	 * @param capacity capacity of the reader
	 * @see #reader
	 */
	public BufferedIteratorDisk(InputStream istream, File tmpFile, Comparator<CharSequence> comparator, int capacity) {
		this.inputStream= istream;
		this.tmpFile= tmpFile;
		this.comparator= comparator;
		this.capacity= capacity;
	}
	
	/**
	 * Creates an instance with BED lines read from a file that is
	 * <i>already sorted</i>.
	 * @param inputFile <b>already sorted</b> input file
	 */
	public BufferedIteratorDisk(File inputFile) {
		this(inputFile, null, null, -1);
	}
	
	/**
	 * Creates an instance with BED lines read from a file that is
	 * <i>already sorted</i>, using the given number of bytes as
	 * reading buffer capacity.
	 * @param inputFile <b>already sorted</b> input file
	 * @param capacity number of bytes used for the reading buffer
	 */
	public BufferedIteratorDisk(File inputFile, int capacity) {
		this(inputFile, null, null, capacity);
	}
	
	/**
	 * Creates an instance with BED lines read from an <i>unsorted
	 * file</i> and sorts it to a temporary file employing the 
	 * given <code>Comparator</code> instance.
	 * @param inputFile <b>unsorted</b> input file
	 * @param tmpFile intermediate <b>sorted</b> file that is created
	 * @param comparator comparison applied for sorting 
	 */
	public BufferedIteratorDisk(File inputFile, File tmpFile, Comparator<CharSequence> comparator) {
		this(inputFile, tmpFile, comparator, -1);
	}
	
	/**
	 * Creates an instance with BED lines read from an <i>unsorted
	 * file</i> and sorts it to a temporary file employing the 
	 * given <code>Comparator</code> instance and the given 
	 * number of bytes as reading buffer capacity.
	 * @param inputFile <b>unsorted</b> input file
	 * @param tmpFile intermediate <b>sorted</b> file that is created
	 * @param comparator comparison applied for sorting 
	 * @param capacity number of bytes used for the reading buffer
	 */
	public BufferedIteratorDisk(File inputFile, File tmpFile, Comparator<CharSequence> comparator, int capacity) {
		this.inputFile= inputFile;
		this.tmpFile= tmpFile;
		this.capacity= capacity;
		this.comparator= comparator;
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
	public File getTmpFile() throws ExecutionException, IOException {
		
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

		if (tmpFile== null&& (!(comparator== null&& inputFile!= null)))
			throw new RuntimeException("Temporary file required for iterating");
		
		if (comparator== null) {	// assume sorted
			// save sorted stream to file
			if (inputFile== null) {				
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
				this.captain= null; // TODO deadlocks: Execute.getExecutor().submit(callme);
			}
			
		} else {	// unsorted
			
			if (inputFile== null) {
				FileOutputStream fos= new FileOutputStream(tmpFile);
				Sorter s= Sorter.create(inputStream, fos, true, null)
					.field(comparator);
				this.captain= s.sortInBackground();

			} else {
				
				Callable<Void> callme= new Callable<Void>() {
					@Override
					public Void call() throws Exception {
						PipedOutputStream pout= new PipedOutputStream();
						PipedInputStream pin= new PipedInputStream(pout);
						Copier copy= new Copier(inputFile, pout);
						Future future1= Execute.getExecutor().submit(copy);
						FileOutputStream fos= new FileOutputStream(tmpFile);
						Sorter s= Sorter.create(pin, fos, true, null)
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
			if (pos> 0)
				reader.currBytes= pos;
		}

		return reader;
	}
	
	/**
	 * Returns the current instance implementing the 
	 * <code>Iterator</code> interface.
	 * @return <code>this</code> instance
	 */
	@Override
	public Iterator<ByteArrayCharSequence> iterator() {		
		return this;
	}

	/**
	 * Determines whether there are more elements by checking 
	 * the number of bytes available for reading in the buffer.
	 * @return <code>true</code> if there are more characters
	 * available, <code>false</code> otherwise
	 * @see barna.io.BufferedBACSReader#available()
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
	public ByteArrayCharSequence next() {
		if (!hasNext())
			return null;
		try {
			reader= getReader(-1);
			cs= reader.readLine(cs);
			// bedObject2 clones byte[]
			return cs;
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
				Log.error(e + "Error reading from file " + tmpFile.getName());
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

	/**
	 * removes the temporary file
	 */
	@Override
	public void clear() {
		boolean b= false;
		try {
			reader.close();
			b= tmpFile.delete();
		} catch (IOException e) {
			if (!b)
				throw new RuntimeException(e);
		}
	}
	
}