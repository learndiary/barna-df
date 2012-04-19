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

import java.io.*;
import java.util.zip.GZIPInputStream;


/**
 * Default implementation of the IOWrapper Interface
 * 
 * @author Micha Sammeth (gmicha@gmail.com)
 *
 */
public abstract class AbstractFileIOWrapper extends AbstractIOWrapper {
	
	/**
	 * The file from which all data is read
	 */
	protected File inputFile= null;
	
	/**
	 * Line separator used in <code>inputFile</code>.
	 */
	protected String lineSeparator= null;	// "\n"

	/**
	 * Number of bytes currently read from the start of the
	 * underlying file.
	 */
	protected long bytesRead= 0;
	
	/**
	 * Size of the underlying file in bytes.
	 */
	protected long size= -1;
	
	
	/**
	 * Creates an instance with initializing the 
	 * mandatory <code>inputFile</code> attribute.
	 * @param inputFile
	 */
	public AbstractFileIOWrapper(File inputFile) {
		setInputFile(inputFile);
	}
	
	/**
	 * Sets an <code>inputFile</code> the wrapper 
	 * is based on.
	 * @param inputFile the file the wrapper is 
	 * reading data from 
	 */
	public void setInputFile(File inputFile) {
		this.inputFile= inputFile;
		this.size= -1;
		this.bytesRead= 0;
		this.lineSeparator= null;
	}
	
	/** 
	 * Returns the underlying file.
	 * @return the file from which is read
	 */
	public File getInputFile() {
		return inputFile;
	}
	
	/**
	 * Returns line separator used in file, if none 
	 * has yet been specified, it is guessed from 
	 * the first line break.
	 * @return an instance of <code>String</code> 
	 * describing the characters used for a line 
	 * break
	 */
	public String getLineSeparator() {
		if (lineSeparator == null) {
			lineSeparator= FileHelper.guessFileSep(inputFile);
		}

		return lineSeparator;
	}	
	
	/**
	 * Returns the size of the file from which is read.
	 * @return size of <code>inputFile</code> in bytes.
	 */
	public long getInputSize() {
		
		if (size< 0) {
			size = inputFile.length();
		}

		return size;
	}
	
	/**
	 * Returns the number of bytes currently read from 
	 * <code>inputFile</code>.
	 * @return the number of bytes currently read from
	 * the underlying file.
	 */
	public long getBytesRead() {
		return bytesRead;
	}
	
	/**
	 * Scans the underlying file and reads corresponding 
	 * statistics.
	 */
	public abstract void scanFile();
	
	/**
	 * Reports the number of invalid lines in the file encountered
	 * from the file start to the current point by calls to 
	 * the <code>read()</code> or <code>scanFile()</code> methods.
	 * @see read, scanFile
	 * @return the number of lines that have been rejected
	 * so far while reading through the file.
	 */
	public abstract int getNrInvalidLines();


    /**
     * Opens a new {@link InputStream} on the input file. This
     * does NOT create any buffered streams, but handles gzip
     * files
     *
     *
     * @return inputStream a new input stream on the input file
     * @throws IOException in case the file is not found or the gzip stream could not be created
     * @since 1.0
     */
    protected InputStream getInputStream() throws IOException {
        if(inputFile == null ) throw new NullPointerException("No input file specified");
        InputStream inputStream = new FileInputStream(inputFile);
        if(inputFile.getName().toLowerCase().endsWith(".gz")){
            inputStream = new GZIPInputStream(inputStream);
        }
        return inputStream;
    }

    /**
     * Opens a new {@link Reader} on the input file. This
     * does NOT create any buffered streams, but handles gzip
     * files
     *
     * @return reader a new reader on the input file
     * @throws IOException in case the file is not found or the gzip reader could not be created
     * @since 1.0
     */
    protected Reader getReader() throws IOException {
        return new InputStreamReader(getInputStream());
    }


}
