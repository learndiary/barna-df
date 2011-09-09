/*
 * This file is part of the Flux Library.
 *
 * The code of the Flux Library may be freely distributed and modified under the terms of the
 * European Union Public Licence (EUPL) published on the web site <http://www.osor.eu/eupl/european-union-public-licence-eupl-v.1.1>.
 * Copyright for the code is held jointly by the individual authors, who should be listed
 * in @author doc comments. According to Article 5 and Article 11 of the EUPL, publications that
 * include results produced by the Flux Library are liable to reference the Work,
 * see the Flux Library homepage <http://flux.sammeth.net> for more information.
 */

package fbi.genome.io;

import java.io.File;



/**
 * Default implementation of the IOWrapper Interface
 * 
 * @author Micha Sammeth (gmicha@gmail.com)
 *
 */
public abstract class AbstractFileIOWrapper implements IOWrapper {
	
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
}
