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

import barna.model.Mapping;
import barna.model.rna.UniversalReadDescriptor;


public interface MappingReader extends Iterable<Mapping>{

    /**
     * Reads mappings from the underlying <code>InputStream</code>
     * overlapping the area specified by <code>chr</code>,
     * <code>start</code> and <code>end</code> and returns them as
     * an <code>MSIterator</code>.
     * @param chromosome chromosome name of the specified area
     * @param start start position of the specified area
     * @param end end position of the specified area
     */
    public MSIterator<Mapping> read(String chromosome, int start, int end);

	/**
	 * Retrieve the number of unique reads in the mapping set.
	 * @return the number of unique reads
	 */
	public int getCountReads();
	
	/**
	 * Retrieve the number of mappings in the mapping set.
	 * @return the number of mappings
	 */
	public int getCountMappings();
	
	/**
	 * Retrieve the number of continuous mappings in the mapping set.
	 * @return the number of continuous mappings
	 */
	public int getCountContinuousMappings();
	
	/**
	 * Retrieve the number of split mappings in the mapping set.
	 * @return the number of split mappings
	 */
	public int getCountSplitMappings();

	/**
	 * Checks if a given read descriptor can be applied to
	 * the data wrapped by the <code>this</code> instance.
	 * @param descriptor a read descriptor
	 * @return <code>true</code> if the mappings wrapped by <code>this</code>
	 * are applicable too the rules of the read descriptor, <code>false</code>
	 * otherwise 
	 */
	public boolean isApplicable(UniversalReadDescriptor descriptor);

    /**
     * Close the reader
     * @return true if the reader was closed without problems
     */
    public boolean close();

    /**
     * Reset the reader to the beginning of the file
     */
    public void reset();

    /**
     * Reset the reader to the start of the specified chromosome
     * @param chr the chromosome
     * @return true if the reader can be reset to the beginning of the chromosome
     */
    public boolean reset(String chr);

    boolean isPaired();
}
