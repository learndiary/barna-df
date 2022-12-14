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

import java.util.Iterator;

/**
 * The interface defines methods for marking a position in 
 * the data, to which the iterator subsequently can be 
 * repositioned.
 * 
 * @author Micha Sammeth (gmicha@gmail.com)
 *
 */
public interface MSIterator<T> extends Iterator<T>, Iterable<T>{

	/**
	 * marks the actual element
	 */
	public void mark();
	
	/**
	 * resets to last marked element (if any)
	 */
	public void reset();

    /**
     * reset to the start position
     */
    public void setAtStart();
	
	/**
	 * frees resources occupied by the iterator
	 */
	public void clear();

    /**
     * In case of paired reads retrieve the mate given a <code>Iterator&lt;Mapping&gt;</code>
     * @return the mate
     */
    public Iterator<Mapping> getMates(Mapping firstMapping);

    /**
     * Get the number of mappings.
     * @return (-1) if the size is unknown, otherwise the number of mappings stored in the underlying set
     * (which might be 0).
     */
    public int size();
}
