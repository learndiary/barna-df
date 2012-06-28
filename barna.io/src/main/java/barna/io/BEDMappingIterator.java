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
import barna.model.bed.BEDMapping;

import java.util.ArrayList;
import java.util.Iterator;

/**
 * A class implementing the <code>MSIterator</code> interface
 * by iterating on an array in memory.
 * 
 * @author Micha Sammeth (gmicha@gmail.com)
 * @see BEDMappingIteratorDisk
 */
public class BEDMappingIterator implements MSIterator<BEDMapping>{
	
	/**
	 * Array of BED lines that are iterated
	 */
	ArrayList<BEDMapping> elements;
	
	/**
	 * Current position of iterator in underlying 
	 * array.
	 * @see #elements
	 */
	int currentIndex;
	
	/**
	 * Position marked for re-positioning later on.
	 * @see #currentIndex
	 * @see #mark()
	 * @see #reset()
	 */
	int markedIndex;
	
	/**
	 * Creates an instance iterating the elements 
	 * provided starting with the first one.
	 * @param elements array of BED lines
	 */
	public BEDMappingIterator(ArrayList<BEDMapping> elements) {
		this.elements= elements;
		currentIndex= 0;
	}

	/**
	 * Returns <code>this</code> instance implementing the
	 * <code>Iterator</code> interface.
	 * @return <code>this</code> iterator instance
	 */
	@Override
	public Iterator<BEDMapping> iterator() {
		return this;
	}

	@Override
	public boolean hasNext() {		
		return (elements!= null&& currentIndex< elements.size());
	}

	
	@Override
	public BEDMapping next() {
		
		if (elements== null|| currentIndex>= elements.size())
			return null;
		
		return elements.get(currentIndex++);
	}
	
	@Override
	public void remove() {
		// TODO Auto-generated method stub
	}
	
	@Override
	public void mark() {
		markedIndex= currentIndex;
	}
	
	@Override
	public void reset() {
		
		if(markedIndex< 0|| markedIndex>= elements.size())
			return;
		currentIndex= markedIndex;
	}

    @Override
    public void setAtStart() {
        currentIndex = 0;
    }

    /**
	 * frees memory occupied by the elements
	 */
	@Override
	public void clear() {
		elements= null;
	}
}
