package fbi.genome.io;

import java.util.Iterator;

import fbi.commons.ByteArrayCharSequence;

/**
 * A class implementing the <code>BufferedBEDiterator</code> interface
 * by iterating on an array in memory.
 * 
 * @author Micha Sammeth (gmicha@gmail.com)
 * @see BufferedIteratorDisk
 */
public class BufferedIteratorRAM implements BufferedIterator{
	
	/**
	 * Array of BED lines that are iterated
	 */
	ByteArrayCharSequence[] elements;
	
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
	public BufferedIteratorRAM(ByteArrayCharSequence[] elements) {
		this.elements= elements;
		currentIndex= 0;
	}

	/**
	 * Returns <code>this</code> instance implementing the
	 * <code>Iterator</code> interface.
	 * @return <code>this</code> iterator instance
	 */
	@Override
	public Iterator<ByteArrayCharSequence> iterator() {
		return this;
	}

	@Override
	public boolean hasNext() {		
		return (elements!= null&& currentIndex< elements.length);
	}

	
	@Override
	public ByteArrayCharSequence next() {
		
		if (elements== null|| currentIndex>= elements.length)
			return null;
		
		return elements[currentIndex++];
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
		
		if(markedIndex< 0|| markedIndex>= elements.length)
			return;
		currentIndex= markedIndex;
	}
}
