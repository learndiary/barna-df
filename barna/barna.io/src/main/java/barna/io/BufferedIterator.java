package barna.io;

import barna.commons.ByteArrayCharSequence;

import java.util.Iterator;

/**
 * The interface defines methods for marking a position in 
 * the data, to which the iterator subsequently can be 
 * repositioned.
 * 
 * @author Micha Sammeth (gmicha@gmail.com)
 *
 */
public interface BufferedIterator extends Iterable<ByteArrayCharSequence>, Iterator<ByteArrayCharSequence>{

	/**
	 * marks the actual element
	 */
	public void mark();
	
	/**
	 * resets to last marked element (if any)
	 */
	public void reset();
}
