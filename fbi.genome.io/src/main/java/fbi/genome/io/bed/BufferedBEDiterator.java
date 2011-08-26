package fbi.genome.io.bed;

import java.util.Iterator;

import fbi.genome.model.bed.BEDobject2;

/**
 * The interface defines methods for marking a position in 
 * the data, to which the iterator subsequently can be 
 * repositioned.
 * 
 * @author Micha Sammeth (gmicha@gmail.com)
 *
 */
public interface BufferedBEDiterator extends Iterable<BEDobject2>, Iterator<BEDobject2>{

	/**
	 * marks the actual element
	 */
	public void mark();
	
	/**
	 * resets to last marked element (if any)
	 */
	public void reset();
}
