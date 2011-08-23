package fbi.genome.io.bed;

import java.util.Iterator;

import fbi.genome.model.bed.BEDobject2;

public interface BufferedBEDiterator extends Iterable<BEDobject2>, Iterator<BEDobject2>{

	/**
	 * marks the actual element
	 */
	public void mark();
	
	/**
	 * resets to last marked element
	 * @return success of reset operation
	 */
	public boolean reset();
}
