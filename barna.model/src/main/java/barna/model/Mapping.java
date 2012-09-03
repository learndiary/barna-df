/**
 * 
 */
package barna.model;

/**
 * @author emilio
 *
 */
public interface Mapping {
	
	/**
	 * 
	 * @return the name of the mapping
	 */
	public CharSequence getName();

	/**
	 * 
	 * @return the name of the reference chromosome	
	 */
	public CharSequence getChromosome();
	
	/**
	 * 
	 * @return mapping start position
	 */
	public int getStart();
	
	/**
	 * 
	 * @return mapping end position
	 */
	public int getEnd();
	
	/**
	 * 
	 * @return mapping length
	 */
	public int getLength();
	
	/**
	 * 
	 * @return mapping score
	 */
	public int getScore();
	
	/**
	 * 
	 * @return mapping strand direction
	 */
	public byte getStrand(); // TODO maybe something like boolean isReverseStrand()?
	
	/**
	 * 
	 * @return number of blocks
	 */
	public int getBlockCount();
	
	/**
	 * 
	 * @return the start position of the next block
	 */
	public int getNextBlockStart();
	
	/**
	 * 
	 * @return the size of the next block
	 */
	public int getNextBlockSize();
}
