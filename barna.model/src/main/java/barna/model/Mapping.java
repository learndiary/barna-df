/**
 * 
 */
package barna.model;

import barna.commons.ByteArrayCharSequence;

/**
 * @author emilio
 *
 */
public interface Mapping {
	
	/**
	 * 
	 * @return the name of the mapping
	 */
	public ByteArrayCharSequence getName();

	/**
	 * 
	 * @return the name of the reference chromosome	
	 */
	public ByteArrayCharSequence getChromosome();
	
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
	
	/**
	 * Sets the name of the mapping
	 * @param name the name to be set
	 */
	public void setName(CharSequence name);
	
	/**
	 * Sets the name of the reference chromosome
	 * @param chr the name to be set
	 */
	public void setChromosome(CharSequence chr);
	
	/**
	 * Sets the start position of the mapping
	 * @param start the position to be set
	 */
	public void setStart(int start);
	
	/**
	 * Sets the end position of the mapping
	 * @param end the position to be set
	 */
	public void setEnd(int end);
	
	/**
	 * Sets the score of the mapping
	 * @param score the score to be set
	 */
	public void setScore(int score);
	
	/**
	 * Sets the strand direction
	 * @param strand the direction to be set
	 */
	public void setStrand(byte strand); // TODO boolean?
	
	/**
	 * Sets the number of blocks of the mappings
	 * @param count the number to be set
	 */
	public void setBlockCount(int count);
	
	/**
	 * Sets the start position of the next block
	 * @param start the position to be set
	 */
	public void setNextBlockStart(int start);
	
	/**
	 * Sets the size of the next block
	 * @param size the size to be set
	 */
	public void setNextBlockSize(int size);	
}
