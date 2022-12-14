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
	 * @param appendMateNumber flag
	 * @return the name of the mapping
	 */
	public CharSequence getName(Boolean appendMateNumber);

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
     * @return mapping mate number
     */
    public byte getMateFlag();
	
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
     *
     * @return the mapped read sequence
     */
    public CharSequence getSequence();

    /**
     *
     * @return the extended Cigar string
     */
    public CharSequence getCigar();

    /**
     *
     * @param weighted whether reporting the count weighted by the number of hits for the read the mapping refers to
     * @return the (weighted) number of mappings
     */
    public double getCount(boolean weighted);

    /**
     * @param readStrand the directionality of the reads
	 * @return the strand
     */
    public byte getReadStrand(String readStrand);
}
