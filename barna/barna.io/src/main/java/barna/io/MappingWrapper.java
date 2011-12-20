package barna.io;


public interface MappingWrapper {

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

}
