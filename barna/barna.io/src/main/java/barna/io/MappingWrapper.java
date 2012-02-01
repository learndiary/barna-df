package barna.io;

import barna.io.rna.ReadDescriptor;
import barna.io.rna.UniversalReadDescriptor;


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

	/**
	 * Checks if a given read descriptor can be applied to
	 * the data wrapped by the <code>this</code> instance.
	 * @param descriptor a read descriptor
	 * @return <code>true</code> if the mappings wrapped by <code>this</code>
	 * are applicable too the rules of the read descriptor, <code>false</code>
	 * otherwise 
	 */
	public boolean isApplicable(UniversalReadDescriptor descriptor);
}
