/**
 * 
 */
package barna.io;

import barna.commons.CharsequenceComparator;
import barna.model.Mapping;
import barna.model.rna.UniversalReadDescriptor;
import barna.model.rna.UniversalReadDescriptor.Attributes;

import java.util.Comparator;

/**
 * @author Emilio
 *
 */
public class MappingComparator implements Comparator<Mapping> {
	
	/**
	 * Wrapped comparator to compare general objects implementing 
	 * the <code>CharSequence</code> interface. 
	 */
	CharsequenceComparator comp;
	/**
	 * The descriptor defining the format of attributes concatenated
	 * in the Mapping name.
	 */
	UniversalReadDescriptor descriptor;
	
	/**
	 * Creates a default instance applying the provided 
	 * <code>UniversalReadDescriptor</code> for subsequent
	 * comparisons.
	 * @param descriptor rules how to extract attributes from
	 * the mapping name
	 */
	/**
	 * 
	 */
	public MappingComparator(UniversalReadDescriptor descriptor) {
		this.descriptor= descriptor;
		this.comp= new CharsequenceComparator();
	}

	@Override
	public int compare(Mapping o1, Mapping o2) {
		CharSequence ss1= o1.getName(),
				ss2= o2.getName();
		if (ss1== null|| ss2== null)
			throw new RuntimeException("failed to get mapping name: "+
					(ss1== null? o1: "")+
					(ss2== null? o2: ""));
		Attributes a1= descriptor.getAttributes(ss1, null),
			a2= descriptor.getAttributes(ss2, null);
		
		if (a1== null) 
			return (a2== null? 0:-1);		
		if (a2== null)
			return 1;
		
		// check ID
		int c= comp.compare(a1.id, a2.id);
		if (c!= 0)
			return c;
		
		// check mates
		if (descriptor.isPaired()) {
			c= a1.flag- a2.flag;
			if (c!= 0)
				return c;
		}
		
		// check strand
		if (descriptor.isStranded()) {
			c= a1.strand- a2.strand;
			if (c!= 0)
				return c;
		}					
		
		return 0;
	}

}
