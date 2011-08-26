package fbi.genome.io.bed;

import java.util.Comparator;

import fbi.commons.ByteArrayCharSequence;
import fbi.commons.CharSequences;
import fbi.commons.CharsequenceComparator;
import fbi.genome.io.rna.UniversalReadDescriptor;
import fbi.genome.io.rna.UniversalReadDescriptor.Attributes;
import fbi.genome.model.bed.BEDobject2;
import fbi.genome.model.commons.Sequence;
import fbi.genome.model.constants.Constants;

/**
 * Compares two <code>CharSequence</code> instances describing  
 * an entire BED line hierarchically by their attributes ID, 
 * mate, strand. Requires an <code>UniversalReadDescriptor</code>
 * instance for parsing instructions. 
 * <code>UniversalReadDescriptor</code>.
 * 
 * @author Micha Sammeth (gmicha@gmail.com)
 *
 */
public class BEDDescriptorComparator implements Comparator<CharSequence> {

	/**
	 * Wrapped comparator to compare general objects implementing 
	 * the <code>CharSequence</code> interface. 
	 */
	CharsequenceComparator comp;
	/**
	 * The descriptor defining the format of attributes concatenated
	 * in the BED name field.
	 */
	UniversalReadDescriptor descriptor;
	
	/**
	 * Creates a default instance applying the provided 
	 * <code>UniversalReadDescriptor</code> for subsequent
	 * comparisons.
	 * @param descriptor rules how to extract attributes from
	 * the BED name
	 */
	public BEDDescriptorComparator(UniversalReadDescriptor descriptor) {
		this.descriptor= descriptor;
		this.comp= new CharsequenceComparator();
	}
	
	@Override
	/**
	 * Compares two BED lines by the attributes in their 
	 * name field (4th column). Criteria are checked 
	 * hierachically: ID, mate, strand.<br>
	 * Example:<br>
	 * XXX:S/1<br>
	 * XXX:A/1<br>
	 * XXX:A/2<br>
	 * YYY:S/2<br>
	 * ZZZ:S/2<br>
	 * The provided <code>CharSequence</code> instances
	 * describe the <u>entire line</u> of a BED of which
	 * the name field is obtained for comparison. 
	 * @param s1 a BED line 
	 * @param s1 another BED line
	 * @return (-1), 0, 1 corresponding to whether the 
	 * first BED line is less, equal to, or greater than
	 * the second BED line argument. 
	 */
	public int compare(CharSequence s1, CharSequence s2) {
		
		CharSequence ss1= CharSequences.getField(s1, '\t', 3),
			ss2= CharSequences.getField(s2, '\t', 3);
		if (ss1== null|| ss2== null)
			throw new RuntimeException("failed to extract BED name: "+
					(ss1== null? s1: "")+
					(ss2== null? s2: ""));
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
