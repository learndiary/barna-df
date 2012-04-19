/*
 * Copyright (c) 2010, Micha Sammeth
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * The names of its contributors may be not used to endorse or promote
 *       products derived from this software without specific prior written
 *       permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL MICHA SAMMETH BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

package barna.io.bed;

import barna.commons.CharsequenceComparator;
import barna.commons.utils.StringUtils;
import barna.io.rna.UniversalReadDescriptor;
import barna.io.rna.UniversalReadDescriptor.Attributes;

import java.util.Comparator;

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
		
		CharSequence ss1= StringUtils.getField(s1, '\t', 3),
			ss2= StringUtils.getField(s2, '\t', 3);
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
