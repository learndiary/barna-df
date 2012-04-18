/*
 * This file is part of the Flux Library.
 *
 * The code of the Flux Library may be freely distributed and modified under the terms of the
 * European Union Public Licence (EUPL) published on the web site <http://www.osor.eu/eupl/european-union-public-licence-eupl-v.1.1>.
 * Copyright for the code is held jointly by the individual authors, who should be listed
 * in @author doc comments. According to Article 5 and Article 11 of the EUPL, publications that
 * include results produced by the Flux Library are liable to reference the Work,
 * see the Flux Library homepage <http://flux.sammeth.net> for more information.
 */

package barna.io.bed;

import barna.model.bed.BEDobject;

public class BEDreadObject extends BEDobject {
	public static char SEP_PE= '_', SEP= ':';
	byte pEnd= 0, mm= -1;
	int tot= -1;
	
	public BEDreadObject(String chromName, byte newStrand, int newStart, int newEnd) {
		super(chromName, newStrand, newStart, newEnd);
	}
	
	@Override
	public void setName(String name) {
		
		int lp= name.length(), p= lp- 1;
		while (p>= 0&& mm< 0) {
			while(p>= 0&& name.charAt(p)!=SEP&& name.charAt(p)!=SEP_PE)
				--p;
			if (p< 0)
				break;
			if (name.charAt(p)== SEP_PE)
				pEnd= Byte.parseByte(name.substring(p+1,lp));
			else if (tot< 0)
				tot= Integer.parseInt(name.substring(p+1,lp));
			else 
				mm= Byte.parseByte(name.substring(p+1,lp));
			lp= p;
		}		
		
		super.setName(name.substring(0,p));
	}

	public byte getPEnd() {
		return pEnd;
	}

	public byte getMm() {
		return mm;
	}

	public int getTot() {
		return tot;
	}
}
