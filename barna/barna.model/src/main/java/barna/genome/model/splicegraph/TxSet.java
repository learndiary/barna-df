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

package barna.genome.model.splicegraph;

public class TxSet {

	long[] transcripts;
	
	public TxSet(long[] trpts) {
		this.transcripts= trpts;
	}
	
	@Override
	public int hashCode() {
		int sum= 0;
		for (int i = 0; transcripts!= null&& i < transcripts.length; i++) {
			sum+= transcripts[i];
		}
		return sum;
	}
	
	@Override
	public boolean equals(Object obj) {
		
		if (!(obj instanceof TxSet))
			return false;
		
		TxSet p= (TxSet) obj;
		if (transcripts== null) {
			if (p.transcripts== null)
				return true;	// ??
			else 
				return false;
		}
		if (transcripts.length!= p.transcripts.length)
			return false;
		
		for (int i = 0; i < transcripts.length; i++) {
			if (transcripts[i]!= p.transcripts[i])
				return false;
		}
		return true;
	}
}
