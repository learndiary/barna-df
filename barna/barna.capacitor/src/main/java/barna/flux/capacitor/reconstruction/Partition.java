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

package barna.flux.capacitor.reconstruction;

import barna.model.Transcript;
import barna.model.splicegraph.SimpleEdge;

public class Partition {

	SimpleEdge e;
	Transcript t;
	String id;
	
	
	public Partition(SimpleEdge e, Transcript t) {
		this.e= e;
		this.t= t;
		this.id= e.toString()+ t.toString();
	}
	
	@Override
	public int hashCode() {		
		return id.hashCode();
	}
	
	@Override
	public boolean equals(Object obj) {
		if (!(obj instanceof Partition))
			return false;
		String id2= ((Partition) obj).id;
		
		return id.equals(id2);
	}
	
}
