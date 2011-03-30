package fbi.genome.sequencing.rnaseq.reconstruction;

import fbi.genome.model.Transcript;
import fbi.genome.model.splicegraph.Edge;

public class Partition {

	Edge e;
	Transcript t;
	String id;
	
	
	public Partition(Edge e, Transcript t) {
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
