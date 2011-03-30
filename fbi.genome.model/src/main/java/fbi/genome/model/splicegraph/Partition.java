package fbi.genome.model.splicegraph;

import java.util.HashMap;
import java.util.Iterator;


public class Partition {
	HashMap<PartitionSet, PartitionSet> parents;
	long[] transcripts;
	public Partition() {
		transcripts= null;
		parents= new HashMap<PartitionSet,PartitionSet>(4,1f);
	}
	public Partition(long[] tx) {
		this();
		this.transcripts= tx;
	}
	
	
	public void addParent(PartitionSet newParent) { 
		newParent.partitions.put(this,this);
		parents.put(newParent,newParent);
	}
	
	@Override
	public String toString() {
		StringBuilder sb= new StringBuilder("[");
		if (transcripts== null)
			sb.append("null]");
		else {
			for (int i = 0; i < transcripts.length; i++) 
				sb.append(Long.toString(transcripts[i])+",");
			sb.deleteCharAt(sb.length()- 1);
			sb.append("]");
		}
		return sb.toString();
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
		
		if (!(obj instanceof Partition))
			return false;
		
		Partition p= (Partition) obj;
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
	
	public Partition clonePartitionWithoutTx() {
		
		Partition p= new Partition();
		p.parents= (HashMap<PartitionSet, PartitionSet>) parents.clone();
		Iterator<PartitionSet> iter= parents.keySet().iterator();
		while (iter.hasNext())
			iter.next().partitions.put(p, p);
		return p;
	}
	public HashMap<PartitionSet, PartitionSet> getParents() {
		return parents;
	}
	public void setParents(HashMap<PartitionSet, PartitionSet> parents) {
		this.parents = parents;
	}
	public long[] getTranscripts() {
		return transcripts;
	}
	public void setTranscripts(long[] transcripts) {
		this.transcripts = transcripts;
	}
}


