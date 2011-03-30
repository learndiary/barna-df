package fbi.genome.model.splicegraph;

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
