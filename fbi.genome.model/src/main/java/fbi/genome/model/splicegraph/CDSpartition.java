package fbi.genome.model.splicegraph;

import fbi.genome.model.Translation;

import java.util.HashMap;

public class CDSpartition {

	Partition parent= null;
	boolean cdsValid53;
	int currFrame5= Translation.FRAME_BYTEVAL[Translation.FRAME_BITNI];
	int frame3= Translation.FRAME_BYTEVAL[Translation.FRAME_BITNI];
	public CDSpartition(Partition parent) {
		this.parent= parent;
		currFrame5= frame3= Translation.FRAME_BYTEVAL[Translation.FRAME_BITNI];
		cdsValid53= false;
	}
	public CDSpartition(Partition parent, int frame3) {
		this(parent);
		this.frame3= frame3;
	}

	public long[] getTranscripts() {
		return parent.getTranscripts();
	}
	public HashMap<PartitionSet, PartitionSet> getParents() {
		return parent.getParents();
	}
	
	@Override
	public boolean equals(Object obj) {
		
		if (!(obj instanceof CDSpartition))
			return false;
		
		if (!super.equals(obj))
			return false;
		
		CDSpartition p= (CDSpartition) obj;
		if (p.currFrame5!= currFrame5
				|| p.frame3!= frame3
				|| p.cdsValid53!= p.cdsValid53)
			return false;
		return true;
	}
	
	public CDSpartition clonePartitionWithoutTx() {
		CDSpartition p= new CDSpartition(parent);
		p.frame3= frame3;
		p.currFrame5= currFrame5;
		p.cdsValid53= cdsValid53;
		return p;
	}
	public int getCurrFrame5() {
		return currFrame5;
	}
	public void setCurrFrame5(int currFrame5) {
		this.currFrame5 = currFrame5;
		int combi= Translation.getCombinedFrame(getCurrFrame5(), getFrame3());
		if (combi== -65528)
			System.currentTimeMillis();
	}
	public int getFrame3() {
		return frame3;
	}
	public void setFrame3(int frame3) {
		this.frame3 = frame3;
		int combi= Translation.getCombinedFrame(getCurrFrame5(), getFrame3());
		if (combi== -65528)
			System.currentTimeMillis();
	}
	
	
}
