package fbi.genome.io.bed;

import fbi.genome.model.bed.BEDobject;

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
