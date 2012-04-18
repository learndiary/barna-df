package barna.flux.capacitor.graph;


public class Mappings {

	
	/**
	 * Coverage of genomic positions.
	 */
	float[] coverage= null;
	float[] coverageRev= null;
	int readNr = 0;
	int revReadNr = 0;

	public void incrReadNr() {
		++readNr;
	}

	public void incrRevReadNr() {
		++revReadNr;
	}

	public void decrReadNr() {
		--readNr;
	}

	public void decrRevReadNr() {
		--revReadNr;
	}

	public int getReadNr() {
		return readNr;
	}

	public void setReadNr(int readNr) {
		this.readNr = readNr;
	}

	public int getRevReadNr() {
		return revReadNr;
	}

	public void setRevReadNr(int revReadNr) {
		this.revReadNr = revReadNr;
	}


}
