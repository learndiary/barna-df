package fbi.genome.io.rna;

public class GingerasDescriptor extends SolexaPairedEndDescriptor {
	@Override
	public boolean allowsStranded() {
		return true;
	}
	public byte getStrand(CharSequence descriptor) {
		// TODO Auto-generated method stub
		return 0;
	}
	public boolean isStranded(CharSequence descriptor) {
		// TODO Auto-generated method stub
		return false;
	}
}
