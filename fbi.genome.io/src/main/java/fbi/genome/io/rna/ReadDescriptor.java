package fbi.genome.io.rna;

public interface ReadDescriptor {

	public CharSequence getUniqueDescriptor(CharSequence descriptor);
	public byte getPairedEndInformation(CharSequence descriptor);
	public byte getStrand(CharSequence descriptor);
	public boolean isPairedEnd(CharSequence descriptor);
	public boolean isStranded(CharSequence descriptor);
	public boolean isApplicable(CharSequence descriptor);
	public boolean allowsPend();
	public boolean allowsStranded();

}
