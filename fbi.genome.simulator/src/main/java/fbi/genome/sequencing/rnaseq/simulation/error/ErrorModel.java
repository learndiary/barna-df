package fbi.genome.sequencing.rnaseq.simulation.error;

public interface ErrorModel {

	public static final byte TYPE_MUTATION= 1, TYPE_INSERTION= 2, TYPE_DELETION= 3;
	
	public void setBaseProbability(double p);
	public double getBaseProbability();
	public void apply(byte[] quals);
	public void apply(byte[] quals, int from, int to);
	
	//public void addChar(char[] sequence, int i, String template, int j, Random rnd);
}
