package fbi.genome.sequencing.rnaseq.simulation;


import commons.ByteArrayCharSequence;

public interface ProgressablePlotter {
	public boolean plot(int start, int end, int segLen, ByteArrayCharSequence ID);
	public void addBase(ByteArrayCharSequence ID, int length, int mol);
	public void reset(String message);
	public void setMolTot(long value);
	public void paint();
}
