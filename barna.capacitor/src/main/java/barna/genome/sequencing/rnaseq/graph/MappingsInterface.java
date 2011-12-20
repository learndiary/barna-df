package barna.genome.sequencing.rnaseq.graph;

import barna.model.Transcript;

public interface MappingsInterface {

	/**
	 * Return the maximum mappable length that can map to the edge. 
	 * @param maxMapLength maximum mapped read length allowed
	 * @return maximum mapping length that can map to the edge.
	 */
	public abstract int getMapLength(Transcript tx, int maxMapLength);
	
	public abstract double getNtCoverage(Transcript tx, byte dir, int mapLenMax);
	
	public abstract float[] getCoverage(boolean sense, int minMapLen, int maxMapLen);

	public abstract Mappings getMappings();

}
