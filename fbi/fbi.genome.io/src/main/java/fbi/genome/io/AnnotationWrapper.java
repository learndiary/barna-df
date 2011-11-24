package fbi.genome.io;

public interface AnnotationWrapper {

	/**
	 * Retrieve the number of genes in the annotation.
	 * @return the number of genes
	 */
	public int getNrGenes();
	
	/**
	 * Retrieve the number of transcripts in the annotation.
	 * @return the number of transcripts
	 */
	public int getNrTranscripts();
	
	/**
	 * Retrieve the number of exons in the annotation.
	 * @return the number of exons
	 */
	public int getNrExons();
}
