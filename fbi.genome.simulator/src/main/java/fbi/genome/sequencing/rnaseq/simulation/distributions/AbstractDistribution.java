package fbi.genome.sequencing.rnaseq.simulation.distributions;

/**
 * Generic wrapper for empiric or analytic distributions.
 * 
 * @author micha
 *
 */
public abstract class AbstractDistribution {

	public abstract double getP(double x);
	
	public abstract double getRelFreq(double x);
	
	public abstract double getP(double x, double mean);
	
	public abstract double getMean();
	
	/**
	 * The weight of the distribution in composites.  
	 */
	double weight= Double.NaN;
	
	/**
	 * The mean of the distribution.
	 */
	double mean= Double.NaN;

	public double getWeight() {
		return weight;
	}

	public void setWeight(double weight) {
		this.weight = weight;		
	}

	
}
