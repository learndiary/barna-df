package fbi.genome.sequencing.rnaseq.simulation.distributions;

public class NormalDistribution extends AbstractDistribution {

	double mean= Double.NaN, sd= Double.NaN;
	double sdSquare= Double.NaN;
	
	public NormalDistribution(double mean) {
		this.mean= mean;
	}
	
	public NormalDistribution(double mean, double sd) {
		this(mean);
		this.sd= sd;
		this.sdSquare= Math.pow(sd, 2d);
	}
	
	public double getP(double x) {
		
		return getP(x, getMean());
	}

	public double getP(double x, double mean) {
		
		double a1= Math.pow(x- mean, 2d);
		double p= Math.exp(-a1/ (2d* sdSquare))/ Math.sqrt(2* Math.PI* sdSquare);
		
		return p;
	}

	public double getRelFreq(double x) {
		return getRelFreq(x, getMean());
	}
	
	public double getRelFreq(double x, double mean) {
		
		return (getP(x, mean)/ getMax(mean));
	}

	private double getMax(double mean) {		
		return getP(mean, mean);
	}

	public double getMean() {
		return mean;
	}

}
