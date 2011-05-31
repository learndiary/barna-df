package fbi.genome.sequencing.rnaseq.simulation.distributions;

public class CompositeDistribution extends AbstractDistribution {

	AbstractDistribution[] d= null;
	
	public CompositeDistribution(AbstractDistribution[] distributions) {
		d= distributions;
	}
	
	@Override
	public double getP(double x) {
		double p= 0;
		for (int i = 0; i < d.length; i++) 
			p+= d[i].getWeight()* d[i].getP(x, mean);
		return p;
	}

	@Override
	public double getRelFreq(double x) {
		double f= 0;
		for (int i = 0; i < d.length; i++) 
			f+= d[i].getWeight()* d[i].getRelFreq(x);
		return f;
	}

	@Override
	public double getP(double x, double mean) {
		double p= 0;
		for (int i = 0; i < d.length; i++) 
			p+= d[i].getWeight()* d[i].getP(x, mean);
		return p;
	}

	@Override
	public double getMean() {
		if (Double.isNaN(mean)) {
			mean= 0d;
			for (int i = 0; i < d.length; i++) 
				mean+= d[i].getWeight()* d[i].getMean();
			mean/= getWeightSum();
		}

		return mean;
	}
	
	public double getWeightSum() {
		double sum= 0d;
		for (int i = 0; i < d.length; i++) 
			sum+= d[i].getWeight();
		return sum;
	}

}
