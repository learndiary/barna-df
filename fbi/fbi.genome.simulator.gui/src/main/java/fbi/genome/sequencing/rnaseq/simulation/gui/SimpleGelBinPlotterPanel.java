package fbi.genome.sequencing.rnaseq.simulation.gui;


import com.jhlabs.image.BoxBlurFilter;
import fbi.commons.gui.SimpleBinPlotterPanel;
import org.apache.commons.math.MathRuntimeException;
import org.apache.commons.math.random.JDKRandomGenerator;
import org.apache.commons.math.random.RandomDataImpl;
import org.apache.commons.math.util.MathUtils;

import java.awt.*;
import java.awt.image.BufferedImage;

public class SimpleGelBinPlotterPanel extends SimpleBinPlotterPanel {

	/*
	 * http://en.wikipedia.org/wiki/Poisson_distribution
	 * http://www.itl.nist.gov/div898/handbook/eda/section3/eda366j.htm
	 * 
	 * For sufficiently large values of lambda, (say lambda>1000), 
	 * the normal distribution with mean lambda and variance 
	 * lambda (standard deviation \sqrt{\lambda}), is an excellent 
	 * approximation to the Poisson distribution. If lambda is greater 
	 * than about 10, then the normal distribution is a good 
	 * approximation if an appropriate continuity correction is 
	 * performed, i.e., P(X <= x), where (lower-case) x is a 
	 * non-negative integer, is replaced by P(X <= x + 0.5).
	 */
	public static JDKRandomGenerator rand= new JDKRandomGenerator();
	public static RandomDataImpl rndGel= new RandomDataImpl(rand);
//	static RandomEngine engine = new DRand();
//	private static Poisson poisson = new Poisson(10d, engine);

	public static int getEffectiveLength(int len) {
		long llen= nextPoisson(len);
		//int poissonObs = poisson.nextInt(len);
		return ((int) llen);
	}
	
    public static long nextPoisson(double mean) {
        if (mean <= 0) {
            throw MathRuntimeException.createIllegalArgumentException(
                  "the Poisson mean must be positive ({0})", mean);
        }

        //RandomGenerator rand = getRan();

        double pivot = 6.0;
        if (mean < pivot) {
            double p = Math.exp(-mean);
            long n = 0;
            double r = 1.0d;
            double rnd = 1.0d;

            while (n < 1000 * mean) {
                rnd = rand.nextDouble();
                r = r * rnd;
                if (r >= p) {
                    n++;
                } else {
                    return n;
                }
            }
            return n;
        } else {
            double mu = (int) mean;	// Math.floor()
            double delta = (int) (pivot + (mu - pivot) / 2.0); // integer, Math.floor
            // between 6
            // and mean
            double mu2delta = 2.0 * mu + delta;
            double muDeltaHalf = mu + delta / 2.0;
            double logMeanMu = Math.log(mean / mu);

            double muFactorialLog = MathUtils.factorialLog((int) mu);

            double c1 = Math.sqrt(Math.PI * mu / 2.0);
            double c2 = c1 +
                        Math.sqrt(Math.PI * muDeltaHalf /
                                  (2.0 * Math.exp(1.0 / mu2delta)));
            double c3 = c2 + 2.0;
            double c4 = c3 + Math.exp(1.0 / 78.0);
            double c = c4 + 2.0 / delta * mu2delta *
                       Math.exp(-delta / mu2delta * (1.0 + delta / 2.0));

            double y = 0.0;
            double x = 0.0;
            double w = Double.POSITIVE_INFINITY;

            boolean accept = false;
            while (!accept) {
                double u = rndGel.nextUniform(0.0, c);
                double e = rndGel.nextExponential(mean);

                if (u <= c1) {
                    double z = rndGel.nextGaussian(0.0, 1.0);
                    y = -Math.abs(z) * Math.sqrt(mu) - 1.0;
                    x = (int) (y);	// Math.floor
                    w = -z * z / 2.0 - e - x * logMeanMu;
                    if (x < -mu) {
                        w = Double.POSITIVE_INFINITY;
                    }
                } else if (c1 < u && u <= c2) {
                    double z = rndGel.nextGaussian(0.0, 1.0);
                    y = 1.0 + Math.abs(z) * Math.sqrt(muDeltaHalf);
                    x = Math.ceil(y);
                    w = (-y * y + 2.0 * y) / mu2delta - e - x * logMeanMu;
                    if (x > delta) {
                        w = Double.POSITIVE_INFINITY;
                    }
                } else if (c2 < u && u <= c3) {
                    x = 0.0;
                    w = -e;
                } else if (c3 < u && u <= c4) {
                    x = 1.0;
                    w = -e - logMeanMu;
                } else if (c4 < u) {
                    double v = rndGel.nextExponential(mean);
                    y = delta + v * 2.0 / delta * mu2delta;
                    x = Math.ceil(y);
                    w = -delta / mu2delta * (1.0 + y / 2.0) - e - x * logMeanMu;
                }
                accept = (w <= x * Math.log(mu) -
                         MathUtils.factorialLog((int) (mu + x)) /
                         muFactorialLog);
            }
            // cast to long is acceptable because both x and mu are whole
            // numbers.
            return (long) (x + mu);
        }
    }
	
	public void paintOSI(Object data) {

		if (offScrImg== null) {
			offScrImg= new BufferedImage(getWidth()- getInsets().left- getInsets().right, 
					getHeight()- getInsets().top- getInsets().bottom, 
					BufferedImage.TYPE_INT_ARGB);
		} 
		Graphics gi= offScrImg.getGraphics();
		gi.setColor(plotBg);
		if ((!(paintMode== MODE_LINE))|| (osiCalled== lineColors.length)) {
			gi.fillRect(0, 0, offScrImg.getWidth(), offScrImg.getHeight());
			if (paintMode== MODE_LINE) {
				osiCalled= 0;
				lastXY= null;
			}
		}
		
		if (data== null) {
			if (paintMode== MODE_HISTO|| paintMode== MODE_BARPLOT|| paintMode== MODE_LINE) {
				paintOSIhisto(gi, globalBins, globalMaxXY);
			} else if (paintMode== MODE_GEL) {
				paintOSIgel(gi, globalBins, globalMaxXY);
			}			
		} else {
			MyArrays a= new MyArrays(data);
			if (paintMode== MODE_HISTO|| paintMode== MODE_BARPLOT|| paintMode== MODE_LINE) {
				double[] maxXY= new double[2];
				double[] bins= getBins(a, maxXY);
				paintOSIhisto(gi, bins, maxXY);
			} else if (paintMode== MODE_GEL) {
				double[] maxXY= new double[3];
				double[] bins= getBins(a, maxXY);
				paintOSIgel(gi, bins, maxXY);
			}
		}		
		++osiCalled;
	}
	
	void paintOSIgel(Graphics gi, double[] bins, double[] maxXY) {
			
			int h= offScrImg.getHeight();
			int xOff= 10, barWidth= 30;
			double minLogX= Math.log(Math.min(maxXY[2], ladder[0]));
			double maxLogX= Math.log(Math.max(maxXY[0], ladder[ladder.length-1]));
			//double m= h/ (minLogX- maxLogX);
			double m= -(h/ (maxLogX- minLogX+ 1));
			//double t= h- (m* minLogX);
			double t= h;
			
			// paint ladder
			for (int i = 0; i< ladder.length; i++) {
				int sc= 255;
				if (!invert)
					sc= 255- sc;	// invert
				Color c= new Color(sc, sc, sc);
				
				gi.setColor(c);
				double logLen= Math.log(ladder[i]);	// log10
				double ypos= (m* (maxLogX- logLen))+ t; 
				//System.out.println(ladder[i]+" "+logLen+" "+ypos);
				gi.fillRect(xOff, (int) ypos, barWidth, 3);			
			}
		
			// paint lane
			xOff+= 3+ barWidth;
			for (int i = 0; i < bins.length; i++) {
				paintOSIValgel(gi, bins, i, maxXY, m, t, maxLogX);
			}
			
			// blur
			BufferedImage tgt= new BufferedImage(getWidth(), getHeight(), BufferedImage.TYPE_INT_ARGB);
			BoxBlurFilter filter= new BoxBlurFilter();
			filter.setRadius(2);
			filter.setIterations(5);
			filter.filter(offScrImg, tgt);
			offScrImg= tgt;
			gi= offScrImg.getGraphics();
			
				// draw legend for ladder
			xOff= 2;
			gi.setColor(plotFg);
			for (int i = 0; i< ladder.length; i++) {
				double logLen= Math.log(ladder[i]);
				double ypos= (m* (maxLogX- logLen))+ t;
				String s= Integer.toString(ladder[i]); 	
				gi.drawString(s, xOff, (int) ypos+ 10);	
				
			}		
		
		}

	private void paintOSIValgel(Graphics gi, double[] bins, int pos,
				double[] maxXY, double m, double t, double maxLogX) {
			int sc= (int) ((bins[pos]/ maxXY[1])* 255);
			if (sc> 255)
				sc= 255;
			double binWidth= maxXY[0]/ bins.length;
			double minLen= pos* binWidth;
			double maxLen= (pos+1)* binWidth;
			if (!invert)
				sc= 255- sc;	// invert
			//double realLen= Math.round(((pos+1)/ (double) bins.length)* maxXY[0]);
			//double logLen= Math.log10(realLen);		
			//double ypos= (m* logLen)+ t;
			double logLenMin= Math.log(minLen);
			double logLenMax= Math.log(maxLen);
			double yposMin= (m* (maxLogX- logLenMin))+ t;
			double yposMax= (m* (maxLogX- logLenMax))+ t;
			double diff= yposMin- yposMax;
			if (diff< 1)
				diff= 1;
	//		sc*= (maxLen+ minLen)/ (2*diff);
			sc/= Math.log10(diff);
			if (sc> 255)
				sc= 255;
			sc= (int) (255/ (maxXY[1]/ bins[pos]));
			Color c= new Color(sc, sc, sc);
			gi.setColor(c);
			//gi.fillRect(50, (int) ypos, 30, 3);
			//gi.fillRect(50, (int) yposMax, 30, (int) diff);
			gi.fillRect(50, (int) yposMin, 30, (int) yposMin+ 1);
		}
}
