/*
 * This file is part of the Flux Library.
 *
 * The code of the Flux Library may be freely distributed and modified under the terms of the
 * European Union Public Licence (EUPL) published on the web site <http://www.osor.eu/eupl/european-union-public-licence-eupl-v.1.1>.
 * Copyright for the code is held jointly by the individual authors, who should be listed
 * in @author doc comments. According to Article 5 and Article 11 of the EUPL, publications that
 * include results produced by the Flux Library are liable to reference the Work,
 * see the Flux Library homepage <http://flux.sammeth.net> for more information.
 */

package fbi.genome.sequencing.rnaseq.simulation.distributions;

import fbi.commons.file.FileHelper;
import fbi.genome.model.commons.DoubleVector;

import java.io.*;

/**
 * An empirical distribution stored in a histogram.
 *
 * @author micha
 */
public class EmpiricalDistribution extends AbstractDistribution {

    /**
     * Minimum, maximum of the underlying distribution.
     */
    double min = Double.NaN, max = Double.NaN;

    /**
     * Integral under histogram stored in <code>a</code>.
     */
    double sum = Double.NaN;

    /**
     * Maximum frequency in histogram <code>a</code>, can be a float
     * when convolving/strething histograms.
     */
    double fmax = Long.MIN_VALUE;

    /**
     * The histogram.
     */
    double[] a = null;


    public static void main(String[] args) {
        double[] a = new double[300];
        double min = 1, max = a.length;
        for (int i = 0; i < a.length; i++) {
            a[i] = (i + 1);
        }
        int nrBins = (int)
                Math.ceil((Math.log(a.length) / Math.log(2)) + 1);    //Sturges formula
        bin(a, min, max, nrBins);
    }

    /**
     * Creates from array <code>a</code> a new array with <code>nrBins</code> bins
     * between <code>min</code> and <code>max</code>.
     *
     * @param a      the original array
     * @param min    minimum value to be included in the new array
     * @param max    maximum value to be included in the new array
     * @param nrBins number of elements in the new array
     * @return
     */
    public static double[] bin(double[] a, double min, double max, int nrBins) {

//		nrBins= (int) 
//			Math.ceil((Math.log(a.length)/ Math.log(2))+ 1);	// Sturges' estimation

        assert (!(Double.isNaN(min) || Double.isNaN(max)));
        double range = max - min; //(Math.ceil(max))- Math.floor(min));

        double[] b = new double[nrBins];
        for (int i = 0; i < b.length; ++i) {
            b[i] = 0;
        }
        double binSize = range / (nrBins - 1);
        for (int i = 0; i < a.length; ++i) {
            if (a[i] < min || a[i] > max) {
                continue;
            }
            // getbin
            //int p= (int) Math.round((b.length- 1)* ((a[i]- min)/ range));
            double pfloat = (a[i] - min) / binSize;
            int p = (int) Math.round(pfloat);

            if (a[i] == 58) {
                System.currentTimeMillis();
            }
            p = (int) ((nrBins - 1) * ((a[i] - min) / (double) (max - min)));
            ++b[p];
        }

        return b;
    }

    public static double[] bin(double[] a, int nrBins) {
        double min = Integer.MAX_VALUE, max = Integer.MIN_VALUE;
        for (int i = 0; i < a.length; i++) {
            if (a[i] < min) {
                min = a[i];
            }
            if (a[i] > max) {
                max = a[i];
            }
        }
        return bin(a, min, max, nrBins);
    }

    /**
     * Divides the values in <code>a</code> by the values in <code>b</code>,
     * i.e., normalizes <code>a</code> (observed distribution) by <code>b</code>
     * (underlying distribution) to obtain a distribution that multiplied with
     * <code>b</code> results in <code>a</code>.<br>
     * <b>Warning:</b> tt is required that <code>a</code> and <code>b</code> are
     * of common lengths. The method is return-by-value, the result is stored in
     * <code>a</code>
     *
     * @param a a distribution (observed distribution), also contains the result
     * @param b another distribution (underlying distribution)
     */
    public void divide(double[] a, double[] b) {
        assert (a.length == b.length);
        for (int i = 0; i < a.length; i++) {
            a[i] = (b[i] == 0 ? 0 : a[i] / b[i]);
        }
    }

    /**
     * Normalizes distribution <code>a</code> to meet the properties of a
     * probability distribution, i.e. Sum(a[i])= 1.<br>
     * <b>Warning:</b> the method is return-by-value, the result is stored in
     * <code>a</code>.
     *
     * @param a distribution
     */
    public static void normalizeToProbabilityDistribution(double[] a) {
        double sum = 0;
        for (int i = 0; i < a.length; i++) {
            sum += a[i];
        }
        for (int i = 0; i < a.length; i++) {
            a[i] /= sum;
        }
    }

    /**
     * Normalizes distribution <code>a</code> to meet the properties of a
     * probability distribution, i.e. Sum(a[i])= 1.<br>
     * <b>Warning:</b> the method is return-by-value, the result is stored in
     * <code>a</code>.
     *
     * @param a distribution
     */
    public static void normalizeToMaximum1(double[] a) {
        double max = Double.MIN_VALUE;
        for (int i = 0; i < a.length; i++) {
            if (a[i] > max) {
                a[i] = max;
            }
        }
        for (int i = 0; i < a.length; i++) {
            a[i] /= max;
        }
    }

    public static EmpiricalDistribution create(File f, boolean fragFile) throws FileNotFoundException, IOException {
        return create(f, Double.NaN, Double.NaN, -1, fragFile);
    }

    public static EmpiricalDistribution create(File f, int nrBins, boolean fragFile) throws FileNotFoundException, IOException {
        return create(f, Double.NaN, Double.NaN, nrBins, fragFile);
    }

    public static EmpiricalDistribution create(File f, double min, double max, int nrBins, boolean fragFile) throws FileNotFoundException, IOException {
        return create(FileHelper.countLines(f.getAbsolutePath()), new FileInputStream(f), min, max, nrBins, fragFile);
    }

    public static EmpiricalDistribution create(int size, InputStream f, double min, double max, int nrBins, boolean fragFile) throws IOException {
        int total = 0;
        String[] ss;
        //int ll= FileHelper.countLines(f.getAbsolutePath());
        DoubleVector v = new DoubleVector(size, 1);

        BufferedReader buffy = new BufferedReader(new InputStreamReader(f));
        for (String s = null; (s = buffy.readLine()) != null; ++total) {
            //Log.progress(current, size);
            ss = s.split("\\s");
            double x = (fragFile ? Double.parseDouble(ss[1]) - Double.parseDouble(ss[0]) + 1
                    : Double.parseDouble(ss[0]));
            // TODO allow multiplier field
            if ((Double.isNaN(min) || x >= min) && (Double.isNaN(max) || x <= max)) {
                v.add(x);
            }
        }
        buffy.close();

        // get attributes
        double[] a = v.toDoubleArray();
        double sum = 0d, lo = Double.MAX_VALUE, hi = Double.MIN_VALUE, mean = Double.NaN;
        for (int i = 0; i < a.length; i++) {
            if (a[i] < lo) {
                lo = a[i];
            }
            if (a[i] > hi) {
                hi = a[i];
            }
            sum += a[i];
        }
        mean = sum / a.length;

        // perform binning
        if (Double.isNaN(min)) {
            min = lo;
        }
        if (Double.isNaN(max)) {
            max = hi;
        }
        if (nrBins < 1) {
            nrBins = (int) ((max - min) + 1);
        }
        double[] h = bin(a, min, max, nrBins);
        System.currentTimeMillis();
        return new EmpiricalDistribution(h, lo, hi, mean);
    }

    /**
     * Creates an instance with the histogram <code>a</code> and
     * <code>min, max, mean</code> attributes of the underlying distribution
     *
     * @param x a value from the empirical distribution
     */
    public EmpiricalDistribution(double[] a, double min, double max, double mean) {
        this.a = a;
        this.min = min;
        this.max = max;
        this.mean = mean;
        updateHistogram();
        weight = 1d;
    }

    /**
     * Returns the probability for the bin of value <code>x</code>.
     *
     * @param x a value from the empirical distribution
     */
    public double getP(double x) {
        int p = getBin(x);
        if (p < 0 || p >= a.length) {
            return 0;
        }
        return (a[p] / sum);
    }

    /**
     * Returns the frequency for the bin of value <code>x</code>
     * relative to the maximum of the distribution.
     *
     * @param x a value from the empirical distribution
     */
    public double getRelFreq(double x) {
        int p = getBin(x);
        if (p < 0 || p >= a.length) {
            return 0;
        }
        return (a[p] / fmax);
    }

    /**
     * Returns the bin (i.e., index of the stored array <code>a</code>)
     * for value <code>x</code>.
     *
     * @param x the value from the empirical distribution
     * @return
     */
    public int getBin(double x) {
        int p = (int) ((a.length - 1) * ((x - min) / (double) (max - min)));
        return p;
    }

    public double[] getBins() {
        return a;
    }

    /**
     * Returns the mid value of bin at index <code>p</code>.
     *
     * @param i index of the bin
     * @return
     */
    public double getBinMid(double i) {
        double lo = min + (i * (max - min) / a.length);
        double hi = lo + (max - min) / a.length;

        return ((lo + hi) / 2d);
    }

    public double getMin() {
        return min;
    }

    public double getMax() {
        return max;
    }

    public double getP(double x, double mean) {
        double diff = x - mean;
        return getP(getMean() + diff);
    }

    public void normalizeToPrior(EmpiricalDistribution d) {

        assert (d.getMin() >= min && d.getMax() <= max);
        double binSize = (max - min) / a.length;
        double offset = min + binSize / 2d;
        double sumD = 0d;
        double[] b = d.a;
        for (int i = 0; i < a.length; i++) {
            //gelProb[i]= (oriProb[i]== 0? 0: gelProb[i]/ oriProb[i]);
//			int lo= d.getBin(min+ i* binSize), hi= d.getBin(min+ (i+ 1)* binSize);
//			double w= 0d;
//			for (int j = lo; j <= hi; j++) 
//				w+= b[j];
//			sum+= w;
//			a[i]= (w== 0? 0: a[i]/ w);
            //double plo= d.getP(min+ i* binSize), phi= d.getP(min+ (i+ 1)* binSize);
            //double pb= (plo+ phi)/ 2d;
            //a[i]= (pb== 0? 0: (a[i]/ sum)/ pb);
            a[i] = (b[i] == 0 ? 0 : a[i] / b[i]);

            //a[i]/= Math.max(d.getBins()[i]/ d.sum, 1);
            System.currentTimeMillis();
        }
        updateHistogram();
    }


    protected void updateHistogram() {
        sum = 0d;
        fmax = 0l;
        for (int i = 0; i < a.length; i++) {
            sum += a[i];
            if (a[i] > fmax) {
                fmax = a[i];
            }
        }
    }


    @Override
    public double getMean() {
        return mean;
    }

    /**
     * @return
     * @deprecated for DEBUG
     */
    public double getPSum() {
        double mySum = 0d;
        for (int i = 0; i < a.length; i++) {
            double p = a[i] / sum;
            mySum += p;
        }
        return mySum;
    }

    public void writeToFile(String path) throws Exception {
        BufferedWriter writer = new BufferedWriter(new FileWriter(path));
        for (int i = 0; i < a.length; i++) {
            writer.write(Double.toString(a[i] / sum));
            writer.write("\n");
        }
        writer.flush();
        writer.close();
    }

    public double getSum() {
        return sum;
    }

    public double getFmax() {
        return fmax;
    }

}
