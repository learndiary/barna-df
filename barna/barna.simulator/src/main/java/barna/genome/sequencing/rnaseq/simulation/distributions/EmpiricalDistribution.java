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

package barna.genome.sequencing.rnaseq.simulation.distributions;

import barna.commons.log.Log;

import java.io.*;

/**
 * An empirical distribution stored in a histogram.
 *
 * @author micha
 */
public class EmpiricalDistribution extends AbstractDistribution {
    /**
     * Default parser that treats a non empty line as a single double value
     */
    public static final LineParser LINE_PARSER_SINGLE_VALUE = new SingleValueLineParser();

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


    /**
     * Adds the value value to the bins if it is within range.
     *
     * @param value the value
     * @param bins the bins
     * @param min minimum value
     * @param max maximum value
     */
    public static void addToBin(double value, double[] bins, double min, double max){
        assert (!(Double.isNaN(min) || Double.isNaN(max)));
        if(value < min || value > max) return;
        int nrBins = bins.length;
        double binSize = (max - min) / (nrBins - 1);
        double bin = (value-min)/binSize;
        /*
        We use bin with values from >= x to < x+range so we
        have to round here
         */
        int binPosition = (int) Math.round(bin);
        ++bins[binPosition];
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

    /**
     * Creates a new distribution from the given file, using the line parser to
     * resolve the values.
     * <p>
     *     NOTE that this iterates the file twice !!!
     *     First iteration counts lines and computes min/max if they are not given (Double.NaN) before
     *     the second interation creates the distribution
     * </p>
     * @param f the file
     * @param min the min value or Doubel.NaN to force computation of min value
     * @param max the max value or Doubel.NaN to force computation of max value
     * @param nrBins number of bins
     * @param lineParser the line parser
     * @return dist the distribution
     * @throws IOException in case of error while reading the data
     */
    public static EmpiricalDistribution create(File f, double min, double max, int nrBins, LineParser lineParser) throws IOException {
        if(f == null) throw new NullPointerException("NULL file not permitted");
        if(!f.exists()) throw new IllegalArgumentException("Unknown file : " + f.getAbsolutePath());
        if(!f.canRead()) throw new IllegalArgumentException("Unable to read from " + f.getAbsolutePath());
        if(lineParser == null) throw new NullPointerException("NULL line parser not permitted");

        /**
         * To enable IO-15 we have to make sure we also get the values if min or max is not set
         * so we can not use the FileUtils.countLines() here
         */
        DistributionAttributes attributes = getAttributes(f, lineParser);
        if(Double.isNaN(max)) max = attributes.getMax();
        if(Double.isNaN(min)) min = attributes.getMin();
        long lineCount = attributes.getSize();
        if(lineCount == 0){
            throw new IllegalArgumentException("Empirical distribution can not be created from an empty file !");
        }
        Log.debug("Creating empirical distribution, counted lines : "+lineCount + " in "+ f.getAbsolutePath() + " with min : " + min + " and max : " + max);
        return create(new FileInputStream(f),lineCount, min, max, nrBins, lineParser);
    }

    /**
     * Creates a new distribution from the data read from the input stream.
     * <p>
     *     NOTE that ALL values have to be specified. If you want to get the attributes,
     *     use {@link #getAttributes(java.io.Reader, barna.genome.sequencing.rnaseq.simulation.distributions.EmpiricalDistribution.LineParser)}
     *     to get the line count and min/max value from your data stream.
     * </p>
     * @param f the input stream
     * @param size number of elements
     * @param min the min value
     * @param max the max value
     * @param nrBins the number of bins
     * @param lineParser the line parser
     * @return dist the distribution
     * @throws IOException in case of IO error while reading data
     */
    public static EmpiricalDistribution create(InputStream f, long size, double min, double max, int nrBins, LineParser lineParser) throws IOException {
        if(size <= 0){
            throw new IllegalArgumentException("Empirical distribution can not be created from an empty file !");
        }
        if(Double.isNaN(min) || Double.isNaN(max)){
            throw new IllegalArgumentException("When reading from an input stream, you have to specify min and max " +
                    "values, because we can not reopen the stream easily and therefore can not compute the values " +
                    "dynamically!");
        }

        if (nrBins < 1) {
            nrBins = (int) ((max - min) + 1);
        }
        double sum = 0d;
        long count = 0l;
        double[] bins = new double[nrBins];
        double mean;
        BufferedReader buffy = new BufferedReader(new InputStreamReader(f));
        for (String s = null; (s = buffy.readLine()) != null;) {
            double x = lineParser.parse(s);
            if(!Double.isNaN(x) && x >= min && x <= max){
                // add to bin
                addToBin(x, bins, min, max);
                // add to stats
                count++;
                sum += x;
            }
        }
        buffy.close();
        // compute the mean
        mean = sum / count;
        return new EmpiricalDistribution(bins, min, max, mean);
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
     * Sample an empirical distribution from the given normal distribution using
     * the given number of bins. The distance parameter is used to compute the range, starting at
     * {@code -(distance*normal.getSd())} going to {@code (distance*normal.getSd())}.
     * The values are normalized to be in the range from 0 to 1.0 reflecting the probabilities
     * from the normal distribution sampled with the given number of bins. Min value is 0 max value is 1.
     *
     * @param normal the normal distribution
     * @param bins the number of bins
     * @param distance the distance to compute the range
     */
    public EmpiricalDistribution(NormalDistribution normal, int bins, final double distance) {
        a = new double[bins+2];
        min = 0;
        max = 1.0;
        double start= -distance *normal.getSd();
        double end = distance*normal.getSd();
        double s = (end-start)/((double)bins);
        double step = start;
        int c = 0;
        sum = 0;
        while(c < a.length-2){
            double p = normal.pdf(step);
            double bb = p / ( 1.0 / s);
            int index = 1 + (c++);
            a[index] = bb;
            sum += bb;
            step += s;
        }
        updateHistogram();
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


    /*
    IO-17 Add custom line parser interface to remove dependency to fragmenter file layout
     */
    /**
     * Implementations parse a single string to a double value
     */
    public static interface LineParser{
        /**
         * Parse the given line and return the double value or return Double.NaN if
         * the line does not resolve to a double value
         *
         * @param line the line
         * @return value the double value or Double.NaN
         */
        double parse(String line);
    }

    /**
     * Parses a line as single double value
     */
    public static class SingleValueLineParser implements LineParser{
        @Override
        public double parse(String line) {
            if(line.trim().isEmpty()) return Double.NaN;
            return Double.parseDouble(line);
        }
    }

    /**
     * Contains information about a distribution. We use this to wrap
     * min/max and size before we create the distribution from a data source.
     */
    public static class DistributionAttributes{
        private double min = Double.NaN;
        private double max = Double.NaN;
        private long size = 0;

        private DistributionAttributes(double min, double max, long size) {
            this.min = min;
            this.max = max;
            this.size = size;
        }

        public double getMin() {
            return min;
        }

        public double getMax() {
            return max;
        }

        public long getSize() {
            return size;
        }
    }

    /**
     * Reads all content from the given file and uses the line parser to create
     * basic attributes such as min and max value and the number of lines.
     *
     * @param file the file
     * @param parser the line parser
     * @return attributes the basic attributes
     * @throws FileNotFoundException in case the file was not found
     */
    public static DistributionAttributes getAttributes(File file, LineParser parser) throws FileNotFoundException {
        if(file == null) throw new NullPointerException("NULL file not permitted");
        if(parser == null) throw new NullPointerException("NULL line parser not permitted");
        if(!file.exists()) throw new FileNotFoundException("Unknown file : " + file.getAbsolutePath());
        if(!file.canRead()) throw new IllegalArgumentException("Unable to read from " + file.getAbsolutePath());
        return getAttributes(new FileReader(file), parser);
    }

    /**
     * Reads all content from the given reader and uses the line parser to create
     * basic attributes such as min and max value and the number of lines.
     * <p>
     * NOTE that the reader is closed after a call to this method.
     * </p>
     *
     * @param reader the content
     * @param parser the line parser
     * @return attributes basic distribution attributes
     */
    public static DistributionAttributes getAttributes(Reader reader, LineParser parser) {
        /**
         * To enable IO-15 we have to make sure we also get the values if min or max is not set
         * so we can not use the FileUtils.countLines() here
         */
        long lineCount = 0;
        BufferedReader buffy = null;
        double min = Double.MAX_VALUE;
        double max = Double.MIN_VALUE;

        try {
            buffy = new BufferedReader(reader);
            for (String s; (s = buffy.readLine()) != null; ++lineCount) {
                double x = parser.parse(s);
                if(!Double.isNaN(x)){
                    min = Math.min(min, x);
                    max = Math.max(max, x);
                }
            }
        } catch (Exception e) {
            throw new RuntimeException("Error while counting lines for empirical distribution creation", e);
        } finally {
            if (buffy != null) {
                try {buffy.close();} catch (IOException ignore) {
                    Log.debug("Error while closing reader while " +
                            "counting lines for empirical distribution : "
                            + ignore.getMessage());
                }
            }
        }
        return new DistributionAttributes(min, max, lineCount);
    }
}
