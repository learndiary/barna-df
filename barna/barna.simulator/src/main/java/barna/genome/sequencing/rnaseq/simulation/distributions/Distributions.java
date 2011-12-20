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

/**
 *
 * Helper class to parse strings to distributions
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class Distributions {

    private static final char FILTER_DISTRIBUTION_NORMAL = 'N';
    private static final char FILTER_DISTRIBUTION_UNIFORM = 'U';
    private static final char FILTER_DISTRIBUTION_WEIBULL = 'W';

    /**
     * Parses a string to a distribution. Currently supports
     *
     * <pre>
     * Normal Distribution as : N(mean, sd) or N(mean)
     * </pre>
     *
     * @param s
     * @return
     */
    public static AbstractDistribution parseDistribution(String s) {

        // find parameter block delimiters
        int p = s.indexOf('('), q = s.indexOf(')');
        assert (p > 0 && q > p);

        // parse parameters
        String[] ss = s.substring(p + 1, q).split(",");
        double[] par = new double[ss.length];
        for (int i = 0; i < par.length; i++) {
            par[i] = Double.parseDouble(ss[i]);
        }

        // nature of function
        AbstractDistribution d = null;
        if (s.charAt(p - 1) == FILTER_DISTRIBUTION_NORMAL) {
            d = (par.length == 1 ? new NormalDistribution(par[0]) : new NormalDistribution(par[0], par[1]));
        } else if (s.charAt(p - 1) == FILTER_DISTRIBUTION_UNIFORM) {
            ; // TODO
            throw new UnsupportedOperationException("Uniform distribution is currently not supported!");
        } else if (s.charAt(p - 1) == FILTER_DISTRIBUTION_WEIBULL) {
            ; // TODO
            throw new UnsupportedOperationException("Weibull distribution is currently not supported!");
        }

        // weight (in sums of functions)
        if (p > 1) {
            double f = Double.parseDouble(s.substring(0, p - 1));
            d.setWeight(f);
        }

        return d;
    }

}
