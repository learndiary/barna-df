/*
 * Copyright (c) 2010, Micha Sammeth
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * The names of its contributors may be not used to endorse or promote
 *       products derived from this software without specific prior written
 *       permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL MICHA SAMMETH BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

package barna.flux.simulator.distributions;

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
