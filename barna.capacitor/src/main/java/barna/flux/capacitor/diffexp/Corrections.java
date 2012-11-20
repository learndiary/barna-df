/*
 * Copyright (c) 2012, Micha Sammeth, Thasso Griebel, Emilio Palumbo
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *      * Redistributions of source code must retain the above copyright
 *        notice, this list of conditions and the following disclaimer.
 *      * Redistributions in binary form must reproduce the above copyright
 *        notice, this list of conditions and the following disclaimer in the
 *        documentation and/or other materials provided with the distribution.
 *      * The names of its contributors may be not used to endorse or promote
 *        products derived from this software without specific prior written
 *        permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 *  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 *  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 *  DISCLAIMED. IN NO EVENT SHALL MICHA SAMMETH BE LIABLE FOR ANY
 *  DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 *  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 *  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 *  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 *  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 *  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

package barna.flux.capacitor.diffexp;

import barna.commons.log.Log;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;

/**
 * Computes corrected p-values
 *
 * @author Thasso Griebel <thasso.griebel@gmail.com>
 */
class Corrections {
    /**
     * Takes the list of expressions and sets the FDR-corrected p-value
     * using Benjamini-Hochberg Correction. This expression list is sorted
     * by p-value and ranked. The corrected value for entry i is then computed as:
     *
     * <pre>
     *
     *   pc_i = min( p_(i+1), (m / i) * p_i)
     * </pre>
     *
     * where m is the number of given expressions.
     *
     * @param expressions the differential expressions
     */
    public static void fdr(List<DifferentialExpression> expressions){
        // sort the input
        Log.progressStart("FDR Correction");
        Collections.sort(expressions, new Comparator<DifferentialExpression>() {
            @Override
            public int compare(DifferentialExpression a, DifferentialExpression b) {
                int c = Double.compare(a.getP(), b.getP());
                if(c == 0){
                    return Double.compare(a.getFoldChange(), b.getFoldChange());
                }
                return c;
            }
        });

        int m = expressions.size();
        double min = 1.0;
        for (int i = m; i > 0; i--) {
            Log.progress(i-m, m);
            double pc = (m / (double)i) * expressions.get(i-1).getP();
            min = Math.min(min, pc);
            expressions.get(i-1).setFdrP(min);
        }
        Log.progressFinish("Done", true);
    }

    /**
     * Comput Bonferroni correction for P values
     *
     * @param expressions the expressions
     * @param n number of samples
     */
    public static void bonferroni(List<DifferentialExpression> expressions, int n){
        // sort the input
        Log.progressStart("Bonferroni Correction");
        for (DifferentialExpression expression : expressions) {
            double pp = expression.getP() * (double) n;
            pp = Math.min(1.0, pp);
            expression.setBonferroniP(pp);
        }
        Log.progressFinish("Done", true);
    }
}
