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

import barna.flux.capacitor.diffexp.math.FishersExactTest;

import java.util.concurrent.Callable;

/**
 * Do fisher exact test in a callable to be used in a separate thread
 *
 * @author Thasso Griebel <thasso.griebel@gmail.com>
 */
class QuantificationEntryComparator implements Callable<DifferentialExpression> {
    private QuantificationEntry source;
    private QuantificationEntry target;
    private double sourceReads;
    private double targetReads;

    public QuantificationEntryComparator(QuantificationEntry source, QuantificationEntry target, double sourceReads, double targetReads) {
        this.source = source;
        this.target = target;
        this.sourceReads = sourceReads;
        this.targetReads = targetReads;
    }

    @Override
    public DifferentialExpression call() throws Exception {
        double p = 1.0;
        double difference = 0;
        double foldChange = 0;
        if(source != null && target != null){
            // fisher test
            p = FishersExactTest.fishersExactTest(
                    (int)Math.round(source.getReadCount()),
                    (int)Math.round(sourceReads - source.getReadCount()),
                    (int)Math.round(target.getReadCount()),
                    (int)Math.round(targetReads - target.getReadCount())
                    )[0];
            difference = target.getRpkm() - source.getRpkm();
            if(target.getRpkm() > 0 && source.getRpkm() == 0){
                foldChange = Double.POSITIVE_INFINITY;
            }else if(target.getRpkm() == 0 && source.getRpkm() > 0){
                foldChange = Double.NEGATIVE_INFINITY;
            }else {
                if(target.getRpkm() == 0 && source.getRpkm() == 0){
                    foldChange = 0;
                }else{
                    foldChange = Math.log(target.getRpkm() / source.getRpkm()) / Math.log(2); // log2
                }
            }
        }else if(source != null){
            difference = -source.getRpkm();
            foldChange = Double.NEGATIVE_INFINITY;
        }else if(target != null){
            difference = target.getRpkm();
            foldChange = Double.POSITIVE_INFINITY;
        }else{
            throw new RuntimeException("No transcript specified");
        }
        return new DifferentialExpression(source, target, p, difference, foldChange);
    }
}
