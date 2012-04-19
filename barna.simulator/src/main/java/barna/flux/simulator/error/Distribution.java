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

package barna.flux.simulator.error;

/**
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
abstract class Distribution {
    /**
     * Number of quality values
     */
    protected int size = 0;
    /**
     * The values
     */
    protected int[] values;

    /**
     * Number of reads
     */
    protected int reads;

    public Distribution(int size) {
        this.size = size;
        this.values = new int[size];
    }

    public abstract void addRead(Read read);

    public double getValue(int position) {
        return values[position] / (double) reads;
    }

    public String toString() {
        StringBuffer b = new StringBuffer();
        for (int i = 0; i < size; i++) {
            b.append(i).append("\t").append(getValue(i)).append("\n");
        }
        return b.toString();
    }

    public double[][] getDistribution() {
        double[][] d = new double[2][size];
        for (int i = 0; i < size; i++) {
            d[0][i] = i;
            d[1][i] = getValue(i);
        }
        return d;
    }


}
