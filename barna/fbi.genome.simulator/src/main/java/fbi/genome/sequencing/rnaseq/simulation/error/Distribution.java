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

package fbi.genome.sequencing.rnaseq.simulation.error;

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
