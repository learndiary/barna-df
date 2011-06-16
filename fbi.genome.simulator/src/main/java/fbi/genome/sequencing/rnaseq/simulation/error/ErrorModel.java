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

public interface ErrorModel {

    public static final byte TYPE_MUTATION = 1, TYPE_INSERTION = 2, TYPE_DELETION = 3;

    public void setBaseProbability(double p);

    public double getBaseProbability();

    public void apply(byte[] quals);

    public void apply(byte[] quals, int from, int to);

    //public void addChar(char[] sequence, int i, String template, int j, Random rnd);
}
