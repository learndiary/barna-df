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

package fbi.genome.sequencing.rnaseq.reconstruction;

public interface Matrix {

	public void add(int p, int len, byte dir);
	public void add(int p, byte dir);
	public int get(int p1, int p2, int readLen, int[] insertMinMax, byte dir);
	public void add(int p1, int p2, int len);
	public int get(int p1, int p2, int p3, int p4, int readLen);
	public int getSum();
	public int getLength();
	public void merge(Matrix n, int readLen);
	public void fill(int[] insertSize, int readLen);
	public byte[] toByteArray();
	public int project(int[][] b);
}
