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

import fbi.genome.model.constants.Constants;

public class CopyOfSimpleMatrix implements Matrix {

	int[][] m;
	int s;
	
	public CopyOfSimpleMatrix(int len, byte dir) {
		if (dir== Constants.DIR_BOTH)
			m= new int[2][];
		else
			m= new int[1][];
		for (int i = 0; i < m.length; i++) {
			m[i]= new int[len];
			for (int j = 0; j < m.length; j++) 
				m[i][j]= 0;
		}
	}
	
	/**
	 * adds a breakpoint (len ignored)
	 */
	public void add(int p, int len, byte dir) {
		add(p, dir);
	}	
	
	
	public void add(int p, byte dir) {
		if (p< 0|| p>= getLength()) {
			System.err.println("Matrix ArrayIndexOutOfBounds "+ p+" in "+ getLength());
			if (p< 0)
				p= 0;
			else
				p= getLength()- 1;
		}
		
		if (dir== Constants.DIR_BACKWARD)
			++m[1][p];
		else
			++m[0][p];
		++s;
	}

	public void add(int p1, int p2, int len) {
		; // :)
	}

	public int getSum() {
		return s;
	}

	public int getLength() {
		return m[0].length;
	}
	
	public void merge(Matrix nn, int readLen) {
		
		CopyOfSimpleMatrix n= (CopyOfSimpleMatrix) nn;
		
		double fac= getLength()/ (double) n.getLength();
		for (int x = 0; x < m.length; x++) {
			for (int i = 0; i < n.getLength(); i++) {
//				int p= (int) Math.round(fac* i);
//				if (p+ readLen- 1>= m[x].length)
//					continue; // falls outside
//				p= (p== m[x].length?p-1:p);
				
				// make blurry, dont respect readlength..
				// ..delegate counted ones to TSuperProfile
				double floatPos= fac* i;
				int minPos= (int) Math.floor(floatPos), maxPos= (int) Math.ceil(floatPos);
				try {assert(floatPos> 0&& minPos> 0&& maxPos> 0);} catch (AssertionError err) {
					if (Constants.verboseLevel> Constants.VERBOSE_NORMAL)
						System.err.println("[ASSERTION] "+getClass().getName()+".merge() matrix index negative:\n\t"
							+"length1= "+ getLength()
							+", length2= "+ n.getLength()
							+", floatPos= "+ floatPos
							+", min/max= "+ minPos+"/"+maxPos
						);
					continue;
				}
				for (int p = minPos; p <= maxPos&& p< m[x].length; p++) {
					m[x][p]+= n.m[x][i];
					s+= n.m[x][i];
				}
			}
		}
	}
	
	public void fill(int[] insertSize, int readLen) {
		for (int i = 0; i < m.length; i++) {
			for (int j = 0; j < m[i].length; j++) {
				++m[i][j];
				++s;
			}
		}
	}

	public int get(int p1, int p2, int readLen, int[] insertMinMax, byte dir) {
		int count= 0;
		try {assert(p1>= 0&& p2<= m[0].length);} catch (AssertionError err){
			if (Constants.verboseLevel> Constants.VERBOSE_NORMAL)
				System.err.println("[ASSERTION] "+getClass().getName()+".get():\n\t" +
						"index pair ("+p1+","+p2+") out of matrix bounds (length "+m[0].length+")!");
			return 0;
		};
		// 090820 < p2 
		// exclusive regions for back-normalization needed
		// otherwise fracs> transcriptcount for gene (and reads also)
		for (int i = p1; i < p2; i++) {
			if (dir== Constants.DIR_FORWARD|| dir== Constants.DIR_BOTH)
				count+= m[0][i];
			if (dir== Constants.DIR_BACKWARD|| dir== Constants.DIR_BOTH)
				count+= m[1][i];
		}
		return count;
	}

	public int get(int p1, int p2, int p3, int p4, int readLen) {
		return 0;	// :)
	}

	private StringBuilder toStringBuilder() {
		StringBuilder sb= new StringBuilder();
		for (int i = 0; i < m[0].length; i++) {
			for (int j = 0; j < m.length; j++) {
				sb.append(Integer.toString(m[j][i]));
				sb.append("\t");
			}
			sb.deleteCharAt(sb.length()- 1);
			sb.append("\n");
		}
		return sb;
	}
	
	@Override
	public String toString() {
		return toStringBuilder().toString();
	}
	
	public byte[] toByteArray() {
		StringBuilder sb= toStringBuilder();
		
		byte[] b= new byte[sb.length()];
		for (int i = 0; i < b.length; i++) 
			b[i]= (byte) sb.charAt(i);
		
		return b;
	}

	public int project(int[][] b) {
		// TODO Auto-generated method stub
		return 0;
	}
}
