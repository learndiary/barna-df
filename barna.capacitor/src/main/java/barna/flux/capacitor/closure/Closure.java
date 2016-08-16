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

package barna.flux.capacitor.closure;


/**
 * @author micha
 */
 public class Closure {

	int seqNbr;
	Sequence[] seq;
	int maxLong;

	PositionSet[] aligSet;
	int nbrAligSets, oldNbrAligSets;

	int[][] predFrontier, succFrontier;

	int[] topolog;
	int[] gauche1, gauche2, droite1, droite2;
	int[][] pos_;
//    private static final long serialVersionUID = -6418424151944216603L;

    /**
	 * Constructor for Closure.
	 */
	public Closure() {
		super();
	}

	public Closure(Closure aClosure) {
		
		seqNbr= aClosure.seqNbr;
		maxLong= aClosure.maxLong;
		nbrAligSets= aClosure.nbrAligSets;
		oldNbrAligSets= aClosure.oldNbrAligSets;
		
		if (aClosure.seq!= null) {
			seq= new Sequence[aClosure.seq.length];
			for (int i = 0; i < seq.length; i++) 
				seq[i]= new Sequence(aClosure.seq[i]);
		}
		if (aClosure.aligSet!= null) {
			aligSet= new PositionSet[aClosure.aligSet.length];
			for (int i = 0; i < aligSet.length; i++) {
				aligSet[i]= new PositionSet(aClosure.aligSet[i]);
			}
		}
		
		topolog= clone(aClosure.topolog);
		gauche1= clone(aClosure.gauche1);
		gauche2= clone(aClosure.gauche2);
		droite1= clone(aClosure.droite1);
		droite2= clone(aClosure.droite2);
		
		predFrontier= clone(aClosure.predFrontier);
		succFrontier= clone(aClosure.succFrontier);
		pos_= clone(aClosure.pos_);
	}
	
	public static int[] clone(int[] a) {
		if (a== null)
			return null;
		int[] b= new int[a.length];
		System.arraycopy(a, 0, b, 0, a.length);
		return b;
	}
	
	public static int[][] clone(int[][] a) {
		if (a== null)
			return null;
		int[][] b= new int[a.length][];
		for (int i = 0; i < b.length; i++) {
			b[i]= new int[a[i].length];
			System.arraycopy(a[i], 0, b[i], 0, a[i].length);
		}
		return b;
	}
}
