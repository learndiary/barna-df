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

package barna.genome.reconstruction.closure;


/**
 * @author micha
 *
 * To change this generated comment edit the template variable "typecomment":
 * Window>Preferences>Java>Templates.
 * To enable and disable the creation of type comments go to
 * Window>Preferences>Java>Code Generation.
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
