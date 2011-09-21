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

package fbi.genome.reconstruction.closure;

import java.io.Serializable;

/**
 * @author micha
 *
 * To change this generated comment edit the template variable "typecomment":
 * Window>Preferences>Java>Templates.
 * To enable and disable the creation of type comments go to
 * Window>Preferences>Java>Code Generation.
 */
 public class Closure implements Serializable{

	int seqNbr;
	Sequence[] seq;
	int maxLong;

	PositionSet[] aligSet;
	int nbrAligSets, oldNbrAligSets;

	int[][] predFrontier, succFrontier;

	int[] topolog;
	int[] gauche1, gauche2, droite1, droite2;
	int[][] pos_;
    private static final long serialVersionUID = -6418424151944216603L;

    /**
	 * Constructor for Closure.
	 */
	public Closure() {
		super();
	}

	public static void main(String[] args) {
	}
}
