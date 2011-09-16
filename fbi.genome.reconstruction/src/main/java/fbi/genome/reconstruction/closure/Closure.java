/*
 * Epos Phylogeny Framework
 * Copyright (C) 2009.  University of Jena
 *
 * This file is part of Epos.
 *
 * Epos is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Epos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with Epos.  If not, see <http://www.gnu.org/licenses/>;.
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
