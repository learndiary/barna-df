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

import java.io.FileReader;

/**
 * @author micha
 *
 * To change this generated comment edit the template variable "typecomment":
 * Window>Preferences>Java>Templates.
 * To enable and disable the creation of type comments go to
 * Window>Preferences>Java>Code Generation.
 */
public class AliGraphClosure {

	public static void main(String[] args) {
		
			// >PIR:GPUGNI Nonlegume hemoglobin I - Swamp oak
			// >PIR:GGZLB Bacterial hemoglobin - Vitreoscilla sp.
			// > LCK       MLCK  RAT SKELETAL MUSCLE MYOSIN LIGHT CHAIN KINASE
		String[] seqs= {
			"ALTEKQEALLKQSWEVLKQNIPAHSLRLFALIIEAAPESKYVFSFLKDSNEIPENNPKLKAHAAVIFKTICESATELRQKGHAVWDNNTLKRLGSIHLKNKITDPHFEVMKGALLGTIKEAIKENWSDEMGQAWTEAYNQLVATIKAEMKE",
			"MLDQQTINIIKATVPVLKEHGVTITTTFYKNLFAKHPEVRPLFDMGRQESLEQPKALAMTVLAAAQNIENLPAILPAVKKIAVKHCQAGVAAAHYPIVGQELLGAIKEVLGDAATDDILDAWGKAYGVIADVFIQVEADLYAQAVEFGNLKDGVNDIKNHK",
			"FSMNSKEALGGGKFGAVCTCTEKSTGLKLAAKVIKKQTPKDKEMVMLEIEVMNQLNHRNLIQLYAAIETPHEIVLFMEYIEGGELFERIVDEDYHLTEVDTMVFVRQICDGILFMHKMRVLHLDLKPENILCVNTTGHLVKIIDFGLARRYNPNEKLKVNFGTPEFLSPEVVNYDQISDKTDMWSLGVITYMLLSGLSPFLGDDDTETLNNVLSGNWYFDEETFEAVSDEAKDFVSNLIVKEQGARMSAAQCLAHPWLNNL"
		};
		int[] lengths= new int[seqs.length];
		for (int i= 0; i< seqs.length; ++i)
			lengths[i]= seqs[i].length();
		
		newAligGraphClosure(seqs.length, lengths, 0, null);
	}

	public static void computeClosure(Closure clos) {
		
		int[][] Succ, Pred;
		int[] NSucc, NPred, npred;
		int nsucc, ni, nj, s, top, bottom, n0, p, n, i, k, pos_n;
		int x;

		Succ= new int[clos.nbrAligSets+ 2][clos.seqNbr];
		Pred= new int[clos.nbrAligSets+ 2][clos.seqNbr];

		NSucc= new int[clos.nbrAligSets + 2];
		NPred= new int[clos.nbrAligSets + 2];
		npred= new int[clos.nbrAligSets + 2];
		clos.topolog= new int[1];

		/* C A L C U L des Succ[n][x] et NPred[n] */

		for (n= 1; n <= clos.nbrAligSets; n++)
			NPred[n]= 0;

		for (n= 1; n <= clos.nbrAligSets; n++) {
			nsucc= 0;
			for (x= 0; x < clos.seqNbr; x++)
				if (clos.aligSet[n].pos[x] > 0) {
					pos_n= clos.aligSet[n].pos[x];
					for (i= pos_n + 1;
						i <= clos.seq[x].longueur
							&& clos.seq[x].aligSetNbr[i] == 0;
						i++)
						clos.seq[x].predAligSetPos[i]= pos_n;
					if (i <= clos.seq[x].longueur) {
						clos.seq[x].predAligSetPos[i]= pos_n;
						if (clos.aligSet[clos.seq[x].aligSetNbr[i]].nbr > 0) {
							n0= Succ[n][nsucc]= clos.seq[x].aligSetNbr[i];
							clos.aligSet[n0].nbr= -clos.aligSet[n0].nbr;
							nsucc++;
						}
					}
					for (i= pos_n - 1;
						i > 0 && clos.seq[x].aligSetNbr[i] == 0;
						i--)
						clos.seq[x].succAligSetPos[i]= pos_n;
					if (i > 0)
						clos.seq[x].succAligSetPos[i]= pos_n;
				}
			for (p= 0; p < nsucc; p++) {
				n0= Succ[n][p];
				Pred[n0][NPred[n0]]= n;
				NPred[n0]++;

				clos.aligSet[n0].nbr= -clos.aligSet[n0].nbr;
			}
			NSucc[n]= nsucc;
		}

		/* C A L C U L de clos.topolog */

		clos.topolog= new int[clos.nbrAligSets + 2];

		bottom= top= 0;

		for (n= 1; n <= clos.nbrAligSets; n++) {
			npred[n]= NPred[n];
			if (npred[n] == 0) {
				top++;
				clos.topolog[top]= n;
			}
		}

		while (bottom != top) {
			bottom++;
			ni= clos.topolog[bottom];
			for (s= 0; s < NSucc[ni]; s++) {
				nj= Succ[ni][s];
				npred[nj]--;
				if (npred[nj] == 0) {
					top++;
					clos.topolog[top]= nj;
				}
			}
		}

		for (x= 0; x < clos.seqNbr; x++) {
			clos.predFrontier[0][x]= 0;
			clos.succFrontier[clos.nbrAligSets + 1][x]=
				clos.seq[x].longueur + 1;
		}

		for (k= 1; k <= clos.nbrAligSets; k++) {
			n0= clos.topolog[k];
			for (x= 0; x < clos.seqNbr; x++) {
				if (clos.aligSet[n0].pos[x] > 0)
					clos.predFrontier[n0][x]= clos.aligSet[n0].pos[x];
				else
					for (p= 0, clos.predFrontier[n0][x]= 0;
						p < NPred[n0];
						p++) {
						n= Pred[n0][p];
						if (clos.predFrontier[n][x] > clos.predFrontier[n0][x])
							clos.predFrontier[n0][x]= clos.predFrontier[n][x];
					}
			}
		}

		for (k= clos.nbrAligSets; k > 0; k--) {
			n0= clos.topolog[k];
			for (x= 0; x < clos.seqNbr; x++) {
				if (clos.aligSet[n0].pos[x] > 0)
					clos.succFrontier[n0][x]= clos.aligSet[n0].pos[x];
				else
					for (p= 0,
						clos.succFrontier[n0][x]= clos.seq[x].longueur + 1;
						p < NSucc[n0];
						p++) {
						n= Succ[n0][p];
						if (clos.succFrontier[n][x] < clos.succFrontier[n0][x])
							clos.succFrontier[n0][x]= clos.succFrontier[n][x];
					}
			}
		}

	}

	public static void moveAligSet(Closure clos, int n1, int n2) {
		
		int x;
		int k;

		for (x= 0; x < clos.seqNbr; x++) {
			
			k= clos.aligSet[n1].pos[x]= clos.aligSet[n2].pos[x];
			if (k > 0)
				clos.seq[x].aligSetNbr[k]= n1;

			clos.predFrontier[n1][x]= clos.predFrontier[n2][x];
			clos.succFrontier[n1][x]= clos.succFrontier[n2][x];
		}

		clos.aligSet[n1].nbr= clos.aligSet[n2].nbr;
	}

	private static void read_closure(Closure clos, int nbreancr, int[][] ancrages) {
		
		FileReader f;
		int x;
		int i, ind, k, n;
		int[][] Succ, Pred;
		int[] NSucc, NPred, npred;

		for (n= 0; n < nbreancr; n++) {
			clos.nbrAligSets++;
			realloc_closure(clos);

			clos.aligSet[clos.nbrAligSets].nbr= 0;
			for (x= 0; x < clos.seqNbr; x++) {
				ind= clos.aligSet[clos.nbrAligSets].pos[x]= ancrages[n][x];
				if (ind > 0) {
					clos.aligSet[clos.nbrAligSets].nbr++;
					clos.seq[x].aligSetNbr[ind]= clos.nbrAligSets;
				}
			}
		}

		computeClosure(clos);
	}

	private static void init_closure(Closure clos, int nbreancr, int[][] ancrages) {
		
		int x;
		int i;
		int[] longsequ;

		longsequ= new int[clos.seqNbr];

		for (x= 0; x < clos.seqNbr; x++) {
			longsequ[x]= clos.seq[x].longueur;
			for (i= 1; i <= clos.seq[x].longueur; i++)
				clos.seq[x].aligSetNbr[i]=
					clos.seq[x].succAligSetPos[i]=
						clos.seq[x].predAligSetPos[i]= 0;
		}

		clos.nbrAligSets= 0;

		if (nbreancr > 0)
			read_closure(clos, nbreancr, ancrages);

		for (x= 0; x < clos.seqNbr; x++)
			clos.seq[x].longueur= longsequ[x];

	}

	public static boolean path(Closure clos, int x, int i, int y, int j) {
		int n2, k;

		if (x == y)
			return (i <= j);

		n2= clos.seq[y].aligSetNbr[j];

		if (n2 == 0) {
			k= clos.seq[y].predAligSetPos[j];
			if (k > 0)
				n2= clos.seq[y].aligSetNbr[k];
		}

		if (n2 == 0)
			return (false);
		else
			return (i <= clos.predFrontier[n2][x]);
	}

	private static void alloc_closure(Closure clos) {

		long nmax, na;
		int x;

		/* sera re'alloue' */
		clos.predFrontier= new int[clos.maxLong + 2][clos.seqNbr + 1];
		clos.succFrontier= new int[clos.maxLong + 2][clos.seqNbr + 1];

		/* sera re'alloue' */
		clos.aligSet= new PositionSet[clos.maxLong + 2];
		for (na= 0; na <= clos.maxLong + 1; na++) {
			
			clos.aligSet[(int) na]= new PositionSet();
			clos.aligSet[(int) na].pos= new int[clos.seqNbr];
		}
		clos.oldNbrAligSets= clos.maxLong;

		for (x= 0; x < clos.seqNbr; x++) {
			clos.seq[x].aligSetNbr= new int[clos.seq[x].longueur + 2];
			clos.seq[x].predAligSetPos= new int[clos.seq[x].longueur + 2];
			clos.seq[x].succAligSetPos= new int[clos.seq[x].longueur + 2];
		}

		clos.gauche1= new int[clos.seqNbr];
		clos.gauche2= new int[clos.seqNbr];
		clos.droite1= new int[clos.seqNbr];
		clos.droite2= new int[clos.seqNbr];
		clos.pos_= new int[clos.seqNbr][clos.seqNbr];
	}

	private static void free_closure(Closure clos) {
		
		System.gc();
/*		long nmax, na;
		int x;

		liberer(clos.gauche1);
		liberer(clos.gauche2);
		liberer(clos.droite1);
		liberer(clos.droite2);
		liberer_mat((void * *) clos.pos_, clos.seqNbr);

		liberer_mat((void * *) clos.succFrontier, clos.oldNbrAligSets + 2);
		liberer_mat((void * *) clos.predFrontier, clos.oldNbrAligSets + 2);

		for (x= 0; x < clos.seqNbr; x++) {
			liberer(clos.seq[x].aligSetNbr);
			liberer(clos.seq[x].predAligSetPos);
			liberer(clos.seq[x].succAligSetPos);
		}

		for (na= 0; na <= clos.oldNbrAligSets + 1; na++) {
			liberer(clos.aligSet[na].pos);
		}
		liberer(clos.aligSet);
*/		
	}

	private static void realloc_closure(Closure clos) {
		
		int na;

		if (clos.nbrAligSets > clos.oldNbrAligSets) {
			
			clos.predFrontier= new int[clos.nbrAligSets+ 2][clos.seqNbr+ 1];
			clos.succFrontier=new int[clos.nbrAligSets+ 2][clos.seqNbr+ 1];
			PositionSet[] newAligSet= new PositionSet[clos.nbrAligSets+ 2];
			for (int i= 0; i< clos.oldNbrAligSets+ 2; ++i) 	// realloc !
				newAligSet[i]= clos.aligSet[i];
			clos.aligSet= newAligSet;
			for (na= clos.oldNbrAligSets + 2;
				na <= clos.nbrAligSets + 1;
				na++) {
				clos.aligSet[na]= new PositionSet();
				clos.aligSet[na].pos= new int[clos.seqNbr];
			}
			clos.oldNbrAligSets= clos.nbrAligSets;
		}
	}

	public static void print_aligSets(Closure clos, int nseq, int i) {
		
		char nouveau_, terminer;
		int n, ng, nd, nn, k;
		int x, y;

		n= ng= nd= clos.seq[nseq].aligSetNbr[i];

		if (ng == 0) {
			k= clos.seq[nseq].predAligSetPos[i];
			if (k > 0)
				ng= clos.seq[nseq].aligSetNbr[k];
			k= clos.seq[nseq].succAligSetPos[i];
			if (k > 0)
				nd= clos.seq[nseq].aligSetNbr[k];
		}

		System.out.print("echelle "+n+": ");
		if (n != 0)
			for (x= 0; x < clos.seqNbr; x++)
				System.out.print(clos.aligSet[n].pos[x]+ " ");

		System.out.print("\nfrontiere clos.gauche "+ng+": ");
		if (ng != 0)
			for (x= 0; x < clos.seqNbr; x++)
				System.out.print(clos.predFrontier[ng][x]+ " ");

		System.out.print("\nfrontiere clos.droite "+nd+": ");
		if (nd != 0)
			for (x= 0; x < clos.seqNbr; x++)
				System.out.print(clos.succFrontier[nd][x]+ " ");

		System.out.print("\n");

	}

	private static void init_seq(Closure clos, int nbreseq, int[] longseq) {
		
		int x;

		clos.seqNbr= nbreseq;

		clos.seq= new Sequence[clos.seqNbr];

		for (x= clos.maxLong= 0; x < clos.seqNbr; x++) {
			
			clos.seq[x]= new Sequence();
			clos.seq[x].longueur= longseq[x];
			if (clos.maxLong < longseq[x])
				clos.maxLong= longseq[x];
		}
	}

	private static void desinit_seq(Closure clos) {
		
		int x;

		clos.seq= null;
	}

	/*********************************************************/
	/************** EXTERN FONCTIONS *************************/
	/*********************************************************/

	public static Closure newAligGraphClosure(
			int nbreseq,
			int[] longseq,
			int nbreancr,
			int[][] ancrages) {

		Closure clos= new Closure();

		init_seq(clos, nbreseq, longseq);
		
		alloc_closure(clos); /* utilise clos.maxLong */

		init_closure(clos, nbreancr, ancrages);

		return clos;
	}

	private static void freeAligGraphClosure(Closure clos) {
		
		free_closure(clos);

		desinit_seq(clos);

		clos= null;
	}

	public static void addAlignedPositions(
		Closure clos,
		int seq1,
		int i,
		int seq2,
		int j) {

		char nouveau_, terminer;
		int n, n1, n2, ng1, ng2, nd1, nd2, nn, k;
		int x, y;

		n1= ng1= nd1= clos.seq[seq1].aligSetNbr[i];
		n2= ng2= nd2= clos.seq[seq2].aligSetNbr[j];

		if (n1 == 0 || n2 == 0 || n1 != n2) {

			if (ng1 == 0) {

				k= clos.seq[seq1].predAligSetPos[i];
				if (k > 0)
					ng1= clos.seq[seq1].aligSetNbr[k];
				k= clos.seq[seq1].succAligSetPos[i];
				if (k > 0)
					nd1= clos.seq[seq1].aligSetNbr[k];
			}

			if (ng2 == 0) {

				k= clos.seq[seq2].predAligSetPos[j];
				if (k > 0)
					ng2= clos.seq[seq2].aligSetNbr[k];
				k= clos.seq[seq2].succAligSetPos[j];
				if (k > 0)
					nd2= clos.seq[seq2].aligSetNbr[k];
			}

			if (ng1 == 0)
				for (x= 0; x < clos.seqNbr; x++)
					clos.gauche1[x]= 0;
			else
				for (x= 0; x < clos.seqNbr; x++)
					clos.gauche1[x]= clos.predFrontier[ng1][x];
			if (nd1 == 0)
				for (x= 0; x < clos.seqNbr; x++)
					clos.droite1[x]= clos.seq[x].longueur + 1;
			else
				for (x= 0; x < clos.seqNbr; x++)
					clos.droite1[x]= clos.succFrontier[nd1][x];
			if (ng2 == 0)
				for (x= 0; x < clos.seqNbr; x++)
					clos.gauche2[x]= 0;
			else
				for (x= 0; x < clos.seqNbr; x++)
					clos.gauche2[x]= clos.predFrontier[ng2][x];
			if (nd2 == 0)
				for (x= 0; x < clos.seqNbr; x++)
					clos.droite2[x]= clos.seq[x].longueur + 1;
			else
				for (x= 0; x < clos.seqNbr; x++)
					clos.droite2[x]= clos.succFrontier[nd2][x];

			clos.gauche1[seq1]= clos.droite1[seq1]= i;
			clos.gauche2[seq2]= clos.droite2[seq2]= j;

			nn= clos.nbrAligSets + 1;

			for (x= 0; x < clos.seqNbr; x++) {
				clos.aligSet[nn].pos[x]= 0;
				if (n1 > 0 && clos.aligSet[n1].pos[x] > 0)
					clos.aligSet[nn].pos[x]= clos.aligSet[n1].pos[x];
				else {
					if (n2 > 0 && clos.aligSet[n2].pos[x] > 0)
						clos.aligSet[nn].pos[x]= clos.aligSet[n2].pos[x];
				}

				if (clos.aligSet[nn].pos[x] == 0) {
					clos.predFrontier[nn][x]=
						Math.max(clos.gauche1[x], clos.gauche2[x]);
					clos.succFrontier[nn][x]=
						Math.min(clos.droite1[x], clos.droite2[x]);
				} else
					clos.predFrontier[nn][x]=
						clos.succFrontier[nn][x]= clos.aligSet[nn].pos[x];
			}
			clos.predFrontier[nn][seq1]=
				clos.succFrontier[nn][seq1]= clos.aligSet[nn].pos[seq1]= i;
			clos.predFrontier[nn][seq2]=
				clos.succFrontier[nn][seq2]= clos.aligSet[nn].pos[seq2]= j;

			for (x= clos.aligSet[nn].nbr= 0; x < clos.seqNbr; x++)
				if (clos.aligSet[nn].pos[x] > 0) {
					k= clos.aligSet[nn].pos[x];
					clos.seq[x].aligSetNbr[k]= nn;
					clos.aligSet[nn].nbr++;
				}

			for (x= 0; x < clos.seqNbr; x++)
				if (clos.droite1[x] != clos.droite2[x])
					/* => la front. clos.gauche peut changer */
					for (y= 0; y < clos.seqNbr; y++) {
						clos.pos_[x][y]= 0;
						k= clos.succFrontier[nn][x];
						if (k == clos.aligSet[nn].pos[x])
							k= clos.seq[x].succAligSetPos[k];
						if (k <= clos.seq[x].longueur)
							while (k > 0) {
								n= clos.seq[x].aligSetNbr[k];
								if (clos.predFrontier[n][y]
									< clos.predFrontier[nn][y]) {
									clos.pos_[x][y]= k;
									k= clos.seq[x].succAligSetPos[k];
								} else
									k= 0;
							}
					}

			for (x= 0; x < clos.seqNbr; x++)
				if (clos.droite1[x] != clos.droite2[x])
					/* => la front. gauche peut changer */
					for (y= 0; y < clos.seqNbr; y++) {
						k= clos.succFrontier[nn][x];
						if (k == clos.aligSet[nn].pos[x])
							k= clos.seq[x].succAligSetPos[k];
						if (clos.pos_[x][y] > 0)
							while (k > 0 && k <= clos.pos_[x][y]) {
								n= clos.seq[x].aligSetNbr[k];
								clos.predFrontier[n][y]=
									clos.predFrontier[nn][y];
								k= clos.seq[x].succAligSetPos[k];
							}
					}

			for (x= 0; x < clos.seqNbr; x++)
				if (clos.gauche1[x] != clos.gauche2[x])
					/* => la front. droite peut changer */
					for (y= 0; y < clos.seqNbr; y++) {
						clos.pos_[x][y]= 0;
						k= clos.predFrontier[nn][x];
						if (k > 0 && k == clos.aligSet[nn].pos[x])
							k= clos.seq[x].predAligSetPos[k];
						while (k > 0) {
							n= clos.seq[x].aligSetNbr[k];
							if (clos.succFrontier[n][y]
								> clos.succFrontier[nn][y]) {
								clos.pos_[x][y]= k;
								k= clos.seq[x].predAligSetPos[k];
							} else
								k= 0;
						}
					}

			for (x= 0; x < clos.seqNbr; x++)
				if (clos.gauche1[x] != clos.gauche2[x])
					/* => la front. clos.droite peut changer */
					for (y= 0; y < clos.seqNbr; y++) {
						k= clos.predFrontier[nn][x];
						if (k > 0 && k == clos.aligSet[nn].pos[x])
							k= clos.seq[x].predAligSetPos[k];
						if (clos.pos_[x][y] > 0)
							while (k >= clos.pos_[x][y]) {
								n= clos.seq[x].aligSetNbr[k];
								clos.succFrontier[n][y]=
									clos.succFrontier[nn][y];
								k= clos.seq[x].predAligSetPos[k];
							}
					}

			if (n1 == 0) {
				for (k= i - 1; k > 0 && clos.seq[seq1].aligSetNbr[k] == 0; k--)
					clos.seq[seq1].succAligSetPos[k]= i;
				if (k > 0)
					clos.seq[seq1].succAligSetPos[k]= i;
				for (k= i + 1;
					k <= clos.seq[seq1].longueur
						&& clos.seq[seq1].aligSetNbr[k] == 0;
					k++)
					clos.seq[seq1].predAligSetPos[k]= i;
				if (k <= clos.seq[seq1].longueur)
					clos.seq[seq1].predAligSetPos[k]= i;
			}

			if (n2 == 0) {
				for (k= j - 1; k > 0 && clos.seq[seq2].aligSetNbr[k] == 0; k--)
					clos.seq[seq2].succAligSetPos[k]= j;
				if (k > 0)
					clos.seq[seq2].succAligSetPos[k]= j;
				for (k= j + 1;
					k <= clos.seq[seq2].longueur
						&& clos.seq[seq2].aligSetNbr[k] == 0;
					k++)
					clos.seq[seq2].predAligSetPos[k]= j;
				if (k <= clos.seq[seq2].longueur)
					clos.seq[seq2].predAligSetPos[k]= j;
			}

			if (n1 > n2) {
				n= n1;
				n1= n2;
				n2= n;
			}

			if (n2 == 0) {
				clos.nbrAligSets++;

				realloc_closure(clos);
			} else {
				if (n1 == 0) {
					moveAligSet(clos, n2, nn);
				} else {
					moveAligSet(clos, n1, nn);

					if (n2 < clos.nbrAligSets)
						moveAligSet(clos, n2, clos.nbrAligSets);
					clos.nbrAligSets--;

					realloc_closure(clos);
				}
			}
		}
	}

	// return was int
	public static void addAlignedSegments(
		Closure clos,
		int x,
		int i,
		int y,
		int j,
		int l) {

		int k;

		for (k= 0; k < l; i++, j++, k++)
			addAlignedPositions(clos, x, i, y, j);

	}

	public static boolean alignablePositions(Closure clos, int x, int i, int y, int j) {

		if (path(clos, x, i, y, j))
			return (path(clos, y, j, x, i));
		else
			return (!path(clos, y, j, x, i));
	}
	public static boolean alignableSegments(
		Closure clos,
		int x,
		int i,
		int y,
		int j,
		int l) {
		int k;

		for (k= 0;
			(k < l) && alignablePositions(clos, x, i, y, j);
			i++, j++, k++);

		return (k == l);
	}

	public static boolean alignedPositions(Closure clos, int x, int i, int y, int j) {

		return (x == y && i == j)
			|| (clos.seq[x].aligSetNbr[i] != 0
				&& clos.seq[x].aligSetNbr[i] == clos.seq[y].aligSetNbr[j]);
	}

	public static boolean alignedSegments(
		Closure clos,
		int x,
		int i,
		int y,
		int j,
		int l) {
		int k;

		for (k= 0; k < l && alignedPositions(clos, x, i, y, j); i++, j++, k++);

		return (k == l);
	}

	public static int predFrontier(
		Closure clos,
		int x,
		int i,
		int y) /* on suppose que x!=y */ {
		int n, k;

		n= clos.seq[x].aligSetNbr[i];

		if (n == 0) {
			k= clos.seq[x].predAligSetPos[i];
			if (k > 0)
				n= clos.seq[x].aligSetNbr[k];
		}

		if (n > 0)
			return (clos.predFrontier[n][y]);
		else
			return (0);
	}
	
	/**
	 * converts:
	 * 		1) the 0-based positions of the params to 1-based
	 * 		2) the 1-based positions of return to 0-based
	 * 		3) the 'transitivity exclusions' to 'boundary inclusives', eg
	 * matches: ub= lb, unaligned areas= lb< ub both included, and
	 * insertions: lb= ub+1 (x, x+1) is converted to (x+1,x)
	 * @param
	 * 		Closure clos	consistency structure for sequences
	 * 		int x			seqNo1 (0-based)
	 * 		int i			pos in seq#1 (0-based)
	 * 		int y			seqNo2 (0-based)
	 * @deprecated don't use, converts e.g. [0,2] to [0,0]
	 */
	public static int getLowerBound(
		Closure clos,
		int x,
		int i,
		int y) {
			
			// get boundaries (ask for 1-based position)
		int lb= predFrontier(clos, x, (i+1), y);
		int ub= succFrontier(clos, x, (i+1), y);
		
			// convert transitivity to boundaries
		if (ub != lb) 
			++lb;
		
			// convert to 0-based return
		--lb;

		return lb;
	}
	
	/**
	 * converts:
	 * 		1) the 0-based positions of the params to 1-based
	 * 		2) the 1-based positions of return to 0-based
	 * 		3) the 'transitivity exclusions' to 'boundary inclusives'
	 * @param
	 * 		Closure clos	consistency structure for sequences
	 * 		int x			seqNo1 (0-based)
	 * 		int i			pos in seq#1 (0-based)
	 * 		int y			seqNo2 (0-based)
	 * @deprecated don't use, converts e.g. [0,2] to [0,0]
	 */
	public static int getUpperBound(
		Closure clos,
		int x,
		int i,
		int y) {
			
			// get boundaries (ask for 1-based position)
		int lb= predFrontier(clos, x, (i+1), y);
		int ub= succFrontier(clos, x, (i+1), y);
		
			// convert transitivity to boundaries
		if (ub != lb) 
			--ub;
			
			// convert to 0-based return
		--ub;

		return ub;
	}
	
	
	
	/**
	 * converts:
	 * 		1) the 0-based positions of the params to 1-based
	 * 		2) the 1-based positions of return to 0-based
	 * 		3) the 'transitivity exclusions' to 'boundary inclusives'
	 * @param
	 * 		Closure clos	consistency structure for sequences
	 * 		int x			seqNo1 (0-based)
	 * 		int i			pos in seq#1 (0-based)
	 * 		int y			seqNo2 (0-based)
	 */
	public static int[] getBoundaries(
		Closure clos,
		int x,
		int i,
		int y) {

			// get boundaries (ask for 1-based position)
		int[] bounds= new int[2];			
		bounds[0]= predFrontier(clos, x, (i+1), y);
		bounds[1]= succFrontier(clos, x, (i+1), y);
		
			// convert transitivity to boundaries
		if (bounds[0] != bounds[1]) {
			++bounds[0];
			--bounds[1];	// excluded ??!
		}
		
	
			// convert to 0-based return
		--bounds[0];
		--bounds[1];
				
		return bounds;
	}
	
	/**
	 * WARNING: return positions are 1-based !!!
	 * (sequence numbers are 0-based),
	 * (parameter-pos corrected)
	 */
	public static int[] getTransitivityBounds(
		Closure clos,
		int x,
		int i,
		int y) {

			// get boundaries (ask for 1-based position)
		int[] bounds= new int[2];			
		bounds[0]= predFrontier(clos, x, (i+1), y);
		bounds[1]= succFrontier(clos, x, (i+1), y);

		return bounds;
	}			
	

	public static int succFrontier(
		Closure clos,
		int x,
		int i,
		int y) /* on suppose que x!=y */ {
		int n, k;

		n= clos.seq[x].aligSetNbr[i];

		if (n == 0) {
			k= clos.seq[x].succAligSetPos[i];
			if (k > 0)
				n= clos.seq[x].aligSetNbr[k];
		}

		if (n > 0)
			return (clos.succFrontier[n][y]);
		else
			return (clos.seq[y].longueur + 1);
	}

}
