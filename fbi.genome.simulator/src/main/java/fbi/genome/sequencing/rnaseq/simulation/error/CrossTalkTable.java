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

import java.util.Arrays;
import java.util.Random;

public class CrossTalkTable {

    public static final String ERR_ID = "CROSSTALK";
    public static final char[] SYMBOLS = new char[]{'A', 'C', 'G', 'N', 'T'};

    Object table; // 5 x qual x 5
    boolean qualities;


    public CrossTalkTable(boolean withQualities) {
        this.qualities = withQualities;
        if (qualities) {
            double[][][] table = new double[SYMBOLS.length][][];
            for (int i = 0; i < table.length; i++) {
                table[i] = new double[ModelPool.qualLevels[1] - ModelPool.qualLevels[0] + 1][];
                for (int j = 0; j < table[i].length; j++) {
                    table[i][j] = new double[SYMBOLS.length];
                    for (int h = 0; h < table[i][j].length; h++) {
                        table[i][j][h] = 0;
                    }
                }
            }
            this.table = table;
        } else {
            double[][] table = new double[SYMBOLS.length][];
            for (int i = 0; i < table.length; i++) {
                table[i] = new double[SYMBOLS.length];
                for (int h = 0; h < table[i].length; h++) {
                    table[i][h] = 0;
                }
            }
            this.table = table;
        }
    }

    public boolean hasQualities() {
        return qualities;
    }


    /**
     * requires upper-case
     *
     * @param a
     * @param b
     */
    public void add(char a, char b, int qual) {
        int p1 = Arrays.binarySearch(SYMBOLS, a);
        int p2 = Arrays.binarySearch(SYMBOLS, b);
        if (p1 < 0 || p2 < 0) {
            System.currentTimeMillis();
        }
        if (hasQualities()) {
            int q = qual - ModelPool.qualLevels[0];
            ++((double[][][]) table)[p1][q][p2];
        } else {
            ++((double[][]) table)[p1][p2];
        }
//		for (int i = 0; i < table[p1].length; i++) 
//			table[p1][i]= ((table[p1][i]* cases[p1])+ ((i==p2)?1:0))/ (cases[p1]+ 1);

//		ModelPool.incr(
//				table[p1][q],
//				p2,
//				cases[p1][q]
//			);

    }

    @Override
    public String toString() {

        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < SYMBOLS.length; i++) {
            if (SYMBOLS[i] == 'N') {
                continue;
            }
            sb.append(ModelPool.ERR_FIELD_TAG + " " + ERR_ID + " " + SYMBOLS[i] + "\n");
            if (hasQualities()) {
                double[][][] table = (double[][][]) this.table;
                for (int j = 0; j < table[i].length; j++) {
                    sb.append(Integer.toString(ModelPool.qualLevels[0] + j) + "\t");
                    for (int h = 0; h < table[i][j].length; h++) {
                        sb.append(Double.toString(table[i][j][h]));
                        sb.append(" ");
                    }
                    sb.deleteCharAt(sb.length() - 1);
                    sb.append("\n");
                }
            } else {
                double[][] table = (double[][]) this.table;
                for (int h = 0; h < table[i].length; h++) {
                    sb.append(Double.toString(table[i][h]));
                    sb.append(" ");
                }
                sb.deleteCharAt(sb.length() - 1);
                sb.append("\n");
            }
            sb.append("\n");
        }

        return sb.toString();
    }

    public Object getTable() {
        return table;
    }

    Random rndMutator = new Random();

    public char mutate(char c, byte b) {
        int p1 = Arrays.binarySearch(SYMBOLS, Character.toUpperCase(c));

        double[] subst = null;
        if (hasQualities()) {
            int q = 0;
            if (b != Byte.MIN_VALUE) {
                b -= (ModelPool.qualLevels[1] == ModelPool.QLEVEL_ILLUMINA[1]) ? 64 : 33;    // bugfix 20101210
                q = b - ModelPool.qualLevels[0];
            }
            subst = ((double[][][]) table)[p1][q];
        } else {
            subst = ((double[][]) table)[p1];
        }

        double r = rndMutator.nextDouble();
        int p = Arrays.binarySearch(subst, r);
        p = (p < 0) ? -p - 1 : p;

        if (p >= SYMBOLS.length) {
            p = SYMBOLS.length - 1;    // TODO find out why
        }
        char c2 = Character.toLowerCase(SYMBOLS[p]);
        return c2;
    }

    /**
     * interpolates x-talk for qualities that did not incur
     */
    public void fill() {

        // wo quals no interpolation
        if (!hasQualities()) {
            return;
        }

        double[] a = new double[SYMBOLS.length];
        for (int j = 0; j < a.length; a[j++] = 0) {
            ;
        }


        double[][][] table = (double[][][]) this.table;
        for (int i = 0; i < table.length; i++) {
            for (int j = 0; j < table[i].length; j++) {
                for (int k = 0; k < a.length; k++) {
                    a[k] += table[i][j][k];
                }
            }
        }


        for (int i = 0; i < table.length; i++) {
            int first = -1, ultim = -1;
            for (int j = 0; j < table[i].length; j++) {
                int x = 0;
                for (; x < table[i][j].length && table[i][j][x] == 0; x++) {
                    ;
                }
                if (x < table[i][j].length) {
                    first = j;
                    break;
                }
            }
            for (int j = table[i].length - 1; j >= 0; --j) {
                int x = 0;
                for (; x < table[i][j].length && table[i][j][x] == 0; x++) {
                    ;
                }
                if (x < table[i][j].length) {
                    ultim = j;
                    break;
                }
            }

            int last = -1, next = -1;
            for (int j = 0; j < table[i].length; j++) {
                int x = 0;
                for (; x < table[i][j].length && table[i][j][x] == 0; ++x) {
                    ;
                }
                if (x < table[i][j].length) {
                    last = j;
                    continue;
                }
                if (last < 0) {
                    if (first >= 0) {
                        System.arraycopy(table[i][first], 0, table[i][j], 0, table[i][j].length);
                    } else {
                        System.arraycopy(a, 0, table[i][j], 0, table[i][j].length);
                    }
                    continue;
                }
                if (next == table[i].length) {
                    if (ultim >= 0) {
                        System.arraycopy(table[i][ultim], 0, table[i][j], 0, table[i][j].length);
                    } else {
                        System.arraycopy(a, 0, table[i][j], 0, table[i][j].length);
                    }
                    continue;
                }
                for (int k = j + 1; (next < j) && k < table[i].length; k++) {
                    for (int h = 0; h < table[i][k].length; h++) {
                        if (table[i][k][h] > 0) {
                            next = k;
                            break;
                        }
                    }
                }
                if (next < j) {
                    next = table[i].length;
                    if (ultim >= 0) {
                        System.arraycopy(table[i][ultim], 0, table[i][j], 0, table[i][j].length);
                    } else {
                        System.arraycopy(a, 0, table[i][j], 0, table[i][j].length);
                    }
                    continue;
                }
                for (int k = 0; k < table[i][j].length; k++) {
                    table[i][j][k] = table[i][last][k] + table[i][next][k];
                }
            }
        }
    }
}
