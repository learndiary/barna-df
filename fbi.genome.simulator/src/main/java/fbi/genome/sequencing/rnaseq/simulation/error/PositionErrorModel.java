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


import fbi.commons.ByteArrayCharSequence;

import java.util.Arrays;
import java.util.Random;

public class PositionErrorModel implements ErrorModel {

    int start = -1, extension = -1;
    double baseProb;
    double[][] quals;
    boolean qualities = false;

    public PositionErrorModel(boolean withQualities, int start, int extension) {
        this.qualities = withQualities;
        this.start = start;
        this.extension = extension;
        if (withQualities) {
            quals = new double[extension][];
            for (int i = 0; i < quals.length; i++) {
                quals[i] = new double[ModelPool.qualLevels[1] - ModelPool.qualLevels[0] + 1];
                for (int j = 0; j < quals[i].length; j++) {
                    quals[i][j] = 0d;
                }
            }
        }
        baseProb = 0l;
    }

    public int getExtension() {
        return extension;
    }

    public double getBaseProbability() {
        return baseProb;
    }

    public boolean add(int start, int extension, ByteArrayCharSequence cs) {
        if (start != this.start || extension != this.extension) {
            return false;
        }

        for (int i = 0; qualities && i < extension; i++) {
//			ModelPool.incr(quals[i], 
//					cs.byteAt(start+ i) 
//					- (ModelPool.qualLevels== ModelPool.QLEVEL_ILLUMINA?64:33)
//					- ModelPool.qualLevels[0],
//				 cases);	// +1-1
            ++quals[i][cs.byteAt(start + i)
                    - (ModelPool.qualLevels[1] == ModelPool.QLEVEL_ILLUMINA[1] ? 64 : 33)
                    - ModelPool.qualLevels[0]];
        }

        ++baseProb;
        return true;
    }

    public void apply(char[] seq) {

    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder(ModelPool.ERR_FIELD_TAG + " " + ERR_ID + " " +
                Integer.toString(start) + " " + Integer.toString(extension) + " " + Double.toString(baseProb) + "\n");
        if (qualities) {
            for (int i = 0; i < quals.length; i++) {
                sb.append(Integer.toString(start + i));
                sb.append("\t");
                for (int j = 0; j < quals[i].length; j++) {
                    sb.append(Double.toString(quals[i][j]));
                    sb.append("\t");
                }
                sb.deleteCharAt(sb.length() - 1);
                sb.append("\n");
            }
        }

        return sb.toString();
    }

    public void setBaseProbability(double baseProb) {
        this.baseProb = baseProb;
    }

    public int getStart() {
        return start;
    }

    public double[][] getQuals() {
        return quals;
    }

    Random rndQualChooser = new Random();
    public static final String ERR_ID = "POSITIONERRORPROFILE";

    public void apply(byte[] quals) {
        for (int i = start; i < start + this.quals.length; i++) {
            double r = rndQualChooser.nextDouble();
            int p = Arrays.binarySearch(this.quals[i - start], r);
            p = (p < 0) ? -p - 1 : p;
//			if (p< 0)
//				System.currentTimeMillis();
            int q = (int) p + ModelPool.qualLevels[0];
            quals[i] = (byte) ((q < quals[i]) ? q : quals[i]);
        }
    }

    public void apply(byte[] quals, int from, int to) {
        for (int i = start; i < start + extension; i++) {
            double r = rndQualChooser.nextDouble();
            int p = Arrays.binarySearch(this.quals[i - start], r);
            p = (p < 0) ? -p - 1 : p;
//			if (p< 0)
//				System.currentTimeMillis();			
            byte q = ModelPool.getAscii((int) p + ModelPool.qualLevels[0]);
            int cpos = i + from;
            quals[cpos] = ((q < quals[cpos]) ? q : quals[cpos]);
        }
    }

}
