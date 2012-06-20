package barna.flux.capacitor.graph;

import barna.model.splicegraph.Node;

import java.util.Arrays;

/**
 * Created with IntelliJ IDEA.
 * User: Emilio Palumbo
 * Date: 6/7/12
 * Time: 11:04 AM
 */
public class SimpleEdgeIntronMappings extends SimpleEdgeMappings {

    final int binCount = 10;
    int[] binLimits = new int[binCount];
    int[] binReads = new int[binCount];

    public SimpleEdgeIntronMappings(Node newTail, Node newHead) {
        super(newTail, newHead);
        binLimits = setLimits(binCount);
    }

    protected int[] setLimits(int count) {
        int[] limits = new int[count];
        for (int i=0;i<count;i++) {
            limits[i] = this.getTail().getSite().getPos();
            limits[i] += i==9?length()+1:(length()/binCount)*(i+1);
        }
        return limits;
    }

    public void incrReadNr(int readStartPos, int readEndPos, boolean count) {

        if (binLimits[0] > 0) {
            int p = Arrays.binarySearch(binLimits,readStartPos);
            if (p<0) {
                p = -(p+1);
            }
            int q = Arrays.binarySearch(binLimits,readEndPos);
            if (q<0) {
                q = -(q+1);
            }
            for (int i = p;i<=q;i++) {
                binReads[i]++;
            }
            if (count)
                mappings.incrReadNr();
        } else
            incrRevReadNr(readStartPos,readEndPos,count);
    }

    public void incrRevReadNr(int readStartPos, int readEndPos, boolean count) {

        int gstart = -readEndPos;
        int gend = -readStartPos;

        int p = Arrays.binarySearch(binLimits,gstart);
        if (p<0) {
            p = -(p+1);
        }
        int q = Arrays.binarySearch(binLimits,gend);
        if (q<0) {
            q = -(q+1);
        }
        for (int i = p;i<=q;i++) {
            binReads[i]++;
        }
        if (count)
            mappings.incrRevReadNr();
    }

    public float getBinCoverage() {
        float count = 0;
        for (int i = 0;i<binReads.length;i++) {
            if (binReads[i]>0)
                count++;
        }
        return count / binReads.length;
    }

    @Override
    public boolean isAllIntronic() {
        return true;
    }

    /*protected double getReadDist(int[] binReads, int n) {
        int count = 0, sum = 0;
        for (int i = 0;i<binReads.length;i++) {
            if (binReads[i]>0) {
                count++;
                sum+=binReads[i];
            }
        }
        return (double)sum/count;
    } */
}
