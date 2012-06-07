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

    public void incrReadNr(int readStartPos, int readEndPos) {
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
        mappings.incrReadNr();
    }

    public int getReadNr() {
        int n = 0;
        /*if (Math.round(getReadDist()*100)>90) {

        } */
        return n;
    }

    protected double getReadDist(int[] binReads, int n) {
        int[] reads = Arrays.copyOf(binReads, binReads.length);
        /*Arrays.sort(reads);
        int max = reads[reads.length-1];
        double value = 0;
        for (int i = reads.length-1;i>=0 && reads[i]==max;i--) {
            value++;
        }*/
        double value = 0;
        int mean = n/binReads.length;//mappings.getReadNr()/binReads.length;
        for (int i =0;i< reads.length;i++) {
            value+=Math.abs(reads[i]-mean);
        }
        return value;
    }
}
