package barna.flux.capacitor.graph;

import barna.model.splicegraph.Node;
import barna.model.splicegraph.SplicingGraph;

import java.util.Arrays;

/**
 * Extension of <code>SimpleEdgeMappings</code> for all-intronic regions.
 *
 * @author  Emilio Palumbo (emiliopalumbo@gmail.com)
 */
public class SimpleEdgeIntronMappings extends SimpleEdgeMappings {

    /**
     * Number of bins.
     */
    final int binCount = 10;

    /**
     * Array with the end position of each bin.
     */
    int[] binLimits = new int[binCount];

    /**
     * Array with the distribution of reads over the edge.
     */
    int[] binReads = new int[binCount];

    public SimpleEdgeIntronMappings(Node newTail, Node newHead) {
        super(newTail, newHead);
        binLimits = setLimits(binCount);
    }

    /**
     * Splits the edge in a certain number of beans.
     *
     * @param count number of bins
     * @return the array whith the limits of the bins
     */
    protected int[] setLimits(int count) {
        int[] limits = new int[count];
        for (int i=0;i<count;i++) {
            limits[i] = this.getTail().getSite().getPos();
            limits[i] += i==9?length()+1:(length()/binCount)*(i+1);
        }
        return limits;
    }

    /**
     * Increase the number of reads and update the bins array.
     *
     * @param readStartPos start postion of the read
     * @param readEndPos end position of the read
     * @param count whether increase the read number or only update the read distribution
     */
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

    /**
     * Increase the number of reads on the reverse direction and update the bins array.
     *
     * @param readStartPos start postion of the read
     * @param readEndPos end position of the read
     * @param count whether increase the read number or only update the read distribution
     */
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

    /**
     * Return the bins coverage.
     *
     * @return fractions of the bin covered out of the total number of bins
     */
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

    @Override
    public boolean isIntronic() {
        return true;
    }

    /*@Override
    public boolean equals(Object obj) {
        if (!obj.getClass().isAssignableFrom(SimpleEdgeIntronMappings.class))
            return false;
        SimpleEdgeIntronMappings e= (SimpleEdgeIntronMappings) obj;
        if (getTail().equals(e.getTail())&& getHead().equals(e.getHead())
                && SplicingGraph.equalSet(getTranscripts(), e.getTranscripts())
                && isExonic()== isExonic()		// multiple edges exonic, intronic for eg intron retention
                && isIntronic()== isIntronic())
            return true;
        return false;
    } */


}
