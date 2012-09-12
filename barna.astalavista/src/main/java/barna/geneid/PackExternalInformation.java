package barna.geneid;

/**
 * Created with IntelliJ IDEA.
 * User: micha
 * Date: 9/11/12
 * Time: 4:21 PM
 * To change this template use File | Settings | File Templates.
 */
public class PackExternalInformation {

    static class HSP {
        long Pos1;
        long Pos2;
        float Score;
    }

    static class PackHSP {
        HSP[][] sPairs;
        long[] nSegments;
        long nTotalSegments;
        int visited;
    }

    static class PackEvidence {
        Site vSites;
        long nvSites;
        ExonGFF vExons;
        long nvExons;
    }


    Dictionary locusNames;

    PackHSP[] homology;
    PackEvidence[] evidence;

    long nSequences;
    long nvExons;
    long nHSPs;

    long i1vExons;
    long i2vExons;
    long ivExons;

    long[] iSegments;
    float[][] sr;
    float[][] readcount;
}
