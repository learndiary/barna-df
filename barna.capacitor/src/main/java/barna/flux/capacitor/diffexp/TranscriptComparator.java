package barna.flux.capacitor.diffexp;

import barna.flux.capacitor.diffexp.math.FishersExactTest;

import java.util.concurrent.Callable;

/**
 * Do fisher exact test
 */
class TranscriptComparator implements Callable<DifferentialExpression> {
    private Transcript source;
    private Transcript target;
    private double sourceReads;
    private double targetReads;

    public TranscriptComparator(Transcript source, Transcript target, double sourceReads, double targetReads) {
        this.source = source;
        this.target = target;
        this.sourceReads = sourceReads;
        this.targetReads = targetReads;
    }

    @Override
    public DifferentialExpression call() throws Exception {
        double p = 1.0;
        double difference = 0;
        double foldChange = 0;
        if(source != null && target != null){
            // fisher test
            p = FishersExactTest.fishersExactTest(
                    (int)Math.round(source.getReadCount()),
                    (int)Math.round(sourceReads - source.getReadCount()),
                    (int)Math.round(target.getReadCount()),
                    (int)Math.round(targetReads - target.getReadCount())
                    )[0];
            difference = target.getRpkm() - source.getRpkm();
            if(target.getRpkm() > 0 && source.getRpkm() == 0){
                foldChange = Double.POSITIVE_INFINITY;
            }else if(target.getRpkm() == 0 && source.getRpkm() > 0){
                foldChange = Double.NEGATIVE_INFINITY;
            }else {
                foldChange = Math.log(target.getRpkm() / source.getRpkm()) / Math.log(2); // log2
            }
        }else if(source != null){
            difference = -source.getRpkm();
            foldChange = Double.NEGATIVE_INFINITY;
        }else if(target != null){
            difference = target.getRpkm();
            foldChange = Double.POSITIVE_INFINITY;
        }else{
            throw new RuntimeException("No transcript specified");
        }
        return new DifferentialExpression(source, target, p, difference, foldChange);
    }
}
