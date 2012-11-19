package barna.flux.capacitor.diffexp;

public class DifferentialExpression {
    private GFFEntry source;
    private GFFEntry target;
    private double p;
    private double difference;
    private double foldChange;

    public DifferentialExpression(GFFEntry source, GFFEntry target, double p, double difference, double foldChange) {
        this.source = source;
        this.target = target;
        this.p = p;
        this.difference = difference;
        this.foldChange = foldChange;
    }

    public GFFEntry getSource() {
        return source;
    }

    public GFFEntry getTarget() {
        return target;
    }

    public double getP() {
        return p;
    }

    public double getDifference() {
        return difference;
    }

    public double getFoldChange() {
        return foldChange;
    }
}
