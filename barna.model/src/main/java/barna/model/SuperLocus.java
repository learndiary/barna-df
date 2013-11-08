package barna.model;

/**
 * Created with IntelliJ IDEA.
 * User: micha
 * Date: 11/7/13
 * Time: 6:19 PM
 * To change this template use File | Settings | File Templates.
 */
public class SuperLocus extends Gene {

    /**
     * Genes joint in this locus
     */
    Gene[] genes= null;

    public SuperLocus(String id, Gene[] genes) {
        super(id);
        this.genes= genes;
        assert(genes!= null&& genes.length> 0);
        this.chromosome= genes[0].getChromosome();
        int min= Integer.MAX_VALUE, max= Integer.MIN_VALUE;
        for (int i = 0; i < genes.length; i++) {
            int s= Math.abs(genes[i].getStart()),
                    e= Math.abs(genes[i].getEnd());
            min= Math.min(min, s);
            max= Math.max(max, e);
        }
        this.start= min;
        this.end= max;
        if (isAntisense())
            this.strand= 2;
        else
            this.strand= genes[0].getStrand();
    }

    /**
     * Determines whether there are genes from both DNA strands
     * contained in <code>this</code> super-locus.
     * @return <code>true</code> if there is at least one gene from either DNA strand
     * contained in <code>this</code> super-locus, <code>false</code> otherwise.
     */
    public boolean isAntisense() {

        assert(genes!= null);

        byte last= genes[0].getStrand();
        for (int i = 1; i < genes.length; i++) {
            if(last!= genes[i].getStrand())
                return true;
        }
        return false;
    }


}
