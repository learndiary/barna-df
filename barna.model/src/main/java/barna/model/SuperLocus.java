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
    }

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
