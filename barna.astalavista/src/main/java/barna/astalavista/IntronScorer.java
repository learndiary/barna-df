package barna.astalavista;

import barna.io.gtf.GTFwrapper;
import barna.model.Gene;
import barna.model.IntronModel;
import barna.model.splicegraph.SplicingGraph;

import java.io.File;

/**
 * Created with IntelliJ IDEA.
 * User: micha
 * Date: 12/3/12
 * Time: 9:11 PM
 * @deprecated to be merged with some tool in a better moment
 */
public class IntronScorer {

    public static void extractSpliceJunctions(int eFlankDon, int eFlankAcc, IntronModel iModel, File inF, File outF) {
        if (outF!= null&& outF.exists())
            outF.delete();
        GTFwrapper reader= new GTFwrapper(inF.getAbsolutePath());
        Gene[] g= null;
        try {
            for (reader.read(); (g= reader.getGenes())!= null; reader.read()) {
                for (int i = 0; i < g.length; i++) {
                    SplicingGraph gr= new SplicingGraph(g[i]);
                    gr.constructGraph();

                    gr.writeSpliceJunctionSeqs(eFlankDon, eFlankAcc, iModel, outF);
                }
            }
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }


}
