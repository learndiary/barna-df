package barna.flux.capacitor.graph;

import barna.flux.capacitor.reconstruction.FluxCapacitorSettings;
import barna.io.MSIterator;
import barna.model.Gene;
import barna.model.Mapping;
import barna.model.SuperLocus;
import barna.model.splicegraph.EdgeSet;

import java.io.File;
import java.util.EnumSet;
import java.util.HashMap;

/**
 * Class to manage the mapping to a set of genes.
 *
 * @author Michael Sammeth
 */
public class SuperAnnotationMapper {

    AnnotationMapper[] annos= null;

    HashMap<EdgeSet, Integer> esetMap;

    public SuperAnnotationMapper(SuperLocus sl, boolean paired, boolean stranded,
                                 boolean weighted, FluxCapacitorSettings.ReadStrand readStrand,
                                 EnumSet<ComplexCounter.CounterType> counterTypes) {

        Gene[] genes= sl.getGenes();
        annos= new AnnotationMapper[genes.length];
        int i= 0;
        for(Gene g:genes)
            annos[i++]= new AnnotationMapper(g, paired, stranded, weighted, readStrand, counterTypes);
    }

    /**
     * Delegates the alignment of <code>mappings</code> to all wrapped <code>AnnotationMapper</code>
     * instances.
     * @param mappings iterator over input mappings
     * @param insertFile file to write insert sizes to
     */
    public void map(MSIterator<Mapping> mappings, File insertFile) {

        esetMap= new HashMap<EdgeSet, Integer>(Math.max(mappings.size(), 1), 1f);

        for (int i = 0; i < annos.length; i++) {
            annos[i].map(mappings, insertFile);
        }

    }


}
