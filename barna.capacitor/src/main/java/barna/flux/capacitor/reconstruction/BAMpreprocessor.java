package barna.flux.capacitor.reconstruction;

import barna.io.sam.SAMReader;
import barna.model.Gene;
import barna.model.SuperLocus;
import barna.model.commons.IntVector;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

import java.util.*;

/**
 * Created with IntelliJ IDEA.
 * User: micha
 * Date: 11/8/13
 * Time: 2:05 PM
 * To change this template use File | Settings | File Templates.
 */
public class BAMpreprocessor {

    /**
     * Atomary gene loci.
     */
    protected Gene[] genes= null;

    /**
     * Reader from which source SAM/BAM annotation is read
     */
    protected SAMReader reader= null;

    /**
     * @return hashtable with gene loci per reference sequence, sorted by start position
     */
    public HashMap<String, Gene[]> getHashStart() {
        return hashStart;
    }

    /**
     * Hashtables to localize a mapping in gene space.
     * The <code>String</code> storing a reference name
     * (= chromosome), maps to an <code>Gene[]</code>
     * storing the gene loci along the reference, sorted
     * by their start position.
     */
    private HashMap<String,Gene[]> hashStart= null;

    /**
     * Comparator for sorting by genomic position and detecting overlaps.
     */
    private static Gene.OverlapComparator ocompi= null;

    /**
     * Returns a comparator for sorting by genomic position and detecting overlaps.
     * @return the comparator
     */
    public static Gene.OverlapComparator getOverlapComparator() {
        if (ocompi== null) {
            ocompi= new Gene.OverlapComparator();
        }
        return ocompi;
    }

    /**
     * A comparator for gene start positions.
     */
    private static Gene.BoundaryComparator startCompi= null;

    /**
     * Returns a comparator for gene start positions.
     * @return the comparator
     */
    public static Gene.BoundaryComparator getStartComparator() {
        if (startCompi== null) {
            startCompi= new Gene.BoundaryComparator(true);
        }
        return startCompi;
    }

    /**
     * Creates an instance with genes and mappings.
     * @param genes gene loci to work on
     * @param reader interface to get mappings
     */
    public BAMpreprocessor(Gene[] genes, SAMReader reader) {
        this.genes= genes;
        this.reader= reader;
        collapse(genes);
        hashStart= new HashMap<String, Gene[]>();
    }

    /**
     * Determines genomic overlap of two genomic regions.
     * @param begin1
     * @param end1
     * @param begin2
     * @param end2
     * @return <code>true</code> if both genes intersect in their genomic coordinates,
     * <code>false</code> otherwise.
     */
    public static boolean overlaps(int begin1, int end1, int begin2, int end2) {
        if ((begin1>= begin2&& begin1<= end2)|| (begin2>= begin1&& begin2<= end1))
            return true;   // overlap
        return false;
    }

    /**
     * Collapses overlapping genes into super-loci. The method does not alternate the original array and performs in
     * O(3N) memory, for an initial array of N elements.
     * @param genes array with gene loci that can mutually overlap
     * @return array with non-overlapping (super-) loci.
     */
    public static Gene[] collapse(Gene[] genes) {

        // sort
        Gene.OverlapComparator ocompi= getOverlapComparator();
        genes= genes.clone();
        Arrays.sort(genes, ocompi);

        // collapse
        int glen= genes.length;
        Vector<Gene> ovlGenes= new Vector<Gene>();
        for (int i = 0; i < glen; i++) {

            int b1= Math.abs(genes[i].getStart()), e1= Math.abs(genes[i].getEnd());

            // lazily collect overlapping genes
            int j = i+ 1;
            for (; j < glen; j++) {
                int b2= Math.abs(genes[j].getStart()), e2= Math.abs(genes[j].getEnd());

                if (overlaps(b1, e1, b2, e2)) {
                    ovlGenes.add(genes[j]);
                    b1= Math.min(b1, b2);
                    e1= Math.max(e1, e2);
                } else
                    break;
            }

            // create super-locus
            if (ovlGenes.size()> 0) {
                Gene[] oGenes= new Gene[ovlGenes.size()+ 1];
                oGenes[0]= genes[i];
                for (int h = 1; h < oGenes.length; h++)
                    oGenes[h]= ovlGenes.elementAt(h- 1);

                SuperLocus sl= new SuperLocus(Gene.getUniqueID(), oGenes);
                genes[i]= sl;
                System.arraycopy(genes, j, genes, i+ 1, glen- j);
                glen-= ovlGenes.size();
                ovlGenes.removeAllElements();
            }
        }

        // shorten array to result size
        Gene[] result= new Gene[glen];
        System.arraycopy(genes, 0, result, 0, glen);
        return result;
    }

    /**
     * Creates a hash with the start and another with the end positions of gene loci per
     * reference sequence.
     * @param genes gene loci to be indexed
     * @param hashGenes build up hash with start values per reference sequence
     * @return the minimum number of genes found on one of the reference sequences
     */
    public static int index(Gene[] genes, HashMap<String,Gene[]> hashGenes) {

        if (hashGenes== null)
            throw new RuntimeException("Hashmap cannot be null.");
        if (genes== null|| genes.length== 0)
            return 0;

        // initial scan for estimates to init dynamic arrays
        int n= 0, min= Integer.MAX_VALUE, nChr= 0;
        String last= null, chr= null;
        for (int i = 0; i < genes.length; i++) {
            chr= genes[i].getChromosome();
            if (chr.equals(last))
                ++n;
            else {
                if (last!= null) {
                    if (n< min)
                        min= n;
                    n= 1;
                }
                last= chr;
                ++nChr;
            }
        }

        HashMap<String, Vector<Gene>> mapGene= new HashMap<String, Vector<Gene>>(nChr,1f);
        Vector<Gene> v= null;
        last= null;
        for (int i = 0; i < genes.length; i++) {

            // sync with hashes
            chr= genes[i].getChromosome();
            if (!chr.equals(last)) {
                v= mapGene.get(chr);
                if (v== null) {
                    v= new Vector<Gene>(min);
                    mapGene.put(chr, v);
                }
            }

            // add gene coordinates
            v.add(genes[i]);
        }

        // copy result
        String[] keys= new String[mapGene.size()];
        keys= mapGene.keySet().toArray(keys);
        for (int i = 0; i < keys.length; i++) {
            String key = keys[i];
            v= mapGene.remove(key);
            Gene[] gg= new Gene[v.size()];
            for (int j = 0; j < gg.length; j++)
                gg[j]= v.elementAt(j);
            Arrays.sort(gg, getOverlapComparator());
            hashGenes.put(key, gg);
        }

        return min;
    }

    /**
     * Finds the gene locus (if one exists) that overlaps with the query sequence.
     * <b>Assumption:</b> there is maximally one gene locus that overlaps with the query region.
     * @param mapGenes index with gene loci sorted by start positions
     *                 per reference sequence
     * @param chr query reference sequence name
     * @param start query start position
     * @param end query end position
     * @param contained flag indicating whether it the region has to be completely contained in the target locus,
     *                  alternatively target loci overlapping the query are considered
     * @return the gene locus was found that harbors the queried region,
     * otherwise <code>null</code>.
     */
    public static Gene getGene(HashMap<String,Gene[]> mapGenes, String chr, int start, int end, boolean contained) {

        Comparator sCompi= getStartComparator();
        Gene[] gg= mapGenes.get(chr);
        int p= Arrays.binarySearch(gg, start, sCompi);

        if (p< 0)
            p= -(p+ 1)- 1; // insertion point- 1


        // check left neighbor
        if (p>= 0) {
            int gs= Math.abs(gg[p].getStart()),
                    ge= Math.abs(gg[p].getEnd());
            if (start>= gs&& start<= ge) {
                if ((!contained)|| (end>= gs&& end<= ge))
                    return gg[p];
            }
        }
        if (!contained) {
            ++p;    // check right neighbor
            if (p< gg.length) {
                int gs= Math.abs(gg[p].getStart()),
                        ge= Math.abs(gg[p].getEnd());
                if (end>= gs&& end<= ge)
                    return gg[p];
            }
        }

        return null;    // no locus found
    }

    /**
     * Iterates mappings ordered by query name (= readID) and joins loci
     * that are bound to each other by multi-mappings.
     * @param inGenes loci that are already collapsed, i.e., that occupy
     *                distinct, unique genomic regions
     * @param reader source of mappings, must be sorted by query name,
     *               i.e., readID
     * @return gene loci collapsed by multi-mapping constraints
     */
    public static Gene[] bind(Gene[] inGenes, SAMFileReader reader) {

        if (!reader.getFileHeader().getSortOrder().equals(SAMFileHeader.SortOrder.queryname))
            throw new RuntimeException("SAM/BAM file must by sorted by query name (readID), but it is "
                    + reader.getFileHeader().getSortOrder());

        SAMRecordIterator iter= reader.iterator();
        SAMRecord map= null;
        int minStratum= -1;
        String name= null, lastName= null;
        while(iter.hasNext()) {
            map= iter.next();
            map.getReadName();
            map.getFirstOfPairFlag();
            int mm= map.getIntegerAttribute("NM");

            if (!name.equals(lastName)) {
                if (lastName!= null) {
                    // process
                }
                lastName= name;
            }
        }

        return null;
    }


}
