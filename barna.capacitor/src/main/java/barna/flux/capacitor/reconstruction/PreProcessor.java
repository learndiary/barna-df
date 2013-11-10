package barna.flux.capacitor.reconstruction;

import barna.commons.Execute;
import barna.commons.log.Log;
import barna.io.FileHelper;
import barna.io.gtf.GTFwrapper;
import barna.io.sam.FilteredSAMRecordSet;
import barna.io.sam.SAMConstants;
import barna.io.sam.SAMReader;
import barna.model.Gene;
import barna.model.SuperLocus;
import barna.model.sam.SAMMapping;
import net.sf.samtools.*;

import java.io.*;
import java.util.*;
import java.util.concurrent.*;

/**
 * Preprocessor for transcriptome annotation and read mappings.
 * User: micha
 */
public class PreProcessor implements Callable<File> {

    /**
     * File with transcriptome annotation in GTF format.
     */
    protected File annotation= null;

    /**
     * Number of different chromosomes in the transcriptome annotation.
     */
    protected int nrChr= -1;

    /**
     * Current set of gene loci.
     */
    protected Gene[] genes= null;

    /**
     * File with mappings in SAM/BAM format.
     */
    protected File mappings= null;

    /**
     * @return hashtable with gene loci per reference sequence, sorted by start position
     */
    public HashMap<String, Gene[]> getHashGenes() {
        if (hashGenes== null) {
            int n= (nrChr< 0? 20: nrChr);
            hashGenes = new HashMap<String, Gene[]>(nrChr, 1f);
        }
        return hashGenes;
    }

    /**
     * Hashtables to localize a mapping in gene space.
     * The <code>String</code> storing a reference name
     * (= chromosome), maps to an <code>Gene[]</code>
     * storing the gene loci along the reference, sorted
     * by their start position.
     */
    private HashMap<String,Gene[]> hashGenes = null;

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
     * Creates an instance with the refernce annotation and corresponding mappings.
     * @param annotation a transcriptome annotation in GTF format
     * @param mappings mappings in SAM/BAM format
     */
    public PreProcessor(File annotation, File mappings) {
        this.annotation= annotation;
        this.mappings= mappings;
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
    public static HashMap<String,Gene[]> index(Gene[] genes, HashMap<String,Gene[]> hashGenes) {

        if (hashGenes== null)
            throw new RuntimeException("Hashmap cannot be null.");
        if (genes== null|| genes.length== 0)
            return hashGenes;

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

        return hashGenes;
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
        if (gg== null)
            return null;    // no locus on this chr
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
     * Finds a gene locus overlapping the queried coordinates.
     * @param chr reference sequence
     * @param start start position
     * @param end end position
     * @param contained <code>true</code> if the query coordinates have to be completely contained in the query region,
     *                  <code>false</code> otherwise (partial overlap)
     * @return
     */
    public Gene getGene(String chr, int start, int end, boolean contained) {

        return getGene(hashGenes, chr, start, end, contained);
    }

    /**
     * Iterates mappings ordered by query name (= readID) and joins loci that are bound to each other by multi-mappings.
     * Works in O(2N+M) for N input Genes and M being the read with most multimaps.
     * @param inGenes loci that are already collapsed, i.e., that occupy
     *                distinct, unique genomic regions
     * @param reader source of mappings, must be sorted by query name,
     *               i.e., readID
     * @return gene loci collapsed by multi-mapping constraints
     */
    public Gene[] bind(Gene[] inGenes, SAMFileReader reader, SAMFileWriter writer) {

        if (!reader.getFileHeader().getSortOrder().equals(SAMFileHeader.SortOrder.queryname))
            throw new RuntimeException("SAM/BAM file must by sorted by query name (readID), but it is "
                    + reader.getFileHeader().getSortOrder());

        long t0= System.currentTimeMillis();
        SAMRecordIterator iter= reader.iterator();
        SAMRecord map= null;
        String name= null, lastName= null, chr= null;
        int start, end;
        FilteredSAMRecordSet set= new FilteredSAMRecordSet();
        HashMap<Gene,Integer> currentBatch= new HashMap<Gene,Integer>(10, 1f);
        HashMap<Gene,Gene> superHash= new HashMap<Gene, Gene>(inGenes.length, 1f);
        HashSet<Gene> hashSet= new HashSet<Gene>(10,1f), hashSet2= new HashSet<Gene>(10, 1f);
        while(iter.hasNext()) {
            map= iter.next();
            name= map.getReadName();

            if (!name.equals(lastName)) {
                // output previous batch
                if (lastName!= null) {
                    set.output(writer);
                    // post-process super-locus
                    if (currentBatch.size()> 1) {   // fair enough to create a new super-locus
                        Gene[] gg= new Gene[currentBatch.size()];
                        gg= currentBatch.keySet().toArray(gg);   // atomary genes that are joined by this batch
                        for (int i = 0; i < gg.length; i++) {    // collect non-redundant list of all members
                            if (superHash.containsKey(gg[i])) {
                                SuperLocus sl= (SuperLocus) superHash.remove(gg[i]);
                                if (!hashSet2.contains(sl)) {
                                    hashSet2.add(sl);
                                    Gene[] g2= sl.getGenes();
                                    for (int j = 0; j < g2.length; j++)
                                        hashSet.add(g2[j]);
                                }
                            } else
                                hashSet.add(gg[i]);
                        }
                        // create SL
                        Gene[] baseG= new Gene[hashSet.size()];
                        baseG= hashSet.toArray(baseG);
                        SuperLocus sl= new SuperLocus(Gene.getUniqueID(), baseG);
                        for (int i = 0; i < baseG.length; i++) {
                            superHash.put(baseG[i], sl);
                        }
                        hashSet.clear();
                        hashSet2.clear();
                    } else {
                        // TODO branch here for learn()
                    }
                }
                set.reset();
                currentBatch.clear();
                lastName= name;
            }

            // add the mapping to the set
            chr= map.getReferenceName();
            start= map.getAlignmentStart();
            end= map.getAlignmentEnd();
            Gene g= getGene(chr, start, end, true);
            if (g== null)
                continue;
            if (set.add(map)) {     // else
                currentBatch.put(g, (currentBatch.containsKey(g)? currentBatch.get(g)+ 1: 1));
            }
        }
        // output last batch
        set.output(writer);

        // build result, hierachy max 3 (Gene->Antisense->SL)
        int[] count= new int[20];
        Arrays.fill(count, 0);
        for (int i = 0; i < inGenes.length; i++) {
            if (superHash.containsKey(inGenes[i])) {
                if (!hashSet.contains(inGenes[i])) {
                    SuperLocus sl= (SuperLocus) superHash.get(inGenes[i]);
                    int v= Math.min(sl.getGenes().length, count.length);
                    ++count[v- 1];
                    hashSet.add(sl);
                }
            } else {
                hashSet.add(inGenes[i]);
                ++count[0];
            }
        }

        // provide summary
        StringBuilder msg= new StringBuilder(), val= new StringBuilder();
        int sum= 0;
        for (int i = 0; i < count.length; i++) {
            sum+= count[i];
            msg.append((i+ 1)+ " ");
            val.append(count[i]+ " ");
            while (msg.length()< val.length())
                msg.append(" ");
            while (val.length()< msg.length())
                val.append(" ");
        }
        Log.info(msg.toString());
        Log.info(val.toString());
        Log.info("Concatenated "+ sum+ " loci in "+ hashSet.size()+ " loci connected by multiple mappings.");
        Log.info("Took "+ (System.currentTimeMillis()- t0)/ 1000+ " sec.");

        // prepare result
        Gene[] outGenes= new Gene[hashSet.size()];
        outGenes= hashSet.toArray(outGenes);
        return outGenes;
    }

    /**
     * Performs the tasks: (1) sort mappings by query name (if not already), (2) load annotation, (3) cluster genes.
     * Writes a file with mappings sorted by genomic positions and indexed.
     * @return file handle of the preprocessed mappings which is already indexed
     */
    @Override
    public File call() {

        SAMFileReader inReader= new SAMFileReader(mappings);    // do NOT use eager decoding
        SAMFileHeader inHeader= inReader.getFileHeader();
        SAMFileWriterFactory factory= SAMConstants.getFactory();
        // streaming
        PipedInputStream pipi= null;
        PipedOutputStream pipo= null;
        Future<Integer> captain= null;

        // if not presorted by query name, do so now (likely to take long)
        if (!inHeader.getSortOrder().equals(SAMFileHeader.SortOrder.queryname)) {
            inReader.close();
            pipi= new PipedInputStream();
            pipo= null;
            try {
                pipo= new PipedOutputStream(pipi);
            } catch (IOException e) {
                throw new RuntimeException(e);
            }
            SAMFileHeader outHeader= inHeader.clone();
            outHeader.setSortOrder(SAMFileHeader.SortOrder.queryname);
            factory.makeSAMWriter(outHeader, false, pipo);
            SAMConstants.SAMSorter sorter= new SAMConstants.SAMSorter(mappings, pipo, true, false);
            captain= Execute.getExecutor().submit(sorter);
            inReader= new SAMFileReader(pipi);
            inHeader= inReader.getFileHeader();
        }

        // load annotation, cluster genes
        GTFwrapper gtfWrapper= new GTFwrapper(annotation);
        gtfWrapper.loadAllGenes();
        genes= gtfWrapper.getGenes();
        nrChr= gtfWrapper.getReadChr().size();
        gtfWrapper.close();
        gtfWrapper= null;

        // process annotation
        int n1= genes.length;
        genes= collapse(genes); // collapse anti/sense loci
        Log.info("Collapsed "+ n1+ " atomary loci to "+ genes.length+ " loci joined by anti-sense transcription.");
        hashGenes= index(genes, getHashGenes());

        // bind and write to position sorted file
        SAMFileHeader outHeader= inHeader.clone();
        outHeader.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        // got to write a file here, for indexing later-on
        File outFile= null;
        try {
            outFile= FileHelper.createTempFile(mappings.getName(),"sam");
        } catch (IOException e) {
            throw new RuntimeException(e);
        }

        factory= new SAMFileWriterFactory().setCreateIndex(true);
        SAMFileWriter writer= factory.makeSAMWriter(outHeader, false, outFile);
        if (captain!= null) {
            int n= 0;
            try {
                n= captain.get();
            } catch (Exception e) { // Interrupted, Execution
                throw new RuntimeException(e);
            }
            Log.info("Sorted SAM/BAM by query name, "+ n+ " lines.");
        }
        int n= genes.length;
        genes= bind(genes, inReader, writer);
        Log.info("Clustered "+ n+" genes into "+genes.length +" sets.");

        return outFile;
    }
}
