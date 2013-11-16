package barna.flux.capacitor.reconstruction;

import barna.commons.Execute;
import barna.commons.log.Log;
import barna.flux.capacitor.matrix.UniversalMatrix;
import barna.flux.capacitor.profile.CoverageStats;
import barna.flux.capacitor.profile.MappingStats;
import barna.flux.capacitor.profile.Profile;
import barna.io.AbstractFileIOWrapper;
import barna.io.FileHelper;
import barna.io.MSIterator;
import barna.io.MappingReader;
import barna.io.gtf.GTFwrapper;
import barna.io.rna.UniversalReadDescriptor;
import barna.io.sam.FilteredSAMRecordSet;
import barna.io.sam.SAMConstants;
import barna.io.sam.SAMReader;
import barna.model.*;
import barna.model.commons.MyFile;
import barna.model.constants.Constants;
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
     * Minimum length recorded for a mapping.
     */
    protected int mapLenMin= -1;


    /**
     * Maximum length recorded for a mapping.
     */
    protected int mapLenMax= -1;

    /**
     * The settings
     */
    protected FluxCapacitorSettings settings= null;

    /**
     * The profile storing biases for different transcript lengths.
     */
    protected Profile profile;

    /**
     * Number of different chromosomes in the transcriptome annotation.
     */
    protected int nrChr= -1;

    /**
     * Current set of gene loci.
     */
    protected Gene[] genes= null;

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
     * @param settings a setting file containing the annotation and mapping file
     */
    public PreProcessor(FluxCapacitorSettings settings) {
        this.settings= settings;
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
     * Finds a gene locus which the mapping hits.
     * @param mapping a mapping
     * @param contained <code>true</code> if the mapping is to be completely included in the locus,
     *                  <code>false</code> otherwise
     * @return <code>null</code> or the gene locus that contains the mapping
     */
    public Gene getGene(SAMRecord mapping, boolean contained) {
        return getGene(mapping.getReferenceName(), mapping.getAlignmentStart(), mapping.getAlignmentEnd(), true);
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
        boolean pairedEnd= false;   // dynamically set for every read
        FilteredSAMRecordSet set= new FilteredSAMRecordSet(), set2= new FilteredSAMRecordSet();
        HashSet<Gene> hashSet= new HashSet<Gene>(10,1f), hashSet2= new HashSet<Gene>(10, 1f);
        ArrayList<Gene> gList= null;
        HashMap<Gene,ArrayList<SAMRecord>> gHash= new HashMap<Gene, ArrayList<SAMRecord>>(10), gHash2= null;
        HashMap<Gene,Gene> superHash= new HashMap<Gene, Gene>(inGenes.length, 1f);
        Iterator<Gene> iGene= null;
        long inCount= 0;
        while(iter.hasNext()) {

            ++inCount;
            map= iter.next();
            name= map.getReadName();

            if (!name.equals(lastName)) {
                // output previous batch
                if (lastName!= null) {
                    gList= outputAnnotationMapped(gList, set, set2, gHash, gHash2, writer);
                    // post-process super-locus
                    createSL(gList, superHash, hashSet, hashSet2);
                }
                set.reset();
                if (pairedEnd)
                    set2.reset();
                lastName= name;
                pairedEnd= map.getReadPairedFlag();
                if(pairedEnd&& gHash2== null)
                    gHash2= new HashMap<Gene, ArrayList<SAMRecord>>(10);
            }

            // add the mapping to the set
            if (pairedEnd) {
                if(map.getFirstOfPairFlag())
                    set.add(map);
                else
                    set2.add(map);
            }

        }
        // output last batch
        outputAnnotationMapped(gList, set, set2, gHash, gHash2, writer);
        createSL(gList, superHash, hashSet, hashSet2);


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
     * Creates a super-locus from genes that are fused by multi-mappings.
     * @param gList list of genes that are fused
     * @param superHash hash that stores super-loci created so far for distinct genes
     * @param hashSet re-use set
     * @param hashSet2 re-use set
     * @return hash that stores super-loci created so far and the new one (if created)
     */
    protected HashMap<Gene, Gene> createSL(ArrayList<Gene> gList, HashMap<Gene,Gene> superHash, HashSet<Gene> hashSet, HashSet<Gene> hashSet2) {

        if (gList== null|| gList.size()== 0)
            return superHash;

        if (gList.size()== 1) {
            // TODO branch here for learn()

        } else {    // fair enough to create a new super-locus

            // safety first
            hashSet.clear();
            hashSet2.clear();

            // iterate atomary genes that are joined by this batch
            Iterator<Gene> iGene= gList.iterator();
            while (iGene.hasNext()) {
                // collect non-redundant list of all members
                Gene g= iGene.next();
                if (superHash.containsKey(g)) {
                    SuperLocus sl= (SuperLocus) superHash.remove(g);
                    if (!hashSet2.contains(sl)) {
                        hashSet2.add(sl);
                        Gene[] g2= sl.getGenes();
                        for (int j = 0; j < g2.length; j++)
                            hashSet.add(g2[j]);
                    }
                } else
                    hashSet.add(g);
            }
            // create SL
            Gene[] baseG= new Gene[hashSet.size()];
            baseG= hashSet.toArray(baseG);
            SuperLocus sl= new SuperLocus(Gene.getUniqueID(), baseG);
            for (int i = 0; i < baseG.length; i++) {
                superHash.put(baseG[i], sl);
            }

            // for GC
            hashSet.clear();
            hashSet2.clear();

        }

        return superHash;
    }


    /**
     * Control gateway for file creation from the main class,
     * adds a hook for delete on exit in case.
     *
     * @param f            the file that has been created
     * @param deleteOnExit flag to mark for deletion on exit
     * @return
     */
    protected File createFile(File f, boolean deleteOnExit) {
        if (deleteOnExit)
            f.deleteOnExit();

        return f;
    }

    /**
     * Creates a temporary file in the location provided, iff write access is
     * available there. Otherwise the file is created in the custom or system
     * temporary directory.
     *
     * @param location     a file in the target directory or the directory itself,
     *                     may be <code>null</code>
     * @param name         prefix of the file to be created, class name is appended
     *                     at the beginning
     * @param extension    (optional) suffix of the temporary file that is created
     * @param deleteOnExit flag for calling the <code>deleteOnExit()</code>
     *                     method for the file
     * @return a temporary file according to the specifications
     */
    protected File createTempFile(File location, String name, String extension, boolean deleteOnExit) {

        // get location
        if (location == null)
            location = settings.get(FluxCapacitorSettings.TMP_DIR);
        else {
            if (!location.isDirectory())
                location = location.getParentFile();
            if (!location.canWrite())
                location = settings.get(FluxCapacitorSettings.TMP_DIR);
        }

        // get name
        if (name == null)
            name = getClass().getSimpleName();
        else
            name = getClass().getSimpleName() + "_" + name;

        File f = null;
        try {
            f = FileHelper.createTempFile(name, extension, location);
        } catch (IOException e) {
            throw new RuntimeException(e);
        }

        return createFile(f, deleteOnExit);
    }



    protected File sortAnotation(GTFwrapper wrapper, File inputFile) {

        File sortedDir = settings.get(FluxCapacitorSettings.KEEP_SORTED);
        File f;
        if (sortedDir!=null)
            f = FileHelper.getSortedFile(new File(sortedDir, inputFile.getName()));
        else
            f = FileHelper.getSortedFile(inputFile);
        File lock = FileHelper.getLockFile(f);

        if (f.exists() && !lock.exists()) {

            Log.warn("Assuming file " + f.getName() + " is a sorted version of " + inputFile.getName());

        } else {    // we have to sort

            boolean lockCreated = false;
            if (sortedDir!=null) {//settings.get(FluxCapacitorSettings.KEEP_SORTED)) {    // try to store in original

                if (lock.exists()) {    // switch to sorting to temp
                    Log.warn("Seems that another process is just sorting file " + inputFile +
                            "\nremove lock file " + lock.getName() + " if dead leftover." +
                            "\nContinuing with sorting to temporary file " +
                            (f = createTempFile(f,     // access to non-Temp
                                    FileHelper.getFileNameWithoutExtension(f),
                                    FileHelper.getExtension(f),
                                    false)).getAbsolutePath());

                } else if (!f.getParentFile().canWrite()) {    // sort to temp, but do not delete (parameter)
                    Log.warn("Cannot write sorted file to " + f.getAbsolutePath() +
                            "\nContinuing with sorting to temporary file " +
                            (f = createTempFile(f, // access to non-Temp
                                    FileHelper.getFileNameWithoutExtension(f),
                                    FileHelper.getExtension(f),
                                    false)).getAbsolutePath());

                } else {    // sort to default sorted file
                    try {
                        lock.createNewFile();
                    } catch (Exception e) {
                        throw new RuntimeException(e);
                    }
                    lockCreated = true;
                }

            } else {    // do not keep sorted files, sort to temp and delete on exit
                f = createTempFile(null,
                        FileHelper.getFileNameWithoutExtension(f),
                        FileHelper.getExtension(f),
                        true);
            }

            // doit
            if (wrapper.getInputFile() != null)
                Log.info("Sorting " + wrapper.getInputFile().getAbsolutePath());
            // TODO make sorter only use "exon" feature
            wrapper.sort(f);

            // if locked
            if (lockCreated)
                lock.delete();
        }
        return f;
    }

    /**
     * File handle with the pre-processed annotation (if necessary).
     */
    protected File annotation= null;

    /**
     * Obtains the file handle for the current mapping file, which may be different from the one provided in the
     * <code>FluxCapacitorSettings</code> after pre-processing has been carried out.
     * @return file handle for the pre-processed mapping file
     */
    public File getMappingFile() {
        if (mappings== null ) {
            return settings.get(FluxCapacitorSettings.MAPPING_FILE.getName());
        }
        return mappings;
    }

    /**
     * Obtains the file handle for the current annotation file. The method is a bit obsolete because currently the
     * pre-processed annotation is loaded into gene models.
     * @return file handle to the pre-processed annotation file, which only differs from the one provided if re-sorting
     * has been required to build up the gene models
     */
    public File getAnnotationFile() {
        if (annotation== null ) {
            return settings.get(FluxCapacitorSettings.ANNOTATION_FILE.getName());
        }
        return annotation;
    }

    /**
     * File handle for the pre-processed mapping file. The returned file is likely different than to the one specified
     * in the <code>FluxCapacitorSettings</code> instance provided.
     */
    protected File mappings= null;


    /**
     * Performs the tasks: (1) sort mappings by query name (if not already), (2) load annotation, (3) cluster genes.
     * Writes a file with mappings sorted by genomic positions and indexed.
     * @return file handle of the preprocessed mappings which is already indexed
     */
    @Override
    public File call() {

        annotation= settings.get(FluxCapacitorSettings.ANNOTATION_FILE.getName());
        mappings= settings.get(FluxCapacitorSettings.MAPPING_FILE.getName());

        SAMFileReader inReader= new SAMFileReader(mappings);    // do NOT use eager decoding
        SAMFileHeader inHeader= inReader.getFileHeader();
        SAMFileWriterFactory factory= SAMConstants.getFactory();
        // streaming
        PipedInputStream pipi= null;
        PipedOutputStream pipo= null;
        SAMConstants.SAMSorter sorter= null;
        Future<Long> captain= null;

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
            sorter= new SAMConstants.SAMSorter(mappings, pipo, true, false);
            sorter.setSkippingNotmapped(true);
            captain= Execute.getExecutor().submit(sorter);
        }

        // load annotation, cluster genes
        GTFwrapper gtfWrapper= new GTFwrapper(annotation);
        if (!gtfWrapper.isApplicable()) {
            // do not have to fork another time, mapping sorting is already async
            File f= sortAnotation(gtfWrapper, annotation);
            gtfWrapper.sort(f);
        }
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

        // prepare output
        SAMFileHeader outHeader= inHeader.clone();
        outHeader.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        // got to write a file here, for indexing later-on
        File outFile= null;
        try {
            outFile= FileHelper.createTempFile(mappings.getName(),"bam");
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
        factory= new SAMFileWriterFactory().setCreateIndex(true);
        // factory.setUseAsyncIo(true); // not sure whether this is a good idea
        SAMFileWriter writer= factory.makeBAMWriter(outHeader, false, outFile);

        // process annotation+ mappings
        int n= genes.length;
        if (captain!= null) {
            // don't wait for captain, output has to be read by this thread otherwise piped-output blocks
            inReader= new SAMFileReader(pipi);
            inHeader= inReader.getFileHeader();
        }
        genes= bind(genes, inReader, writer);
        Log.info("Clustered "+ n+" genes into "+genes.length +" sets.");

        // close I/O
        long t0= System.currentTimeMillis();
        writer.close();
        Log.info("Wrote BAM file in "+ (System.currentTimeMillis()- t0)/ 1000+ " sec: "+outFile.getAbsolutePath());

        return outFile;
    }

    protected void outputMapping(SAMRecord map, SAMFileWriter writer) {

        writer.addAlignment(map);

        // update stats about min/max mapping length
        Iterator<AlignmentBlock> blocks= map.getAlignmentBlocks().iterator();
        int sum= 0;
        while(blocks.hasNext())
            sum+= blocks.next().getLength();
        if (mapLenMin < 0 || sum < mapLenMin)
            mapLenMin = sum;
        if (sum > mapLenMax)
            mapLenMax = sum;

    }


    /**
     * Maps single mappings (of one of the mates in the case of paired-end mappings) to genes, stores results in a
     * hash if provided.
     * @param gList re-using <code>ArrayList</code> for passing result, only employed in the case of single mappings
     * @param set set of filtered mappings (for one of the mates if paired-end)
     * @param gHash re-use gene hash for mate 1, unused (<code>null</code>) for single end
     * @param writer I/O interface to write the annotation filtered reads to
     * @return if single ended, a list of genes that are fused by multi-mappings of the read is returned.
     *         For paired-end mappings the re-use list is returned unchanged.
     */
    protected ArrayList<Gene> outputMapped(ArrayList<Gene> gList, FilteredSAMRecordSet set, HashMap<Gene, ArrayList<SAMRecord>> gHash, SAMFileWriter writer) {

        if (gHash== null) {    // single end, need list
            if (gList== null)
                gList= new ArrayList<Gene>(Math.min(set.getSet().size(), 10));  // limit to 10 connected loci initially
            else
                gList.clear();
        }

        Iterator<SAMRecord> i2= set.iterator();
        SAMRecord map= null;
        ArrayList<SAMRecord> hits= null;

        while (i2.hasNext()) {
            map= i2.next();
            Gene g= getGene(map, true);
            if (g== null)
                continue;
            if (gHash== null) {
                gList.add(g);
                outputMapping(map, writer); // single end
                continue;
            }
            // else: paired-end
            hits= gHash.get(g);
            if (hits== null) {
                hits= new ArrayList<SAMRecord>(2);
                gHash.put(g, hits);
            }
            hits.add(map);
        }

        return gList;
    }

    /**
     * Maps mappings to the annotation and determines the loci that are connected by multi-mappings. In the case of
     * paired-end reads the method returns:
     * <ol>
     * <li>(1) only loci to which both mates map, if any</li>
     * <li>(2) all loci to which mate1 or mate2 are mapping iff the other mate has no mapping respectively</li>
     * <li>(3) neigboring genes that are fused by paired-end mappings of mate1 and mate2</li>
     * </ol>
     * @param gList re-using <code>ArrayList</code> for passing result
     * @param set set of filtered mappings (mate1 if paired-end)
     * @param set2 set of filtered mappings for mate2, or <code>null</code> if single end reads
     * @param gHash re-use gene hash for mate 1, unused (<code>null</code>) for single end
     * @param gHash2 re-use gene hash for mate 2, unused (<code>null</code>) for single end
     * @param writer I/O interface to write the annotation filtered reads to
     * @return a list of genes that are fused by multi-mappings of the read
     */
    protected ArrayList<Gene> outputAnnotationMapped(ArrayList<Gene> gList, FilteredSAMRecordSet set, FilteredSAMRecordSet set2,
                                          HashMap<Gene, ArrayList<SAMRecord>> gHash, HashMap<Gene, ArrayList<SAMRecord>> gHash2, SAMFileWriter writer) {

        boolean pairedEnd= (set2!= null);
        gHash.clear();
        if (gHash2!= null)
            gHash2.clear();

        // first round
        gList= outputMapped(gList, set, (pairedEnd? gHash: null), writer);
        if (!pairedEnd)
            return gList; // single end done

        // second round for paired-end
        outputMapped(gList, set2, gHash2, writer);

        // intersect paired-end
        Iterator<Gene> i2= gHash.keySet().iterator();
        int isectMax= Math.min(gHash.size(), gHash2.size());
        if (gList== null)
            gList= new ArrayList<Gene>(isectMax);
        else
            gList.clear();
        Gene g= null;
        if (isectMax> 0) {
            while(i2.hasNext()) {
                g= i2.next();
                if (gHash2.containsKey(g))
                    gList.add(g);
            }
        }

        Iterator<SAMRecord> iRec= null;
        if (gList.size()> 0) {
            // case 1: there are loci with both mappings
            // output mappings for both mates for the corresponding loci
            i2= gList.iterator();
            while(i2.hasNext()) {
                g= i2.next();
                iRec= gHash.get(g).iterator();
                while (iRec.hasNext())
                    outputMapping(iRec.next(), writer);
                iRec= gHash2.get(g).iterator();
                while (iRec.hasNext())
                    outputMapping(iRec.next(), writer);
            }

        } else {
            // case 2: both mates have no common locus
            if (gHash.size()== 0^ gHash2.size()== 0) {
                // 2a: exclusively one mate has mappings (polyA, base caller messed up..)
                // output unpaired mappings
                if (gHash.size()> 0) {
                    i2= gHash.keySet().iterator();
                    while (i2.hasNext()) {
                        g= i2.next();
                        gList.add(g);
                        iRec= gHash.get(g).iterator();
                        while (iRec.hasNext())
                            outputMapping(iRec.next(), writer);
                    }
                } else {
                    i2= gHash2.keySet().iterator();
                    while (i2.hasNext()) {
                        g= i2.next();
                        gList.add(g);
                        iRec= gHash2.get(g).iterator();
                        while (iRec.hasNext())
                            outputMapping(iRec.next(), writer);
                    }
                }
            } else {
                // 2b: both mates have mappings, but not in a common locus
                // TODO check for fusion of neighboring loci
                System.currentTimeMillis();
            }
        }

        return gList;
    }

    /**
     * Learns systematic biases along a transcript
     *
     * @param tx   the Transcript
     * @param mappings the mappings
     */
    private void learn(Transcript tx, MSIterator<Mapping> mappings, boolean paired) {

        if (mappings== null)
            return;
        if (!mappings.hasNext())
            return;

        boolean weighted= settings.get(FluxCapacitorSettings.WEIGHTED_COUNT);
        byte strand= FluxCapacitorConstants.STRAND_ENABLED; // TODO check this, but it is the same in the capacitor
        UniversalReadDescriptor descriptor= settings.get(FluxCapacitorSettings.READ_DESCRIPTOR);

        Mapping mapping, otherMapping;
        UniversalReadDescriptor.Attributes
                attributes = descriptor.createAttributes(),
                attributes2 = descriptor.createAttributes();
        int elen = tx.getExonicLength();    // this is the "effective" length, modify by extensions
//				if (elen< readLenMin)
//					return;	// discards reads

        UniversalMatrix m = profile.getMatrix(elen);
        MappingStats stats = profile.getMappingStats();

        if (settings.get(FluxCapacitorSettings.COVERAGE_FILE) != null) {
            if (profile.getCoverageStats() == null)
                profile.setCoverageStats(new CoverageStats(elen));
            else
                profile.getCoverageStats().getCoverage().reset(elen);
        }

        while (mappings.hasNext()) {
            mapping= mappings.next();

            CharSequence tag = mapping.getName();
            attributes = descriptor.getAttributes(tag, attributes);
            if (paired) {
                if (attributes.flag < 1)
                    Log.warn("Read ignored, error in readID: " + tag);
                if (attributes.flag == 2)    // don't iterate second read
                    continue;
            }
            stats.incrReadsSingleTxLoci(1);

            // use reliable info
            if (mapping instanceof SAMMapping) {
                SAMMapping smap= (SAMMapping) mapping;
                if (!smap.isPrimary())
                    continue;
                if (paired&& (!smap.isProperlyPaired()))
                    continue;
            }

            if (strand == 1) {
                if ((tx.getStrand() == mapping.getStrand() && attributes.strand == 2)
                        || (tx.getStrand() != mapping.getStrand() && attributes.strand == 1)) {
                    stats.incrMappingsWrongStrand(1);
                    continue;
                }
            }

            int bpoint1 = getBpoint(tx, mapping);
            if (bpoint1 < 0 || bpoint1 >= elen) {    // outside tx area, or intron (Int.MIN_VALUE)
                stats.incrMappingsSingleTxLociNoAnn(1);
                continue;
            }

            stats.incrMappingsSingleTxLoci(mapping.getCount(weighted)); // the (first) read maps

            if (paired) {

//                    mappings.mark();
                Iterator<Mapping> mates = mappings.getMates(mapping,descriptor);
                while(mates.hasNext()) {
                    otherMapping= mates.next();
//                        attributes2 = settings.get(FluxCapacitorSettings.READ_DESCRIPTOR).getAttributes(bed2.getName(), attributes2);
//                        if (attributes2 == null)
//                            continue;
//                        if (!attributes.id.equals(attributes2.id))
//                            break;
//                        if (attributes2.flag == 1)    // not before break, inefficient
//                            continue;

                    // use reliable info, independent of Mate_only
                    if (otherMapping instanceof SAMMapping) {
                        SAMMapping oMap= (SAMMapping) otherMapping;
                        if (!oMap.isMateOf((SAMMapping) mapping))
                            continue;
                    }

                    int bpoint2 = getBpoint(tx, otherMapping);
                    if (bpoint2 < 0 || bpoint2 >= elen) {
                        stats.incrMappingsSingleTxLociNoAnn(1);
                        continue;
                    }

                    // check again strand in case one strand-info had been lost
                    if (strand == 1) {
                        if ((tx.getStrand() == otherMapping.getStrand() && attributes2.strand == 2)
                                || (tx.getStrand() != otherMapping.getStrand() && attributes2.strand == 1)) {
                            stats.incrMappingsWrongStrand(1);
                            continue;
                        }
                    }

                    // check directionality (sequencing-by-synthesis)
                    if ((mapping.getStrand() == otherMapping.getStrand())
                            || ((mapping.getStart() < otherMapping.getStart()) && (mapping.getStrand() != DirectedRegion.STRAND_POS))
                            || ((otherMapping.getStart() < mapping.getStart()) && (otherMapping.getStrand() != DirectedRegion.STRAND_POS))) {
                        stats.incrPairsWrongOrientation(2);
                        continue;
                    }

                    m.add(bpoint1, bpoint2, -1, -1, elen);    // 5TODO rlen currently not used
                    // update coverage
                    if (settings.get(FluxCapacitorSettings.COVERAGE_FILE) != null) {
                        if (bpoint1 < bpoint2) {
                            for (int i = bpoint1; i < bpoint1 + mapping.getLength(); i++)
                                profile.getCoverageStats().getCoverage().increment(i);
                            for (int i = bpoint2 - otherMapping.getLength() + 1; i <= bpoint2; i++)
                                profile.getCoverageStats().getCoverage().increment(i);
                        } else {
                            for (int i = bpoint2; i < bpoint2 + otherMapping.getLength(); i++)
                                profile.getCoverageStats().getCoverage().increment(i);
                            for (int i = bpoint1 - mapping.getLength() + 1; i <= bpoint1; i++)
                                profile.getCoverageStats().getCoverage().increment(i);
                        }
                    }
                    //addInsertSize(Math.abs(bpoint2- bpoint1)+ 1);	// TODO write out insert size distribution

                    //nrReadsSingleLociPairsMapped += 2;
                    stats.incrMappingPairsSingleTxLoci(mapping.getCount(weighted)+mapping.getCount(weighted));
                }
//                    mappings.reset();

            } else {    // single reads
                m.add(bpoint1, -1, elen,
                        mapping.getStrand() == tx.getStrand() ? Constants.DIR_FORWARD : Constants.DIR_BACKWARD);
                // update coverage
                if (settings.get(FluxCapacitorSettings.COVERAGE_FILE) != null) {
                    if (mapping.getStrand() == tx.getStrand()) {
                        for (int i = bpoint1; i < bpoint1 + mapping.getLength(); i++)
                            profile.getCoverageStats().getCoverage().increment(i);
                    } else {
                        for (int i = bpoint1 - mapping.getLength() + 1; i <= bpoint1; i++)
                            profile.getCoverageStats().getCoverage().increment(i);
                    }
                }
            }

        } // iterate bed objects

        // output coverage stats
        if (settings.get(FluxCapacitorSettings.COVERAGE_FILE) != null) {


            if (!profile.getCoverageStats().writeCoverageStats(getCoverageWriter(),
                    tx.getGene().getLocusID(),
                    tx.getTranscriptID(),
                    tx.isCoding(),
                    tx.getExonicLength(),
                    paired ? stats.getMappingPairsSingleTxLoci() : stats.getMappingsSingleTxLoci())){
                Log.warn("Failed to write coverage statistics to " +
                        settings.get(FluxCapacitorSettings.COVERAGE_FILE).getAbsolutePath() + barna.commons.system.OSChecker.NEW_LINE
                );
            }
        }
    }

    /**
     * Returns the breakpoint indicated by a mapping within a transcript.
     *
     * @param tx  transcript to which a read maps
     * @param bed genomic mappping
     * @return transcript coordinate of the breakpoint indicated by the mapping
     */
    private int getBpoint(Transcript tx, Mapping bed) {

        // TODO add check whether complete read is contained in transcript

        // just depends on genomic position, not on sense/antisense!
        int gpos = bed.getStrand() >= 0 ? bed.getStart() + 1 : bed.getEnd();
        int epos = tx.getExonicPosition(gpos);

        // security check, get distance between both exonic coordinates
        int epos2 = tx.getExonicPosition(bed.getStrand() >= 0 ? bed.getEnd() : bed.getStart() + 1);
        int len = bed.getLength();

        if (len != Math.abs(epos - epos2) + 1)
            return Integer.MIN_VALUE;
        return epos;
    }

    /**
     * Copied from class <code>BiasProfiler</code> to make compiler happy, marked for deletion.
     * @deprecated to be refactored or removed
     */
    private BufferedWriter coverageWriter= null;

    /**
     * Copied from class <code>BiasProfiler</code> to make compiler happy, marked for deletion.
     * @deprecated to be refactored or removed
     */
    private BufferedWriter getCoverageWriter() {
        if (coverageWriter == null) {
            File coverageFile = null;
            try {
                if (settings.get(FluxCapacitorSettings.COVERAGE_FILE) == null) {
                    coverageFile = FileHelper.createTempFile("coverage", "pro");
                    Log.warn("Temporary file for coverage stats created as " + coverageFile);
                } else {
                    coverageFile = settings.get(FluxCapacitorSettings.COVERAGE_FILE);
                }
                coverageWriter =  new BufferedWriter(new FileWriter(coverageFile));
            } catch (IOException ex) {
                throw new RuntimeException(ex);
            }
        }
        return coverageWriter;
    }

}
