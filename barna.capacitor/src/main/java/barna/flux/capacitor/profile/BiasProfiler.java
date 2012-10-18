package barna.flux.capacitor.profile;

import barna.commons.log.Log;
import barna.commons.utils.StringUtils;
import barna.flux.capacitor.matrix.UniversalMatrix;
import barna.flux.capacitor.reconstruction.FluxCapacitorSettings;
import barna.io.MSIterator;
import barna.io.MappingReader;
import barna.io.gtf.GTFwrapper;
import barna.io.rna.UniversalReadDescriptor;
import barna.model.DirectedRegion;
import barna.model.Gene;
import barna.model.Mapping;
import barna.model.Transcript;
import barna.model.commons.Coverage;
import barna.model.constants.Constants;

import java.util.Iterator;
import java.util.concurrent.Callable;

/**
 * A class for profiling read biases on mappings
 *
 * @author Emilio Palumbo <emiliopalumbo@gmail.com>
 */
public class BiasProfiler implements Callable<BiasProfile> {

    /**
     * Settings
     */
    private FluxCapacitorSettings settings;

    /**
     * Mappings stats
     */
    private MappingStats stats;

    /**
     * The profile
     */
    private BiasProfile profile;

    /**
     * Coverage
     */
    Coverage coverage;

    /**
     * Counters
     */
    private int nrSingleTranscriptLoci;
    private int nrReadsLoci;
    private int nrReadsMapped;
    private int nrReadsWrongLength;
    private int nrMappingsWrongStrand;
    private int nrReadsSingleLoci;
    private int nrReadsSingleLociMapped;
    private int nrReadsSingleLociPairsMapped;

    private double checkGTFscanExons;
    private double checkBEDscanMappings;

    private byte strand;
    private boolean paired;

    private long nrReadsSingleLociNoAnnotation;
    private int nrPairsWrongOrientation;
    private double readLenMin;
    private double readLenMax;

    private GTFwrapper gtfReader;
    private MappingReader mappingReader;

    public BiasProfiler(FluxCapacitorSettings settings, byte strand, boolean paired, GTFwrapper gtfReader, MappingReader mappingReader) {
        if (settings == null) {
            throw new NullPointerException("You have to specify settings! NULL not permitted.");
        }
        this.settings = settings;
        this.strand = strand;
        this.paired = paired;

        this.gtfReader = gtfReader;
        this.mappingReader = mappingReader;

        stats = new MappingStats();
        profile = new BiasProfile();
    }

    @Override
    public BiasProfile call() throws Exception {
        profile();
        return profile;
    }

    private void profile() {
        nrSingleTranscriptLoci = 0;
        nrReadsLoci = 0;
        nrReadsMapped = 0;
        nrReadsWrongLength = 0;
        nrMappingsWrongStrand = 0;
        nrReadsSingleLoci = 0;
        nrReadsSingleLociMapped = 0;

        //System.out.println(System.getProperty("java.library.path"));
        long t0 = System.currentTimeMillis();
        try {

            Transcript.setEdgeConfidenceLevel(Transcript.ID_SRC_MOST_INCONFIDENT);
            gtfReader.reset();
            mappingReader.reset();

            if (Constants.verboseLevel > Constants.VERBOSE_SHUTUP) {
                Log.info("", "");
                Log.info("PROFILE", "Scanning the input and getting the attributes.");
            }
            Log.progressStart("profiling");

            gtfReader.read();
            Gene[] gene = null, geneNext = gtfReader.getGenes();

            long tlast = System.currentTimeMillis();
            boolean output = false;

            String lastChr = null;
            byte lastStr = 0;
            int lastEnd = -1;
            int tol = 1000;

            if (geneNext != null) {
                lastChr = geneNext[0].getChromosome();
                lastStr = geneNext[0].getStrand();
            }

            Thread readerThread = null;
            int readObjects = 0;
            while (lastChr != null) {    // MAIN LOOP


                if ((gene = geneNext) == null)
                    break;

                if (readerThread == null)
                    readerThread = new Thread() {

                        @Override
                        public void run() {
                            try {
                                gtfReader.read();
                            } catch (Exception e) {
                                e.printStackTrace();
                            }
                        }
                    };
                //readerThread.start();
                readerThread.run();
                geneNext = gtfReader.getGenes();

                for (int i = 0; (gene != null) && i < gene.length; i++) {


                    // flop strand
                    if (lastChr.equals(gene[i].getChromosome())) {
                        if (lastStr != gene[i].getStrand()) {
                            //System.err.println(lastChr+" "+lastStr+ " "+ readObjects+ " wrote "+ dbgCntWriteMap +" not "+ dbgCntWriteNonmap);
                            readObjects = 0;
                            // jump back
                            mappingReader.reset(gene[i].getChromosome());
                            lastStr = gene[i].getStrand();
                            lastEnd = -1;
                        }
                    } else {                        // flop chr
                        //System.err.println(lastChr+" "+lastStr+ " "+ readObjects+ " wrote "+ dbgCntWriteMap +" not "+ dbgCntWriteNonmap);
                        readObjects = 0;
                        lastChr = gene[i].getChromosome();
                        lastStr = gene[i].getStrand();
                        lastEnd = -1;
                    }

                    if (gene[i].getTranscriptCount() == 1) {
                        ++nrSingleTranscriptLoci;

                        // boundaries
                        int start = gene[i].getStart();
                        int end = gene[i].getEnd();
                        assert (geneNext == null || geneNext.length == 1);

                        if (gene[i].getStrand() < 0) {
                            start = -start;
                            end = -end;
                        }
                        tol = 0;
                        start = Math.max(1, start - tol);
                        end = end + tol;

                        MSIterator<Mapping> mappings= mappingReader.read(gene[i].getChromosome(), start, end);

                        learn(gene[i].getTranscripts()[0], mappings);

                        if (mappings != null)
                            mappings.clear();

                        if (output) {
                            System.out.println(gene[i].getChromosome() + " " +
                                    gene[i].getStrand() +
                                    " cluster " + gene[i].getLocusID());
                            // TODO beds.size() no longer available
                            //", "+beds.size()+" reads.");
                            if ((lastStr != gene[i].getStrand()
                                    || !(lastChr.equals(gene[i].getChromosome())))) {
                                long t = System.currentTimeMillis();
                                if (lastStr != 0 && (!lastChr.equals("")))
                                    System.out.println(lastChr + " " + lastStr +
                                            " " + ((t - tlast) / 1000) + " sec.");
                                tlast = t;
                                lastStr = gene[i].getStrand();
                                lastChr = gene[i].getChromosome();
                            }
                        }
                    }
                }
                //getWriter().flush();


            }    // end iterate GTF

            //mappingReader.finish(); //TODO check

            Log.progressFinish(StringUtils.OK, true);
            if (checkGTFscanExons > 0 && checkGTFscanExons != gtfReader.getNrExons())
                System.err.println("[ERROR] consistency check failed in GTF reader: " + checkGTFscanExons + "<>" + gtfReader.getNrExons());
            checkGTFscanExons = gtfReader.getNrExons();
            if (checkBEDscanMappings > 0&& checkBEDscanMappings != mappingReader.getCountMappings())
                System.err.println("[ERROR] consistency check failed in BED reader "+ checkBEDscanMappings + "<>"+ mappingReader.getCountMappings());
            //checkBEDscanMappings= getBedReader().getNrLines();

        } catch (Exception e1) {
            Log.error("Error while iterating loci:", e1);
            throw new RuntimeException(e1);
        }
    }

    /**
     * Learns systematic biases along a transcript
     *
     * @param tx   the Transcript
     * @param mappings the mappings
     */
    private void learn(Transcript tx, MSIterator<Mapping> mappings) {

        if (mappings== null)
            return;
        if (!mappings.hasNext())
            return;

        Mapping bed1, bed2;
        UniversalReadDescriptor.Attributes
                attributes = settings.get(FluxCapacitorSettings.READ_DESCRIPTOR).createAttributes(),
                attributes2 = settings.get(FluxCapacitorSettings.READ_DESCRIPTOR).createAttributes();
        int elen = tx.getExonicLength();    // this is the "effective" length, modify by extensions
//				if (elen< readLenMin)
//					return;	// discards reads

        UniversalMatrix m = profile.getMatrix(elen);
        if (settings.get(FluxCapacitorSettings.COVERAGE_STATS)) {
            if (coverage == null)
                coverage = new Coverage(elen);
            else
                coverage.reset(elen);
        }

        while (mappings.hasNext()) {
            ++nrReadsSingleLoci;
            bed1= mappings.next();
            CharSequence tag = bed1.getName();
            attributes = settings.get(FluxCapacitorSettings.READ_DESCRIPTOR).getAttributes(tag, attributes);
            if (paired) {
                if (attributes.flag < 1)
                    Log.warn("Read ignored, error in readID: " + tag);
                if (attributes.flag == 2)    // don't iterate second read
                    continue;
            }

            if (strand == 1) {
                if ((tx.getStrand() == bed1.getStrand() && attributes.strand == 2)
                        || (tx.getStrand() != bed1.getStrand() && attributes.strand == 1)) {
                    ++nrMappingsWrongStrand;
                    continue;
                }
            }

            int bpoint1 = getBpoint(tx, bed1);
            if (bpoint1 < 0 || bpoint1 >= elen) {    // outside tx area, or intron (Int.MIN_VALUE)
                ++nrReadsSingleLociNoAnnotation;
                continue;
            }

            ++nrReadsSingleLociMapped;    // the (first) read maps

            if (paired) {

//                    mappings.mark();
                Iterator<Mapping> mates = mappings.getMates(bed1,settings.get(FluxCapacitorSettings.READ_DESCRIPTOR));
                while(mates.hasNext()) {
                    bed2= mates.next();
//                        attributes2 = settings.get(FluxCapacitorSettings.READ_DESCRIPTOR).getAttributes(bed2.getName(), attributes2);
//                        if (attributes2 == null)
//                            continue;
//                        if (!attributes.id.equals(attributes2.id))
//                            break;
//                        if (attributes2.flag == 1)    // not before break, inefficient
//                            continue;

                    int bpoint2 = getBpoint(tx, bed2);
                    if (bpoint2 < 0 || bpoint2 >= elen) {
                        ++nrReadsSingleLociNoAnnotation;
                        continue;
                    }

                    // check again strand in case one strand-info had been lost
                    if (strand == 1) {
                        if ((tx.getStrand() == bed2.getStrand() && attributes2.strand == 2)
                                || (tx.getStrand() != bed2.getStrand() && attributes2.strand == 1)) {
                            ++nrMappingsWrongStrand;
                            continue;
                        }
                    }

                    // check directionality (sequencing-by-synthesis)
                    if ((bed1.getStrand() == bed2.getStrand())
                            || ((bed1.getStart() < bed2.getStart()) && (bed1.getStrand() != DirectedRegion.STRAND_POS))
                            || ((bed2.getStart() < bed1.getStart()) && (bed2.getStrand() != DirectedRegion.STRAND_POS))) {
                        nrPairsWrongOrientation += 2;
                        continue;
                    }

                    m.add(bpoint1, bpoint2, -1, -1, elen);    // 5TODO rlen currently not used
                    // update coverage
                    if (settings.get(FluxCapacitorSettings.COVERAGE_STATS)) {
                        if (bpoint1 < bpoint2) {
                            for (int i = bpoint1; i < bpoint1 + bed1.getLength(); i++)
                                coverage.increment(i);
                            for (int i = bpoint2 - bed2.getLength() + 1; i <= bpoint2; i++)
                                coverage.increment(i);
                        } else {
                            for (int i = bpoint2; i < bpoint2 + bed2.getLength(); i++)
                                coverage.increment(i);
                            for (int i = bpoint1 - bed1.getLength() + 1; i <= bpoint1; i++)
                                coverage.increment(i);
                        }
                    }
                    //addInsertSize(Math.abs(bpoint2- bpoint1)+ 1);	// TODO write out insert size distribution

                    nrReadsSingleLociPairsMapped += 2;

                }
//                    mappings.reset();

            } else {    // single reads
                m.add(bpoint1, -1, elen,
                        bed1.getStrand() == tx.getStrand() ? Constants.DIR_FORWARD : Constants.DIR_BACKWARD);
                // update coverage
                if (settings.get(FluxCapacitorSettings.COVERAGE_STATS)) {
                    if (bed1.getStrand() == tx.getStrand()) {
                        for (int i = bpoint1; i < bpoint1 + bed1.getLength(); i++)
                            coverage.increment(i);
                    } else {
                        for (int i = bpoint1 - bed1.getLength() + 1; i <= bpoint1; i++)
                            coverage.increment(i);
                    }
                }
            }

        } // iterate bed objects
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
        if (readLenMin < 0 || len < readLenMin)
            readLenMin = len;
        if (len > readLenMax)
            readLenMax = len;

        if (len != Math.abs(epos - epos2) + 1)
            return Integer.MIN_VALUE;
        return epos;
    }
}
