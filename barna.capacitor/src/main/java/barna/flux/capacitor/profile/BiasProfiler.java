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
import com.google.gson.Gson;
import com.google.gson.GsonBuilder;

import java.io.*;
import java.util.Arrays;
import java.util.Enumeration;
import java.util.Iterator;
import java.util.Vector;
import java.util.concurrent.Callable;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;

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

    private double checkGTFscanExons;
    private double checkBEDscanMappings;

    private byte strand;
    private boolean paired;

    private double readLenMin;
    private double readLenMax;

    private GTFwrapper gtfReader;
    private MappingReader mappingReader;

    public BiasProfiler(FluxCapacitorSettings settings, MappingStats stats, byte strand, boolean paired, GTFwrapper gtfReader, MappingReader mappingReader) {
        if (settings == null) {
            throw new NullPointerException("You have to specify settings! NULL not permitted.");
        }
        this.settings = settings;
        this.strand = strand;
        this.paired = paired;

        this.gtfReader = gtfReader;
        this.mappingReader = mappingReader;

        this.stats = stats;
        profile = new BiasProfile();
    }

    @Override
    public BiasProfile call() throws Exception {
        profile();
        writeProfiles(settings.get(FluxCapacitorSettings.PROFILE_FILE),true);
        stats.writeStats(settings.get(FluxCapacitorSettings.STATS_FILE),settings.get(FluxCapacitorSettings.STATS_FILE_APPEND));
        return profile;
    }

    private void profile() {

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
                        stats.incrSingleTxLoci(1);

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
            stats.incrReadsSingleTxLoci(1);
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
                    stats.incrMappingsWrongStrand(1);
                    continue;
                }
            }

            int bpoint1 = getBpoint(tx, bed1);
            if (bpoint1 < 0 || bpoint1 >= elen) {    // outside tx area, or intron (Int.MIN_VALUE)
                stats.incrMappingsSingleTxLociNoAnn(1);
                continue;
            }

            stats.incrMappingsSingleTxLoci(1); // the (first) read maps

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
                        stats.incrMappingsSingleTxLociNoAnn(1);
                        continue;
                    }

                    // check again strand in case one strand-info had been lost
                    if (strand == 1) {
                        if ((tx.getStrand() == bed2.getStrand() && attributes2.strand == 2)
                                || (tx.getStrand() != bed2.getStrand() && attributes2.strand == 1)) {
                            stats.incrMappingsWrongStrand(1);
                            continue;
                        }
                    }

                    // check directionality (sequencing-by-synthesis)
                    if ((bed1.getStrand() == bed2.getStrand())
                            || ((bed1.getStart() < bed2.getStart()) && (bed1.getStrand() != DirectedRegion.STRAND_POS))
                            || ((bed2.getStart() < bed1.getStart()) && (bed2.getStrand() != DirectedRegion.STRAND_POS))) {
                        stats.incrPairsWrongOrientation(2);
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

                    //nrReadsSingleLociPairsMapped += 2;
                    stats.incrMappingPairsSingleTxLoci(2);
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

    public MappingStats getStats() {
        return stats;
    }

    /**
     * Writes bias profiles to disk.
     */
    public void writeProfiles(File fileProfile, boolean json) {
        try {
            final String MSG_WRITING_PROFILES = "writing profiles";

            Log.progressStart(MSG_WRITING_PROFILES);

            BufferedWriter buffy = new BufferedWriter(new FileWriter(fileProfile));
            UniversalMatrix[] mm = profile.getMasters();

            if (json) {
                Gson gson =  new GsonBuilder().serializeSpecialFloatingPointValues().create();
                String jstring = gson.toJson(mm);
                buffy.write(jstring);
            } else {
                for (int i = 0; i < mm.length; i++) {
                    String lenString = Integer.toString(mm[i].getLength());
                    buffy.write(lenString);
                    buffy.write(barna.commons.system.OSChecker.NEW_LINE);
                    buffy.write(mm[i].toStringBuilder().toString());
                }
            }
            buffy.close();
            Log.progressFinish(StringUtils.OK, true);

        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    /**
     * Read bias profiles from disk.
     */
    public static BiasProfile readProfile(File fileProfile, boolean json, MappingStats stats) {
        try {
            final String MSG_WRITING_PROFILES = "reading profiles";

            Log.progressStart(MSG_WRITING_PROFILES);

            BufferedReader buffy = new BufferedReader(new FileReader(fileProfile));
            BiasProfile profile = new BiasProfile();
            if (json) {
              Gson gson = new GsonBuilder().serializeSpecialFloatingPointValues().create();
              profile.masters = gson.fromJson(buffy,UniversalMatrix[].class);
            } else {
                UniversalMatrix[] mm = profile.getMasters();
                for (int i = 0; i < mm.length; i++) {
                    int length = Integer.parseInt(buffy.readLine());
                    if (length==mm[i].getLength()) {
                        for (int k = 0; k<3;k++) {
                            String[] row = buffy.readLine().split(",");
                            if (k==0)
                                mm[i] = new UniversalMatrix(row.length);
                            for (int j=0;j<row.length;j++) {
                                switch (k) {
                                    case 0:
                                        mm[i].sense[j] = Integer.parseInt(row[j]);
                                        mm[i].sums += mm[i].sense[j];
                                        break;
                                    case 1:
                                        mm[i].asense[j] = Integer.parseInt(row[j]);
                                        mm[1].suma += mm[i].asense[j];
                                        break;
                                }
                            }
                        }
                    } else {
                        throw new RuntimeException("Wrong profile file format");
                    }
                }
            }
            buffy.close();
            Log.progressFinish(StringUtils.OK, true);
            return profile;
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(-1);
        }
        return null;
    }

    /**
     * Reads bias profiles from the provided source file.
     *
     * @return an instance describing the bias distribution
     */
    public static BiasProfile readProfiles(File fileProfile) {

        try {
            BiasProfile profile = new BiasProfile();

            ZipFile zf = new ZipFile(fileProfile);
            Enumeration entries = zf.entries();
            String line;
            Vector<Integer> v = new Vector<Integer>();
            Vector<UniversalMatrix> w = new Vector<UniversalMatrix>();
            System.err.println("[LOAD] getting profiles");
            while (entries.hasMoreElements()) {
                ZipEntry ze = (ZipEntry) entries.nextElement();
                BufferedReader buffy = new BufferedReader(
                        new InputStreamReader(zf.getInputStream(ze)));
                int lcount = 0;
                while ((line = buffy.readLine()) != null)
                    ++lcount;
                buffy.close();
                v.add(lcount);
                UniversalMatrix m = new UniversalMatrix(lcount);
                buffy = new BufferedReader(
                        new InputStreamReader(zf.getInputStream(ze)));
                lcount = 0;
                while ((line = buffy.readLine()) != null) {
                    String[] ss = line.split("\t");
                    assert (ss.length == 2);
                    m.sense[lcount] = Integer.parseInt(ss[0]);
                    m.sums += m.sense[lcount];
                    m.asense[lcount] = Integer.parseInt(ss[1]);
                    m.suma += m.asense[lcount];
                    ++lcount;
                }
                buffy.close();
                assert (lcount == m.sense.length);
                w.add(m);
            }
            zf.close();

            int[] len = new int[v.size()];
            for (int i = 0; i < len.length; i++)
                len[i] = v.elementAt(i);
            Arrays.sort(len);
            profile.masters = new UniversalMatrix[w.size()];
            for (int i = 0; i < len.length; i++) {
                for (int j = 0; j < len.length; j++) {
                    if (len[i] == v.elementAt(j)) {
                        profile.masters[i] = w.elementAt(j);
                        // check
                        for (int n = 0; n < profile.masters[i].getLength(); n++) {
                            if (profile.masters[i].asense[n] == 0 || profile.masters[i].sense[n] == 0) {
                                if (Constants.verboseLevel > Constants.VERBOSE_SHUTUP)
                                    System.err.println("\tprofile with 0-count positions");
                                return null;
                            }
                        }
                    }
                }
            }
            System.err.println("\tfound " + profile.masters.length + " profiles.");

            return profile;

        } catch (Exception e) {
            e.printStackTrace();
            return null;
        }
    }
}
