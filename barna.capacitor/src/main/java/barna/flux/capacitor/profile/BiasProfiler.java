package barna.flux.capacitor.profile;

import barna.commons.log.Log;
import barna.commons.utils.StringUtils;
import barna.flux.capacitor.matrix.UniversalMatrix;
import barna.flux.capacitor.reconstruction.FluxCapacitor;
import barna.flux.capacitor.reconstruction.FluxCapacitorSettings;
import barna.io.FileHelper;
import barna.io.MSIterator;
import barna.io.MappingReader;
import barna.io.gtf.GTFwrapper;
import barna.io.rna.UniversalReadDescriptor;
import barna.model.DirectedRegion;
import barna.model.Gene;
import barna.model.Mapping;
import barna.model.Transcript;
import barna.model.constants.Constants;
import barna.model.sam.SAMMapping;
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
public class BiasProfiler implements Callable<Profile> {

    private final byte DEFAULT_CONFIDENCE = Transcript.ID_SRC_MOST_INCONFIDENT;

    /**
     * Settings
     */
    private FluxCapacitorSettings settings;

    /**
     * The profile
     */
    private Profile profile;

    /**
     * Transcript confidence level
     */
    private byte edgeConfidence = -1;

    private double checkGTFscanExons;
    private double checkBEDscanMappings;

    private byte strand;
    private boolean paired;
    private boolean weighted;

    private int readLenMin = Integer.MAX_VALUE;
    private int readLenMax = -1;

    private GTFwrapper gtfReader;
    private MappingReader mappingReader;
    private BufferedWriter coverageWriter;

    /**
     * Temporary file for coverage statistics of the 5' to 3' read distribution.
     */
    File fileTmpCovStats = null;

    /**
     * Writer of the coverage statistics of the 5' to 3' read distribution.
     */
    private BufferedWriter writerTmpCovStats = null;
    /**
     * The capacitor instance
     */
    private final FluxCapacitor capacitor;

    public BiasProfiler(FluxCapacitor capacitor, byte strand, boolean paired, boolean weighted, GTFwrapper gtfReader, MappingReader mappingReader) {
        if (capacitor == null) {
            throw new NullPointerException("You have to specify settings! NULL not permitted.");
        }
        this.settings = capacitor.getSettings();
        this.capacitor = capacitor;
        this.strand = strand;
        this.paired = paired;
        this.weighted = weighted;

        this.gtfReader = gtfReader;
        this.mappingReader = mappingReader;

        profile = new Profile();
    }

    @Override
    public Profile call() throws Exception {
        profile();
        profile.getMappingStats().setReadLenMin(readLenMin);
        profile.getMappingStats().setReadLenMax(readLenMax);
        if (settings.get(FluxCapacitorSettings.PROFILE_FILE)!=null)
            writeProfiles(settings.get(FluxCapacitorSettings.PROFILE_FILE),true);
        return profile;
    }

    private void profile() {

        //System.out.println(System.getProperty("java.library.path"));
        long t0 = System.currentTimeMillis();
        try {

            Transcript.setEdgeConfidenceLevel(edgeConfidence == -1 ? DEFAULT_CONFIDENCE : edgeConfidence);
            gtfReader.reset();
            mappingReader.reset();

            if (Constants.verboseLevel > Constants.VERBOSE_SHUTUP) {
                Log.info("", "");
                Log.info("PROFILE", "Scanning the input and getting the attributes.");
                if (settings.get(FluxCapacitorSettings.COVERAGE_FILE) != null){
                    Log.info("PROFILE", "Coverage statistics in " + settings.get(FluxCapacitorSettings.COVERAGE_FILE).getAbsolutePath());
                }
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
                        profile.getMappingStats().incrSingleTxLoci(1);

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
            if (coverageWriter !=null) {
                coverageWriter.close();
                coverageWriter = null;
            }

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

        Mapping mapping, otherMapping;
        UniversalReadDescriptor.Attributes
                attributes = settings.getReadDescriptor().createAttributes(),
                attributes2 = settings.getReadDescriptor().createAttributes();
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
            attributes = settings.getReadDescriptor().getAttributes(tag, attributes);
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
                Iterator<Mapping> mates = mappings.getMates(mapping,settings.getReadDescriptor());
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
                String jstring = gson.toJson(profile);
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
    public static Profile readProfile(File fileProfile, boolean json) {
        try {
            final String MSG_WRITING_PROFILES = "reading profiles";

            Log.progressStart(MSG_WRITING_PROFILES);

            BufferedReader buffy = new BufferedReader(new FileReader(fileProfile));
            Profile profile = new Profile();
            if (json) {
              Gson gson = new GsonBuilder().serializeSpecialFloatingPointValues().create();
              profile = gson.fromJson(buffy,Profile.class);
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
            profile.setMasters(new UniversalMatrix[w.size()]);
            for (int i = 0; i < len.length; i++) {
                for (int j = 0; j < len.length; j++) {
                    if (len[i] == v.elementAt(j)) {
                        profile.getMasters()[i] = w.elementAt(j);
                        // check
                        for (int n = 0; n < profile.getMasters()[i].getLength(); n++) {
                            if (profile.getMasters()[i].asense[n] == 0 || profile.getMasters()[i].sense[n] == 0) {
                                if (Constants.verboseLevel > Constants.VERBOSE_SHUTUP)
                                    System.err.println("\tprofile with 0-count positions");
                                return null;
                            }
                        }
                    }
                }
            }
            System.err.println("\tfound " + profile.getMasters().length + " profiles.");

            return profile;

        } catch (Exception e) {
            e.printStackTrace();
            return null;
        }
    }

    public byte getEdgeConfidence() {
        return edgeConfidence;
    }

    public void setEdgeConfidence(byte edgeConfidence) {
        this.edgeConfidence = edgeConfidence;
    }

}
