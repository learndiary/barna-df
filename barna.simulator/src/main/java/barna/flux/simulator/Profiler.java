/*
 * This file is part of the Flux Library.
 *
 * The code of the Flux Library may be freely distributed and modified under the terms of the
 * European Union Public Licence (EUPL) published on the web site <http://www.osor.eu/eupl/european-union-public-licence-eupl-v.1.1>.
 * Copyright for the code is held jointly by the individual authors, who should be listed
 * in @author doc comments. According to Article 5 and Article 11 of the EUPL, publications that
 * include results produced by the Flux Library are liable to reference the Work,
 * see the Flux Library homepage <http://flux.sammeth.net> for more information.
 */

package barna.flux.simulator;

import barna.commons.ByteArrayCharSequence;
import barna.commons.io.ReverseFileReader;
import barna.commons.log.Log;
import barna.commons.utils.StringUtils;
import barna.io.FileHelper;
import barna.io.ThreadedBufferedByteArrayStream;
import barna.io.gtf.GTFwrapper;
import barna.model.Gene;
import barna.model.commons.Distribution;
import barna.model.commons.IntVector;

import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.Callable;


/**
 * Simulator profiler. Manages and creates the expression profile
 */
public class Profiler implements Callable<Void> {

    /**
     * Status: no profile available
     */
    private static final byte STAT_NONE = 0;
    /**
     * Status: annotations available
     */
    private static final byte STAT_ANN = 4;
    /**
     * Status: relative frequencies available
     */
    private static final byte STAT_RELFREQ = 5;

    /**
     * The setting
     */
    private FluxSimulatorSettings settings;
    /**
     * Current state
     */
    private byte status = -1;
    /**
     * Transcript ids, i.e. NM_009826
     */
    private ByteArrayCharSequence[] ids = null;
    /**
     * Locus ids based on chromosome and position, i.e., chr1:6204743-6266185W
     */
    private ByteArrayCharSequence[] locIDs;
    /**
     * Lengths
     */
    private int[] len = null;
    /**
     * Molecules
     */
    private long[] molecules = null;
    /**
     * Cds
     */
    private boolean[] cds = null;

    /**
     * Stores a mapping from the global ID (chromosome+position+transcriptID) to an array
     * where a[0] contains the length and int[1] contains the absolute number of molecules
     */
    private Map<ByteArrayCharSequence, int[]> mapLenExp;

    /**
     * Create a new new profiler
     *
     * @param settings the settings
     */
    public Profiler(FluxSimulatorSettings settings) {
        if (settings == null) {
            throw new NullPointerException("You have to specify settings! NULL not permitted.");
        }
        this.settings = settings;
    }

    /**
     * Start profiling
     *
     * @return null always returns null
     * @throws Exception in case of any errors
     */
    public Void call() throws Exception {
        Log.info("PROFILING", "I am assigning the expression profile");
        status = readStatus();

        if (status == STAT_NONE || !isFinishedReadAnnotation()) {
            // read annotation and write initial profile file without expression
            readAnnotation();
            status = STAT_ANN;
            ProfilerFile.writeProfile(this, settings.get(FluxSimulatorSettings.PRO_FILE));
        }

        // write some info
        Log.info("PROFILING", "Parameters");
        Log.message("\t" + settings.toString(FluxSimulatorSettings.NB_MOLECULES));
        Log.message("\t" + settings.toString(FluxSimulatorSettings.EXPRESSION_K));
        Log.message("\t" + settings.toString(FluxSimulatorSettings.EXPRESSION_X0));
        Log.message("\t" + settings.toString(FluxSimulatorSettings.EXPRESSION_X1));
        Log.message("\t" + settings.toString(FluxSimulatorSettings.PRO_FILE));
        Log.message("");

        long sumMol = profile();
        Log.message("\tmolecules\t" + sumMol);
        Log.message("");
        return null;
    }

    /**
     * Returns true if the profile contains transcripts with length information, loci ID and transcript id
     *
     * @return annotations true if profile contains transcripts
     */
    public boolean isFinishedReadAnnotation() {
        return (locIDs != null && ids != null && ids.length > 0 && cds != null && len != null &&
                ids.length == len.length && ids.length == locIDs.length && ids.length == cds.length);
    }

    /**
     * Returns true if transcripts and an expression profile exists
     *
     * @return expression true if expression profile exists
     */
    public boolean isFinishedExpression() {
        return (isFinishedReadAnnotation() && molecules != null && ids.length == molecules.length);
    }

    /**
     * Read the profiler file and find out the status based on the number of columns
     *
     * @return status the current status
     */
    byte readStatus() {

        File profileFile = settings.get(FluxSimulatorSettings.PRO_FILE);
        if (profileFile == null || !profileFile.exists()) {
            return STAT_NONE;
        }

        try {
            ReverseFileReader rreader = new ReverseFileReader(profileFile.getCanonicalPath());
            String s = rreader.readLine();
            rreader.close();

            if (s == null) {
                return STAT_NONE;
            }

            String[] tokens = s.split("\\s");
            return (byte) tokens.length;
        } catch (Exception e) {
            // ignore the exception
            return STAT_NONE;
        }

    }

    /**
     * Read transcript information from GTF annotations
     *
     * @throws Exception in case of any errors
     */
    private void readAnnotation() throws Exception {
        GTFwrapper reader = createGTFReader();

        Log.progressStart("Reading reference annotation");

        reader.read();

        List<ByteArrayCharSequence> v = new ArrayList<ByteArrayCharSequence>(30000);
        List<ByteArrayCharSequence> vLoc = new ArrayList<ByteArrayCharSequence>(30000);
        List<Boolean> vBoo = new ArrayList<Boolean>(30000);
        IntVector lenV = new IntVector(30000); /*, txLoc= new IntVector(20000);*/

        mapLenExp = new Hashtable<ByteArrayCharSequence, int[]>();
        boolean loadCoding = settings.get(FluxSimulatorSettings.LOAD_CODING);
        boolean loadNonCoding = settings.get(FluxSimulatorSettings.LOAD_NONCODING);
        long totalBytes = settings.get(FluxSimulatorSettings.REF_FILE).length();
        for (Gene[] g; (g = reader.getGenes()) != null; reader.read()) {
            Log.progress(reader.getBytesRead(), totalBytes);
            for (Gene aG : g) {
                for (int j = 0; j < aG.getTranscripts().length; j++) {


                    if (aG.getTranscripts()[j].isCoding() && (!loadCoding)) {
                        continue;
                    }

                    if ((!aG.getTranscripts()[j].isCoding()) && (!loadNonCoding)) {
                        continue;
                    }
                    if (aG.getTranscripts()[j].isCoding()) {
                        vBoo.add(true);
                    } else {
                        vBoo.add(false);
                    }
                    String transcriptID = aG.getTranscripts()[j].getTranscriptID();
                    v.add(new ByteArrayCharSequence(transcriptID));
                    ByteArrayCharSequence locName = new ByteArrayCharSequence(aG.getGeneID());
                    vLoc.add(locName);
                    int[] a = new int[2];
                    a[0] = aG.getTranscripts()[j].getExonicLength();
                    lenV.add(a[0]);
                }
            }
        }

        ids = new ByteArrayCharSequence[v.size()];
        for (int i = 0; i < v.size(); i++) {
            ids[i] = v.get(i);
        }

        locIDs = new ByteArrayCharSequence[vLoc.size()];
        for (int i = 0; i < vLoc.size(); i++) {
            locIDs[i] = vLoc.get(i);
        }

        cds = new boolean[vBoo.size()];
        for (int i = 0; i < vBoo.size(); i++) {
            cds[i] = vBoo.get(i);
        }


        len = lenV.toIntArray();

        Log.progressFinish(StringUtils.OK, true);
        Log.message("\tfound " + ids.length + " transcripts\n");

        v = null;
        vLoc = null;
        vBoo = null;
        lenV = null;
        System.gc();

    }


    /**
     * Get the transcript length distribution
     *
     * @return lengthDistribution length distribution
     */
    public Distribution getLengthDistribution() {
        return new Distribution(len.clone());
    }

    /**
     * Get the distributions of transcripts per loci
     *
     * @return transcriptDistribution
     */
    public Distribution getTranscriptDistribution() {
        if (locIDs == null) {
            throw new RuntimeException("No profile loaded!");
        }
        CharSequence lastLocID = null;
        int asCtr = 1;
        IntVector asV = new IntVector();
        for (ByteArrayCharSequence locID : locIDs) {
            if ((!locID.equals(lastLocID))) {
                if (lastLocID != null) {
                    asV.add(asCtr);
                    asCtr = 1;
                }
                lastLocID = locID;
            } else {
                ++asCtr;
            }
        }
        return new Distribution(asV.toIntArray());
    }


    /**
     * Do the profiling. This assumes that annotations are read and exist
     *
     * @return molecules number of expressed molecules
     */
    protected long profile() {
        Log.progressStart("profiling");
        try {
            if (ids == null) {
                long lines = FileHelper.countLines(settings.get(FluxSimulatorSettings.PRO_FILE).getCanonicalPath());
                if(lines > Integer.MAX_VALUE) throw new RuntimeException("Unable to cache " + lines + " elements... value > Integer.MAX_VALUE");
                ids = new ByteArrayCharSequence[(int) lines];
                len = new int[ids.length];
            }
            double sumRF = 0;
            long sumMol = 0;
            molecules = new long[ids.length];
            double[] relFreq = new double[ids.length];

            if (status < STAT_RELFREQ) {
                double expressionK = settings.get(FluxSimulatorSettings.EXPRESSION_K);
                double expression_x0 = settings.get(FluxSimulatorSettings.EXPRESSION_X0);
                double expression_x1 = settings.get(FluxSimulatorSettings.EXPRESSION_X1);
                long nb_molecules = settings.get(FluxSimulatorSettings.NB_MOLECULES);
                // generate random permutation of ranks
                Random r = new Random();
                for (int i = 0; i < molecules.length; i++) {
                    molecules[i] = 1 + r.nextInt(molecules.length - 1);
                }
                /*
                Expressions
                 */
                for (int i = 0; i < relFreq.length; i++) {
                    double par = pareto(molecules[i], expressionK, expression_x0);
                    double exp = exponential(molecules[i], expression_x1);
                    sumRF += (relFreq[i] = par * exp);
                }
                /*
                Normalize
                 */
                for (int i = 0; i < relFreq.length; ++i) {
                    relFreq[i] /= sumRF;
                    sumMol += (molecules[i] = Math.round(relFreq[i] * nb_molecules));
                    if (molecules[i] > 0) {
                        mapLenExp.put(getCombinedID(i), new int[]{len[i], (int) molecules[i]});
                    }
                    Log.progress(i, relFreq.length);
                }
            }
            Log.progressFinish(StringUtils.OK, false);


            /*
            Append expression information to profile file
             */
            Map<CharSequence, Number> map = new HashMap<CharSequence, Number>(size(), 1f);
            for (int i = 0; i < size(); i++) {
                if (molecules[i] != 0) {
                    ByteArrayCharSequence locNtid = locIDs[i].cloneCurrentSeq();
                    locNtid.append(Character.toString(FluxSimulatorSettings.SEP_LOC_TID));
                    locNtid.append(ids[i]);
                    map.put(locNtid, molecules[i]);
                }
            }
            if (!ProfilerFile.appendProfile(settings.get(FluxSimulatorSettings.PRO_FILE), ProfilerFile.PRO_COL_NR_MOL, map, true)) {
                throw new RuntimeException("Unable to append data to profile file!"); 
            }
            return sumMol;
        } catch (Exception e) {
            Log.progressFailed("FAILED");
            throw new RuntimeException("Error while profiling :" + e.getMessage(), e);
        }
    }

    /**
     * Returns the number of transcripts in the profiler. You
     * can use this, for example, to iterate over all entries.
     *
     * @return size the number of transcripts
     */
    public int size() {
        return ids != null ? ids.length : 0;
    }

    /**
     * Get the transcript IDs at the i'the position
     *
     * @param i the entry position
     * @return id transcript id at the i'th position
     */
    public ByteArrayCharSequence getId(int i) {
        return ids[i];
    }

    /**
     * The the transcript lengths of the i'th entry
     *
     * @param i entry index
     * @return lengths of the i'th entry
     */
    public int getLength(int i) {
        return len[i];
    }

    /**
     * Returns the length of the transcript based on the given global ID
     *
     * @param id the global ID
     * @return length transcript length
     */
    public int getLength(ByteArrayCharSequence id) {
        return mapLenExp.get(id)[0];
    }


    /**
     * Get the number of molecules for the i'th transcript
     *
     * @param i the entry
     * @return molecules number of molecules for the i'th transcript
     */
    public long getMolecules(int i) {
        return molecules[i];
    }

    /**
     * Returns true if the entry at the i'th position is of type CDS
     *
     * @param i the entry
     * @return cds true if CDS
     */
    public boolean isCds(int i) {
        return cds[i];
    }

    /**
     * Load profile from existing profiler file
     *
     * @param profileFile the profiler file
     * @return success true if profile loaded successfully, false otherwise
     */
    public boolean initializeProfiler(File profileFile) {
        if (profileFile == null || (!profileFile.exists()) || !profileFile.canRead()) {
            return false;
        }

        int lim = Integer.MAX_VALUE;    // last working token
        BufferedInputStream istream = null;
        ThreadedBufferedByteArrayStream buffy = null;
        try {
            Log.progressStart("initializing profiler ");
            long ll = FileHelper.countLines(profileFile);
            if(ll > Integer.MAX_VALUE) throw new RuntimeException("Unable to cache " + ll + " elements... value > Integer.MAX_VALUE");
            int lines = (int) ll;
            int separatorLength = Math.min(1, FileHelper.guessFileSep(profileFile).length());

            /*
            Initialize data structures
             */
            ids = new ByteArrayCharSequence[lines];
            locIDs = new ByteArrayCharSequence[lines];
            len = new int[lines];
            molecules = new long[lines];
            cds = new boolean[lines];
            /**
             * Cache to identify duplicated IDs and always use the
             * same object for them
             */
            Map<ByteArrayCharSequence, ByteArrayCharSequence> locIDset = new HashMap<ByteArrayCharSequence, ByteArrayCharSequence>();
            int ptr = -1;
            long bytesRead = 0, bytesTot = profileFile.length();
            mapLenExp = new Hashtable<ByteArrayCharSequence, int[]>();

            // read
            istream = new BufferedInputStream(new FileInputStream(profileFile));
            buffy = new ThreadedBufferedByteArrayStream(10 * 1024, istream, true, false);

            // cache
            ByteArrayCharSequence cs = new ByteArrayCharSequence(1024);
            for (buffy.readLine(cs); cs.end > 0; buffy.readLine(cs)) {
                // reset cache find
                cs.resetFind();
                bytesRead += cs.length() + separatorLength;

                // report progress every 1000 lines
                if (ptr % 1000 == 0) {
                    Log.progress(bytesRead, bytesTot);
                }
                // increase the line counter / array pointer
                ++ptr;


                /*
                 Read Token 0:
                   chromosome and position, i.e.,   chr1:4797974-4836816W
                */
                ByteArrayCharSequence x = cs.getToken(0);
                if (x == null) {
                    // no ID found
                    lim = -1;
                    break;
                } else {
                    // check if we have found this already
                    // if so, use existing id instance
                    // else add a new one
                    ByteArrayCharSequence id = locIDset.get(x);
                    if (id != null) {
                        locIDs[ptr] = id;
                    } else {
                        // create a new object, note that x is only a view on
                        // the char[], we have to clone here to get a copy that
                        // stays and only represents the id
                        locIDs[ptr] = x.cloneCurrentSeq();
                        locIDset.put(locIDs[ptr], locIDs[ptr]);
                    }
                }
                /*
                Read Token 1:
                    transcript ID, i.e., YAL056W
                 */
                x = cs.getToken(1);
                if (x == null) {
                    // no transcript ID found !
                    lim = 0;
                    break;
                } else {
                    // store the transcript ID
                    ids[ptr] = x.cloneCurrentSeq();
                }

                /*
                 Read Token 3
                    Type, either CDS or NC
                 */
                x = cs.getToken(2);
                if (x == null) {
                    // no type defined
                    lim = 1;
                    break;
                } else {
                    // check type
                    if (x.equals(ProfilerFile.PRO_FILE_CDS)) {
                        cds[ptr] = true;
                    } else if (x.equals(ProfilerFile.PRO_FILE_NC)) {
                        cds[ptr] = false;
                    } else {
                        lim = -1;
                        break;
                    }
                }

                /*
                read Token 3
                    Length
                 */
                int y = cs.getTokenInt(3);
                if (y == Integer.MIN_VALUE) {
                    lim = 2;
                    continue;
                } else {
                    len[ptr] = y;
                }

                if (lim < 3) {
                    continue;
                }

                // skip the relative number and
                // read the absolute number of molecules
                y = cs.getTokenInt(5);
                if (y == Integer.MIN_VALUE) {
                    // not found
                    lim = 3;
                } else {
                    molecules[ptr] = y;
                    if (y > 0) {
                        mapLenExp.put(getCombinedID(ptr),    // Issue32: ids[ptr]
                                new int[]{len[ptr], y});
                    }
                }
            }
        } catch (Exception e) {
            lim = -1;
            Log.progressFailed("ERROR");
            Log.error("Error while loading stats: " + e.getMessage(), e);
            return false;
        } finally {
            if (istream != null) {
                try {
                    istream.close();
                } catch (IOException ignore) {
                }
            }
            if (buffy != null) {
                buffy.close();
            }

        }

        // forget everything that
        // we do not have for all the entries
        if (lim < 4) {
            molecules = null;
            mapLenExp.clear();
            mapLenExp = null;
        }
        if (lim < 3) {
            len = null;
        }
        if (lim < 2) {
            ids = null;
            locIDs = null;
        }
        Log.progressFinish();
        return true;
    }

    /**
     * Access the loci ID ath the i'th position
     *
     * @param i the entry
     * @return lociId loci id at the i'th position
     */
    public ByteArrayCharSequence getLociId(int i) {
        return locIDs[i];
    }

    /**
     * Create the global ID for entry i. Global ID consists of
     * chromosome, position, and the transcript ID
     *
     * @param i the entry index
     * @return globalID the global ID
     */
    public ByteArrayCharSequence getCombinedID(int i) {
        if (i < 0 || i >= size()) {
            return null;
        }
        ByteArrayCharSequence locID = locIDs[i];
        ByteArrayCharSequence tID = ids[i];

        ByteArrayCharSequence cs = new ByteArrayCharSequence(locID.length() + tID.length() + 1);
        System.arraycopy(locID.chars, locID.start, cs.chars, cs.end, locID.length());
        cs.end += locID.length();
        cs.chars[cs.end++] = FluxSimulatorSettings.SEP_LOC_TID;
        System.arraycopy(tID.chars, tID.start, cs.chars, cs.end, tID.length());
        cs.end += tID.length();

        return cs;
    }

    /**
     * Create the global ID for entry i. Global ID consists of
     * chromosome, position, and the transcript ID.
     *
     * @param i the entry index
     * @return globalID the global ID
     */
    public String getGlobalID(int i) {
        if (i < 0 || i >= size()) {
            return null;
        }

        StringBuilder bb = new StringBuilder(getLociId(i));
        bb.append(FluxSimulatorSettings.SEP_LOC_TID);
        bb.append(getId(i));
        return bb.toString();
    }


    /**
     * Compute the maximal molecule length. Note that this is NOT cached
     * and iterates over all molecules
     *
     * @return maxLength max molecule length
     */
    public int getMaxMoleculeLength() {
        int maxLen = -1;
        for (int i = 0; i < len.length; i++) {
            if (molecules[i] > 0 && len[i] > maxLen) {
                maxLen = len[i];
            }
        }
        return maxLen;
    }

    /**
     * Compute the median molecule length. Note that the value is not cached
     * and computed each time.
     *
     * @return medianLength median molecule length
     */
    public double getMedMoleculeLength() {
        int sumMol = 0;
        for (long molecule : molecules) {
            sumMol += molecule;
        }
        IntVector v = new IntVector(sumMol);
        for (int i = 0; i < molecules.length; i++) {
            for (int j = 0; j < molecules[i]; j++) {
                v.add(len[i]);
            }
        }

        Distribution dist = new Distribution(v.toIntArray());
        return dist.getMedian();
    }

    /**
     * Delete existing profiles and reset status
     */
    public void resetProfile() {
        ids = null;
        molecules = null;
        locIDs = null;
        cds = null;
        if (mapLenExp != null) {
            mapLenExp.clear();
        }
        mapLenExp = null;
        len = null;

        File profilerFile = settings.get(FluxSimulatorSettings.PRO_FILE);
        if (profilerFile.exists()) {
            boolean b = profilerFile.delete();
            if (!b) {
                Log.error("PROFILER", "Unable to delete profile!");
            }
            status = STAT_NONE;
        }
    }

    /**
     * Create a GTF reader and validate the input file, eventually sort it
     *
     * @return reader ready to use reader on a sorted valid GTF file
     */
    private GTFwrapper createGTFReader() {
        File currentRefFile = settings.get(FluxSimulatorSettings.REF_FILE);
        GTFwrapper gffReader = new GTFwrapper(currentRefFile.getAbsolutePath());
        // make sure the gtf is valid and sorted
        if (!gffReader.isApplicable()) {
            gffReader.close();
            throw new RuntimeException("The reference annotation GTF is not sorted!");
        }
        gffReader.setSilent(true);
        gffReader.setStars(true);
        return gffReader;
    }


    private static double exponential(double rank, double par1) {
        return Math.exp(-(Math.pow(rank / par1, 2)) - (rank / par1));
    }

    private static double pareto(double rank, double par1, double par2) {
        return Math.pow(rank / par2, par1);
    }

}
