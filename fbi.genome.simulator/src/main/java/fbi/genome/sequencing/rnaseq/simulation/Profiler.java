package fbi.genome.sequencing.rnaseq.simulation;

import fbi.commons.ByteArrayCharSequence;
import fbi.commons.Log;
import fbi.commons.StringUtils;
import fbi.commons.file.FileHelper;
import fbi.commons.file.ReverseFileReader;
import fbi.genome.io.ThreadedBufferedByteArrayStream;
import fbi.genome.io.gff.GFFReader;
import fbi.genome.model.Gene;
import fbi.genome.model.commons.Distribution;
import fbi.genome.model.commons.IntVector;

import java.io.*;
import java.util.*;
import java.util.concurrent.Callable;

public class Profiler implements Callable<Void> {

    public static final byte STAT_NONE = 0, STAT_ANN = 4, STAT_RELFREQ = 5, STAT_MOL = 6;
    FluxSimulatorSettings settings;
    byte status = -1;
    ByteArrayCharSequence[] ids = null, locIDs;
    int[] len = null;
    long[] molecules = null;
    boolean[] cds = null;
    int cntLoci = -1;
    float txLocAvg = -1, txLocMed = -1, txLoc1Q = -1, txLoc3Q = -1, txLocSTD = -1, lenMed = -1, lenAvg = -1, lenSTD = -1, len1Q = -1, len3Q = -1, lenMin = -1, lenMax = -1;
    HashSet<CharSequence> sfHi, sfMed, sfLo;
    Hashtable<ByteArrayCharSequence, int[]> mapLenExp;
    long sumMol = 0;

    /**
     * Create a new new profiler
     *
     * @param settings the settings
     */
    public Profiler(FluxSimulatorSettings settings) {
        if(settings == null ) throw new NullPointerException("You have to specify settings! NULL not permitted.");
        this.settings = settings;
    }


    public Void call() throws Exception {
        Log.info("PROFILING", "I am assigning the expression profile");
        status = getStatus();
        if (status == STAT_NONE) {
            readAnnotation();
            status = STAT_ANN;
            writeProfile();
        }

        // write some info
        Log.info("PROFILING", "Parameters");
        Log.message("\t" + settings.toString(FluxSimulatorSettings.NB_MOLECULES));
        Log.message("\t" + settings.toString(FluxSimulatorSettings.EXPRESSION_K));
        Log.message("\t" + settings.toString(FluxSimulatorSettings.EXPRESSION_X0));
        Log.message("\t" + settings.toString(FluxSimulatorSettings.EXPRESSION_X1));
        Log.message("\t" + settings.toString(FluxSimulatorSettings.PRO_FILE));
        Log.message("");

        profile();
        Log.message("\tmolecules\t" + sumMol);
        Log.message("");
        return null;
    }

    public boolean isFinishedReadAnnotation() {
        return (locIDs != null && ids != null && ids.length > 0 && cds != null && len != null &&
                ids.length == len.length && ids.length == locIDs.length && ids.length == cds.length);
    }

    public boolean isFinishedExpression() {
        return (isFinishedReadAnnotation() && molecules != null && ids.length == molecules.length);
    }

    byte getStatus() {

        if (settings.get(FluxSimulatorSettings.PRO_FILE) == null || !settings.get(FluxSimulatorSettings.PRO_FILE).exists()) {
            return STAT_NONE;
        }

        try {
            ReverseFileReader rreader = new ReverseFileReader(settings.get(FluxSimulatorSettings.PRO_FILE).getCanonicalPath());
            String s = rreader.readLine();
            if (s == null) {
                return STAT_NONE;
            }

            String[] tokens = s.split("\\s");
            return (byte) tokens.length;

        } catch (Exception e) {
            //e.printStackTrace();
            return STAT_NONE;
        }

    }


    private void readAnnotation() throws Exception {
        GFFReader reader = createGTFReader();

        Log.progressStart("Reading reference annotation");

        reader.read();

        List<ByteArrayCharSequence> v = new ArrayList<ByteArrayCharSequence>(30000);
        List<ByteArrayCharSequence> vLoc = new ArrayList<ByteArrayCharSequence>(30000);
        List<Boolean> vBoo = new ArrayList<Boolean>(30000);
        IntVector lenV = new IntVector(30000); /*, txLoc= new IntVector(20000);*/

        mapLenExp = new Hashtable<ByteArrayCharSequence, int[]>();
        boolean loadCoding = settings.get(FluxSimulatorSettings.LOAD_CODING);
        boolean loadNonCoding = settings.get(FluxSimulatorSettings.LOAD_NONCODING);
        for (Gene[] g; (g = reader.getGenes()) != null; reader.read()) {
            for (int i = 0; i < g.length; i++) {
                ++cntLoci;
                for (int j = 0; j < g[i].getTranscripts().length; j++) {


                    if (g[i].getTranscripts()[j].isCoding() && (!loadCoding)) {
                        continue;
                    }

                    if ((!g[i].getTranscripts()[j].isCoding()) && (!loadNonCoding)) {
                        continue;
                    }
                    if (g[i].getTranscripts()[j].isCoding()) {
                        vBoo.add(true);
                    } else {
                        vBoo.add(false);
                    }
                    v.add(new ByteArrayCharSequence(g[i].getTranscripts()[j].getTranscriptID()));
                    ByteArrayCharSequence locName = new ByteArrayCharSequence(g[i].getGeneID());
                    vLoc.add(locName);
                    int[] a = new int[2];
                    a[0] = g[i].getTranscripts()[j].getExonicLength();
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
        calcStats();

        Log.message("\tfound " + ids.length + " transcripts\n");

        v = null;
        vLoc = null;
        vBoo = null;
        lenV = null;
        System.gc();

    }

    private void writeProfile() throws IOException {
        BufferedWriter writer = null;
        try {
            writer = new BufferedWriter(new FileWriter(settings.get(FluxSimulatorSettings.PRO_FILE)));
            for (int i = 0; i < ids.length; i++) {
                writer.write(locIDs[i] + ProfilerFile.PRO_FILE_SEP + ids[i] + ProfilerFile.PRO_FILE_SEP + (cds[i] ? ProfilerFile.PRO_FILE_CDS : ProfilerFile.PRO_FILE_NC)
                        + ProfilerFile.PRO_FILE_SEP + Integer.toString(len[i]) + "\n");
            }
            writer.flush();
            writer.close();
        }finally {
            if(writer != null){
                try {writer.close();} catch (IOException e) {}
            }
        }
    }

    public int getLength(ByteArrayCharSequence id) {
        return mapLenExp.get(id)[0];
    }


    private void calcStats() {

        CharSequence lastLocID = null;
        int asCtr = 1;
        IntVector asV = new IntVector();
        cntLoci = 0;
        for (int i = 0; i < locIDs.length; i++) {
            if ((!locIDs[i].equals(lastLocID))) {
                ++cntLoci;
                if (lastLocID != null) {
                    asV.add(asCtr);
                    asCtr = 1;
                }
                lastLocID = locIDs[i];
            } else {
                ++asCtr;
            }
        }

        int[] as = asV.toIntArray();
        Arrays.sort(as);
        Distribution dist = new Distribution(as);
        txLocMed = (float) dist.getMedian();
        txLocAvg = (float) dist.getMean();
        txLocSTD = (float) dist.getStandardDeviation();
        txLoc1Q = (float) dist.get1stQuart();
        txLoc3Q = (float) dist.get3rdQuart();
        as = null;
        dist = new Distribution(len.clone());
        System.gc();

        lenMin = (float) dist.getMin();
        lenMax = (float) dist.getMax();
        lenAvg = (float) dist.getMean();
        lenMed = (float) dist.getMedian();
        len1Q = (float) dist.get1stQuart();
        len3Q = (float) dist.get3rdQuart();
        lenSTD = (float) dist.getStandardDeviation();
    }

    /**
     * Create a GTF reader and validate the input file, eventually sort it
     *
     * @return reader ready to use reader on a sorted valid GTF file
     */
    private GFFReader createGTFReader() {
        File currentRefFile = settings.get(FluxSimulatorSettings.REF_FILE);
        GFFReader gffReader = new GFFReader(currentRefFile.getAbsolutePath());
        // make sure the gtf is valid and sorted
        if (!gffReader.isApplicable()) {
            File refFile = gffReader.createSortedFile();
            File target = new File(settings.get(FluxSimulatorSettings.PRO_FILE).getParent() + File.separator + refFile.getName());
            settings.setRefFile(target);
            if (!refFile.equals(currentRefFile)) {
                if (!FileHelper.move(refFile, target)) {
                    settings.setRefFile(refFile);
                }
            }
            currentRefFile = settings.get(FluxSimulatorSettings.REF_FILE);
            gffReader = new GFFReader(currentRefFile.getAbsolutePath());
        }
        gffReader.setSilent(true);
        gffReader.setStars(true);
        return gffReader;
    }



    public boolean profile() {

        Log.progressStart("profiling");

        try {
            if (ids == null) {
                ids = new ByteArrayCharSequence[FileHelper.countLines(settings.get(FluxSimulatorSettings.PRO_FILE).getCanonicalPath())];
                len = new int[ids.length];
            }
            double sumRF = 0;
            sumMol = 0;
            sfHi = null;
            sfMed = null;
            sfLo = null;
            molecules = new long[ids.length];
            double[] relFreq = new double[ids.length];

            if (status < STAT_RELFREQ) {    // generate ranks

                // generate random permutation of ranks
                Random r = new Random();
                for (int i = 0; i < molecules.length; i++) {
                    molecules[i] = 1 + r.nextInt(molecules.length - 1); //i+1;	// ranks
                }
//				for (int k = molecules.length - 1; k > 0; k--) {
//				    int w = (int) Math.floor(r.nextDouble() * (k+1));
//				    long temp = molecules[w];
//				    molecules[w] = molecules[k];
//				    molecules[k] = temp;
//				}

                relFreq = new double[ids.length];
                for (int i = 0; i < relFreq.length; i++) {
                    double par = pareto(molecules[i], settings.get(FluxSimulatorSettings.EXPRESSION_K), settings.get(FluxSimulatorSettings.EXPRESSION_X0));
                    double exp = exponential(molecules[i], settings.get(FluxSimulatorSettings.EXPRESSION_X1));
                    sumRF += (relFreq[i] = par * exp);
                }
                for (int i = 0; i < relFreq.length; ++i) {    // normalize
                    relFreq[i] /= sumRF;
                    sumMol += (molecules[i] = Math.round(relFreq[i] * settings.get(FluxSimulatorSettings.NB_MOLECULES)));
                    if (molecules[i] > 0) {
                        mapLenExp.put(getCombinedID(i), new int[]{len[i], (int) molecules[i]});
                    }
                    Log.progress(i, relFreq.length);
                }
            }


            Log.progressFinish(StringUtils.OK, false);
            Hashtable<CharSequence, Long> map = new Hashtable<CharSequence, Long>(getMolecules().length);
            for (int i = 0; i < getMolecules().length; i++) {
                if (getMolecules()[i] != 0) {
                    ByteArrayCharSequence locNtid = getLocIDs()[i].cloneCurrentSeq();
                    locNtid.append(Character.toString(FluxSimulatorSettings.SEP_LOC_TID));
                    locNtid.append(getIds()[i]);
                    map.put(locNtid, getMolecules()[i]);
                }
            }
            if (!ProfilerFile.appendProfile(settings, ProfilerFile.PRO_COL_NR_MOL, map)) {
                return false;
            }
            return true;

        } catch (Exception e) {
            Log.progressFailed("FAILED");
            Log.error("Error while profiling :" + e.getMessage(), e);
            return false;
        }
    }

    public static double exponential(double rank, double par1) {
        //		double val= Math.exp(- (Math.pow(rank, 2)/ Math.pow(par1, 2))
        //				- (rank/par1));
        double val = Math.exp(-(Math.pow(rank / par1, 2))    // / 122000000
                - (rank / par1));    // 7000

        return val;
    }

    public static double pareto(double rank, double par1, double par2) {
        //		double val= par2/ Math.pow(rank, par1)
        //			* Math.exp(- (Math.pow(rank, 2)/ 122000000)
        //					- (rank/7000));
        double val = Math.pow(rank / par2, par1)/* 2731598d*/;    // par1= 0,6  par2= (41627d/ 2731598d)
        return val;
    }


    public ByteArrayCharSequence[] getIds() {
        return ids;
    }

    public int[] getLen() {
        return len;
    }

    public long[] getMolecules() {
        return molecules;
    }

    public long getSumMol() {
        return sumMol;
    }


    public float getTxLocAvg() {
        return txLocAvg;
    }


    public void setTxLocAvg(float txLocAvg) {
        this.txLocAvg = txLocAvg;
    }


    public float getTxLocMed() {
        return txLocMed;
    }


    public void setTxLocMed(float txLocMed) {
        this.txLocMed = txLocMed;
    }


    public float getTxLoc1Q() {
        return txLoc1Q;
    }


    public void setTxLoc1Q(float txLoc1Q) {
        this.txLoc1Q = txLoc1Q;
    }


    public float getTxLoc3Q() {
        return txLoc3Q;
    }


    public void setTxLoc3Q(float txLoc3Q) {
        this.txLoc3Q = txLoc3Q;
    }


    public float getTxLocSTD() {
        return txLocSTD;
    }


    public void setTxLocSTD(float txLocSTD) {
        this.txLocSTD = txLocSTD;
    }


    public float getLenMed() {
        return lenMed;
    }


    public float getLenAvg() {
        return lenAvg;
    }


    public float getLenSTD() {
        return lenSTD;
    }


    public float getLen1Q() {
        return len1Q;
    }


    public float getLen3Q() {
        return len3Q;
    }


    public float getLenMin() {
        return lenMin;
    }


    public float getLenMax() {
        return lenMax;
    }


    public int getCntLoci() {
        return cntLoci;
    }


    public boolean loadStats() {
        if (settings.get(FluxSimulatorSettings.PRO_FILE) == null || (!settings.get(FluxSimulatorSettings.PRO_FILE).exists())) {
            return false;
        }

        int lim = Integer.MAX_VALUE;    // last working token
        try {
            Log.progressStart("initializing profiler ");
            int lines = FileHelper.countLines(settings.get(FluxSimulatorSettings.PRO_FILE).getAbsolutePath());
            //String lineSep= FileHelper.getLineSeparator() // TODO
            ids = new ByteArrayCharSequence[lines];
            locIDs = new ByteArrayCharSequence[lines];
            len = new int[lines];
            molecules = new long[lines];
            cds = new boolean[lines];
            HashMap<ByteArrayCharSequence, ByteArrayCharSequence> locIDset =
                    new HashMap<ByteArrayCharSequence, ByteArrayCharSequence>();    // TODO MyHashSet.get(Object o)
            int ptr = -1, perc = 0;
            long bytesRead = 0, bytesTot = settings.get(FluxSimulatorSettings.PRO_FILE).length();
            mapLenExp = new Hashtable<ByteArrayCharSequence, int[]>();

            BufferedInputStream istream = new BufferedInputStream(new FileInputStream(settings.get(FluxSimulatorSettings.PRO_FILE)));
            ThreadedBufferedByteArrayStream buffy =
                    new ThreadedBufferedByteArrayStream(10 * 1024, istream, true, false);
            ByteArrayCharSequence cs = new ByteArrayCharSequence(1024);
            for (buffy.readLine(cs); cs.end > 0; buffy.readLine(cs)) {

                cs.resetFind();
                bytesRead += cs.length() + 1;    // TODO fs
                if (ptr % 1000 == 0) {
                    Log.progress(bytesRead, bytesTot);
                }
                if (lim < 3) {
                    break;    // give up
                }

                ++ptr;
                int tok = 0;
                ByteArrayCharSequence x = cs.getToken(tok++);
                if (x == null) {
                    lim = tok - 2;
                    continue;
                } else {
                    if (locIDset.containsKey(x)) {
                        locIDs[ptr] = locIDset.get(x);
                    } else {
                        locIDs[ptr] = x.cloneCurrentSeq();
                        locIDset.put(locIDs[ptr], locIDs[ptr]);
                    }
                }

                if (lim < tok) {
                    continue;
                }
                x = cs.getToken(tok++);
                if (x == null) {
                    lim = tok - 2;
                    continue;
                } else {
                    ids[ptr] = x.cloneCurrentSeq();
                }

                if (lim < tok) {
                    continue;
                }
                x = cs.getToken(tok++);
                if (x == null) {
                    lim = tok - 2;
                    continue;
                } else {
                    if (x.equals(ProfilerFile.PRO_FILE_CDS)) {
                        cds[ptr] = true;
                    } else if (x.equals(ProfilerFile.PRO_FILE_NC)) {
                        cds[ptr] = false;
                    } else {
                        lim = 1;
                        continue;
                    }
                }

                if (lim < tok) {
                    continue;
                }
                int y = cs.getTokenInt(tok++);
                if (y == Integer.MIN_VALUE) {
                    lim = tok - 2;
                    continue;
                } else {
                    len[ptr] = y;
                }

                tok++;    // perc

                if (lim < tok) {
                    continue;
                }
                y = cs.getTokenInt(tok++);
                if (y == Integer.MIN_VALUE) {
                    lim = tok - 3;
                } else {
                    molecules[ptr] = y;
                    if (molecules[ptr] > 0) {
                        mapLenExp.put(getCombinedID(ptr),    // Issue32: ids[ptr]
                                new int[]{len[ptr], (int) molecules[ptr]});
                    }
                    lim = tok - 1;
                }
            }
            istream.close();
            buffy.close();
        } catch (Exception e) {
            lim = -1; // :)
            Log.progressFailed("ERROR");
            Log.error("Error while loading stats: " + e.getMessage(), e);
            return false;
        }

        if (lim < 4) {
            molecules = null;
        }
        if (lim < 3) {
            len = null;
        }
        if (lim < 2) {
            ids = null;
            locIDs = null;
        } else {
            calcStats();
        }


        Log.progressFinish();
        return true;
    }


    public ByteArrayCharSequence[] getLocIDs() {
        return locIDs;
    }


    public void setLocIDs(ByteArrayCharSequence[] locIDs) {
        this.locIDs = locIDs;
    }


    public boolean[] getCds() {
        return cds;
    }


    public void setCds(boolean[] cds) {
        this.cds = cds;
    }


    public void setMolecules(long[] molecules) {
        this.molecules = molecules;
    }


    public ByteArrayCharSequence getCombinedID(int i) {
        if (i < 0 || i >= getIds().length) {
            return null;
        }
        ByteArrayCharSequence locID = getLocIDs()[i];
        ByteArrayCharSequence tID = getIds()[i];

        ByteArrayCharSequence cs = new ByteArrayCharSequence(locID.length() + tID.length() + 1);
        System.arraycopy(locID.a, locID.start, cs.a, cs.end, locID.length());
        cs.end += locID.length();
        cs.a[cs.end++] = FluxSimulatorSettings.SEP_LOC_TID;
        System.arraycopy(tID.a, tID.start, cs.a, cs.end, tID.length());
        cs.end += tID.length();

        return cs;
    }


    public int getMaxMoleculeLength() {
        int maxLen = -1;
        for (int i = 0; i < len.length; i++) {
            if (molecules[i] > 0 && len[i] > maxLen) {
                maxLen = len[i];
            }
        }
        return maxLen;
    }


    public double getMedMoleculeLength() {

        int sumMol = 0;
        for (int i = 0; i < molecules.length; i++) {
            sumMol += molecules[i];
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
}
