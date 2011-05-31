package fbi.genome.sequencing.rnaseq.simulation.fragmentation;

import fbi.commons.ByteArrayCharSequence;
import fbi.commons.Execute;
import fbi.commons.Log;
import fbi.commons.StringUtils;
import fbi.commons.file.FileHelper;
import fbi.commons.io.IOHandler;
import fbi.commons.io.IOHandlerFactory;
import fbi.genome.io.gff.GFFReader;
import fbi.genome.model.Gene;
import fbi.genome.model.Graph;
import fbi.genome.model.Transcript;
import fbi.genome.sequencing.rnaseq.simulation.FluxSimulatorSettings;
import fbi.genome.sequencing.rnaseq.simulation.PWM;
import fbi.genome.sequencing.rnaseq.simulation.Profiler;
import fbi.genome.sequencing.rnaseq.simulation.ProfilerFile;
import fbi.genome.sequencing.rnaseq.simulation.distributions.AbstractDistribution;
import fbi.genome.sequencing.rnaseq.simulation.distributions.EmpiricalDistribution;
import fbi.genome.sequencing.rnaseq.simulation.distributions.NormalDistribution;
import org.apache.commons.math.random.RandomDataImpl;

import java.io.*;
import java.net.URL;
import java.util.*;
import java.util.concurrent.Callable;
import java.util.concurrent.Future;


public class Fragmenter implements Callable<Void> {
    public static final byte MODE_NOT_INITED = -1;
    public static final byte MODE_NONE = 0;
    public static final byte MODE_FRAG_EZ = 1;
    public static final byte MODE_NEBU = 2;
    public static final byte MODE_FRAG = 3;
    public static final byte MODE_RT = 4;
    public static final byte MODE_FILT_REJ = 5;
    public static final byte MODE_FILT_ACC = 6;
    public static final byte MODE_FILT_MH = 7;

    private static int GEL_NB_BINS_LENGTH = 100;

    private static final char FILTER_DISTRIBUTION_NORMAL = 'N';
    private static final char FILTER_DISTRIBUTION_UNIFORM = 'U';
    private static final char FILTER_DISTRIBUTION_WEIBULL = 'W';

    /**
     * The profiler
     */
    private Profiler profiler;
    private FluxSimulatorSettings settings;
    /**
     * Map from global ID to molecule count. This is later used to add mol counts to
     * the profile
     */
    private Map<ByteArrayCharSequence, Long> mapFrags;
    private Map<CharSequence, CharSequence> mapTxSeq = null;

    AbstractDistribution originalDist = null;
    AbstractDistribution[] filterDist = null;

    RandomDataImpl rndTSS;
    Random rndPA, rndPlusMinus;
    /**
     * The maximum start offset to the left. This is used to compute the
     * weight matrices if TSS MEAN is set and the start could be shifted to the left.
     */
    private int startOffset = 100;

    /**
     * If POLYA scale or shape are set, this is the maximum length of the tail
     */
    private int endOffset = 300;
    /**
     * Cache the length information
     */
    private long maxLength;


    public Fragmenter(FluxSimulatorSettings settings, Profiler profiler) {
        this.settings = settings;
        this.profiler = profiler;
    }


    public Void call() throws Exception {
        Log.message("[LIBRARY] creating the cDNA libary");

        File libraryFile = settings.get(FluxSimulatorSettings.LIB_FILE);
        if (libraryFile.exists()) {
            Log.message("[LIBRARY] Library file exists, skipping...");
            return null;
        }


        // if filtering is enabled, load the size distribution in background
        // just a minor improvement, but still
        // size selection
        Future<AbstractDistribution[]> sizeDistFuture = null;
        if (settings.get(FluxSimulatorSettings.FILTERING)) {
            File sizeDist = settings.get(FluxSimulatorSettings.SIZE_DISTRIBUTION);
            String distName = null;
            if(sizeDist != null){
                distName = sizeDist.getAbsolutePath();
            }

            final String finalDistName = distName;
            sizeDistFuture = Execute.getExecutor().submit(new Callable<AbstractDistribution[]>() {
                @Override
                public AbstractDistribution[] call() throws Exception {
                    return parseFilterDistribution(finalDistName, Double.NaN, Double.NaN, GEL_NB_BINS_LENGTH, false);
                }
            });
        }



        // now write the initial file
        // initialize random sampler
        // for initial file
        rndTSS = new RandomDataImpl();
        rndPA = new Random();
        rndPlusMinus = new Random();
        File tmpFile = writeInitialFile();
        if (tmpFile == null) {
            return null;
        }


        // do it
        // count transcripts and initialize the fragments map
        File profilerFile = settings.get(FluxSimulatorSettings.PRO_FILE);
        int nbTx = FileHelper.countLines(profilerFile);
        mapFrags = new Hashtable<ByteArrayCharSequence, Long>(nbTx, 1f);

        FluxSimulatorSettings.Substrate substrate = settings.get(FluxSimulatorSettings.FRAG_SUBSTRATE);
        FluxSimulatorSettings.FragmentationMethod fragMode = settings.get(FluxSimulatorSettings.FRAG_METHOD);

        byte mode;
        switch (fragMode){
            case NB: mode = MODE_NEBU; break;
            case EZ: mode = MODE_FRAG_EZ; break;
            case NONE: mode = MODE_NONE; break;
            default: mode = MODE_FRAG; break;
        }


        if (substrate == FluxSimulatorSettings.Substrate.RNA) {
            if (settings.get(FluxSimulatorSettings.FRAGMENTATION)) {
                if (!process(mode, tmpFile)) {
                    return null;
                }
            }
            if (settings.get(FluxSimulatorSettings.RTRANSCRIPTION)) {
                if (!process(MODE_RT, tmpFile)) {
                    return null;
                }
            }
            // start with RT
        } else {
            if (settings.get(FluxSimulatorSettings.RTRANSCRIPTION)) {
                if (!process(MODE_RT, tmpFile)) {
                    return null;
                }
            }
            if (settings.get(FluxSimulatorSettings.FRAGMENTATION)) {
                if (!process(mode, tmpFile)) {
                    return null;
                }
            }
        }


        // size selection
        if (settings.get(FluxSimulatorSettings.FILTERING)) {
            Log.message("\t\tinitializing Selected Size distribution");
            filterDist = sizeDistFuture.get();
            EmpiricalDistribution eDist = ((EmpiricalDistribution) filterDist[0]);
            if (filterDist[0] instanceof EmpiricalDistribution) {
                originalDist = parseFilterDistribution(tmpFile.getAbsolutePath(), eDist.getMin(), eDist.getMax(),eDist.getBins().length, true)[0];
                eDist.normalizeToPrior((EmpiricalDistribution) originalDist);
            }
            mode = parseFilterSampling(settings.get(FluxSimulatorSettings.SIZE_SAMPLING));
            if ((!process(mode, tmpFile))) {
                return null;
            }
        }


        if(!FileHelper.move(tmpFile, libraryFile, null)){
            throw new RuntimeException();
        }
        Log.message("\tCopied results to " + libraryFile.getAbsolutePath());

        // todo : do this in every step
        if (!ProfilerFile.appendProfile(profilerFile, ProfilerFile.PRO_COL_NR_FRG, mapFrags)){
            throw new RuntimeException();
        }
        return null;
    }


    /**
     * Reader on the GTF annotation
     *
     * @return reader the GFF reader
     */
    private GFFReader createGFFReader() {
        File ref_file = settings.get(FluxSimulatorSettings.REF_FILE);
        GFFReader gffReader = new GFFReader(ref_file.getAbsolutePath());
        gffReader.setSilent(true);
        gffReader.setStars(false);
//        if (!gffReader.isApplicable()) {
//            gffReader.close();
//            throw new RuntimeException("The reference annotation GTF is not sorted!");
//        }
        return gffReader;

    }

    /**
     * Iterate over the reads in the map and apply the pwm
     *
     * @param mapSeq the reads
     * @param pwm the pwm
     * @return map maps from read to position weights
     */
    static Map<CharSequence, double[]> getMapWeight(Map<CharSequence, CharSequence> mapTxSeq, Map<CharSequence, CharSequence> mapSeq, PWM pwm) {
        HashMap<CharSequence, double[]> map = new HashMap<CharSequence, double[]>(mapSeq.size(), 1f);
        Iterator<CharSequence> iter = mapSeq.keySet().iterator();
        while (iter.hasNext()) {
            CharSequence id = iter.next();
            CharSequence seq = mapTxSeq.get(id);
            double[] a = new double[seq.length()];
            for (int p = 0; p < a.length; ++p) {
                double pb = pwm.apply(seq, p);
                a[p] = pb;
                assert (!(Double.isNaN(pb) || Double.isInfinite(pb)));
            }
            map.put(id, a);
        }
        return map;
    }

    private Map<CharSequence, CharSequence> getMapTxSeq() {

        if (mapTxSeq == null) {
            mapTxSeq = new HashMap<CharSequence, CharSequence>(10000);
            GFFReader reader = createGFFReader();
            if (settings.get(FluxSimulatorSettings.GEN_DIR) != null) {
                Graph.overrideSequenceDirPath = settings.get(FluxSimulatorSettings.GEN_DIR).getAbsolutePath();
            }

            double tssMean = settings.get(FluxSimulatorSettings.TSS_MEAN);
            double polyaShape = settings.get(FluxSimulatorSettings.POLYA_SHAPE);
            double polyaScale = settings.get(FluxSimulatorSettings.POLYA_SCALE);

            try {
                int total = FileHelper.countLines(settings.get(FluxSimulatorSettings.REF_FILE));
                reader.read();
                Log.progressStart("preparing reads");

                int leftFlank =  Double.isNaN(tssMean) ? 0 : this.startOffset;
                int rightFlank = Double.isNaN(polyaScale) || Double.isNaN(polyaShape) ? 0 : this.endOffset;


                StringBuffer polyA = null;
                if(rightFlank >0){
                    polyA = new StringBuffer();
                    for (int i = 0; i < rightFlank; i++) {
                        polyA.append("A");
                    }
                }
                for (Gene[] g; (g = reader.getGenes()) != null; reader.read()) {
                    Log.progress(reader.getNrLinesRead(), total);
                    for (int i = 0; i < g.length; i++) {
                        for (int j = 0; j < g[i].getTranscripts().length; j++) {
                            Transcript t = g[i].getTranscripts()[j];
                            String s = t.getSplicedSequence(leftFlank, 0, "N", "A");    // TODO check chr pre-loading
                            int sourceLength = s.length();
                            if(rightFlank > 0){
                                if(s.length() < sourceLength+rightFlank){
                                    s = s+ polyA.substring(0, (sourceLength+rightFlank)-s.length());
                                }
                            }
                            ByteArrayCharSequence combID = new ByteArrayCharSequence(g[i].getGeneID());
                            combID.append((byte) FluxSimulatorSettings.SEP_LOC_TID);
                            combID.append(t.getTranscriptID());
                            mapTxSeq.put(combID, s);
                        }
                    }
                }
                Log.progressFinish(StringUtils.OK, true);
            } catch (Exception e) {
                Log.progressFailed("ERROR");
                Log.error("Error while preparing sequences: " + e.getMessage(), e);
            }
        }

        return mapTxSeq;
    }


    /**
     * Perform the fragmentation
     *
     * @param mode the fragmentation mode
     * @param tmpFile the current library file
     * @return success true if success
     */
    boolean process(byte mode, File tmpFile) {
        if(mode == MODE_NONE) return true;
        if (tmpFile == null) {
            throw new NullPointerException("No library file given !");
        }

        // reset
        long maxLen = 0;
        long cumuLen = 0;
        long currMols = 0;
        long newMols = 0;
        long totalWritten = 0;

        if (mapFrags != null) {
            mapFrags.clear();    // 20101215 re-init
        }


        // IO
        FileInputStream fis = null;
        BufferedWriter fos = null;
        IOHandler rw  = null;

        double tssMean = settings.get(FluxSimulatorSettings.TSS_MEAN);
        double polyaShape = settings.get(FluxSimulatorSettings.POLYA_SHAPE);
        double polyaScale = settings.get(FluxSimulatorSettings.POLYA_SCALE);
        //int lengthOffset = (!Double.isNaN(tssMean) || !Double.isNaN(polyaScale) || !Double.isNaN(polyaShape) ) ? 100 : 0;
        int leftFlank =  Double.isNaN(tssMean) ? 0 : this.startOffset;
        int rightFlank = Double.isNaN(polyaScale) || Double.isNaN(polyaShape) ? 0 : this.endOffset;



        try {
            for (int roundCtr = 0; roundCtr < 1; ++roundCtr) {

                fis = new FileInputStream(tmpFile);
                rw = IOHandlerFactory.createDefaultHandler();
                rw.addStream(fis);

                File tmpWriteFile = File.createTempFile("Framgmenter-write", ".tmp");
                fos = new BufferedWriter(new FileWriter(tmpWriteFile));




                ByteArrayCharSequence cs = new ByteArrayCharSequence(100);
                long byteTot = tmpFile.length();
                long byteNow = 0l;

                /// setup filter
                FragmentProcessor processor = null;
                switch(mode){
                    case MODE_FILT_REJ: processor = new FragmentFilterRejection(filterDist, true); break;
                    case MODE_FILT_ACC: processor = new FragmentFilterRejection(filterDist, false); break;
                    case MODE_FILT_MH: processor = new FragmentFilterMCMC(originalDist, filterDist); break;
                    case MODE_NEBU:
                        // determine C
                        // C~ f(lambda,M), adjust that 1.5 lambda => pb= 0.5
                        // C= lambda(1.5- (-ln(0.5))^(1/M))
                        // e.g., C=486 f. M=9, lambda=900
                        double nb_lambda = settings.get(FluxSimulatorSettings.FRAG_NB_LAMBDA);
                        Double M = settings.get(FluxSimulatorSettings.FRAG_NB_M);
                        Double thold = settings.get(FluxSimulatorSettings.FRAG_NB_THOLD);

                        processor = new FragmentNebulization(nb_lambda, M, thold, this.maxLength);
                        break;
                    case MODE_FRAG:
                        double d0 = settings.get(FluxSimulatorSettings.FRAG_UR_D0);
                        double delta = settings.get(FluxSimulatorSettings.FRAG_UR_DELTA);
                        double eta = settings.get(FluxSimulatorSettings.FRAG_UR_ETA);
                        boolean filtering = settings.get(FluxSimulatorSettings.FILTERING);
                        processor = new FragmentUniformRandom(d0, delta, eta, profiler.getMedMoleculeLength(), filtering);
                        break;
                    case MODE_FRAG_EZ:
                        File motif = settings.get(FluxSimulatorSettings.FRAG_EZ_MOTIF);
                        processor = new FragmentEnzymatic(motif, getMapTxSeq(), leftFlank, rightFlank);
                        break;
                    case MODE_RT:
                        processor = new FragmentReverseTranscription(
                                settings.get(FluxSimulatorSettings.RT_PRIMER),
                                settings.get(FluxSimulatorSettings.RT_MOTIF),
                                settings.get(FluxSimulatorSettings.RT_MIN),
                                settings.get(FluxSimulatorSettings.RT_MAX),
                                settings.get(FluxSimulatorSettings.RT_GC_LO),
                                getMapTxSeq(),
                                profiler,
                                leftFlank, rightFlank,
                                settings.get(FluxSimulatorSettings.RT_LOSSLESS)
                        );
                        break;
                }

                String name = processor.getName();
                String config = processor.getConfiguration();

                if(name != null){
                    Log.info("LIBRARY", name);
                }
                if(config != null){
                    Log.info("LIBRARY", "Configuration");
                    Log.message(config);
                }
                Log.progressStart("Processing Fragments");


                while (rw.readLine(fis, cs) > 0) {
                    // progress
                    byteNow += cs.length() + 1;
                    Log.progress(byteNow, byteTot);
                    ++currMols;

                    cs.resetFind();
                    int start = cs.getTokenInt(0);
                    int end = cs.getTokenInt(1);
                    assert (start <= end);
                    int len = end - start + 1;


                    ByteArrayCharSequence id = cs.getToken(2);
                    ByteArrayCharSequence ccs = cs.cloneCurrentSeq(); // TODO check whether needed


                    List<Fragment> fragments = processor.process(id, ccs, start, end, len);
                    if(fragments != null && fragments.size() > 0){
                        addFragCount(id, 1l);
                        newMols += fragments.size() > 1 ? fragments.size()-1 : 0;
                        for (Fragment frag : fragments) {
                            cumuLen += frag.length();
                            maxLen = Math.max(maxLen, frag.length());
                            totalWritten++;
                            // write the fragment
                            fos.write(frag.toString());
                            fos.write("\n");
                        }
                    }

                }



                if (!tmpFile.delete()) {
                    throw new IOException("Couldn't delete source");
                }
                if (!tmpWriteFile.renameTo(tmpFile)) {
                    throw new IOException("Couldn't move file");
                }


                Log.progressFinish(StringUtils.OK, true);

                String status = processor.done();
                if(status != null){
                    Log.message(status);
                }

                // sum up counts and prepare stats
                long total = currMols + newMols;
                Log.message("\t\t" + total + " mol: in " + currMols + ", new " + newMols + ", out " + totalWritten);
                Log.message("\t\tavg Len " + (cumuLen / (float) total) + ", maxLen " + maxLen);
            }
            System.gc();
            return true;
        } catch (Exception e) {
            Log.progressFailed("FAILED");
            Log.error("Error while fragmenting : " + e.getMessage(), e);
        }finally {

            this.maxLength = maxLen;
            if(fis != null){
                try {fis.close();} catch (IOException ignore) {}
            }
            if(fos != null){
                try {fos.close();} catch (IOException ignore) {}
            }

            rw.close();
        }
        return false;
    }

    private synchronized void addFragCount(ByteArrayCharSequence string, Long i) {
        ByteArrayCharSequence clone = string.cloneCurrentSeq();
        if (mapFrags.containsKey(clone)) {
            // 20101215: whyever, dont use in put clause
            // see tstCurMol and tstNewMol divergence
            long otherVal = mapFrags.get(string) + i;
            mapFrags.put(clone, otherVal);
        } else {
            mapFrags.put(clone, i);
        }

    }

    /**
     * Write and return the initial fragmentation file
     *
     * @return fragmenterFile the initial file
     */
    private File writeInitialFile() {
        if(profiler == null) throw new NullPointerException("Null Profiler not permitted in Fragmenter");
        if(profiler.size() == 0) throw new IllegalArgumentException("Profiler size is 0. Are you sure Profiling was done ?");

        int moleculesInitilized = 0;
        BufferedWriter fos = null;
        File tmpFile = null;

        double tssMean = settings.get(FluxSimulatorSettings.TSS_MEAN);
        double polyaShape = settings.get(FluxSimulatorSettings.POLYA_SHAPE);
        double polyaScale = settings.get(FluxSimulatorSettings.POLYA_SCALE);

        try {
            Log.progressStart("Initializing Fragmentation File");
            tmpFile = File.createTempFile("Framgmenter-tmp", ".tmp");
            fos = new BufferedWriter(new FileWriter(tmpFile));
            //molInit = 0;
            int profileSize = profiler.size();
            for (int i = 0; i < profileSize; i++) {
                Log.progress(i, profileSize);
                int origLen = profiler.getLength(i);
                String compIDstring = profiler.getGlobalID(i);
                long molecules = profiler.getMolecules(i);
                for (int x = 0; x < molecules; x++) {
                    ++moleculesInitilized;
                    int[] ints = processInitial(origLen, tssMean, polyaShape, polyaScale);
                    fos.write(ints[0] + "\t" + ints[1] + "\t" + compIDstring + "\n");
                }
            }
            Log.progressFinish("OK", true);
            Log.message("\t" + moleculesInitilized + " mol initialized");


            //medLen = lenSum / (double) lenNb;
            //lastMinLen = minLen;
            return tmpFile;
        } catch (Exception e) {
            Log.progressFailed(" FAILED");
            Log.error("Error while creating initial fragmentation: " + e.getMessage(), e);
            if (tmpFile != null) {
                tmpFile.delete();
            }
            return null;
        } finally {
            if (fos != null) {
                try {
                    fos.close();
                } catch (IOException ignore) {
                }
            }
        }

    }

    /**
     * Returns an error message if something is broken or missing and null if everything is fine
     *
     * @return message error message or null
     */
    public String isReady() {
        if (settings == null) {
            return "No Setting specified!";
        }
        if (profiler == null) {
            return "Profiler is not initialized!";
        }
        if (profiler.size() == 0) {
            return "Profiler has no molecules!";
        }
        if (settings.get(FluxSimulatorSettings.LIB_FILE) == null) {
            return "No fragmentation file specified!";
        }
        return null;
    }

    /**
     * Load library
     *
     * @return status true if successfully loaded
     */
    public boolean loadStats(File file) {
        if (file == null) {
            throw new NullPointerException("Null library file not permitted!");
        }
        if (!file.canRead()) {
            throw new IllegalArgumentException("Unable to read library file " + file.getAbsolutePath());
        }
        return true;
    }



    /**
     * @return a value representing the sum of the values that have been turned into a CDF
     */
    public static double toCDF(double[] a, int start, int end, int leftFlank, int rightFlank) {
        start = leftFlank + start;
        end = leftFlank + end;

        end = Math.min(a.length - 1, end);
        for (int i = start + 1; i <= end; i++) {
            a[i] += a[i - 1];
        }
        double sum = a[end];
        if (sum == 0d) {
            return sum;
        }
        for (int i = start; i <= end; i++) {
            a[i] /= sum;
            assert (!(Double.isNaN(a[i]) || Double.isInfinite(a[i]) || a[i] < 0 || a[i] > 1));
        }

        return sum;
    }

    public static void toPDF(double[] a, int start, int end, double sum, int leftflank, int rightFlank) {
        start = leftflank + start;
        end = leftflank + end;

        end = Math.min(a.length - 1, end);
        if (sum == 0) {
            return;
        }
        for (int i = start; i <= end; i++) {
            a[i] *= sum;
        }
        for (int i = end; i > start; --i) {
            a[i] -= a[i - 1];
            assert (!(Double.isNaN(a[i]) || Double.isInfinite(a[i]) || a[i] < 0/*|| a[i]> 1*/));
        }
    }

    /**
     * Initial processing of a read with given length.
     * If tssMean is a number, the start coordinate is shifted.
     * If polya shape and scale are set, the end is shifted.
     *
     *
     * @param origLen the reads length
     * @param tssMean the tss mean or NaN
     * @param polyaShape the polyA shape or NaN
     * @param polyaScale the polyA scale or NaN
     * @return startEnd array containing the start end end coordinates
     */
    int[] processInitial(int origLen, double tssMean, double polyaShape, double polyaScale) {

        // transcript variation
        int start = 0;
        int end = origLen - 1;
        if (!Double.isNaN(tssMean)) {
            start = origLen;


            while (start >= Math.min(startOffset, origLen)) {
                if(origLen <= 0) throw new RuntimeException("Transcript length <= 0?");
                if(tssMean <= 0) throw new RuntimeException("TSS mean <= 0");
                start = (int) Math.round(rndTSS.nextExponential(Math.min(tssMean, origLen / 4)));    // exp mean= 25: exceeds bounds, nextGaussian(1,100))-100;
            }
            double r = rndPlusMinus.nextDouble();
            if (r < 0.5) {
                start = -start;
            }
        }
        if (!(Double.isNaN(polyaShape) || Double.isNaN(polyaScale))) {
            int pAtail = endOffset+1; // was 301
            while (pAtail > endOffset) { // was 300
                pAtail = (int) Math.round(sampleWeibull(rndPA, polyaScale, polyaShape));    // 300d, 2d
            }
            end = origLen + pAtail;
        }
        if (end < origLen-1) {
            Log.error("end < length in Fragmenter!");
        }
        assert (start < end);
        return new int[]{start, end};
    }

    /**
     * Parses command line string and returns corresponding subsampling mode identifier.
     *
     * @param s the sampling method
     * @return sampling mode identifier
     */
    public static byte parseFilterSampling(FluxSimulatorSettings.SizeSamplingModes s) {
        switch (s) {
            case RJ:
                return MODE_FILT_REJ;
            case AC:
                return MODE_FILT_ACC;
            case MH:
                return MODE_FILT_MH;
            default:
                return MODE_NOT_INITED;
        }
    }

    /**
     * Parses the distribution argument, either a path to a file with an empirical
     * description of the function, or a set of analytically described functions.
     *
     * <p>
     *     Currently function parsing is disabled and this will
     * </p>
     *
     * @param s string describing the distribution
     * @return a (set of) distributions as initialized from the argument
     */
    public static AbstractDistribution[] parseFilterDistribution(String s, double min, double max, int nrBins, boolean fragFile) {

        if(s != null){
            File f = new File(s);
            if (f.exists()) {
                try {
                    return new AbstractDistribution[]{EmpiricalDistribution.create(f, min, max, nrBins, fragFile)};
                } catch (IOException e) {
                    throw new RuntimeException("Unable to read size filter distribution : " + e.getMessage(), e);
                }
            }else{
                throw new RuntimeException("Unable to find size distribution " + s);
            }
        }else{
            // load the default
            URL url = Fragmenter.class.getResource("/expAll.isizes");
            try {
                return new AbstractDistribution[]{EmpiricalDistribution.create(209208, url.openStream(), min, max, nrBins, fragFile)};
            } catch (IOException e) {
                throw new RuntimeException("Unable to read size filter distribution : " + e.getMessage(), e);
            }
        }

        // no file
//        String[] ss = s.split("\\+");
//        AbstractDistribution[] d = new AbstractDistribution[ss.length];
//        for (int i = 0; i < ss.length; i++) {
//            d[i] = parseFilterDistributionFunction(ss[i]);
//        }
//
//        return d;
    }


    private static AbstractDistribution parseFilterDistributionFunction(String s) {

        // find parameter block delimiters
        int p = s.indexOf('('), q = s.indexOf(')');
        assert (p > 0 && q > p);

        // parse parameters
        String[] ss = s.substring(p + 1, q).split(",");
        double[] par = new double[ss.length];
        for (int i = 0; i < par.length; i++) {
            par[i] = Double.parseDouble(ss[i]);
        }

        // nature of function
        AbstractDistribution d = null;
        if (s.charAt(p - 1) == FILTER_DISTRIBUTION_NORMAL) {
            d = (par.length == 1 ? new NormalDistribution(par[0]) : new NormalDistribution(par[0], par[1]));
        } else if (s.charAt(p - 1) == FILTER_DISTRIBUTION_UNIFORM) {
            ; // TODO
        } else if (s.charAt(p - 1) == FILTER_DISTRIBUTION_WEIBULL) {
            ; // TODO
        }

        // weight (in sums of functions)
        if (p > 1) {
            double f = Double.parseDouble(s.substring(0, p - 1));
            d.setWeight(f);
        }

        return d;
    }


    /**
     * @param rnd
     * @param lamda scale
     * @param k     shape
     * @return
     */
    public static double sampleWeibull(Random rnd, double lamda, double k) {
        double ret = lamda * Math.pow((-1d * Math.log(1d - rnd.nextDouble())), (1d / k));
        return ret;
    }

}
