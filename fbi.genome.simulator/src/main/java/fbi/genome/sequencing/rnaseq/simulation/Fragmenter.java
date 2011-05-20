package fbi.genome.sequencing.rnaseq.simulation;

import fbi.commons.ByteArrayCharSequence;
import fbi.commons.Log;
import fbi.commons.file.FileHelper;
import fbi.commons.thread.SyncIOHandler2;
import fbi.genome.io.BufferedByteArrayReader;
import fbi.genome.io.gff.GFFReader;
import fbi.genome.model.Gene;
import fbi.genome.model.Graph;
import fbi.genome.model.Transcript;
import fbi.genome.model.constants.Constants;
import genome.sequencing.rnaseq.simulation.distributions.AbstractDistribution;
import genome.sequencing.rnaseq.simulation.distributions.EmpiricalDistribution;
import genome.sequencing.rnaseq.simulation.distributions.NormalDistribution;
import org.apache.commons.math.random.RandomDataImpl;
import org.apache.commons.math.special.Gamma;

import java.io.*;
import java.util.*;
import java.util.concurrent.Callable;


public class Fragmenter implements Callable<Void> {
    public static final byte MODE_NOT_INITED = -1;
    public static final byte MODE_FRAG_EZ = 1;
    public static final byte MODE_NEBU = 2;
    public static final byte MODE_FRAG = 3;
    public static final byte MODE_RT = 4;
    public static final byte MODE_FILT_REJ = 5;
    public static final byte MODE_FILT_ACC = 6;
    public static final byte MODE_FILT_MH = 7;
    private static final String MODE_RT_MESSAGE = "Reverse Transcription";
    private static final String MODE_NEBU_MESSAGE = "Nebulization";
    private static final String MODE_FILT_MESSAGE = "Segregating cDNA";
    private static final String MODE_FRAG_MESSAGE = "Fragmentation UR";
    private static final String MODE_FRAG_EZ_MESSAGE = "Enzymatic Digestion";
    /**
     * Default median size after fragmentation, according to 2010 Illumina protocol.
     */
    private static final double DEFAULT_MED_SIZE = 200;


    // 2.85 - <0.5%
    private static final double CUT_OFF_GAUSSIAN_VAL = 2.85f;
    // gel bins
    private static int GEL_NB_BINS_LENGTH = 100;

    private static final char FILTER_DISTRIBUTION_NORMAL = 'N';
    private static final char FILTER_DISTRIBUTION_UNIFORM = 'U';
    private static final char FILTER_DISTRIBUTION_WEIBULL = 'W';

    /**
     * The profiler
     */
    private Profiler profiler;
    private FluxSimulatorSettings settings;
    private Map<ByteArrayCharSequence, Long> mapFrags;
    private Map<CharSequence, double[]> mapWeightSense = null;
    private Map<CharSequence, double[]> mapWeightAsense = null;

    private PWM pwmSense;
    private PWM pwmAsense;
    private Map<CharSequence, CharSequence> mapTxSeq = null;



    Random rndBreak;
    Random rndBP = new Random();
    Random rndGel = new Random();
    Random rnd1 = new Random();
    Random rnd2 = new Random();
    Random rnd3 = new Random();

    //GaussianRndThread rndNebu;
    RandomDataImpl rdiNebuBP;
    long newMols = 0;
    long currMols = 0;
    long tgtMols = 0;
    long tstNewMols = 0;
    long tstCurrMols = 0;    // un-synchronized, but should be ok -- are estimates anyway

    /**
     * Parameter to adapt breaking probability distribution to 0.5 for (1.5*lambda).
     */
    double nebuC = -1;

    /**
     * Number of Iterations for recursive nebulization.
     */
    int nebuRecursionDepth = -1;
    // tmporary variables re-used (nebulization, RT)
    int[] index1 = null;
    Random rtRndWhere = new Random();

    double rtC = Double.NaN;

    SyncIOHandler2 rw;
    long cumuLen;
    int[] med = new int[3];
    AbstractDistribution originalDist = null;
    AbstractDistribution[] filterDist = null;

    FileOutputStream fos;

    private int maxLen;
    RandomDataImpl rndTSS;
    Random rndPA, rndPlusMinus;
    int cntMolInit = 0;

    int lastLen = -1; //MCMC
    double lastP = -1;




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
     * @param s string describing the distribution
     * @return a (set of) distributions as initialized from the argument
     */
    public static AbstractDistribution[] parseFilterDistribution(String s, double min, double max, int nrBins, boolean fragFile) {

        File f = new File(s);
        if (f.exists()) {
            try {
                return new AbstractDistribution[]{EmpiricalDistribution.create(f, min, max, nrBins, fragFile)};
            } catch (Exception e) {
                e.printStackTrace();
                return null;
            }
        }

        // no file
        String[] ss = s.split("\\+");
        AbstractDistribution[] d = new AbstractDistribution[ss.length];
        for (int i = 0; i < ss.length; i++) {
            d[i] = parseFilterDistributionFunction(ss[i]);
        }

        return d;
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
     * min<= r <= max
     *
     * @param random
     * @param min
     * @param max
     * @return
     */
    //public static RandomDataImpl rndINebu= new RandomDataImpl();
    public static strictfp double nextGaussianDouble(Random random, double min, double max) {
        double rdm = 3d;    // gaussian value, stddev 1
        while (rdm < -CUT_OFF_GAUSSIAN_VAL || rdm > CUT_OFF_GAUSSIAN_VAL) {
            rdm = random.nextGaussian();
        }
        //rdm= rndINebu.nextGaussian(0d, 1d);
        double mid = ((double) min) + (max - min) / 2f;
        double realValue = mid + (rdm * (max - mid) / CUT_OFF_GAUSSIAN_VAL);    // 0..CUT_OFF_GAUSSIAN = mid..max
        assert (realValue >= min && realValue <= max);
        return realValue;
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


    public Fragmenter(FluxSimulatorSettings settings, Profiler profiler) {
        this.settings = settings;
        this.profiler = profiler;
    }


    public Void call() throws Exception {
        Log.message("[LIBRARY] creating the cDNA libary");

        int nbTx = FileHelper.countLines(settings.get(FluxSimulatorSettings.PRO_FILE).getCanonicalPath());
        mapFrags = new Hashtable<ByteArrayCharSequence, Long>(nbTx, 1f);
        if (settings.get(FluxSimulatorSettings.LIB_FILE).exists()) {
            Log.message("[LIBRARY] Library file exists, skipping...");
            return null;
        }

        File tmpFile = writeInitialFile();
        if (tmpFile == null) {
            return null;
        }


        // do it
        double gc_lo = settings.get(FluxSimulatorSettings.RT_GC_LO);
        double gc_high = settings.get(FluxSimulatorSettings.RT_GC_HI);
        FluxSimulatorSettings.Substrate substrate = settings.get(FluxSimulatorSettings.FRAG_SUBSTRATE);


        FluxSimulatorSettings.FragmentationMethod fragMode = settings.get(FluxSimulatorSettings.FRAG_METHOD);
        if (substrate == FluxSimulatorSettings.Substrate.RNA) {
            rndBreak = new Random();
            if ((boolean) settings.get(FluxSimulatorSettings.FRAGMENTATION)) {
                byte mode = MODE_NOT_INITED;
                if (fragMode == FluxSimulatorSettings.FragmentationMethod.NB) {
                    mode = MODE_NEBU;
                    rdiNebuBP = new RandomDataImpl();
                } else if (fragMode == FluxSimulatorSettings.FragmentationMethod.UR) {
                    mode = MODE_FRAG;
                    rdiNebuBP = new RandomDataImpl();
                    rndBP = new Random();
                }
                if (!process(mode, tmpFile)) {
                    return null;
                }
            }
            if ((boolean) settings.get(FluxSimulatorSettings.RTRANSCRIPTION)) {
                initRTpar(gc_lo, gc_high);
                if (settings.get(FluxSimulatorSettings.RT_MOTIF) != null) {
                    getMapTxSeq();
                    //getWeights(settings.getRTmotif());
                    try {
                        pwmSense = PWM.create(settings.get(FluxSimulatorSettings.RT_MOTIF));
                        for (int i = 0; i < 100; i++) {
                            pwmSense.multiply();
                        }
                        pwmSense.makePDF();
                        mapWeightSense = getMapWeight(getMapTxSeq(), pwmSense);
                        pwmSense.invert();
                        pwmAsense = pwmSense;
                        mapWeightAsense = getMapWeight(getMapTxSeq(), pwmAsense);
                    } catch (Exception e) {
                        e.printStackTrace();
                    }

                }
                if (!process(MODE_RT, tmpFile)) {
                    return null;
                }
            }

            // start with RT
        } else {
            if (settings.get(FluxSimulatorSettings.RTRANSCRIPTION)) {
                initRTpar(gc_lo, gc_high);
                if (settings.get(FluxSimulatorSettings.RT_MOTIF) != null) {
                    getMapTxSeq();
                    getWeights(settings.get(FluxSimulatorSettings.RT_MOTIF));
                }
                if (!process(MODE_RT, tmpFile)) {
                    return null;
                }
            }
            rndBreak = new Random();
            if ((boolean) settings.get(FluxSimulatorSettings.FRAGMENTATION)) {
                byte mode = MODE_NOT_INITED;
                if (fragMode == FluxSimulatorSettings.FragmentationMethod.NB) {
                    mode = MODE_NEBU;
                    rdiNebuBP = new RandomDataImpl();
                } else if (fragMode == FluxSimulatorSettings.FragmentationMethod.UR) {
                    mode = MODE_FRAG;
                    rndBP = new Random();
                } else if (fragMode == FluxSimulatorSettings.FragmentationMethod.EZ) {
                    mode = MODE_FRAG_EZ;
                    rnd1 = new Random();
                    rnd2 = new Random();
                    rnd3 = new Random();
                    try {
                        getMapTxSeq();
                        pwmSense = PWM.create(settings.get(FluxSimulatorSettings.FRAG_EZ_MOTIF));
                        mapWeightSense = getMapWeight(getMapTxSeq(), pwmSense);
                        pwmSense.invert();
                        pwmAsense = pwmSense;
                        mapWeightAsense = getMapWeight(getMapTxSeq(), pwmAsense);
                    } catch (Exception e) {
                        e.printStackTrace();
                    }
                }
                if (!process(mode, tmpFile)) {
                    return null;
                }
            }
        }


        // size selection
        if ((boolean) settings.get(FluxSimulatorSettings.FILTERING)) {

            //if (settings.getFileFilterDistr()!= null) {
            //	initGelBins(new File(settings.getFilterDistribution()));
            //}
            System.err.println("Initializing Selected Size distribution");
            //Log.progressStart("Initializing Selected Size distribution");
            filterDist = parseFilterDistribution(settings.get(FluxSimulatorSettings.SIZE_DISTRIBUTION).getAbsolutePath(), Double.NaN, Double.NaN, GEL_NB_BINS_LENGTH, false);
            //Log.progressFinish();
            EmpiricalDistribution eDist = ((EmpiricalDistribution) filterDist[0]);

            if (filterDist[0] instanceof EmpiricalDistribution) {

                System.err.println("Initializing Current Size distribution");
                //Log.progressStart("Initializing Current Size distribution");	// TODO to be done during earlier steps
                originalDist = parseFilterDistribution(tmpFile.getAbsolutePath(), eDist.getMin(), eDist.getMax(),
                        eDist.getBins().length, true)[0];
                // Log.progressFinish();

                eDist.normalizeToPrior((EmpiricalDistribution) originalDist);
            }

            byte mode = parseFilterSampling(settings.get(FluxSimulatorSettings.SIZE_SAMPLING));
            if ((!process(mode, tmpFile))) {
                return null;
            }
        }

        if (!writeFinalFile(tmpFile)) {
            throw new RuntimeException();
        }
        File proFile = settings.get(FluxSimulatorSettings.PRO_FILE);
        if (!ProfilerFile.appendProfile(proFile, ProfilerFile.PRO_COL_NR_FRG, mapFrags))    // writeProFile()
        {
            throw new RuntimeException();
        }
        return null;
    }

    private GFFReader getGFFReader() {
        File ref_file = settings.get(FluxSimulatorSettings.REF_FILE);
        GFFReader gffReader = new GFFReader(ref_file.getAbsolutePath());
        try {
            if (!gffReader.isApplicable()) {
                File refFile = gffReader.createSortedFile();
                if (refFile == null) {
                    return null;
                }
                settings.setRefFile(new File(settings.get(FluxSimulatorSettings.PRO_FILE).getParent() + File.separator + refFile.getName()));
                if (!refFile.equals(ref_file)) {
                    if (!FileHelper.move(refFile, ref_file, null)) {
                        settings.setRefFile(refFile);
                    }
                }
                gffReader = new GFFReader(ref_file.getAbsolutePath());
            }
            gffReader.setSilent(true);
            gffReader.setStars(true);

        } catch (Exception e) {
            return null;
        }

        return gffReader;

    }


    private Map<CharSequence, double[]> getMapWeight(Map<CharSequence, CharSequence> mapSeq, PWM pwm) {
        HashMap<CharSequence, double[]> map = new HashMap<CharSequence, double[]>(mapSeq.size(), 1f);
        Iterator<CharSequence> iter = mapSeq.keySet().iterator();
        while (iter.hasNext()) {
            CharSequence id = iter.next();
            CharSequence seq = getMapTxSeq().get(id);
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
            GFFReader reader = getGFFReader();
            if (settings.get(FluxSimulatorSettings.GEN_DIR) != null) {
                Graph.overrideSequenceDirPath = settings.get(FluxSimulatorSettings.GEN_DIR).getAbsolutePath();
            }
            try {
                reader.read();
                for (Gene[] g; (g = reader.getGenes()) != null; reader.read()) {
                    for (int i = 0; i < g.length; i++) {
                        for (int j = 0; j < g[i].getTranscripts().length; j++) {
                            Transcript t = g[i].getTranscripts()[j];
                            String s = t.getSplicedSequence();    // TODO check chr pre-loading
                            ByteArrayCharSequence combID = new ByteArrayCharSequence(g[i].getGeneID());
                            combID.append((byte) FluxSimulatorSettings.SEP_LOC_TID);
                            combID.append(t.getTranscriptID());
                            mapTxSeq.put(combID, s);
                        }
                    }
                }
            } catch (Exception e) {
                e.printStackTrace();
            }
        }

        return mapTxSeq;
    }

    private void getWeights(File filePWM) {
        try {
            pwmSense = PWM.create(filePWM);
            mapWeightSense = getMapWeight(getMapTxSeq(), pwmSense);
            pwmSense.invert();
            pwmAsense = pwmSense;
            mapWeightAsense = getMapWeight(getMapTxSeq(), pwmAsense);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }


    private boolean writeFinalFile(File tmpFile) {
        Log.message("Copying results");
        if (FileHelper.move(tmpFile, settings.get(FluxSimulatorSettings.LIB_FILE), null)) {
            return true;
        }
        return false;
    }


    boolean process(byte mode, File tmpFile) {
        if (tmpFile == null) {
            throw new NullPointerException();
        }
        try {
            long t0 = System.currentTimeMillis();
            double breakRatio = 1;
            Long longNull = new Long(0);
            double thrTgt = 0.95;
            double tgtFrac = 0d;
            cumuLen = 0;
            for (int i = 0; i < profiler.size(); i++) {
                long molecules = profiler.getMolecules(i);
                cumuLen += molecules * profiler.getLength(i);
            }

            if (mode == MODE_NEBU) {
                initNebulizationParameters(
                        settings.get(FluxSimulatorSettings.FRAG_NB_LAMBDA),
                        settings.get(FluxSimulatorSettings.FRAG_NB_M),
                        settings.get(FluxSimulatorSettings.FRAG_NB_THOLD),
                        profiler.getMaxMoleculeLength()
                );
            }

            // TODO estimate required disk space
            for (int roundCtr = 0; roundCtr < 1; ++roundCtr) {


                //lengthV= new IntVector();
                maxLen = 0;
                cumuLen = 0;
                maxLen = 0;
                if (mapFrags != null) {
                    mapFrags.clear();    // 20101215 re-init
                }

                String msg = null;
                if (mode == MODE_RT) {
                    msg = MODE_RT_MESSAGE;
                } else if (mode == MODE_NEBU) {
                    msg = MODE_NEBU_MESSAGE + "-" + Integer.toString(roundCtr + 1);
                } else if (mode == MODE_FILT_REJ || mode == MODE_FILT_ACC || mode == MODE_FILT_MH) {
                    msg = MODE_FILT_MESSAGE;
                } else if (mode == MODE_FRAG) {
                    msg = MODE_FRAG_MESSAGE + "-" + Integer.toString(roundCtr + 1);
                } else if (mode == MODE_FRAG_EZ) {
                    msg = MODE_FRAG_EZ_MESSAGE;
                }

                Log.progressStart(msg);

                Object[] keys = mapFrags.keySet().toArray();
                for (int i = 0; i < keys.length; i++) {
                    mapFrags.put((ByteArrayCharSequence) keys[i], longNull);    // 20101205: changed from cast to string, ClassCastException for ByteArrayCharSequence
                }

                currMols = 0;
                newMols = 0;
                tgtMols = 0;
                tstNewMols = 0;
                tstCurrMols = 0;

                // reset medians
                Arrays.fill(med, 0);

                // IO 
//				BufferedReader buffy= new BufferedReader(new FileReader(tmpFile));
//				ThreadedQWriter qwriter= getQWriter();
//				qwriter.init();
//				qwriter.start();
//				ThreadedBufferedByteArrayStream buffy= 
//					new ThreadedBufferedByteArrayStream(10* 1024* 1024, fis, true, false);
                //BufferedOutputStream outStream
                //writer= new BufferedOutputStream(new FileOutputStream(tmpWriteFile), 1024* 1024);
                //ThreadedBufferedByteArrayStream inBacs= null;
                FileInputStream fis = new FileInputStream(tmpFile);
                File tmpWriteFile = File.createTempFile("Framgmenter-write", ".tmp");
                fos = new FileOutputStream(tmpWriteFile);
                rw = new SyncIOHandler2(2);
                rw.addStream(fis, 10 * 1024);
                rw.addStream(fos, 10 * 1024);

                ByteArrayCharSequence cs = new ByteArrayCharSequence(100);

                int perc = 0;
                long byteTot = tmpFile.length(), byteNow = 0l;
                //for (String s; (!isStop())&& (s= buffy.readLine())!= null;/*++currMols*/) {
                //for (buffy.readLine(cs); cs.end> 0; buffy.readLine(cs)) {
                while (rw.readLine(fis, cs) > 0) {
                    byteNow += cs.length() + 1;    // TODO fs length
                    Log.progress(byteNow, byteTot);

                    ++currMols;

                    cs.resetFind();
                    int start = cs.getTokenInt(0);
                    int end = cs.getTokenInt(1);
                    assert (start <= end);
                    int len = end - start + 1;
//					if (processFragDiss(len))
//						return;
                    if (len <= settings.get(FluxSimulatorSettings.FRAG_NB_LAMBDA)) {
                        ++tgtMols;
                    }

                    ByteArrayCharSequence id = cs.getToken(2);
                    addFragCount(id, 1l);
                    tstCurrMols += 1;

                    ByteArrayCharSequence ccs = cs.cloneCurrentSeq(); // TODO check whether needed
                    if (mode == Fragmenter.MODE_FILT_REJ) {
                        processFilterRejection(ccs, start, end, len, id, filterDist, true);
                    } else if (mode == Fragmenter.MODE_FILT_ACC) {
                        processFilterRejection(ccs, start, end, len, id, filterDist, false);
                    } else if (mode == Fragmenter.MODE_FILT_MH) {
                        processFilterMCMC(ccs, start, end, len, id, originalDist, filterDist);
                    } else if (mode == Fragmenter.MODE_NEBU) {
                        processNebu(ccs, true, start, end, len, id);
                    } else if (mode == Fragmenter.MODE_FRAG) {
                        processFrag(ccs, false, start, end, len, id);
                    }
                    //processFragUniform(false, start, end, len, id);
                    else if (mode == Fragmenter.MODE_FRAG_EZ) {
                        processFragEnzyme(ccs, start, end, len, id);
                    } else if (mode == Fragmenter.MODE_RT) {
                        //rw.writeLine(cs, 1);
                        processRT(ccs, start, end, len, id);
                    }


                }

                rw.close();

                //fis.close();
                //writer.flush();
                //writer.close();

                //buffy.close();
/*				if (isStop())
					qwriter.close();
				qwriter.flush();
				qwriter.close();
*/
/*				if (mode== MODE_NEBU) {					
					rndNebu.setStop();
					rndNebu.interrupt();
				}
*/
/*				while (qwriter.isAlive()) // while deadlocks
					try {
						qwriter.join();
					} catch (InterruptedException e) {
						qwriter.close(); // ?? Thread.currentThread().interrupt();
					}
*/
                //Distribution dist= new Distribution(lengthV.toIntArray());
                //medLen = lenSum / (double) lenNb;
                //lastMinLen = minLen;
                //lengthV= null;

                breakRatio = newMols / (2d * currMols);
                tgtFrac = ((double) tgtMols) / (currMols + newMols);
/*				if (mode== MODE_NEBU|| mode== MODE_FRAG)
					System.out.println("bratio "+Double.toString(breakRatio)
							+ ", all "+ (currMols+newMols)
							+ ", in "+ currMols
							+ ", new "+ newMols
							+ ", tgt "+ tgtMols
							+ ", trat "+ tgtFrac
					);
*/
                //DEBUG
/*								File save= new File("N:\\tmp\\round_"+roundCtr);
								if (save.exists())
									save.delete();									
								FileHelper.copy(tmpFile, save);
*/

                boolean b = tmpFile.delete();
                if (!b) {
                    throw new IOException("Couldn't delete source");
                }
                b = tmpWriteFile.renameTo(tmpFile);
                if (!b) {
                    throw new IOException("Couldn't move file");
                }
                ;

                Log.progressFinish(null, true);
                long total = 0;
                Iterator iter = mapFrags.keySet().iterator();
                while (iter.hasNext()) {
                    Long val = mapFrags.get(iter.next());
                    if (val != null) {
                        total += val;
                    }
                }
                if (Constants.verboseLevel >= Constants.VERBOSE_NORMAL) {
                    System.err.println("\t\t" + (currMols + newMols) + " mol, " + total + ": " + currMols + "," + tstCurrMols + " " + newMols + "," + tstNewMols);
                }
                if (Constants.verboseLevel >= Constants.VERBOSE_NORMAL) {
                    //System.err.println("\t\t" + (currMols + newMols) + " mol: in " + currMols + ", new " + newMols + ", out " + cntLinesWr);
                    // removed "out" ... number of written lines
                    System.err.println("\t\t" + (currMols + newMols) + " mol: in " + currMols + ", new " + newMols + ", out ");
                    System.err.println("\t\tavg Len " + (cumuLen / (float) currMols) + ", maxLen " + maxLen);
                    System.err.println("tgt frac " + tgtFrac);
                }
            }

            if (Constants.verboseLevel >= Constants.VERBOSE_ERRORS) {
                System.err.println(" OK");
            }
            System.gc();
            return true;
        } catch (Exception e) {
            if (Constants.verboseLevel >= Constants.VERBOSE_SHUTUP) {
                e.printStackTrace();
            }
        }
        if (Constants.verboseLevel > Constants.VERBOSE_SHUTUP) {
            System.err.print(" FAILED");
        }
        return false;

    }


    public void initRTpar(double minGC, double maxGC) {
        rtC = minGC * (1.5d + Math.log(0.5d));

    }


    private void initNebulizationParameters(double lambda, double M, double thold, int maxLen) {
        // determine C
        // C~ f(lambda,M), adjust that 1.5 lambda => pb= 0.5
        // C= lambda(1.5- (-ln(0.5))^(1/M))
        // e.g., C=486 f. M=9, lambda=900
        nebuC = lambda * (1.5d - Math.pow(-Math.log(0.5d), 1d / M));

        // expected recursion depth
        // p_t<= 0.1
        // tmax= ceil( log2((maxlen-C)/(lambda*(-ln(0.9)^(1/M))) )
        // e.g. len= 1500 -> ceil(0.53), len= 15k -> ceil(4.37)
        // maxLen= 10000; lambda= 900; nebuC= 486;
        nebuRecursionDepth =
                (int) Math.ceil(Math.log10((maxLen - nebuC) / (lambda * Math.pow(-Math.log(1d - thold), 1d / M))) / Math.log10(2));

    }

    private double getFragURdelta(double len) {

        if (Double.isNaN(settings.get(FluxSimulatorSettings.FRAG_UR_DELTA))) {
            return Math.max(Math.log10(len), 1);
        }
        return settings.get(FluxSimulatorSettings.FRAG_UR_DELTA);
    }

    /**
     * Provides eta ("intensity of fragmentation") of the uniform random fragmentation process;
     * if no eta has been provided as input parameter, eta is optimized to provide the median
     * molecule length a value that corresponds to the median of the subsequently filtered
     * values, or constant <code>DEFAULT_MED_SIZE</code>.
     *
     * @return
     */
    private double getFragUReta() {

        if (Double.isNaN(settings.get(FluxSimulatorSettings.FRAG_UR_ETA))) {
            double medLen = profiler.getMedMoleculeLength();
            double medDelta = getFragURdelta(medLen);
            double d0 = settings.get(FluxSimulatorSettings.FRAG_UR_D0);
            double medFilt = DEFAULT_MED_SIZE;
            if ((boolean) settings.get(FluxSimulatorSettings.FILTERING)) {
                medFilt = 170; //getFiltMedian();
            }

            settings.setFragUReta((medFilt - d0) / Math.exp(Gamma.logGamma(1d + 1d / medDelta)));
        }

        return settings.get(FluxSimulatorSettings.FRAG_UR_ETA);
    }


    private synchronized void addFragCount(ByteArrayCharSequence string, Long i) {
        if (mapFrags.containsKey(string)) {
            long otherVal = mapFrags.get(string) + i;    // 20101215: whyever, dont use in put clause
            // see tstCurMol and tstNewMol divergence
            mapFrags.put(string.cloneCurrentSeq(), otherVal);
        } else {
            mapFrags.put(string.cloneCurrentSeq(), i);
        }

    }


    private File writeInitialFile() {
        rndTSS = new RandomDataImpl();
        rndPA = new Random();
        rndPlusMinus = new Random();
        cntMolInit = 0;
        BufferedWriter fos = null;
        File tmpFile = null;
        try {
            Log.progressStart("Initializing Fragmentation File");
            tmpFile = File.createTempFile("Framgmenter-tmp", ".tmp");
            fos = new BufferedWriter(new FileWriter(tmpFile));
            //molInit = 0;
            maxLen = Integer.MIN_VALUE;
            int profileSize = profiler.size();
            for (int i = 0; i < profileSize; i++) {
                Log.progress(i, profileSize);
                int origLen = profiler.getLength(i);
                String compIDstring = profiler.getGlobalID(i);
                long molecules = profiler.getMolecules(i);
                for (int x = 0; x < molecules; x++) {
                    ++cntMolInit;
                    int[] ints = processInitial(origLen);
                    fos.write(ints[0] + "\t" + ints[1] + "\t" + compIDstring);
                }
            }
            Log.progressFinish("OK", true);
            Log.message("\t" + cntMolInit + " mol inited");


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

        FileInputStream fos = null;
        BufferedByteArrayReader buffy = new BufferedByteArrayReader();
        try {
            Log.progressStart("initializing library");
            fos = new FileInputStream(settings.get(FluxSimulatorSettings.LIB_FILE));
            long totBytes = file.length();
            long bytesRead = 0l;
            ByteArrayCharSequence cs = new ByteArrayCharSequence(50);
            currMols = 0;
            newMols = 0;
            for (; (buffy.readLine(fos, cs).length() > 0); ++currMols) {
                cs.resetFind();
                bytesRead += cs.length() + 1;
                Log.progress(bytesRead, totBytes);
            }
            Log.progressFinish();
            return true;
        } catch (Exception e) {
            Log.progressFailed("ERROR");
            Log.error("Error while initializing library: " + e.getMessage(), e);
            return false;
        } finally {
            if (fos != null) {
                try {
                    fos.close();
                } catch (IOException ignore) {
                }
            }
        }
    }


    void processRT(ByteArrayCharSequence cs, int start, int end, int len, ByteArrayCharSequence id) {

        if (index1 == null) {
            index1 = new int[getRTeventNr(profiler.getMaxMoleculeLength())];
        }
        double[] wSense = null, wAsense = null;
        double n1 = Double.NaN, n2 = Double.NaN;
        int txLen = -1, howmany = -1;
        if (settings.get(FluxSimulatorSettings.RT_MOTIF) == null) {
            if (settings.get(FluxSimulatorSettings.RT_PRIMER) == FluxSimulatorSettings.RtranscriptionMode.PDT) {
                txLen = profiler.getLength(id);
                if (end < txLen - 1)    // 0-based coordinates
                {
                    return;
                } else {
                    howmany = getRTeventNr(txLen - end);
                }
            } else {
                howmany = getRTeventNr(len);
            }
            howmany = Math.min(howmany, index1.length);
        } else {
            howmany = getRTeventNr(len);
            howmany = Math.min(howmany, index1.length);
            wAsense = mapWeightAsense.get(id);
            n2 = toCDF(wAsense, start, end);
            wSense = mapWeightSense.get(id);
            n1 = toCDF(wSense, start, end);
        }
        if (howmany == 0) {
            return;
        }


        //processFragNot(start, end, len, id);
        // choose new 3' end(s)
        for (int i = 0; i < howmany; i++) {
            int p;
            if (wAsense == null) {
                if (settings.get(FluxSimulatorSettings.RT_PRIMER) == FluxSimulatorSettings.RtranscriptionMode.PDT) {
                    p = end + (int) Math.floor(rtRndWhere.nextDouble() * (txLen - end));
                } else {
                    assert (settings.get(FluxSimulatorSettings.RT_PRIMER) == FluxSimulatorSettings.RtranscriptionMode.RH);
                    p = start + (int) Math.round(rtRndWhere.nextDouble() * (len - 1));
                }
            } else {
                p = Arrays.binarySearch(wAsense, start, end, rtRndWhere.nextDouble());
                p = (p >= 0 ? p : -(p + 1));
                // TODO put a minimum threshold on weight (thermodynamics affinity)?!
                ++p; // anti-sense matrix, cut 3' of 0-position
            }
            //new3Prime= Math.max(start+ p, new3Prime);
            index1[i] = p;
        }

        // extend first strand
        // resolve conflicts
        // choose new 5' (second strand synthesis) 
        // output
        Arrays.sort(index1, 0, howmany);
        for (int i = howmany - 1; i >= 0; --i) {
            if (index1[i] < 0) {
                continue; // got displaced
            }
            int ext = settings.get(FluxSimulatorSettings.RT_MIN) + (int) rnd2.nextDouble() * ((int) settings.get(FluxSimulatorSettings.RT_MAX) - settings.get(FluxSimulatorSettings.RT_MIN));
            int to = Math.max(index1[i] - ext, start);
            // check GC
            double gc = getGCcontent(id, to, index1[i] - 1);    // bp is first nt of next fragment
            double gc_lo = settings.get(FluxSimulatorSettings.RT_GC_LO);
            double pg = gc < rtC ? 0d : 1d - Math.exp(-(gc - rtC) / gc_lo);
            Math.exp((-1d) * (gc - gc_lo) / gc_lo);
            if (pg > 1 || rnd1.nextDouble() > pg) {
                index1[i] = -1;
                continue;
            }

            // resolve displacements
            boolean displaced = false;
            for (int j = i - 1; j >= 0; --j) {
                if (index1[j] < 0) {
                    continue;
                }
                if (index1[j] < to) {
                    break;
                }
                double f;
                if (wAsense == null) {    // displacement function of distance
                    int dist = Math.min(index1[i] - index1[j], index1[j] - start);
                    f = Math.exp((-1d) * dist / settings.get(FluxSimulatorSettings.RT_MIN));
                } else {    // displacement is a function of motifs
                    double mi = wAsense[index1[i] - 1] - (index1[i] == 1 ? 0 : wAsense[index1[i] - 2]);    // (-1) as cut after asense 0-pos
                    double mj = wAsense[index1[j] - 1] - (index1[j] == 1 ? 0 : wAsense[index1[j] - 2]);
                    f = (mj == 0 ? 2 : mi / mj); // smaller -> more likely displaced
                }
                if (f > 1 || rnd3.nextDouble() <= f) {
                    index1[j] = -1;     // displaced other pol
                } else {
                    displaced = true;
                    break;
                }
            }
            if (displaced) {
                continue;
            }

            // choose 5'-end for second strand synthesis
            int howmany2 = 5;    // roll 5 values
            int to2 = Math.min(to + 50, index1[i] - 1);    // within 50nt closest to 5'
            int new5Prime = index1[i];
            for (int j = 0; j < howmany2; j++) {
                int p;
                double r = rtRndWhere.nextDouble();
                if (wAsense == null) {
                    p = to + (int) Math.floor(r * (to2 - to));
                } else {
                    r = wSense[to] + (r * (wSense[to2] - wSense[to]));
                    p = Arrays.binarySearch(wSense, to, to2, r);
                    p = (p >= 0 ? p : -(p + 1));
                }
                new5Prime = Math.min(p, new5Prime);
            }

            // write it
            int nuLen = index1[i] - new5Prime + 1;
            // updateMedian(nuLen);
            cs.replace(0, new5Prime);
            cs.replace(1, index1[i]);
            cumuLen += nuLen;
            ++newMols;
            rw.writeLine(cs, fos);    // id is invalid now
            if (nuLen <= settings.get(FluxSimulatorSettings.FRAG_NB_LAMBDA)) {
                ++tgtMols;
            }
            if (len <= settings.get(FluxSimulatorSettings.FRAG_NB_LAMBDA)) {
                --tgtMols;
            }

        }

        // re-normalize
        if (settings.get(FluxSimulatorSettings.RT_MOTIF) != null) {
            toPDF(wAsense, start, end, n2);
            toPDF(wSense, start, end, n1);
        }

    }

    /**
     * @param id ID of the sequence
     * @param i  start index (included)
     * @param j  end index (included)
     * @return
     */
    private double getGCcontent(ByteArrayCharSequence id, int i, int j) {

        CharSequence seq = getMapTxSeq().get(id);
        int g = 0, n = 0;
        for (int k = i; k <= j; ++k) {
            // todo: assume gc content before transcription start
            if (k >= seq.length()) {
                n++;
            } else if (k >= 0) {
                char c = seq.charAt(k);
                if (c == 'G' || c == 'C') {
                    ++g;
                } else if (c != 'N') {
                    ++n;
                }
            }
        }
        if (n + g == 0) {
            return 0;
        }
        return (g / (double) (n + g));
    }

    private int getRTeventNr(int len) {

        int howmany = (int) Math.ceil(Math.abs(len) / (double) 100);    //  settings.getRTPrimerPerNucleotides()
//		howmany= (int) Math.round(rndHowMany.nextPoisson(((new3Prime- start)+1)/ 50d));	// TODO Poisson..;
        return howmany;
    }


    int[] tmpFragEnzyme = null;

    /**
     * Implements enzymatic cleavage of the fragment according to a provided motif (PWM).
     *
     * @param start
     * @param end
     * @param len
     * @param id
     */
    void processFragEnzyme(ByteArrayCharSequence cs, int start, int end, int len, ByteArrayCharSequence id) {

        assert (mapWeightSense != null && mapWeightAsense != null);
        double[] wsense = mapWeightSense.get(id), wasense = mapWeightAsense.get(id);

        // short-hand value to get number of enzyme attacks as function of length
        // could be a parameter
        int n = (int) Math.sqrt(len);
        int k = 0;
        int[] pos = (tmpFragEnzyme == null ? new int[n] : tmpFragEnzyme);

        // get breakpoints
        double norm = Fragmenter.toCDF(wsense, start, end);
        for (int i = 0; i < n; i++) {

            // localize potential breakpoint
            boolean sense = rnd1.nextFloat() < 0.5;
            double[] a = (sense ? wsense : wasense);
            int bp = Arrays.binarySearch(a, start, end + 1, rnd2.nextDouble());
            bp = (bp >= 0 ? bp : -(bp + 1));
            if (bp >= end) {
                continue;
            }

            // lookup in array with bp positions
            double pb = (bp == 0 ? a[bp] : a[bp] - a[bp - 1]);    // breaking probability
            // ok, the pwm recognizes the base where it is cut
            // but has no concept about 5'/3'. If cutting always
            // consistently 5' of 0-position, the bp in sense
            // marks the start of the next fragment, when motif
            // is recognized in anti-sense, we got to add 1.
            // (e.g., palindrome ACGT cut at G is found at 
            // position 1166 in sense, and position 1165 in 
            // anti-sense)
            bp = (sense ? bp : bp + 1);    // set bp consistently to 1st nt of next fragment
            int idx = Arrays.binarySearch(pos, 0, k, bp);
            if (idx >= 0) {
                continue;
            } else {
                idx = -(idx + 1);
            }

            // breaking probability
            pb /= (sense ? pwmSense.getMaximumP() : pwmSense.getMaximumP()); // should be the same if only transposed
            if (rnd3.nextDouble() <= pb) {
                System.arraycopy(pos, idx, pos, idx + 1, (k - idx) + 1);
                pos[idx] = bp;
                ++k;
            }
        }
        Fragmenter.toPDF(wsense, start, end, norm);

        // break
        int last = start;
        for (int i = 0; i < k; i++) {
            int now = start + pos[i] - 1;    // end
            int nuLen = now - last;
            if (nuLen == 0) {
                continue;
            }
            //updateMedian(nuLen);
            ++newMols;
            cumuLen += nuLen;
            cs.replace(0, last);
            cs.replace(1, now);
            rw.writeLine(cs, fos);
            last = now + 1;    // start
        }
        int now = end;
        int nuLen = now - last + 1;
        //updateMedian(nuLen);
        cumuLen += nuLen;
        cs.replace(0, last);
        cs.replace(1, now);
        rw.writeLine(cs, fos);
    }

    /**
     * @return a value representing the sum of the values that have been turned into a CDF
     */
    public static double toCDF(double[] a, int start, int end) {
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

    public static void toPDF(double[] a, int start, int end, double sum) {
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
     * Implements a uniform random fragmentation model according to
     * Tenchov and Yanev.
     *
     * @param nebu
     * @param start
     * @param end
     * @param len
     * @param id
     */
    void processFrag(ByteArrayCharSequence cs, boolean nebu, int start, int end, int len, ByteArrayCharSequence id) {

        // get parameters
        double d0 = settings.get(FluxSimulatorSettings.FRAG_UR_D0);
        assert (d0 >= 1); // problem with integer breakpoints, when fragment size << 1 !
        double delta = getFragURdelta(len);
        double eta = getFragUReta();
        //eta= 170;

        // [2a] dmax= eta((delta- 1)/delta)^(1/ delta)
        double dmax = eta * Math.pow((delta - 1) / delta, 1 / delta);
        // expectancy E of the Weibull distribution W(delta,eta)
        // [3a] E= d0+ eta*gamma(1/delta+ 1)
        double E = d0 + eta * Math.exp(Gamma.logGamma(1d + 1d / delta));
        double D2 = Math.pow(eta, 2) *
                (Math.exp(Gamma.logGamma((2d / delta) + 1d)) -
                        Math.pow(Math.exp(Gamma.logGamma((1d / delta) + 1d)), 2d));

        // DEBUG
/*		if (out1&& len< 1000) {
			System.err.println(len+"\t"+delta+"\t"+E+"\t"+D2+"\t"+dmax);
			out1= false;
			System.currentTimeMillis();
		} else if (out2&& len> 1200&& len< 1500) {
			System.err.println(len+"\t"+delta+"\t"+E+"\t"+D2+"\t"+dmax);
			out2= false;
			System.currentTimeMillis();
		} else if (out3&& len> 2000) {
			System.err.println(len+"\t"+delta+"\t"+E+"\t"+D2+"\t"+dmax);
			out3= false;
			System.currentTimeMillis();
		}
*/

        // determine n, the number of fragments (i.e. (n-1) breakpoints)
        double nn = ((double) len) / E;
        int n = (int) nn;
        // (int) Math.round(nn+ (rndBreak.nextGaussian()));	 
        // (int) rndImpl.nextPoisson(nn);		// too variable
        double nr = nn - n;
        double r = rndBreak.nextDouble();
        if (r <= nr) {
            ++n;
        }

        // molecule does not break
        if (n <= 1 || len <= 1 || (n * d0) >= len) {
            updateMedian(len);
            cumuLen += len;
            cs.replace(0, start);
            cs.replace(1, end);
            maxLen = Math.max(maxLen, len);
            rw.writeLine(cs, fos);
            //increaseFragCount(id);	// 20101215 deactivated, counts for pro-file
            return;
        }

        // uniformly cut (n-1) times unit space
        double[] x = new double[n];
        for (int i = 0; i < (n - 1); ++i) {
            x[i] = rndBreak.nextDouble();
        }
        x[x.length - 1] = 1;    // last breakpoint is end

        // get fragment lengths (in unit space)
        Arrays.sort(x);
        for (int i = (n - 1); i > 0; --i) {
            x[i] -= x[i - 1];
        }

        // compute c, transform to molecule space
        float sum = 0;
        for (int i = 0; i < x.length; i++) {
            sum += Math.pow(x[i], 1 / delta);
        }
        double c = Math.pow((len - n * d0) / sum, -delta);
        for (int i = 0; i < n; i++) {
            x[i] = d0 + Math.pow(x[i] / c, (1 / delta));
        }

        double dsum = 0;
        for (int i = 0; i < n; i++) {
            int nuStart = start + (int) Math.round(dsum);
            dsum += x[i];
            int nuEnd = start + (int) Math.round(dsum) - 1;
            //double frac= dsum/ len;
            int nuLen = (nuEnd - nuStart) + 1;
            if (nuLen < 0) {
                System.currentTimeMillis();
            }
            cs.replace(0, nuStart);
            cs.replace(1, nuEnd);
            cumuLen += nuLen;
            updateMedian(nuLen);
            rw.writeLine(cs, fos);    // id is invalid now
        }
        assert (Math.round(dsum) == len);
    }

    void processNebu(ByteArrayCharSequence cs, boolean nebu, int start, int end, int len, ByteArrayCharSequence id) {

        if (nebuRecursionDepth < 1) {
            rw.writeLine(cs, fos);    // no breakage
            return;
        }

        // parameters:
        // pb: M, lambda, len
        // bp: Sigma, length
        double lambda = this.settings.get(FluxSimulatorSettings.FRAG_NB_LAMBDA);
        double M = this.settings.get(FluxSimulatorSettings.FRAG_NB_M);
        int recDepth = Fragmenter.this.nebuRecursionDepth;

        if (index1 == null) {
            index1 = new int[(int) Math.pow(2, recDepth)];
        }
        Arrays.fill(index1, -1);
        index1[0] = len;
        int fragmentNb = 1;

        // now break it!
        for (int i = 0; i < recDepth; ++i) {
            for (int j = 0; index1[j] > 0; ++j) {

                // breakpoint location
                // N(length/2,sigma)= (N(0,1)*sigma*(length/2)+ length/2
                int L = index1[j];
                //int bp= (int) rdiNebuBP.nextGaussian(len/ 2d, len/ 4d);
                //int bp= (int) ((rndBP.nextGaussian()* sigma
                //		* ((L-1)/2d))+ (L-1)/2d);	// bp index [0;L[						
                int bp = (int) nextGaussianDouble(rndBP, 0, L - 1);

                // breaking probability (pb)
                // pb= 1- exp^(-((x-C)/lambda)^M)
                int minL = (bp < (L - bp) ? bp + 1 : L - bp - 1);
                double pb = minL < nebuC ? 0d : 1d - Math.exp(-Math.pow((minL - nebuC) / lambda, M));
                double r = rndBreak.nextDouble();
                if (r > pb) {
                    continue;
                }

                // fragment j breaks
                int rest = fragmentNb - j - 1;
                if (rest > 0) {
                    assert (j + 2 + rest <= index1.length);
                    System.arraycopy(index1, j + 1, index1, j + 2, rest);
                }
                assert (bp + 1 > 0 && L - bp - 1 > 0);
                index1[j] = bp + 1;
                index1[j + 1] = L - bp - 1;    // L- (bp+1)
                ++fragmentNb;
                ++j;
            }

        }

        // write result to disk
        for (int j = 0; index1[j] > 0; ++j) {
            cumuLen += index1[j];
            cs.replace(0, start);
            start += index1[j];
            cs.replace(1, start - 1);
            rw.writeLine(cs, fos);    // id is invalid now

            if (index1[j] <= settings.get(FluxSimulatorSettings.FRAG_NB_LAMBDA)) {
                ++tgtMols;
            }
            if (len <= settings.get(FluxSimulatorSettings.FRAG_NB_LAMBDA)) {
                --tgtMols;
            }
        }
        assert (start == len);

    }


    /**
     * see <code>http://personal.strath.ac.uk/gary.koop/extra_material_on_metropolis.pdf</code>,
     * <code>http://www.ps.uci.edu/~markm/numerical_methods/Metropolis%96Hastings%20algorithm%20-%20Wikipedia,%20the%20free%20encyclopedia.pdf</code>
     *
     * @param start
     * @param end
     * @param len
     * @param id
     */
    void processFilterMCMC(ByteArrayCharSequence cs, int start, int end, int len, ByteArrayCharSequence id,
                           AbstractDistribution dGenerate, AbstractDistribution[] dProposal) {

        // first value always accepted to init algorithm (but not output)
        double p = 0d;
        for (int i = 0; i < dProposal.length; i++) {
            p += dProposal[i].getP(len);
        }
        if (lastLen < 0) {
            lastLen = len;
            lastP = p;
            return;
        }

        // Metropolis/Hastings/Ema
        double a1 = p / lastP;
        double a2 = dGenerate.getP(lastLen, len) / dGenerate.getP(len, lastLen);
        double a = a1 * a2;

        // accept 
        if (a >= 1 || rndGel.nextDouble() <= a) {
            lastLen = len;
            lastP = p;
            rw.writeLine(cs, fos);
            addFragCount(id, 1l);
        } else {
            --newMols;
        }

    }

    void processFilterRejection(ByteArrayCharSequence cs, int start, int end, int len, ByteArrayCharSequence id,
                                AbstractDistribution[] d, boolean probDistr) {

        // get (possibly cumulative) probability for length being in result
        double plen = 0d;
        for (int i = 0; i < d.length; i++) {
            double p = (probDistr ? d[i].getP(len) : d[i].getRelFreq(len));
            plen += d[i].getWeight() * p;
        }

        // Bernoulli trial
        if (plen > 1 || rndGel.nextDouble() < plen) {
            rw.writeLine(cs, fos);
            addFragCount(id, 1l);
        } else {
            --newMols;
        }
    }

    private void updateMedian(int nuLen) {


        if (med[0] == 0) {
            med[0] = nuLen;
            for (int i = 0; i < med.length - 1; i++) {
                if (med[i] > med[i + 1]) {
                    int h = med[i];
                    med[i] = med[i + 1];
                    med[i + 1] = h;
                }
            }
            return;
        }

        int p = Arrays.binarySearch(med, nuLen);
        if (p >= 0) {
            return;
        }

        p = -(p + 1);
        if (p == 0 || p >= med.length) {
            return;
        }
        if (p + 1 < med.length) {
            System.arraycopy(med, p, med, p + 1, med.length - p - 1);
        }
        med[p] = nuLen;
    }

    int[] processInitial(int origLen) {

        // transcript variation
        int start = 0;
        int end = origLen - 1;
        Double tssMean = settings.get(FluxSimulatorSettings.TSS_MEAN); // todo : move this up so we call it just once
        if (!Double.isNaN(tssMean)) {
            start = origLen;
            while (start >= Math.min(100, origLen)) {
                start = (int) Math.round(rndTSS.nextExponential(Math.min(tssMean, origLen / 4)));    // exp mean= 25: exceeds bounds, nextGaussian(1,100))-100;
            }
            double r = rndPlusMinus.nextDouble();
            if (r < 0.5) {
                start = -start;
            }
        }
        if (!(Double.isNaN(settings.get(FluxSimulatorSettings.POLYA_SHAPE)) || Double.isNaN(settings.get(FluxSimulatorSettings.POLYA_SCALE)))) {
            int pAtail = 301;
            while (pAtail > 300) {
                pAtail = (int) Math.round(sampleWeibull(rndPA, settings.get(FluxSimulatorSettings.POLYA_SCALE), settings.get(FluxSimulatorSettings.POLYA_SHAPE)));    // 300d, 2d
            }
            end = origLen + pAtail;
        }
        if (end < origLen) {
            Log.error("end < length in Fragmenter!");
        }

        int newLen = end - start + 1;
        assert (start < end);
        cumuLen += newLen;
        maxLen = Math.max(maxLen, newLen);
        return new int[]{start, end};
    }


}
