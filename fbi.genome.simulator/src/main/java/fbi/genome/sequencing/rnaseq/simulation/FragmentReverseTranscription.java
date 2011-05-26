package fbi.genome.sequencing.rnaseq.simulation;

import fbi.commons.ByteArrayCharSequence;
import fbi.commons.Log;

import java.io.File;
import java.util.*;

/**
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class FragmentReverseTranscription implements FragmentProcessor{
    private FluxSimulatorSettings.RtranscriptionMode mode;
    private File pwmFile;
    private int leftFlank;
    private int rightFlank;
    private boolean customMotif;
    private Map<CharSequence, double[]> mapWeightSense = null;
    private Map<CharSequence, double[]> mapWeightAsense = null;
    private int[] index1;
    private Random rnd1 = new Random();
    private Random rnd2 = new Random();
    private Random rnd3 = new Random();
    private Random rtRndWhere = new Random();
    private Profiler profiler;
    private int rtMin;
    private int rtMax;
    private double gc_lo;// = settings.get(FluxSimulatorSettings.RT_GC_LO);
    private double rtC;// = gc_lo * (1.5d + Math.log(0.5d));
    private Map<CharSequence, CharSequence> mapTxSeq;

    public FragmentReverseTranscription(FluxSimulatorSettings.RtranscriptionMode mode,
                                        File pwmFile,
                                        final int rtMin,
                                        final int rtMax,
                                        final double gc_lo,
                                        final Map<CharSequence, CharSequence> mapTxSeq,
                                        final Profiler profiler,
                                        final int leftFlank, final int rightFlank) {
        this.mode = mode;
        this.pwmFile = pwmFile;
        this.leftFlank = leftFlank;
        this.rightFlank = rightFlank;
        this.customMotif = pwmFile != null;

        this.rtMin = rtMin;
        this.rtMax = rtMax;
        this.gc_lo = gc_lo;
        this.mapTxSeq = mapTxSeq;
        this.profiler = profiler;
        this.rtC = gc_lo * (1.5d + Math.log(0.5d));

        if (pwmFile != null) {
            // read the sequence annotations
            try {
                PWM pwmSense = PWM.create(pwmFile);
                for (int i = 0; i < 100; i++) {
                    pwmSense.multiply();
                }
                pwmSense.makePDF();
                mapWeightSense = Fragmenter.getMapWeight(mapTxSeq, mapTxSeq, pwmSense);
                pwmSense.invert();
                PWM pwmAsense = pwmSense;
                mapWeightAsense = Fragmenter.getMapWeight(mapTxSeq,mapTxSeq, pwmAsense);
            } catch (Exception e) {
                Log.error("Error while initializing PWM : " + e.getMessage(), e);
            }
        }



    }

    @Override
    public List<Fragment> process(final ByteArrayCharSequence id, final ByteArrayCharSequence cs, final int start, final int end, final int len) {
        if (index1 == null) {
            index1 = new int[getRTeventNr(profiler.getMaxMoleculeLength())];
        }
        double[] wSense = null, wAsense = null;
        double n1 = Double.NaN, n2 = Double.NaN;
        int txLen = -1, howmany = -1;

        //File motifFile = null; // FluxSimulatorSettings.RT_MOTIF
        if (!customMotif) {
            if (mode == FluxSimulatorSettings.RtranscriptionMode.PDT) {
                txLen = profiler.getLength(id);
                if (end < txLen - 1){    // 0-based coordinates
                    return null;
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
            n2 = Fragmenter.toCDF(wAsense, start, end, leftFlank, rightFlank);

            wSense = mapWeightSense.get(id);
            n1 = Fragmenter.toCDF(wSense, start, end, leftFlank, rightFlank);
        }

        if (howmany <= 0) {
            throw new RuntimeException("No RT events !");
        }

        //processFragNot(start, end, len, id);
        // choose new 3' end(s)
//        for (int i = 0; i < howmany; i++) {
//            int p;
//            if (wAsense == null) {
//                if (mode == FluxSimulatorSettings.RtranscriptionMode.PDT) {
//                    p = end + (int) Math.floor(rtRndWhere.nextDouble() * (txLen - end));
//                } else {
//                    assert (mode == FluxSimulatorSettings.RtranscriptionMode.RH);
//                    p = start + (int) Math.round(rtRndWhere.nextDouble() * (len - 1));
//                }
//            } else {
//                p = Arrays.binarySearch(wAsense, start, end, rtRndWhere.nextDouble());
//                p = (p >= 0 ? p : -(p + 1));
//                // TODO put a minimum threshold on weight (thermodynamics affinity)?!
//                ++p; // anti-sense matrix, cut 3' of 0-position
//            }
//            //new3Prime= Math.max(start+ p, new3Prime);
//            index1[i] = p;
//        }

        for (int i = 0; i < howmany; i++) {
            int p = end;
            if (wAsense == null) {
                if (mode == FluxSimulatorSettings.RtranscriptionMode.PDT) {
                    p = end + (int) Math.floor(rtRndWhere.nextDouble() * (txLen - end));
                } else {

                    int ext = rtMin + (int) (rnd2.nextDouble() * (rtMax - rtMin));
                    int maxEnd = start + (Math.min(end-start, ext)-1);

                    for (int j = 0; j < maxEnd; j++) {
                        if(rtRndWhere.nextBoolean()){
                            // found end
                            p = maxEnd-i;
                            break;
                        }
                    }

                    //p = start + (int) Math.round(rtRndWhere.nextDouble() * (len - 1));
                }
            } else {

                int ext = rtMin + (int) (rnd2.nextDouble() * (rtMax - rtMin));
                int maxEnd = start + (Math.min(end-start, ext)-1);

                //p = Arrays.binarySearch(wAsense, start, end, rtRndWhere.nextDouble());
                p = Arrays.binarySearch(wAsense, leftFlank+start, leftFlank+rightFlank+maxEnd, rtRndWhere.nextDouble());
                p = (p >= 0 ? p : -(p + 1));
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

        List<Fragment> fragments = new ArrayList<Fragment>();
        for (int i = howmany - 1; i >= 0; --i) {
            if (index1[i] < 0) {
                continue; // got displaced
            }
            int ext = rtMin + (int) (rnd2.nextDouble() * (rtMax - rtMin));
            int from = Math.max(index1[i] - ext, start);
            // check GC
            double gc = getGCcontent(id, from, index1[i] - 1);    // bp is first nt of next fragment

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
                if (index1[j] < from) {
                    break;
                }
                double f;
                if (wAsense == null) {    // displacement function of distance
                    int dist = Math.min(index1[i] - index1[j], index1[j] - start);
                    f = Math.exp((-1d) * dist / rtMin);
                } else {    // displacement is a function of motifs
                    double mi = wAsense[rightFlank+index1[i] - 1] - (index1[i] == 1 ? 0 : wAsense[rightFlank+index1[i] - 2]);    // (-1) as cut after asense 0-pos
                    double mj = wAsense[rightFlank+index1[j] - 1] - (index1[j] == 1 ? 0 : wAsense[rightFlank+index1[j] - 2]);
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


            int new5Prime = index1[i];
            for (int j = 0; j < howmany2; j++) {
                int p = start;
                double r = rtRndWhere.nextDouble();
                if (wAsense == null) {
                    //p = from + (int) Math.floor(r * (to2 - from));
                    for (int k = start; k < index1[i]; k++) {
                        if(rtRndWhere.nextBoolean()){
                            p = k;
                            break;
                        }
                    }
                } else {
                    int to2 = Math.min(leftFlank+from + 50, index1[i] - 1);    // within 50nt closest to 5'
                    r = wSense[leftFlank+from] + (r * (wSense[to2] - wSense[leftFlank+from]));
                    p = Arrays.binarySearch(wSense, leftFlank + from, to2, r);
                    p = (p >= 0 ? p : -(p + 1));
                }
                new5Prime = Math.min(p, new5Prime);
            }

            // write it
            cs.replace(0, new5Prime);
            cs.replace(1, index1[i]);
            fragments.add(new Fragment(id, new5Prime, index1[i]));
        }

        // re-normalize
        if (customMotif) {
            Fragmenter.toPDF(wAsense, start, end, n2, leftFlank, rightFlank);
            Fragmenter.toPDF(wSense, start, end, n1, leftFlank, rightFlank);
        }
        return fragments;
    }

    private int getRTeventNr(int len) {
        return (int) Math.ceil(Math.abs(len) / (double) 100);
    }


    /**
     * Compute the relative GC content
     *
     * @param id ID of the sequence
     * @param i  start index (included)
     * @param j  end index (included)
     * @return gc gc content
     */
    private double getGCcontent(ByteArrayCharSequence id, int i, int j) {

        CharSequence seq = mapTxSeq.get(id);
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

    @Override
    public String getName() {
        return "Reverse Transcription";
    }

    @Override
    public String getConfiguration() {
        StringBuffer b = new StringBuffer();
        b.append("\t\t").append("Mode: ").append(mode).append("\n");
        b.append("\t\t").append("Custom PWM: ").append(pwmFile == null ? "No" : pwmFile.getName()).append("\n");
        b.append("\t\t").append("RT MIN: ").append(rtMin).append("\n");
        b.append("\t\t").append("RT MAX: ").append(rtMax).append("\n");
        b.append("\t\t").append("GC LO: ").append(gc_lo).append("\n");
        b.append("\t\t").append("RT C: ").append(rtC).append("\n");
        return b.toString();
    }
}


