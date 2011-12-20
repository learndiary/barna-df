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

package barna.genome.sequencing.rnaseq.simulation.fragmentation;

import barna.commons.ByteArrayCharSequence;
import barna.commons.log.Log;
import barna.genome.sequencing.rnaseq.simulation.FluxSimulatorSettings;
import barna.genome.sequencing.rnaseq.simulation.PWM;
import barna.genome.sequencing.rnaseq.simulation.Profiler;

import java.io.File;
import java.io.InputStreamReader;
import java.util.*;

/**
 * Apply reverse transcription
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class FragmentReverseTranscription implements FragmentProcessor {
    private FluxSimulatorSettings.RtranscriptionMode mode;
    private File pwmFile;
    private int leftFlank;
    private int rightFlank;
    private boolean lossless;
    private boolean customMotif;
    //private Map<CharSequence, double[]> mapWeightSense = null;
    //private Map<CharSequence, double[]> mapWeightAsense = null;
    private int[] index1;
    private Random rnd1 = new Random();
    private Random rnd2 = new Random();
    private Random rnd3 = new Random();
    private Random rtRndWhere = new Random();
    private Profiler profiler;
    private int rtMin;
    private int rtMax;
    private Map<CharSequence, CharSequence> mapTxSeq;

    private PWM pwmSense;
    private PWM pwmASense;
    private Map<CharSequence, double[]> mapWeightAsense;
    private Map<CharSequence, double[]> mapWeightSense;

    public FragmentReverseTranscription(FluxSimulatorSettings.RtranscriptionMode mode,
                                        File pwmFile,
                                        final int rtMin,
                                        final int rtMax,
                                        final Map<CharSequence, CharSequence> mapTxSeq,
                                        final Profiler profiler,
                                        final int leftFlank,
                                        final int rightFlank,
                                        final boolean lossless) {
        this.mode = mode;
        this.pwmFile = pwmFile;
        this.leftFlank = leftFlank;
        this.rightFlank = rightFlank;
        this.lossless = lossless;
        this.customMotif = pwmFile != null;

        this.rtMin = rtMin;
        this.rtMax = rtMax;
        this.mapTxSeq = mapTxSeq;
        this.profiler = profiler;


        if (pwmFile != null) {
            // read the sequence annotations
            try {
                pwmSense = null;

                if(!pwmFile.exists() && pwmFile.getName().equalsIgnoreCase("default")){
                    // load default motif_1mer_0-5
                    pwmSense = PWM.create(new InputStreamReader(getClass().getResource("/motif_1mer_0-5.pwm").openStream()));
                    pwmASense = PWM.create(new InputStreamReader(getClass().getResource("/motif_1mer_0-5.pwm").openStream()));
                }else{
                    pwmSense = PWM.create(pwmFile);
                    pwmASense = PWM.create(pwmFile);
                }
                for (int i = 0; i < 100; i++) {
                    pwmSense.multiply();
                    pwmASense.multiply();
                }
                pwmSense.makePDF();
                pwmASense.makePDF();
                pwmASense.invert();
                //pwmSense.makePDF();
                //mapWeightSense = Fragmenter.getMapWeight(mapTxSeq, mapTxSeq, pwmSense);
                //pwmSense.invert();
                //PWM pwmAsense = pwmSense;
                //mapWeightAsense = Fragmenter.getMapWeight(mapTxSeq, mapTxSeq, pwmAsense);
            } catch (Exception e) {
                Log.error("Error while initializing PWM : " + e.getMessage(), e);
            }
        }
    }


    public void initPWMMap(){
        if(pwmASense != null){
            Log.info("Initializing PWM cache");
            mapWeightAsense = Fragmenter.getMapWeight(mapTxSeq,null, pwmASense);
            mapWeightSense = Fragmenter.getMapWeight(mapTxSeq,null, pwmSense);
            Log.info("Done");
        }
    }

    @Override
    public List<Fragment> process(final ByteArrayCharSequence id, final ByteArrayCharSequence cs, final int start, final int end, final int len) {
        if (index1 == null) {
            index1 = new int[getRTeventNr(profiler.getMaxMoleculeLength())];
        }
        double[] wSense = null;
        double[] wAsense = null;
        double n1 = Double.NaN;
        double n2 = Double.NaN;
        int txLen = -1;
        int howmany = -1;

        if (!customMotif) {
            if (mode == FluxSimulatorSettings.RtranscriptionMode.PDT) {
                txLen = profiler.getLength(id);
                howmany = getRTeventNr(end - txLen);
            } else {
                howmany = getRTeventNr(len);
            }
            howmany = Math.min(howmany, index1.length);
        } else {
            howmany = getRTeventNr(len);
            howmany = Math.min(howmany, index1.length);

            wAsense = mapWeightAsense.get(id);
            //wAsense = Fragmenter.applyPWM(mapTxSeq.get(id), pwmASense);
            n2 = Fragmenter.toCDF(wAsense, start, end, leftFlank, rightFlank);

            wSense = mapWeightSense.get(id);
            //wSense = Fragmenter.applyPWM(mapTxSeq.get(id), pwmSense);
            n1 = Fragmenter.toCDF(wSense, start, end, leftFlank, rightFlank);
        }

        if (lossless) {
            howmany = Math.max(howmany, 1);
        }
        if (howmany <= 0) {
            return null;
        }

        for (int i = 0; i < howmany; i++) {
            int p = end;
            if (wAsense == null) {
                if (mode == FluxSimulatorSettings.RtranscriptionMode.PDT) {
                    int extension = (int) Math.floor(rtRndWhere.nextDouble() * (end - txLen));
                    p = end + extension;
                } else {
                    int ext = rtMin + (int) (rnd2.nextDouble() * (rtMax - rtMin));
                    int maxEnd = start + (Math.min(end - start, ext) - 1);
                    for (int j = 0; j < maxEnd; j++) {
                        if (rtRndWhere.nextBoolean()) {
                            // found end
                            p = maxEnd - i;
                            break;
                        }
                    }
                }
            } else {
                int ext = rtMin + (int) (rnd2.nextDouble() * (rtMax - rtMin));
                int maxEnd = start + (Math.min(end - start, ext) - 1);

                int s = Math.min(wAsense.length-1, leftFlank + start);
                int e = Math.min(wAsense.length, leftFlank + rightFlank + maxEnd);
                p = Arrays.binarySearch(wAsense, s, e, rtRndWhere.nextDouble());
                p = (p >= 0 ? p : -(p + 1));
                ++p; // anti-sense matrix, cut 3' of 0-position
            }
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
//            if (gc_lo > 0) {
//                double gc = getGCcontent(id, from, index1[i] - 1);
//                double pg = gcScore(gc, 0.04);
//                if (rnd1.nextDouble() > pg) {
//                    index1[i] = -1;
//                    gcOut[((int) (100 * gc))]++;
//                    continue;
//                } else {
//                    gcIn[((int) (100 * gc))]++;
//                }
//            }

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
                    int ri1 = Math.min(wAsense.length-1, rightFlank + index1[i] - 1);
                    int ri2 = Math.min(wAsense.length-1, rightFlank + index1[i] - 2);
                    int rj1 = Math.min(wAsense.length-1, rightFlank + index1[j] - 1);
                    int rj2 = Math.min(wAsense.length-1, rightFlank + index1[j] - 2);

                    double mi = wAsense[ri1] - (index1[i] == 1 ? 0 : wAsense[ri2]);    // (-1) as cut after asense 0-pos
                    double mj = wAsense[rj1] - (index1[j] == 1 ? 0 : wAsense[rj2]);
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
            int to2 = Math.min(from + 50, index1[i] - 1);    // within 50nt closest to 5'
            for (int j = 0; j < howmany2; j++) {
                int p = start;
                double r = rtRndWhere.nextDouble();
                if (wAsense == null) {
                    p = from + (int) Math.floor(r * (to2 - from));
                } else {
                    to2 = Math.min(wSense.length-1, Math.min(leftFlank + from + 50, index1[i] - 1));    // within 50nt closest to 5'
                    int from2 = Math.min(leftFlank + from, wSense.length-1);

                    r = wSense[from2] + (r * (wSense[to2] - wSense[from2]));
                    p = Arrays.binarySearch(wSense, from2, to2, r);
                    p = (p >= 0 ? p : -(p + 1));
                }
                new5Prime = Math.min(p, new5Prime);
            }

            // write it
            cs.replace(0, new5Prime);
            cs.replace(1, index1[i]);

            Fragment fragment = new Fragment(id, new5Prime, index1[i]);
            fragments.add(fragment);
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


    @Override
    public String getName() {
        return "Reverse Transcription";
    }

    @Override
    public String getConfiguration() {
        StringBuffer b = new StringBuffer();
        b.append("\t\t").append("Mode: ").append(mode).append("\n");
        b.append("\t\t").append("PWM: ").append(pwmFile == null ? "No" : pwmFile.getName()).append("\n");
        b.append("\t\t").append("RT MIN: ").append(rtMin).append("\n");
        b.append("\t\t").append("RT MAX: ").append(rtMax).append("\n");
        return b.toString();
    }

    @Override
    public String done() {
        return null;
    }

}

