/*
 * Copyright (c) 2010, Micha Sammeth
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * The names of its contributors may be not used to endorse or promote
 *       products derived from this software without specific prior written
 *       permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL MICHA SAMMETH BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

package barna.flux.simulator.fragmentation;

import barna.commons.ByteArrayCharSequence;
import barna.commons.RandomFactory;
import barna.flux.simulator.PWM;

import java.io.File;
import java.io.InputStreamReader;
import java.util.*;

/**
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class FragmentEnzymatic implements FragmentProcessor {
    private Map<CharSequence, double[]> mapWeightSense = null;
    private Map<CharSequence, double[]> mapWeightAsense = null;
    private PWM pwmSense;
    private PWM pwmAsense;
    Random rnd1 = RandomFactory.get();
    Random rnd2 = RandomFactory.get();
    Random rnd3 = RandomFactory.get();
    private File ezMotif;
    private int leftFlank;
    private int rightFlank;


    public FragmentEnzymatic(File ezMotif, Map<CharSequence, CharSequence> mapTxSeq, int leftFlank, final int rightFlank) throws Exception {
        this.ezMotif = ezMotif;
        this.leftFlank = leftFlank;
        this.rightFlank = rightFlank;
        if (ezMotif == null) {
            throw new NullPointerException("Null Motif file for enzymatic fragmentation not permitted !");
        }

        if(!ezMotif.exists() && (ezMotif.getName().equalsIgnoreCase("NlaIII") || ezMotif.getName().equalsIgnoreCase("DpnII"))){
            if(ezMotif.getName().equalsIgnoreCase("NlaIII")){
                pwmSense = PWM.create(new InputStreamReader(getClass().getResource("/NlaIII.pwm").openStream()));
            }else{
                pwmSense = PWM.create(new InputStreamReader(getClass().getResource("/DpnII.pwm").openStream()));
            }

        }else{
            pwmSense = PWM.create(ezMotif);
        }
        mapWeightSense = Fragmenter.getMapWeight(mapTxSeq, mapTxSeq, pwmSense);
        pwmSense.invert();
        pwmAsense = pwmSense;
        mapWeightAsense = Fragmenter.getMapWeight(mapTxSeq, mapTxSeq, pwmAsense);
    }

    @Override
    public List<Fragment> process(final ByteArrayCharSequence id, final ByteArrayCharSequence cs, final int start, final int end, final int len) {
        assert (mapWeightSense != null && mapWeightAsense != null);
        double[] wsense = mapWeightSense.get(id), wasense = mapWeightAsense.get(id);

        // short-hand value to get number of enzyme attacks as function of length
        // could be a parameter
        int n = (int) Math.sqrt(len);
        int k = 0;
        int[] pos = new int[n];

        // get breakpoints
        double norm = Fragmenter.toCDF(wsense, start, end, leftFlank, rightFlank);
        for (int i = 0; i < n; i++) {

            // localize potential breakpoint
            boolean sense = rnd1.nextFloat() < 0.5;
            double[] a = (sense ? wsense : wasense);
            int bp = Arrays.binarySearch(a, Math.max(start,0), Math.min(end + 1, a.length- 1), rnd2.nextDouble());
            bp = (bp >= 0 ? bp : -(bp + 1));
            if (bp >= len- 1) {
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
                System.arraycopy(pos, idx, pos, idx + 1, (k - idx));
                pos[idx] = bp;
                ++k;
            }
        }
        Fragmenter.toPDF(wsense, start, end, norm, leftFlank, rightFlank);

        // break
        int last = start;
        List<Fragment> fragments = new ArrayList<Fragment>();
        for (int i = 0; i < k; i++) {
            int now = start + pos[i] - 1;    // end
            int nuLen = now - last;
            if (nuLen == 0) {
                continue;
            }
            cs.replace(0, last);
            cs.replace(1, now);
            fragments.add(new Fragment(id, last, now));
            //fragments.add(cs.toString());
            //rw.writeLine(cs, fos);
            last = now + 1;    // start
        }
        int now = end;
        int nuLen = now - last + 1;
        //updateMedian(nuLen);

        cs.replace(0, last);
        cs.replace(1, now);
        //rw.writeLine(cs, fos);
        //fragments.add(cs.toString());
        fragments.add(new Fragment(id, last, now));
        return fragments;
    }

    @Override
    public String getName() {
        return "Enzymatic Digestion";
    }

    @Override
    public String getConfiguration() {
        return "Left Flank : " + leftFlank + barna.commons.system.OSChecker.NEW_LINE +
               "Right Flank : " + rightFlank+ barna.commons.system.OSChecker.NEW_LINE+
               "Motif: " + ezMotif.getName();
    }

    @Override
    public String done() {
        return null;
    }
}
