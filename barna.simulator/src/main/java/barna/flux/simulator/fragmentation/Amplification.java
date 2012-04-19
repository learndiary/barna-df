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
import barna.commons.log.Log;
import barna.flux.simulator.PWM;
import barna.flux.simulator.distributions.AbstractDistribution;
import barna.flux.simulator.distributions.EmpiricalDistribution;
import barna.flux.simulator.distributions.GCPCRDistribution;
import barna.flux.simulator.distributions.NormalDistribution;

import java.io.File;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Random;

/**
 *
 * GC filtering and amplification
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class Amplification implements FragmentProcessor{
    private GCPCRDistribution pcrDistribution;
    /**
     * PCR duplication prob in case GC is disabled
     */
    private double pcrProbability = 0.7;
    /**
     * The mean
     */
    private double mean = 0.5;
    /**
     * The standard deviation
     */
    private double sigma = 0.1;

    /**
     * Number of processed fragments
     */
    private long in;
    /**
     * Number of returned fragments
     */
    private long out;
    /**
     * Map IDs to sequences
     */
    private Map<CharSequence, CharSequence> mapTxSeq;

    /**
     * The GC distribution
     */
    private AbstractDistribution distribution;
    /**
     * Maximal number of duplicated fragments
     */
    private long maxFragments = 0;

    /**
     * Radom sampler
     */
    private Random random;
    
    Random randomDELME= new Random();
	Map<CharSequence, CharSequence> mapTx;
	private Map<CharSequence, double[]> mapWeightAsense;
	private Map<CharSequence, double[]> mapWeightSense;
	private PWM pwmASense;
	private PWM pwmSense;



    /**
     * Create a new instance
     *
     * @param pcrDistribution the pcr distribution
     * @param mean the mean
     * @param sigma the standard deviation
     */
    public Amplification(GCPCRDistribution pcrDistribution, final double pcrProbability, final double mean, final double sigma, Map<CharSequence, CharSequence> mapTxSeq) {
        this.pcrDistribution = pcrDistribution;
        this.pcrProbability = pcrProbability;
        this.mean = mean;
        this.sigma = sigma;
        this.mapTxSeq = mapTxSeq;
        if(!Double.isNaN(mean)){
            this.distribution = new EmpiricalDistribution(new NormalDistribution(mean, sigma), 1000, 4);
        }
        this.maxFragments = Math.max(1,(long) Math.pow(2, pcrDistribution.getGenerations())-1);
        this.random = new Random();
    }

    /**
	 * Create a new instance
	 *
	 * @param pcrDistribution the pcr distribution
	 * @param mean the mean
	 * @param sigma the standard deviation
	 */
	public Amplification(GCPCRDistribution pcrDistribution, final double pcrProbability, final double mean, final double sigma, Map<CharSequence, CharSequence> mapTxSeq,
			Map<CharSequence, CharSequence> mapTx, File pwmFile) {
	    this(pcrDistribution, pcrProbability, mean, sigma, mapTxSeq);
    	this.mapTx= mapTx;
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
/*                for (int i = 0; i < 100; i++) {
                    pwmSense.multiply();
                    pwmASense.multiply();
                }
                pwmSense.makePDF();
                pwmASense.makePDF();
*/             
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

	@Override
    public List<Fragment> process(final ByteArrayCharSequence id, final ByteArrayCharSequence cs, final int start, final int end, final int len) {

        double[] wSense = null;
        double[] wAsense = null;
        in++;
    	if (mapWeightAsense!= null) {	// customMotif
            wAsense = mapWeightAsense.get(id);
            wSense = mapWeightSense.get(id);
        }

		// get the gc content
        double gcp = -1, gc= 0;
        if(distribution != null){
            gc = getGCcontent(id, start, end);
            gcp = distribution.getP(gc);
        }

        List<Fragment> fragments = new ArrayList<Fragment>();
        Fragment fragment = new Fragment(id, start, end);

        int nfragments = 1;

        // test DEBUG
//        double p1= wSense[start];
//        double p2= wAsense[end];
//        double pp= p1* p2;
        
//        double q= randomDELME.nextDouble();
//        if (q> p1)
//        	return fragments;
//        q= randomDELME.nextDouble();
//        if (q> p2)
//        	return fragments;
        
        
//      double pgc= Math.abs(mean- gc);
//      double ppp= 10* pp* ((!Double.isNaN(mean))? pgc: pcrProbability);
        
        if(!Double.isNaN(mean)){
            nfragments = (int) (Math.max(0, 25* pcrDistribution.getNext(random.nextDouble(), gcp)));
        }else{
            // use pcr prob
            nfragments = (int) (Math.max(0, 25* pcrDistribution.getNext(random.nextDouble(), pcrProbability)));
        }
        
        // test DEBUG
//        nfragments= (int) Math.pow(nfragments, pp);

        if (nfragments> 0)
            fragments.add(fragment);
        else
        	System.currentTimeMillis();
        
        out+=nfragments;
        fragment.setDuplicates(nfragments);
        return fragments;
    }

    /**
     * Compute the relative GC content
     *
     * @param id ID of the sequence
     * @param start  start index (included)
     * @param stop  end index (included)
     * @return gc gc content
     */
    private double getGCcontent(ByteArrayCharSequence id, int start, int stop) {
        CharSequence seq = mapTxSeq.get(id);
        int g = 0, n = 0;
        for (int k =  start; k <= stop; ++k) {
            if(k >=0 && k< seq.length()){
                char c = seq.charAt(k);
                if (c == 'G' || c == 'C' || c == 'g' || c == 'c' || c == 's' || c == 'S') {
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
        return "Amplification";
    }

    @Override
    public String getConfiguration() {
        StringBuffer b = new StringBuffer();
        b.append("\t\tRounds: " + pcrDistribution.getGenerations()+" \n");
        if(Double.isNaN(mean)){
            b.append("\t\tPCR Probability: " + pcrProbability+" \n");
        }else{
            b.append("\t\tMean: " + mean+" \n");
            b.append("\t\tStandard Deviation: " + sigma+" \n");
        }
        b.append("\n");
        return b.toString();
    }

    @Override
    public String done() {
        return "\tAmplification done.\n\tIn: " + in + " Out: " + out+"\n\n";
    }

	public void initPWMMap(){
	    if(pwmASense != null){
	        Log.info("Initializing PWM cache");
	        mapWeightSense = Fragmenter.getMapWeight(mapTx,null, pwmSense);
	        mapWeightAsense = Fragmenter.getMapWeight(mapTx,null, pwmASense);
	        Log.info("Done");
	    }
	}

	//@Override
	public List<Fragment> process_original(final ByteArrayCharSequence id, final ByteArrayCharSequence cs, final int start, final int end, final int len) {
	    // get the gc content
	    double gcp = -1;
	    if(distribution != null){
	        double gc = getGCcontent(id, start, end);
	        gcp = distribution.getP(gc);
	    }
	
	    List<Fragment> fragments = new ArrayList<Fragment>();
	    Fragment fragment = new Fragment(id, start, end);
	    fragments.add(fragment);
	
	    int nfragments = 1;
	
	    if(!Double.isNaN(mean)){
	        nfragments = (int) Math.max(1, pcrDistribution.getNext(random.nextDouble(), gcp));
	    }else{
	        // use pcr prob
	        nfragments = (int) Math.max(1, pcrDistribution.getNext(random.nextDouble(), pcrProbability));
	    }
	    
	    in++;
	    out+=nfragments;
	    fragment.setDuplicates(nfragments);
	    return fragments;
	}
}
