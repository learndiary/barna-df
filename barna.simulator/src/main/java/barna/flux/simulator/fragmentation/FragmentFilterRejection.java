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

package barna.flux.simulator.fragmentation;

import barna.commons.ByteArrayCharSequence;
import barna.commons.log.Log;
import barna.flux.simulator.PWM;
import barna.flux.simulator.distributions.AbstractDistribution;

import java.io.File;
import java.io.InputStreamReader;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.Random;

/**
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class FragmentFilterRejection implements FragmentProcessor {
    private AbstractDistribution[] d;
    private boolean probDistr;
    private Random rndGel = new Random();
    private Random rndDELME= new Random();
	private PWM pwmASense;
	private PWM pwmSense;
	private Map<CharSequence, double[]> mapWeightAsense;
	private Map<CharSequence, double[]> mapWeightSense;
	Map<CharSequence, CharSequence> mapTx;

    public FragmentFilterRejection(AbstractDistribution[] d, boolean probDistr) {
        this.d = d;
        this.probDistr = probDistr;
    }
    
    public FragmentFilterRejection(AbstractDistribution[] d, boolean probDistr, Map<CharSequence, CharSequence> mapTx, File pwmFile) {
    	this(d, probDistr);
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
    	if (false) {	// customMotif
            wAsense = mapWeightAsense.get(id);
            wSense = mapWeightSense.get(id);
        }
    	// get (possibly cumulative) probability for length being in result
        double plen = 0d;
        for (int i = 0; i < d.length; i++) {
            double p = (probDistr ? d[i].getP(len) : d[i].getRelFreq(len));
            plen += d[i].getWeight() * p;
        }

        // Bernoulli trial
        double b= rndGel.nextDouble();
        if (plen > 1 || b < plen) {
//        	double rnd= rndDELME.nextDouble();
//            double rnd2= rndDELME.nextDouble();
//            
//            double p= start>= 0&& start< wSense.length? wSense[start]: 0;
//            double q= end>= 0&& end< wAsense.length? wAsense[end]: 0;
//            if (rnd<= p&& rnd2<= q)	
            	return Arrays.asList(new Fragment(id, start, end));
//            else
//            	return null;
        } else {
            return null;
        }
    }


    @Override
    public String getName() {
        return "Segregating cDNA ("+(!probDistr ? "Acceptance" : "Rejection")+")";
    }

    @Override
    public String getConfiguration() {
        return null;
    }

    @Override
    public String done() {
        return null;
    }

	public void initPWMMap(){
	    if(pwmASense != null){
	        Log.info("Initializing PWM cache");
	        mapWeightSense = Fragmenter.getMapWeight(mapTx,null, pwmSense);
	        mapWeightAsense = Fragmenter.getMapWeight(mapTx,null, pwmASense);
	        Log.info("Done");
	    }
	}

}
