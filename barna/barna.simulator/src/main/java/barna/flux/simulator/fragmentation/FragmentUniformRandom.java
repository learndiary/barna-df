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

import org.apache.commons.math.special.Gamma;

import java.io.File;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.Random;

/**
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class FragmentUniformRandom implements FragmentProcessor {
    /**
     * Default median size after fragmentation, according to 2010 Illumina protocol.
     */
    private static final double DEFAULT_MED_SIZE = 200;
    private double d0;
    private double urDelta;
    private double urEta;
    private Random rndBreak = new Random();
    private double medMoleculeLength;
    private boolean filtering;

    Map<CharSequence, CharSequence> mapTx;
	private PWM pwmASense;
	private PWM pwmSense;
	private Map<CharSequence, double[]> mapWeightAsense;
	private Map<CharSequence, double[]> mapWeightSense;
    
    public FragmentUniformRandom(double d0, double delta, double eta, double medMoleculeLength, boolean filtering,  
    		Map<CharSequence, CharSequence> mapTx, File pwmFile) {
    	this(d0, delta, eta, medMoleculeLength, filtering);
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
    public FragmentUniformRandom(double d0, double delta, double eta, double medMoleculeLength, boolean filtering) {
        this.filtering = filtering;
        if (d0 < 1.0) {
            throw new IllegalArgumentException("D0 must be >= 1.0!");
        }
        this.d0 = d0;
        this.urDelta = delta;
        this.urEta = eta;
        this.medMoleculeLength = medMoleculeLength;
    }

    Random rndDELME= new Random();
    @Override
    public List<Fragment> process(final ByteArrayCharSequence id, final ByteArrayCharSequence cs, final int start, final int end, final int len) {
        assert (d0 >= 1); // problem with integer breakpoints, when fragment size << 1 !
        double delta = getFragURdelta(len);
        double eta = getFragUReta();
        double[] wSense = null;
        double[] wAsense = null;
        if (true) {	// customMotif
            wAsense = mapWeightAsense.get(id);
            wSense = mapWeightSense.get(id);
        }
        
        double E = d0 + eta * Math.exp(Gamma.logGamma(1d + 1d / delta));
        // determine n, the number of fragments (i.e. (n-1) breakpoints)
        double nn = ((double) len) / E;
        int n = (int) nn;
        double nr = nn - n;
        double r = rndBreak.nextDouble();
        if (r <= nr) {
            ++n;
        }

        // molecule does not break
        List<Fragment> fragments = new ArrayList<Fragment>();
        if (n <= 1 || len <= 1 || (n * d0) >= len) {
            cs.replace(0, start);
            cs.replace(1, end);
            //rw.writeLine(cs, fos);
            //fragments.add(cs.toString());
            fragments.add(new Fragment(id, start, end));
            return fragments;
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
                throw new RuntimeException("Fragment with length < 0! " + nuLen);
            }

            cs.replace(0, nuStart);
            cs.replace(1, nuEnd);
            //rw.writeLine(cs, fos);    // id is invalid now
            //fragments.add(cs.toString());
            
            
            double p= nuStart>= 0&& nuStart<= wSense.length? wSense[nuStart]: 1;
            double q= nuEnd>= 0&& nuEnd<= wAsense.length? wAsense[nuEnd]: 1;
            double rnd= rndDELME.nextDouble();
            double rnd2= rndDELME.nextDouble();
            
            if ((rnd<= p&& rnd2<= q))
            	fragments.add(new Fragment(id, nuStart, nuEnd));

        }
        assert (Math.round(dsum) == len);
        return fragments;
    }

    private double getFragURdelta(double len) {

        if (Double.isNaN(urDelta)) {
            return Math.max(Math.log10(len), 1);
        }
        return urDelta;
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

        if (Double.isNaN(urEta)) {

            double medDelta = getFragURdelta(medMoleculeLength);
            //double d0 = settings.get(FluxSimulatorSettings.FRAG_UR_D0);
            double medFilt = DEFAULT_MED_SIZE;
            if (filtering) {
                medFilt = 170; //getFiltMedian();
            }
            urEta = ((medFilt - d0) / Math.exp(Gamma.logGamma(1d + 1d / medDelta)));
        }

        return urEta;
    }

    @Override
    public String getName() {
        return "Fragmentation UR";
    }

    @Override
    public String getConfiguration() {
        StringBuffer b = new StringBuffer();
        b.append("\t\t").append("D0: ").append(d0).append("\n");
        if(!Double.isNaN(urDelta)){
            b.append("\t\t").append("Delta: ").append(urDelta).append("\n");
        }else{
            b.append("\t\t").append("Delta: ").append(" Not specified, depends on sequence length\n");
        }
        b.append("\t\t").append("Eta: ").append(getFragUReta()).append("\n");
        return b.toString();
    }

    @Override
    public String done() {
        return null;
    }
	public void initPWMMap(){
	    if(pwmASense != null){
	        Log.info("Initializing PWM cache");
	        mapWeightAsense = Fragmenter.getMapWeight(mapTx,null, pwmASense);
	        mapWeightSense = Fragmenter.getMapWeight(mapTx,null, pwmSense);
	        Log.info("Done");
	    }
	}
	
	

}
