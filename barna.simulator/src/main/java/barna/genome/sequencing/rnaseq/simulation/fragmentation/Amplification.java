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
import barna.genome.sequencing.rnaseq.simulation.distributions.AbstractDistribution;
import barna.genome.sequencing.rnaseq.simulation.distributions.EmpiricalDistribution;
import barna.genome.sequencing.rnaseq.simulation.distributions.GCPCRDistribution;
import barna.genome.sequencing.rnaseq.simulation.distributions.NormalDistribution;

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

    @Override
    public List<Fragment> process(final ByteArrayCharSequence id, final ByteArrayCharSequence cs, final int start, final int end, final int len) {
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
}
