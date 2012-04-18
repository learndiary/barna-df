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

package barna.flux.simulator.distributions;

public class CompositeDistribution extends AbstractDistribution {

    AbstractDistribution[] d = null;

    public CompositeDistribution(AbstractDistribution[] distributions) {
        d = distributions;
    }

    @Override
    public double getP(double x) {
        double p = 0;
        for (int i = 0; i < d.length; i++) {
            p += d[i].getWeight() * d[i].getP(x, mean);
        }
        return p;
    }

    @Override
    public double getRelFreq(double x) {
        double f = 0;
        for (int i = 0; i < d.length; i++) {
            f += d[i].getWeight() * d[i].getRelFreq(x);
        }
        return f;
    }

    @Override
    public double getP(double x, double mean) {
        double p = 0;
        for (int i = 0; i < d.length; i++) {
            p += d[i].getWeight() * d[i].getP(x, mean);
        }
        return p;
    }

    @Override
    public double getMean() {
        if (Double.isNaN(mean)) {
            mean = 0d;
            for (int i = 0; i < d.length; i++) {
                mean += d[i].getWeight() * d[i].getMean();
            }
            mean /= getWeightSum();
        }

        return mean;
    }

    public double getWeightSum() {
        double sum = 0d;
        for (int i = 0; i < d.length; i++) {
            sum += d[i].getWeight();
        }
        return sum;
    }

}
