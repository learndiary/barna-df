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

package fbi.genome.errormodel;

/**
 * Wraps around the markov error model and the crosstalk table
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class QualityErrorModel {
    /**
     * the quality model
     */
    private QualityTransitions qualityModel;
    /**
     * the crosstalk table
     */
    private CrossTalkModel crossTalk;

    /**
     * INTERNAL : Empty constructor
     */
    QualityErrorModel() {
    }

    /**
     * Create a new error model
     *
     * @param qualityModel the quality model
     * @param crossTalk the crosstalk table
     */
    public QualityErrorModel(final QualityTransitions qualityModel, final CrossTalkModel crossTalk) {
        this.qualityModel = qualityModel;
        this.crossTalk = crossTalk;
    }

    /**
     * Access the quality model
     *
     * @return transition model
     */
    public QualityTransitions getQualityModel() {
        return qualityModel;
    }

    /**
     * The crosstalk model
     *
     * @return crosstalk crosstalk table
     */
    public CrossTalkModel getCrossTalk() {
        return crossTalk;
    }
}
