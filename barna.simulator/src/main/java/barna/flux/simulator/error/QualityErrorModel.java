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

package barna.flux.simulator.error;

import barna.model.Qualities;

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
     * The technology
     */
    private Qualities.Technology technology;
    /**
     * The read length
     */
    private int readLength;

    /**
     * INTERNAL : Empty constructor
     */
    QualityErrorModel() {
    }

    /**
     * Create a new error model
     * @param technology the sequencing technology
     * @param readLength the read length
     * @param qualityModel the quality model
     * @param crossTalk the crosstalk table
     */
    public QualityErrorModel(final Qualities.Technology technology, final int readLength, final QualityTransitions qualityModel, final CrossTalkModel crossTalk) {
        this.technology = technology;
        this.readLength = readLength;
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

    /**
     * Get the technology used to create this model
     * @return tech the technology used to create this model
     */
    public Qualities.Technology getTechnology() {
        return technology;
    }

    /**
     * Get the readlength of this model
     *
     * @return length read length
     */
    public int getReadLength() {
        return readLength;
    }
}
