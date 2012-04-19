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

package barna.commons.parameters;

import java.math.BigDecimal;

/**
 * Double parameter implementation
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
class DoubleParameter extends NumberParameter<Double> {
    public DoubleParameter(String name) {
        this(name, "");
    }

    public DoubleParameter(String name, String description) {
        this(name, description, 0d, null, null);
    }

    public DoubleParameter(String name, String description, Double defaultValue) {
        this(name, description, defaultValue, null, null);
    }

    public DoubleParameter(String name, String description, Double defaultValue, Double minimumValue, Double maximumValue) {
        this(name, description, defaultValue, minimumValue, maximumValue, null);
    }

    public DoubleParameter(String name, String description, Double defaultValue, Double minimumValue, Double maximumValue, ParameterValidator validator) {
        super(name, description, defaultValue, minimumValue, maximumValue, Double.class, validator);
    }

    @Override
    protected void parse(String value) throws ParameterException {
        if (value.equalsIgnoreCase("nan")) {
            this.value = Double.NaN;
            return;
        }
        if (value.equalsIgnoreCase("inf") || value.equalsIgnoreCase("infinity")) {
            this.value = Double.POSITIVE_INFINITY;
            return;
        }
        try {
            this.value = new BigDecimal(value).doubleValue();
            return;
        } catch (Exception e) {
        }
        throw new ParameterException(this, value, "Unable to parse parameter " + this + " with value " + value);
    }

    @Override
    public Parameter copy() {
        DoubleParameter doubleParameter = new DoubleParameter(getName(), getDescription(), getDefault(), minimumValue, maximumValue, getValidator());
        doubleParameter.set(get());
        return doubleParameter;
    }
}
