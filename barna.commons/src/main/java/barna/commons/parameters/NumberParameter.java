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
 * General number parameters implementation
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
abstract class NumberParameter<T extends Number> extends Parameter<T> {
    T value;
    protected T minimumValue;
    protected T maximumValue;

    protected NumberParameter(String name, String description, T defaultValue, T minimumValue, T maximumValue, Class<T> type, ParameterValidator validator) {
        super(name, description, defaultValue, type, validator);
        this.minimumValue = minimumValue;
        this.maximumValue = maximumValue;
        if (defaultValue == null) {
            throw new NullPointerException();
        }
    }

    @Override
    protected void set(T value) {
        this.value = value;
    }

    @Override
    protected T get() {
        if (value == null) {
            return getDefault();
        }
        return value;
    }

    @Override
    protected void validate(ParameterSchema schema) throws ParameterException {
        if (minimumValue != null || maximumValue != null) {
            if (value == null) {
                value = getDefault();
            }
            if (value != null && !value.equals(Double.NaN)) {
                if (minimumValue != null) {
                    if (new BigDecimal(value.toString()).compareTo(new BigDecimal(minimumValue.toString())) < 0) {
                        throw new ParameterException(this, value.toString(), "Parameter " + this + " value must be >= " + minimumValue);
                    }
                }

                if (maximumValue != null) {
                    if (new BigDecimal(value.toString()).compareTo(new BigDecimal(maximumValue.toString())) > 0) {
                        throw new ParameterException(this, value.toString(), "Parameter " + this + " value must be <= " + maximumValue);
                    }
                }
            }
        }
        super.validate(schema);
    }

    @Override
    public String getValuesString() {
        if (minimumValue != null && maximumValue != null) {
            return minimumValue + "<= number <= " + maximumValue;
        } else if (minimumValue != null) {
            return minimumValue + "<= number";
        } else if (maximumValue != null) {
            return "number <= " + maximumValue;
        }
        return "number";
    }
}
