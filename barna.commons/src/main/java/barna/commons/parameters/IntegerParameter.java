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
 * Integer parameter implementation
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
class IntegerParameter extends NumberParameter<Integer> {
    public IntegerParameter(String name) {
        this(name, "");
    }

    public IntegerParameter(String name, String description) {
        this(name, description, 0, null, null);
    }

    public IntegerParameter(String name, String description, Integer defaultValue) {
        this(name, description, defaultValue, null, null);
    }

    public IntegerParameter(String name, String description, Integer defaultValue, Integer minimumValue, Integer maximumValue) {
        this(name, description, defaultValue, minimumValue, maximumValue, null);
    }

    public IntegerParameter(String name, String description, Integer defaultValue, Integer minimumValue, Integer maximumValue, ParameterValidator validator) {
        super(name, description, defaultValue, minimumValue, maximumValue, Integer.class, validator);
    }

    @Override
    public Integer parse(String value) throws ParameterException {
        try {
            this.value = new BigDecimal(value).intValue();
            return this.value;
        } catch (Exception e) {
        }
        throw new ParameterException(this, value, "Unable to parse parameter " + this + " with value " + value);
    }

    @Override
    public Parameter copy() {
        IntegerParameter intparameter = new IntegerParameter(getName(), getDescription(), getDefault(), minimumValue, maximumValue, getValidator());
        intparameter.longOption(getLongOption()).shortOption(getShortOption());
        intparameter.set(get());
        return intparameter;
    }

}
