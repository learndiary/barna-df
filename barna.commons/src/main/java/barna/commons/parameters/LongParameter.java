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
 * Long parameter implementation
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
class LongParameter extends NumberParameter<Long> {
    public LongParameter(String name) {
        this(name, "");
    }

    public LongParameter(String name, String description) {
        this(name, description, 0l, null, null);
    }

    public LongParameter(String name, String description, Long defaultValue) {
        this(name, description, defaultValue, null, null);
    }

    public LongParameter(String name, String description, Long defaultValue, Long minimumValue, Long maximumValue) {
        this(name, description, defaultValue, minimumValue, maximumValue, null);
    }

    public LongParameter(String name, String description, Long defaultValue, Long minimumValue, Long maximumValue, ParameterValidator validator) {
        super(name, description, defaultValue, minimumValue, maximumValue, Long.class, validator);
    }

    @Override
    public Long parse(String value) throws ParameterException {
        try {
            this.value = new BigDecimal(value).longValue();
            return this.value;
        } catch (Exception e) {
        }
        throw new ParameterException(this, value, "Unable to parse parameter " + this + " with value " + value);
    }

    @Override
    public Parameter copy() {
        LongParameter longParameter = new LongParameter(getName(), getDescription(), getDefault(), minimumValue, maximumValue, getValidator());
        longParameter.longOption(getLongOption()).shortOption(getShortOption());
        longParameter.set(get());
        return longParameter;
    }


}
