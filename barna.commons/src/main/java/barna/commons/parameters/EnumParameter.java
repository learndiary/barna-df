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

import java.util.Arrays;

/**
 * Enum parameter implementation
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
class EnumParameter<E extends Enum<E>> extends Parameter<E> {

    private E[] values;
    private E value;


    public EnumParameter(String name, String description, E defaultValue) {
        this(name, description, defaultValue, null);
    }

    public EnumParameter(String name, String description, E defaultValue, ParameterValidator validator) {
        this(name, description, defaultValue, defaultValue.getDeclaringClass(), validator);
    }

    public EnumParameter(String name, String description, E defaultValue, Class<E> values, ParameterValidator validator) {
        super(name, description, defaultValue, values, validator);
        this.values = values.getEnumConstants();
    }


    @Override
    protected void set(E value) {
        this.value = value;
    }

    protected E get() {
        return value == null ? getDefault() : value;
    }

    public E parse(String value) throws ParameterException {
        for (E e : values) {
            if (e.name().equalsIgnoreCase(value)) {
                this.value = e;
                return this.value;
            }
        }
        throw new ParameterException(this, value, "Unable to parse parameter " + this + " with value " + value);
    }

    @Override
    public String getValuesString() {
        return Arrays.toString(values);
    }

    @Override
    public Parameter copy() {
        EnumParameter enumParameter = new EnumParameter(getName(), getDescription(), getDefault(), getType(), getValidator());
        enumParameter.longOption(getLongOption()).shortOption(getShortOption());
        enumParameter.set(get());
        return enumParameter;
    }
}
