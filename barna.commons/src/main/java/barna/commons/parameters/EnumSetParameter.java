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
import java.util.EnumSet;

/**
 * EnumSet parameter implementation
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
class EnumSetParameter<E extends Enum<E>> extends Parameter<EnumSet<E>> {

    private E[] values;
    private EnumSet<E> value;

    public EnumSetParameter(String name, String description, EnumSet<E> defaultValue, Class<E> values) {
        this(name, description, defaultValue, values, null);
    }

    public EnumSetParameter(String name, String description, EnumSet<E> defaultValue, Class<E> values, ParameterValidator validator) {
        super(name, description, defaultValue, (Class<EnumSet<E>>) defaultValue.getClass(), validator);
        this.values = values.getEnumConstants();
        this.value = getDefault().clone();
    }


    @Override
    protected void set(EnumSet value) {
        this.value = value;
    }

    protected EnumSet<E> get() {
        return value == null ? getDefault().clone() : value;
    }

    /**
     * Parses String representation and sets <code>this.value</code>
     *
     * @param value String representaion of the value(s)
     * @throws ParameterException if something went wrong
     */
    public EnumSet<E> parse(String value) throws ParameterException {
        // 2013-01-31 Micha was here:
        // shouldn't it be that the given String representation
        // *overwrites* the default behaviour?!
        // I moved the init of default to constructor, and here
        // I start from scratch
        //if (this.value == null)
        //    this.value = getDefault().clone();
        this.value.removeAll(this.value);
        String[] vals = value.replaceAll("[\\[\\]\\s]", "").split(",");
        try {
            for (String val : vals) {
                if (!val.isEmpty()) {
                    boolean found = false;
                    for (E e : values) {
                        if (e.name().equalsIgnoreCase(val)) {
                            this.value.add(e);
                            found = true;
                            break;
                        }
                    }
                    if(!found){
                        throw new ParameterException(this, value, "Unable to parse parameter " + this + " with value " + val);
                    }
                }
            }
        } catch (Exception e) {
            throw new ParameterException(this, value, "Unable to parse parameter " + this + " with value " + value);
        }
        return this.value;
    }

    @Override
    public String getValuesString() {
        return Arrays.toString(values);
    }

    @Override
    public Parameter copy() {
        EnumSetParameter enumParameter = new EnumSetParameter(getName(), getDescription(), getDefault(),values.getClass().getComponentType(), getValidator());
        enumParameter.longOption(getLongOption()).shortOption(getShortOption());
        enumParameter.set(get());
        return enumParameter;
    }
}
