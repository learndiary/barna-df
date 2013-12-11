/*
 * Copyright (c) 2012, Micha Sammeth, Thasso Griebel, Emilio Palumbo
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *      * Redistributions of source code must retain the above copyright
 *        notice, this list of conditions and the following disclaimer.
 *      * Redistributions in binary form must reproduce the above copyright
 *        notice, this list of conditions and the following disclaimer in the
 *        documentation and/or other materials provided with the distribution.
 *      * The names of its contributors may be not used to endorse or promote
 *        products derived from this software without specific prior written
 *        permission.
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

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * String parameter implementation
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
class ListParameter extends Parameter<List<String>> {

    private List<List<String>> values;
    private List<String> value;

    public ListParameter(String name) {
        this(name, "", null);
    }

    public ListParameter(String name, String description) {
        this(name, description, new ArrayList<String>(), null);
    }

    public ListParameter(String name, String description, List<String> defaultValue) {
        this(name, description, defaultValue, null);
    }

    public ListParameter(String name, String description, List<String> defaultValue, List<List<String>> values) {
        this(name, description, defaultValue, values, null);
    }

    public ListParameter(String name, String description, List<String> defaultValue, List<List<String>> values, ParameterValidator validator) {
        super(name, description, defaultValue, (Class)List.class, validator);
        this.value = defaultValue;
        this.values = values;
    }


    public List<List<String>> getValues() {
        return Collections.unmodifiableList(values);
    }

    @Override
    protected void set(List<String> value) {
        this.value = value;
    }

    protected List<String> get() {
        return value == null ? getDefault() : value;
    }

    public List<String> parse(String value) throws ParameterException {
        this.value = new ArrayList<String>();
        if (value != null){
            value = value.replaceAll("[\\[\\]\\s]", "");
            Collections.addAll(this.value, value.split(","));
        }
        return this.value;
    }

    @Override
    public String getValuesString() {
        if (values != null) {
            return values.toString();
        }
        return "text";
    }

    @Override
    public Parameter copy() {
        ListParameter stringParameter = new ListParameter(getName(), getDescription(), getDefault(), getValues(), getValidator());
        stringParameter.longOption(getLongOption()).shortOption(getShortOption());
        stringParameter.set(get());
        return stringParameter;
    }
}
