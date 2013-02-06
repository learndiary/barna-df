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

import java.util.Comparator;

/**
 * Represents a parameter parsed from string values
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public abstract class Parameter<T> {

    /**
     * Compare two parameter instances lexicographically by their name.
     */
    public static class ParameterByNameComparator implements Comparator<Parameter> {
        @Override
        public int compare(Parameter tParameter, Parameter tParameter2) {
            return tParameter.getName().compareTo(tParameter2.getName());
        }
    }

    private String name;
    private String description;
    private T defaultValue;
    private Class<T> type;
    private ParameterValidator validator;
    /**
     * Short version of the CLI name
     */
     private char shortOption;

    public String getLongOption() {
        return longOption;
    }

    public char getShortOption() {
        return shortOption;
    }

    /**
      * Long version of the CLI name
      */
     private String longOption;

    protected Parameter(String name, String description, T defaultValue, Class<T> type, ParameterValidator validator) {
        this.name = name;
        this.description = description;
        this.defaultValue = defaultValue;
        this.type = type;
        this.validator = validator;
    }


    
    /**
     * Lazy clone constructor by micha.
     * @param otherParameter
     */
    protected Parameter(Parameter<T> otherParameter) {
    	this.name= otherParameter.name;
        this.description = otherParameter.description;
        this.defaultValue = otherParameter.defaultValue;
        this.type = otherParameter.type;
        this.validator = otherParameter.validator;
    }

    /**
     * Long version of the CLI name
     * @param longOption
     * @return <code>this</code>
     */
    public Parameter<T> longOption(String longOption){
        this.longOption = longOption;
        return this;
    }

    /**
     * Short version of the CLI name
     * @param shortOption
     * @return <code>this</code>
     */
    public Parameter<T> shortOption(char shortOption){
        this.shortOption = shortOption;
        return this;
    }

    public String getName() {
        return name;
    }

    public String getDescription() {
        return description;
    }

    public Class<T> getType() {
        return type;
    }

    public T getDefault() {
        return defaultValue;
    }

    ParameterValidator getValidator() {
        return validator;
    }

    protected abstract T get();

    protected abstract void set(T value);

    public abstract void parse(String value) throws ParameterException;

    protected void validate(ParameterSchema schema) throws ParameterException {
        if (validator != null) {
            validator.validate(schema, this);
        }
    }

    @Override
    public String toString() {
        return name;
    }

    public String getValuesString() {
        return null;
    }

    /**
     * Clone this parameter
     *
     * @return parameter copy of this parameter
     */
    public abstract Parameter copy();
}
