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

package barna.commons.parameters;

/**
 * Represents a parameter parsed from string values
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public abstract class Parameter<T> {
    private String name;
    private String description;
    private T defaultValue;
    private Class<T> type;
    private ParameterValidator validator;

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

    protected abstract void parse(String value) throws ParameterException;

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