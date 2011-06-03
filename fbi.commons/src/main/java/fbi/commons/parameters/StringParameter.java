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

package fbi.commons.parameters;

import java.util.Collections;
import java.util.List;

/**
 * String parameter implementation
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
class StringParameter extends Parameter<String> {

    private List<String> values;
    private String value;

    public StringParameter(String name) {
        this(name, "", null);
    }

    public StringParameter(String name, String description) {
        this(name, description, "", null);
    }

    public StringParameter(String name, String description, String defaultValue) {
        this(name, description, defaultValue, null);
    }

    public StringParameter(String name, String description, String defaultValue, List<String> values) {
        this(name, description, defaultValue, values, null);
    }

    public StringParameter(String name, String description, String defaultValue, List<String> values, ParameterValidator validator) {
        super(name, description, defaultValue, String.class, validator);
        this.value = defaultValue;
        this.values = values;
    }


    public List<String> getValues() {
        return Collections.unmodifiableList(values);
    }

    @Override
    void set(String value) {
        if (values != null && !values.contains(value)) {
            throw new IllegalArgumentException("Unknown value " + value + ". Must be one of " + values);
        }
        this.value = value;
    }

    String get() {
        return value == null ? getDefault() : value;
    }

    void parse(String value) throws ParameterException {
        this.value = value;
    }

    @Override
    void validate(ParameterSchema schema) throws ParameterException {
        if (getValidator() == null) {
            // if values are available, check that this is valid
            if (values != null && values.size() > 0) {
                if (!values.contains(value)) {
                    throw new ParameterException(this, value, "Parameter " + this + " must be one of " + values);
                }
            }
        } else {
            super.validate(schema);
        }
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
        StringParameter stringParameter = new StringParameter(getName(), getDescription(), getDefault(), getValues(), getValidator());
        stringParameter.set(get());
        return stringParameter;
    }
}
