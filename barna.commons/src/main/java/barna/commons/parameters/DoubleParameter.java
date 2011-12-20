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
