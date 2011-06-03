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
    void parse(String value) throws ParameterException {
        try {
            this.value = new BigDecimal(value).longValue();
            return;
        } catch (Exception e) {
        }
        throw new ParameterException(this, value, "Unable to parse parameter " + this + " with value " + value);
    }

    @Override
    public Parameter copy() {
        LongParameter longParameter = new LongParameter(getName(), getDescription(), getDefault(), minimumValue, maximumValue, getValidator());
        longParameter.set(get());
        return longParameter;
    }


}
