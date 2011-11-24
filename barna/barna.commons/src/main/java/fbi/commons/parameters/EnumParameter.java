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

    protected void parse(String value) throws ParameterException {
        for (E e : values) {
            if (e.name().equalsIgnoreCase(value)) {
                this.value = e;
                return;
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
        enumParameter.set(get());
        return enumParameter;
    }
}
