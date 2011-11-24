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

/**
 * Boolean parameter implementation
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
class BooleanParameter extends Parameter<Boolean> {
    private static final String[] TRUE = {"yes", "1", "true"};
    private static final String[] FALSE = {"no", "0", "false"};
    private boolean value;

    protected BooleanParameter(String name) {
        this(name, "");
    }

    protected BooleanParameter(String name, String description) {
        this(name, description, false);
    }

    protected BooleanParameter(String name, String description, boolean defaultValue) {
        this(name, description, defaultValue, null);
    }

    protected BooleanParameter(String name, String description, boolean defaultValue, ParameterValidator validator) {
        super(name, description, defaultValue, Boolean.class, validator);
        value = defaultValue;
    }

    @Override
    protected void set(Boolean value) {
        this.value = value;
    }

    @Override
    protected Boolean get() {
        return value;
    }

    @Override
    protected void parse(String value) throws ParameterException {
        String l = value.toLowerCase();
        for (String s : TRUE) {
            if (l.equals(s)) {
                this.value = true;
                return;
            }
        }

        for (String s : FALSE) {
            if (l.equals(s)) {
                this.value = false;
                return;
            }
        }
        throw new ParameterException(this, value, "Unable to parse parameter " + this + " with value '" + value + "'. Possible values are: [YES,NO]");
    }

    @Override
    public String getValuesString() {
        return "true|false or yes|no";
    }

    @Override
    public Parameter copy() {
        BooleanParameter booleanParameter = new BooleanParameter(getName(), getDescription(), getDefault(), getValidator());
        booleanParameter.set(get());
        return booleanParameter;
    }
}
