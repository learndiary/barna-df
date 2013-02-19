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
    public Boolean parse(String value) throws ParameterException {
        String l = value.toLowerCase();
        for (String s : TRUE) {
            if (l.equals(s)) {
                this.value = true;
                return this.value;
            }
        }

        for (String s : FALSE) {
            if (l.equals(s)) {
                this.value = false;
                return this.value;
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
        booleanParameter.longOption(getLongOption()).shortOption(getShortOption());
        booleanParameter.set(get());
        return booleanParameter;
    }
}
