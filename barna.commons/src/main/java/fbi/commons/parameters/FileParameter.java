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

import java.io.File;

/**
 * File parameter implementation
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
class FileParameter extends Parameter<File> {

    private File value;
    private FileNameParser nameParser;

    public FileParameter(String name) {
        this(name, "", null);
    }

    public FileParameter(String name, String description) {
        this(name, description, new File("."), null);
    }

    public FileParameter(String name, String description, File defaultValue) {
        this(name, description, defaultValue, null);
    }

    public FileParameter(String name, String description, File defaultValue, ParameterValidator validator) {
        this(name, description, defaultValue, validator, null);
    }

    public FileParameter(String name, String description, File defaultValue, ParameterValidator validator, FileNameParser nameParser) {
        super(name, description, defaultValue, File.class, validator);
        this.nameParser = nameParser;
    }


    @Override
    protected void set(File value) {
        this.value = value;
    }

    protected File get() {
        return value == null ? getDefault() : value;
    }

    protected void parse(String value) throws ParameterException {
        if (nameParser != null) {
            this.value = nameParser.parse(value);
        } else {
            this.value = new File(value);
        }
    }

    @Override
    public Parameter copy() {
        FileParameter fileParameter = new FileParameter(getName(), getDescription(), getDefault(), getValidator(), nameParser);
        fileParameter.set(get());
        return fileParameter;
    }
}
