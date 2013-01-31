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

    public void parse(String value) throws ParameterException {
        if (nameParser != null) {
            this.value = nameParser.parse(value);
        } else {
            this.value = new File(value);
        }
    }

    @Override
    public Parameter copy() {
        FileParameter fileParameter = new FileParameter(getName(), getDescription(), getDefault(), getValidator(), nameParser);
        fileParameter.longOption(getLongOption()).shortOption(getShortOption());
        fileParameter.set(get());
        return fileParameter;
    }
}
