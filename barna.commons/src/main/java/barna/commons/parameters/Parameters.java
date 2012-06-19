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
import java.util.Arrays;
import java.util.Collections;
import java.util.EnumSet;

/**
 * Helper class to create parameters
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class Parameters {

    public static Parameter<String> stringParameter(String name) {
        return stringParameter(name, "");
    }

    public static Parameter<String> stringParameter(String name, String description) {
        return stringParameter(name, description, null, null, (String[]) null);
    }

    public static Parameter<String> stringParameter(String name, String description, String defaultValue, ParameterValidator validator) {
        return stringParameter(name, description, defaultValue, validator, (String[]) null);
    }

    public static Parameter<String> stringParameter(String name, String description, String defaultValue, String... values) {
        return stringParameter(name, description, defaultValue, null, values);
    }

    public static Parameter<String> stringParameter(String name, String description, String defaultValue, ParameterValidator validator, String... values) {
        return new StringParameter(name, description, defaultValue, values == null ? Collections.EMPTY_LIST : Arrays.asList(values), validator);
    }

    public static Parameter<Boolean> booleanParameter(String name) {
        return booleanParameter(name, "");
    }

    public static Parameter<Boolean> booleanParameter(String name, String description) {
        return booleanParameter(name, description, false, null);
    }

    public static Parameter<Boolean> booleanParameter(String name, String description, boolean defaultValue) {
        return booleanParameter(name, description, defaultValue, null);
    }

    public static Parameter<Boolean> booleanParameter(String name, String description, boolean defautlValue, ParameterValidator validator) {
        return new BooleanParameter(name, description, defautlValue, validator);
    }

    public static Parameter<Double> doubleParameter(String name) {
        return doubleParameter(name, "");
    }

    public static Parameter<Double> doubleParameter(String name, String description) {
        return doubleParameter(name, description, 0d, null);
    }

    public static Parameter<Double> doubleParameter(String name, String description, double defaultValue) {
        return doubleParameter(name, description, defaultValue, null);
    }

    public static Parameter<Double> doubleParameter(String name, String description, double defaultValue, ParameterValidator validator) {
        return new DoubleParameter(name, description, defaultValue, null, null, validator);
    }

    public static Parameter<Double> doubleParameter(String name, String description, double defaultValue, double min, double max, ParameterValidator validator) {
        return new DoubleParameter(name, description, defaultValue, min, max, validator);
    }


    public static Parameter<Integer> intParameter(String name) {
        return intParameter(name, "");
    }

    public static Parameter<Integer> intParameter(String name, String description) {
        return intParameter(name, description, 0, null);
    }

    public static Parameter<Integer> intParameter(String name, String description, int defaultValue) {
        return intParameter(name, description, defaultValue, null);
    }


    public static Parameter<Integer> intParameter(String name, String description, int defaultValue, ParameterValidator validator) {
        return new IntegerParameter(name, description, defaultValue, null, null, validator);
    }

    public static Parameter<Integer> intParameter(String name, String description, int defaultValue, int min, int max, ParameterValidator validator) {
        return new IntegerParameter(name, description, defaultValue, min, max, validator);
    }


    public static Parameter<Long> longParameter(String name) {
        return longParameter(name, "");
    }

    public static Parameter<Long> longParameter(String name, String description) {
        return longParameter(name, description, 0l, null);
    }

    public static Parameter<Long> longParameter(String name, String description, long defaultValue) {
        return longParameter(name, description, defaultValue, null);
    }


    public static Parameter<Long> longParameter(String name, String description, long defaultValue, ParameterValidator validator) {
        return new LongParameter(name, description, defaultValue, null, null, validator);
    }

    public static Parameter<Long> longParameter(String name, String description, long defaultValue, long min, long max, ParameterValidator validator) {
        return new LongParameter(name, description, defaultValue, min, max, validator);
    }


    public static Parameter<File> fileParameter(String name) {
        return fileParameter(name, "");
    }

    public static Parameter<File> fileParameter(String name, String description) {
        return fileParameter(name, description, new File("."));
    }

    public static Parameter<File> fileParameter(String name, String description, FileNameParser parser) {
        return fileParameter(name, description, null, null, parser);
    }

    public static Parameter<File> fileParameter(String name, String description, File file) {
        return fileParameter(name, description, file, null, null);
    }


    public static Parameter<File> fileParameter(String name, String description, File file, FileNameParser parser) {
        return fileParameter(name, description, file, null, parser);
    }


    public static Parameter<File> fileParameter(String name, String description, File file, ParameterValidator validator) {
        return new FileParameter(name, description, file, validator);
    }

    public static Parameter<File> fileParameter(String name, String description, File file, ParameterValidator validator, FileNameParser parser) {
        return new FileParameter(name, description, file, validator, parser);
    }


    public static <E extends Enum<E>> Parameter<E> enumParameter(String name, String description, E value) {
        return new EnumParameter<E>(name, description, value);
    }

    public static <E extends Enum<E>> Parameter<E> enumParameter(String name, String description, E value, final ParameterValidator validator) {
        return new EnumParameter<E>(name, description, value, validator);
    }

    public static <E extends Enum<E>> Parameter<EnumSet<E>> enumSetParameter(String name, String description, EnumSet<E> value) {
        return new EnumSetParameter<E>(name, description, value);
    }

    public static <E extends Enum<E>> Parameter<EnumSet<E>> enumSetParameter(String name, String description, EnumSet<E> value, final ParameterValidator validator) {
        return new EnumSetParameter<E>(name, description, value, validator);
    }

}
