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
import java.util.Arrays;

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
        return new StringParameter(name, description, defaultValue, Arrays.asList(values), validator);
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

}
