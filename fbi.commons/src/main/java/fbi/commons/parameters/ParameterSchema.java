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

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.lang.reflect.Field;
import java.util.HashMap;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import fbi.commons.Log;

/**
 * A parameter schema is a set of parameters that can be parsed and validated.
 * The class checks itself for static parameter fields and initializes the parameter map
 * accordingly.
 *
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */

public abstract class ParameterSchema {
    /**
     * Parser pattern
     */
    private static final Pattern PROPERTY_PATTERN = Pattern.compile("(.*)\\s+(.*)");
    /**
     * The parameters
     */
    private Map<String, Parameter> parameters;

    /**
     * Initialize a new schema
     */
    protected ParameterSchema() {
        this.parameters = new HashMap<String, Parameter>();

        /*
        do the reflection magic
         */

        Class<? extends ParameterSchema> clazz = getClass();
        Field[] fields = clazz.getFields();
        for (Field field : fields) {
            if (Parameter.class.isAssignableFrom(field.getType())) {
                try {
                    field.setAccessible(true);
                    Parameter p = (Parameter) field.get(null);
                    if (p != null) {
                        register(p.copy());
                    }
                } catch (IllegalAccessException e) {
                    throw new RuntimeException("Unable to access parameter field " + field.getName() + " in " + clazz.getName());
                }
            }
        }
    }

    /**
     * Returns the number of parameters registered in this set
     *
     * @return size the number of registered parameters
     */
    public int size() {
        return parameters.size();
    }

    /**
     * Manually register a parameter
     *
     * @param parameter the parameter
     */
    public void register(Parameter parameter) {
        if (parameters.containsKey(parameter.getName())) {
            throw new IllegalArgumentException("Parameter " + parameter.getName() + " already exists !");
        }
        this.parameters.put(parameter.getName().toUpperCase(), parameter);
    }

    /**
     * Validate the parameters
     *
     * @throws ParameterException in case a parameter could not be parsed
     */
    public void validate() throws ParameterException {
        for (Map.Entry<String, Parameter> p : parameters.entrySet()) {
            p.getValue().validate(this);
        }
    }

    /**
     * Access a parameter value
     *
     * @param parameter the parameter
     * @param <T> the type
     * @return value the parameter value
     */
    public <T> T get(Parameter<T> parameter) {
        // find the parameter
        Parameter local = parameters.get(parameter.getName().toUpperCase());
        if (local == null) {
            throw new IllegalArgumentException("Unknown parameter '" + parameter.getName() + "'");
        }
        return (T) local.get();
    }

    /**
     * Set a parameter value
     *
     * @param parameter the parameter
     * @param value the value
     * @param <T> the type
     */
    public <T> void set(Parameter<T> parameter, T value) {
        Parameter local = parameters.get(parameter.getName().toUpperCase());
        if (local == null) {
            throw new IllegalArgumentException("Unknown parameter '" + parameter.getName() + "'");
        }
        local.set(value);
    }

    /**
     * Return the {@code key-> value} string for the given parameter. This prints the current value.
     *
     * @param parameter the parameter
     * @return string string consisting of {@code PARAMETERNAME VALUE}
     */
    public String toString(Parameter parameter){
        return parameter.getName() +"\t" +get(parameter);
    }

    /**
     * Write the parameter set to a string
     *
     * @param out the target stream
     */
    @Override
    public String toString() {
    	StringBuilder sb= new StringBuilder();
        for (Map.Entry<String, Parameter> entry : parameters.entrySet()) {
            Parameter p = entry.getValue();
            String name = entry.getKey();
            Object o = get(p);
            if (o== null) 
            	continue;
            sb.append(name);
            sb.append(" ");
            sb.append(o.toString());
            sb.append("\n");
        }
        return sb.toString();
    }

    /**
     * Parse parameters from the input stream
     *
     * @param input the input stream
     * @throws ParameterException in case a parameter could not be parsed
     */
    public void parse(InputStream input) throws ParameterException {
        BufferedReader reader = null;
        try {
            reader = new BufferedReader(new InputStreamReader(input));

            int lineCounter = 0;
            String line = null;

            while ((line = reader.readLine()) != null) {
                line = line.trim();
                lineCounter++;
                if (line.isEmpty()) {
                    continue;
                }
                if (line.startsWith("#")) {
                    continue;
                }


                Matcher matcher = PROPERTY_PATTERN.matcher(line);
                if (!matcher.find() || matcher.groupCount() != 2) {
                    throw new ParameterException("Error while parsing line " + lineCounter + ": " + line);
                }

                String name = matcher.group(1).toUpperCase().trim();
                String value = matcher.group(2).trim();

                if (!parameters.containsKey(name)) {
                    throw new ParameterException("Error while parsing line " + lineCounter + ". Parameter " + name + " not found. Check the spelling!");
                }
                try {
                    parameters.get(name).parse(value);
                } catch (ParameterException p) {
                    throw new ParameterException("Error while parsing line " + lineCounter + ". " + p.getMessage());
                }
            }

        } catch (IOException io) {
            Log.error("Error while reading parameter file : " + io.getMessage(), io);
        } finally {
            if (reader != null) {
                try {reader.close();} catch (IOException ignore) {}
            }
        }
    }

    /**
	 * Write the parameter set to the given stream
	 *
	 * @param out the target stream
	 */
	public void write(OutputStream out) {
	    BufferedWriter writer = null;
	    try {
	        writer = new BufferedWriter(new OutputStreamWriter(out));
	
	
	        for (Map.Entry<String, Parameter> entry : parameters.entrySet()) {
	            Parameter p = entry.getValue();
	            String name = entry.getKey();
	            if (p.getDescription() != null) {
	                writer.write("# " + cleanDescription(p.getDescription()) + "\n");
	            }
	            writer.write("#\n");
	            String valuesString = p.getValuesString();
	            if (valuesString != null && !valuesString.isEmpty()) {
	                writer.write("# " + cleanDescription(valuesString) + " default: " + (p.getDefault() != null ? p.getDefault() : "") + "\n");
	            }
	            Object o = get(p);
	            writer.write(name + "\t" + (o != null ? o : "") + "\n");
	            writer.write("\n");
	        }
	    } catch (IOException e) {
	        try {writer.close();} catch (IOException ignore) {}
	    }finally {
	        if(writer != null){
	            try {
	                writer.close();
	            } catch (IOException e) {
	                e.printStackTrace();
	            }
	        }
	    }
	}

	/**
     * Replace all newlines with \n# to stay in comment mode
     *
     * @param s the source
     * @return clean cleaned description
     */
    private static String cleanDescription(String s) {
        return s.replaceAll("\n", "\n# ");
    }
}
