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

import barna.commons.log.Log;
import barna.commons.system.OSChecker;

import java.io.*;
import java.lang.reflect.Field;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

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
    private static final Pattern PROPERTY_PATTERN = Pattern.compile("([^\\s]+)\\s+(.*)");
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
     * @see ParameterSchema()
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
     * Access a parameter value
     *
     * @param paramName the name of the parameter
     * @param <T> the type
     * @return value the parameter value
     */
    public <T> T get(String paramName) {
        // find the parameter
        Parameter local = parameters.get(paramName.toUpperCase());
        if (local == null) {
            throw new IllegalArgumentException("Unknown parameter '" + paramName + "'");
        }
        return (T) local.get();
    }

    /**
     * Retrieve all parameters
     * @return a map with the identifier x parameter tuples
     */
    public Map<String, Parameter> getParameters() {
        return parameters;
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
     * Set a parameter value
     *
     * @param paramName a <code>String</code> identifying the parameter
     * @param value the value
     * @param <T> the type
     */
    public <T> void set(String paramName, T value) {
        Parameter local = parameters.get(paramName.toUpperCase());
        if (local == null) {
            throw new IllegalArgumentException("Unknown parameter '" + paramName + "'");
        }
        local.set(value);
    }

    /**
     * Set a parameter value
     *
     * @param paramName a <code>String</code> identifying the parameter
     * @param value the value
     * @param <T> the type
     */
    public <T> void set(String paramName, String value) throws ParameterException {
        Parameter local = parameters.get(paramName.toUpperCase());
        if (local == null) {
            throw new IllegalArgumentException("Unknown parameter '" + paramName + "'");
        }
        local.parse(value);
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
            sb.append(OSChecker.NEW_LINE);
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
	

            //SIMULATOR-23 sort the parameter list
            ArrayList<Parameter> params = new ArrayList<Parameter>(parameters.values());
            Collections.sort(params, new Comparator<Parameter>() {
                @Override
                public int compare(Parameter o1, Parameter o2) {
                    return o1.getName().compareTo(o2.getName());
                }
            });

            for (Parameter p : params) {
	            String name = p.getName();
	            if (p.getDescription() != null) {
	                writer.write("# " + cleanDescription(p.getDescription()) + OSChecker.NEW_LINE);
	            }
	            writer.write("#"+OSChecker.NEW_LINE);
	            String valuesString = p.getValuesString();
	            if (valuesString != null && !valuesString.isEmpty()) {
	                writer.write("# " + cleanDescription(valuesString) + " default: " + (p.getDefault() != null ? p.getDefault() : "") + OSChecker.NEW_LINE);
	            }
	            Object o = get(p);
	            writer.write(name + "\t" + (o != null ? o : "") + OSChecker.NEW_LINE);
	            writer.write(OSChecker.NEW_LINE);
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
        return s.replaceAll(OSChecker.NEW_LINE, OSChecker.NEW_LINE+"# ");
    }
}
