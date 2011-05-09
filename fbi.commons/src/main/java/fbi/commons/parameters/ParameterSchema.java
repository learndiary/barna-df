package fbi.commons.parameters;

import fbi.commons.Log;

import java.io.*;
import java.lang.reflect.Field;
import java.util.HashMap;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */

public abstract class ParameterSchema {
    private static final Pattern PROPERTY_PATTERN = Pattern.compile("(.*)\\s+(.*)");
    private Map<String, Parameter> parameters;

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

    public void register(Parameter parameter) {
        if (parameters.containsKey(parameter.getName())) {
            throw new IllegalArgumentException("Paramter " + parameter.getName() + " already exists !");
        }
        this.parameters.put(parameter.getName().toUpperCase(), parameter);
    }

    public void validate() throws ParameterException {
        for (Map.Entry<String, Parameter> p : parameters.entrySet()) {
            p.getValue().validate(this);
        }
    }


    public <T> T get(Parameter<T> parameter) {
        // find the parameter
        Parameter local = parameters.get(parameter.getName().toUpperCase());
        if (local == null) {
            throw new IllegalArgumentException("Unknown parameter '" + parameter.getName() + "'");
        }
        return (T) local.get();
    }

    public <T> void set(Parameter<T> parameter, T value) {
        Parameter local = parameters.get(parameter.getName().toUpperCase());
        if (local == null) {
            throw new IllegalArgumentException("Unknown parameter '" + parameter.getName() + "'");
        }
        local.set(value);
    }

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
            }
        } catch (IOException e) {
            try {
                writer.close();
            } catch (IOException e1) {
            }
        }
    }

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
                try {
                    reader.close();
                } catch (IOException e) {
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
