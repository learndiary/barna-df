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

package barna.commons.launcher;

import barna.commons.log.Log;
import barna.commons.utils.TableFormatter;

import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.io.StringWriter;
import java.lang.reflect.Method;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

/**
 * Command line options class.
 * <p/>
 * Command line options can have an arbitrary number of names and must have a description.
 * They are mapped to method names for a given object. The methods are then called using reflection API.
 * <p/>
 * You can use @link{#addOption} to add a boolean parameter or @link{#addParameter} to map a method
 * that evaluates the given parameter string.
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class Options {
    /**
     * The long option prefix
     */
    private static final String CLI_PAR_LONG = "--";
    /**
     * The short option prefix
     */
    private static final String CLI_PAR_SHORT = "-";
    /**
     * The target class
     */
    private Class<?> targetClass;
    /**
     * The target object
     */
    private Object target;
    /**
     * The short options mapping to a method
     */
    private Map<Character, Method> shortOptions;

    /**
     * The long options mapping to a method
     */
    private Map<String, Method> longOptions;

    /**
     * Maps methods to descriptions
     */
    private Map<Method, String> descriptions;

    /**
     * Create a new Options instance for the given target object
     *
     * @param target the target
     */
    public Options(Object target) {
        this.target = target;

        this.targetClass = target.getClass();
        this.shortOptions = new HashMap<Character, Method>();
        this.longOptions = new HashMap<String, Method>();
        this.descriptions = new HashMap<Method, String>();
    }

    /**
     * Add a parameter option. The method that is mapped by the given options must take a single String as parameter!
     *
     * @param method      the method to map to
     * @param description the description
     * @param shortOption the short option (null permitted)
     * @param longOptions the long options (null permitted)
     * @throws RuntimeException in case the method can not be found or a mapping for the method already exists
     */
    public void addParameter(String method, String description, Character shortOption, String... longOptions) {
        if (description == null) {
            throw new NullPointerException("You have to specify a description for option " + shortOption + " " + Arrays.toString(longOptions));
        }
        if (shortOption == null && (longOptions == null || longOptions.length == 0)) {
            throw new IllegalArgumentException("You have to specify at least one option");
        }

        // find the method
        try {
            Method m = targetClass.getDeclaredMethod(method, new Class[]{String.class});
            addMapping(m, description, shortOption, longOptions);
        } catch (NoSuchMethodException e) {
            throw new RuntimeException("Parameter method " + method + " not found in " + targetClass.getName() + ". Make sure the method takes exactly one String parameter!");
        }

    }

    /**
     * Add an option. An option is a boolean switch that takes no parameters.
     * The method that is mapped by the given options must either tak eno parameters or a single boolean parameter.
     *
     * @param method      the method to map to
     * @param description the description
     * @param shortOption the short option (null permitted)
     * @param longOptions the long options (null permitted)
     * @throws RuntimeException in case the method can not be found or a mapping for the method already exists
     */
    public void addOption(String method, String description, Character shortOption, String... longOptions) {
        if (description == null) {
            throw new NullPointerException("You have to specify a description for option " + shortOption + " " + Arrays.toString(longOptions));
        }
        if (shortOption == null && (longOptions == null || longOptions.length == 0)) {
            throw new IllegalArgumentException("You have to specify at least one option");
        }

        // find the method
        try {
            Method m = targetClass.getDeclaredMethod(method, (Class[]) null);
            addMapping(m, description, shortOption, longOptions);
        } catch (NoSuchMethodException e) {
            // check if we find one that takes a boolean
            Method m = null;
            try {
                m = targetClass.getDeclaredMethod(method, new Class[]{boolean.class});
                addMapping(m, description, shortOption, longOptions);
            } catch (NoSuchMethodException e1) {
                throw new RuntimeException("Option method " + method + " not found in " + targetClass.getName() + ". Make sure the method takes exactly one String parameter!", e);
            }
        }
    }


    /**
     * Checks that no mapping exists for the given method and then adds the parameters
     *
     * @param method      the method
     * @param description the description
     * @param shortOption the short parameter
     * @param longOptions the long parameter
     */
    protected void addMapping(Method method, String description, Character shortOption, String[] longOptions) {
        if (method == null) {
            throw new NullPointerException();
        }

        // check that no mapping exists !
        if (descriptions.containsKey(method)) {
            throw new IllegalArgumentException("A parameter mapping for the method " + method + " already exists!");
        }

        descriptions.put(method, description);
        if (shortOption != null) {
            shortOptions.put(shortOption, method);
        }
        if (longOptions != null) {
            for (String longOption : longOptions) {
                if (longOption != null) {
                    this.longOptions.put(longOption, method);
                }
            }
        }
    }

    /**
     * Parse the given parameters
     *
     * @param parameters the parameters
     * @return valid returns true if all parameters could be parsed
     */
    public boolean parse(String[] parameters) throws Exception {
        for (int i = 0; parameters != null && i < parameters.length; i++) {
            if (parameters[i].startsWith(CLI_PAR_LONG)) {
                Method method = longOptions.get(parameters[i].substring(CLI_PAR_LONG.length()));
                if (method == null) {
                    return false;
                }
                i = setParameter(method, parameters, i);
            } else if (parameters[i].startsWith(CLI_PAR_SHORT)) {
                Method m = shortOptions.get(parameters[i].substring(CLI_PAR_SHORT.length()).charAt(0));
                if (m == null) {
                    // see if we find it in long options
                    m = longOptions.get(parameters[i].substring(CLI_PAR_SHORT.length()));
                    if (m == null) {
                        return false;
                    }
                }
                i = setParameter(m, parameters, i);
            } else {
                Log.error("What do you mean by " + parameters[i] + "?\nRunaway argument or bad monday?");
                return false;
            }
        }
        return true;
    }

    /**
     * Apply the parameter
     *
     * @param m    the method
     * @param args all arguments
     * @param i    the index of the parameters
     * @return index next index to use
     * @throws Exception in case the method could not me called
     */
    protected int setParameter(Method m, String[] args, int i) throws Exception {
        int l = m.getParameterTypes().length;
        boolean isBoolean = false;
        if (l == 1 && m.getParameterTypes()[0] == boolean.class) {
            l = 0;
            isBoolean = true;
        }

        Object[] cc = new Object[l];
        if (cc.length + i >= args.length) {
            throw new RuntimeException("Missing arguments for parameter " + args[i] + "!");
        }
        for (int j = 0; j < cc.length; j++) {
            cc[j] = args[i + 1 + j];
        }
        if (!isBoolean) {
            m.invoke(target, cc);
        } else {
            m.invoke(target, true);
        }
        return (i + cc.length);
    }


    /**
     * Print the options usage to the given writer
     *
     * @param out the output writer
     */
    public void printUsage(PrintStream out) {
        TableFormatter tf = new TableFormatter(3);
        tf.addRow("Parameter", "Argument", "Description");
        Iterator<Method> it = descriptions.keySet().iterator();
        while (it.hasNext()) {
            Method m = it.next();
            StringBuilder sb = new StringBuilder("[");
            Object[] oo = shortOptions.entrySet().toArray();
            for (int i = 0; i < oo.length; i++) {
                Map.Entry<Character, Method> en = (Map.Entry<Character, Method>) oo[i];
                if (en.getValue().equals(m)) {
                    sb.append(CLI_PAR_SHORT);
                    sb.append(en.getKey());
                    break;
                }
            }
            oo = longOptions.entrySet().toArray();
            for (int i = 0; i < oo.length; i++) {
                Map.Entry<String, Method> en = (Map.Entry<String, Method>) oo[i];
                if (en.getValue().equals(m)) {
                    if (sb.length() > 1) {
                        sb.append("|");
                    }
                    sb.append(CLI_PAR_SHORT);
                    sb.append(en.getKey());
                    break;
                }
            }
            sb.append("]");
            String s1 = sb.toString();

            sb = new StringBuilder();
            if (m.getParameterTypes() != null && m.getParameterTypes().length > 0) {
                for (int i = 0; i < m.getParameterTypes().length; i++) {
                    sb.append(m.getParameterTypes()[i].getSimpleName().toString());
                    sb.append(",");
                }
                sb.deleteCharAt(sb.length() - 1);
            }
            String s2 = sb.toString();
            String s3 = descriptions.get(m);

            tf.addRow(s1, s2, s3);
        }
        out.println(tf.toString());
    }

    @Override
    public String toString() {
        final StringWriter w = new StringWriter();
        PrintStream o = new PrintStream(new OutputStream() {
            @Override
            public void write(int b) throws IOException {
                w.write(b);
            }
        });
        printUsage(o);
        o.close();
        return w.toString();
    }

}
