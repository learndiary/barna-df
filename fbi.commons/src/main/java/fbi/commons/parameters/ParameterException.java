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
 * Exceptions triggered by the parameters while parsing
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class ParameterException extends Exception {
    /**
     * The parameter
     */
    private Parameter parameter;

    /**
     * Create a new instance
     *
     * @param message the message
     */
    public ParameterException(String message) {
        super(message);
    }

    /**
     * Create a new instance for a parameter that can not parse its value
     *
     * @param causedBy the <code>Exception</code> this <code>ParameterException</code>
     * is deferred from.
     */
    public ParameterException(Exception causedBy) {
    	super(causedBy);
    }

    /**
     * Create a new instance for a parameter that can not parse its value
     *
     * @param parameter the parameter that has a problem
     * @param value the parameter value that caused the problem
     */
    public ParameterException(Parameter parameter, String value) {
        super("Invalid Parameter " + parameter.getName() + " with value " + value);
        this.parameter = parameter;
    }

    /**
     * Create a new instance for a parameter that can not parse its value
     *
     * @param parameter the parameter that has a problem
     * @param value the parameter value that caused the problem
     * @param causedBy the <code>Exception</code> this <code>ParameterException</code>
     * is deferred from.
     */
    public ParameterException(Parameter parameter, String value, Exception causedBy) {
        super("Invalid Parameter " + parameter.getName() + " with value " + value, causedBy);
        this.parameter = parameter;
    }

    /**
     * Create a new instance for a parameter that can not parse its value with a custom message
     *
     * @param parameter the parameter tha has a problem
     * @param value the parameter value that caused the problem
     * @param message the message
     */
    public ParameterException(Parameter parameter, String value, String message) {
        super(message);
        this.parameter = parameter;
    }

}
