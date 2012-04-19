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
