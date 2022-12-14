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

import barna.commons.parameters.ParameterSchema;
import com.martiansoftware.jsap.JSAPResult;
import com.martiansoftware.jsap.Parameter;

import java.util.List;
import java.util.concurrent.Callable;

/**
 * Base interface for flux tools. Implement this interface to add new flux tool.
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public interface Tool<T> extends Callable<T> {

    /**
     * Provides an unique name for the tool
     * @return the unique name of the tool
     */
    String getName();

    /**
     * Provides a brief description of the tool's
     * functionality
     * @return a brief description of the tool
     */
    String getDescription();

    /**
     * Provides a verbose description of the tool's
     * functionality
     * @return a verbose description of the tool
     */
    String getLongDescription();

    /**
     * List of parameters for the tool
     * @return parameter list
     */
    List<Parameter> getParameter();

    /**
     * Checks CLI parameters with respect to their validity
     * @param args result from parsing the command line
     * @return <code>true</code> if everything is ok with the parameters,
     * <code>false</code> otherwise
     */
    boolean validateParameter(JSAPResult args);

}
