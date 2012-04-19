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

package barna.commons.launcher;

import com.martiansoftware.jsap.JSAPResult;
import com.martiansoftware.jsap.Parameter;

import java.util.List;
import java.util.concurrent.Callable;

/**
 * Base interface for flux tools. Implement this interface to add new flux tool.
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public interface FluxTool<T> extends Callable<T> {
    String getName();
    String getDescription();
    List<Parameter> getParameter();
    boolean validateParameter(JSAPResult args);
}
