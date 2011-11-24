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

/**
 * Implementations parse strings to file names.
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public interface FileNameParser {
    /**
     * Parse the given string and return a file or null
     *
     * @param string the string to be parsed
     * @return file the file or null
     */
    File parse(String string);
}
