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

package barna.commons.io;

/**
 * Use this to access the IOHandler implementations
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class IOHandlerFactory {

    /**
     * No instances allowed
     */
    private IOHandlerFactory() {
    }

    /**
     * THe default handler
     */
    private static IOHandler defaultHandler;

    /**
     * Get the default handler
     *
     * @return ioHandler the default handler
     */
    public static IOHandler getDefaultHandler() {
        if (defaultHandler == null) {
            defaultHandler = new SimpleIOHandler();
        }
        return defaultHandler;
    }

    /**
     * Create a new instance of the default handler
     *
     * @return handler new instance of the default handler
     */
    public static IOHandler createDefaultHandler() {
        return new SimpleIOHandler();
    }
}
