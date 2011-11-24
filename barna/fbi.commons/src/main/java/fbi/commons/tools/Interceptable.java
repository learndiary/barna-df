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

package fbi.commons.tools;


/**
 * Implementations intercepting the object at specific points. For example,
 * the {@link StreamSorter} implementation can implement this interface to allow the user to intercept
 * the sorted lines without dealing with a separate thread and piped streams.
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public interface Interceptable<T> {
    /**
     * Add an interceptor. Depending on the implementation, interceptors are called
     * to process an object
     *
     * @param interceptor the interceptor
     */
    void addInterceptor(Interceptor<T> interceptor);

    /**
     * Interceptors catch object from the Interceptable. They can return a changed version of the
     * line, but it is up to the Intercaptable to accept the change.
     */
    public static interface Interceptor<T> {
        /**
         * Intercepts an object. Implementations can return a changed version but it is up to the
         * interceptable to accept the change
         *
         * @param object the intercepted object
         * @return object the mutated object
         */
        T intercept(T object);
    }
}
