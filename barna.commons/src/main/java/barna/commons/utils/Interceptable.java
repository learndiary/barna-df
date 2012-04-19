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

package barna.commons.utils;


/**
 * Implementations intercepting the object at specific points. For example,
 * the Sorter implementation can implement this interface to allow the user to intercept
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
