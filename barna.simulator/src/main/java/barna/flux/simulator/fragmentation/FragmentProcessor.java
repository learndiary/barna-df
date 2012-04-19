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

package barna.flux.simulator.fragmentation;

import barna.commons.ByteArrayCharSequence;

import java.util.List;

/**
 * Fragmentation processor
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public interface FragmentProcessor {
    /**
     * Process the given read
     *
     * @param id    the id
     * @param cs    the read
     * @param start the start
     * @param end   the end
     * @param len   the length
     * @return fragments list of fragments or null
     */
    List<Fragment> process(ByteArrayCharSequence id, ByteArrayCharSequence cs, int start, int end, int len);

    /**
     * Return the name of this processor
     *
     * @return name the name
     */
    String getName();

    /**
     * Return the current configuration
     *
     * @return config the configuration
     */
    String getConfiguration();

    /**
     * Called after all fragments are processed. The returned message is printed
     *
     * @return status status or null
     */
    String done();
}
