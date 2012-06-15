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

package barna.flux.simulator.error;

import org.junit.Test;

import static junit.framework.Assert.assertEquals;

/**
 * @author Thasso Griebel <thasso.griebel@gmail.com>
 */
public class ModelPoolTest {

    @Test
    public void testScalingPosition(){
        assertEquals(0, ModelPool.scalePosition(0, 100, 100));
        for (int i = 0; i < 100; i++) {
            assertEquals(i, ModelPool.scalePosition(i, 100, 100));
        }
        assertEquals(49, ModelPool.scalePosition(100, 50, 100));
        assertEquals(24, ModelPool.scalePosition(50, 50, 100));
        assertEquals(0, ModelPool.scalePosition(2, 50, 100));
        assertEquals(48, ModelPool.scalePosition(98, 50, 100));


        assertEquals(99, ModelPool.scalePosition(50, 100, 50));
        assertEquals(49, ModelPool.scalePosition(25, 100, 50));
        assertEquals(3, ModelPool.scalePosition(2, 100, 50));
        assertEquals(95, ModelPool.scalePosition(48, 100, 50));

    }
}
