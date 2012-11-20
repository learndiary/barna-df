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

package barna.flux.simulator.distributions;

import org.junit.Test;

import static junit.framework.Assert.assertEquals;

/**
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class NormalDistributionTest {

    @Test
    public void testName() throws Exception {
        NormalDistribution d = new NormalDistribution(0, 1, false);
        EmpiricalDistribution e = new EmpiricalDistribution(d, 1000, 4.0);
        assertEquals(0.0, e.getP(0), 0.000001);
        assertEquals(0.0, e.getP(1), 0.000001);
        assertEquals(0.0031, e.getP(0.5), 0.0001);
        assertEquals(0.999, e.getSum(), 0.001);
    }

    @Test
    public void testConversion(){
        NormalDistribution n = new NormalDistribution(300, 30);
        EmpiricalDistribution emp = n.toEmpirical(1000);
        assertEquals(n.getP(300), emp.getP(300), 0.0001);
    }


}
