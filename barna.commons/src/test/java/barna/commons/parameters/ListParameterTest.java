/*
 * Copyright (c) 2012, Micha Sammeth, Thasso Griebel, Emilio Palumbo
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *      * Redistributions of source code must retain the above copyright
 *        notice, this list of conditions and the following disclaimer.
 *      * Redistributions in binary form must reproduce the above copyright
 *        notice, this list of conditions and the following disclaimer in the
 *        documentation and/or other materials provided with the distribution.
 *      * The names of its contributors may be not used to endorse or promote
 *        products derived from this software without specific prior written
 *        permission.
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

import org.junit.Test;

import java.util.Arrays;

import static junit.framework.Assert.*;

/**
 * Test string parameter
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class ListParameterTest {

    @Test
    public void testParse(){
        ListParameter test = new ListParameter("TEST");
        try{
            test.parse("VALUE");
            test.validate(null);
            assertTrue(test.get().size() == 1);
            assertTrue(test.get().get(0) == "VALUE");
        }catch (Exception e){fail();}
        try{
            test.parse("");
            test.validate(null);
            assertEquals("[]", test.get().toString());
        }catch (Exception e){fail();}
        try{
            test.parse(null);
            test.validate(null);
            assertNotNull(test.get());
            assertEquals(test.get().size(), 0);
        }catch (Exception e){e.printStackTrace();fail();}
    }
}
