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

package barna.flux.capacitor.diffexp;

import org.junit.Test;

import java.io.StringReader;

import static junit.framework.Assert.assertEquals;
import static junit.framework.Assert.assertNotNull;
import static junit.framework.Assert.assertNull;

/**
 * @author Thasso Griebel <thasso.griebel@gmail.com>
 */
public class QuantificationModelTest {
    @Test
    public void testTableParsing() throws Exception {
        StringReader reader = new StringReader(
                "ID\tTYPE\tgene_id\n" +
                "11\tgene\tABC\n" +
                "12\tgene\tDEF\n" +
                "13\ttran\t123");
        QuantificationModel model = QuantificationModel.readTable(reader);
        assertNotNull(model);
        assertNotNull(model.getFeatures("gene"));
        assertNotNull(model.getFeatures("tran"));
        assertNotNull(model.getFeatures("xyz"));

        assertEquals(2, model.getFeatures("gene").size());
        assertEquals(1, model.getFeatures("tran").size());

        assertEquals("ABC", model.get("gene", "11", "gene_id"));
        assertEquals("DEF", model.get("gene", "12", "gene_id"));
        assertEquals("123", model.get("tran", "13", "gene_id"));
        assertNull(model.get("tran", "12", "gene_id-xxx"));
        assertNull(model.get("tran", "1", "gene_id"));
        assertNull(model.get("tran-xxx", "12", "gene_id"));
    }




}
