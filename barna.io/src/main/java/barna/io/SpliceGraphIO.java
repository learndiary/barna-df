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

package barna.io;

import barna.io.gtf.GTFwrapper;
import barna.model.*;
import barna.model.commons.MyFile;
import barna.model.constants.Constants;
import barna.model.splicegraph.SplicingGraph;

import java.io.BufferedWriter;
import java.io.File;
import java.io.Writer;
import java.util.Date;

/**
 * IO methods to create a spliece graph
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class SpliceGraphIO {

    public static void extractSpliceJunctions(int eFlankDon, int eFlankAcc, IntronModel iModel, File inF, File outF) {
        if (outF!= null&& outF.exists())
            outF.delete();
        GTFwrapper reader= new GTFwrapper(inF.getAbsolutePath());
        Gene[] g= null;
        try {
            for (reader.read(); (g= reader.getGenes())!= null; reader.read()) {
                for (int i = 0; i < g.length; i++) {
                    SplicingGraph gr= new SplicingGraph(g[i]);
                    gr.constructGraph();
                    gr.writeSpliceJunctionSeqs(eFlankDon, eFlankAcc, iModel, outF);
                }
            }
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }



}
