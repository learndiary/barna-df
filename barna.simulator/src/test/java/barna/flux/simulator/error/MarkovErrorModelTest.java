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

import barna.model.Qualities;
import com.thoughtworks.xstream.XStream;
import org.junit.Test;
import org.jzy3d.chart.Chart;
import org.jzy3d.colors.Color;
import org.jzy3d.colors.ColorMapper;
import org.jzy3d.colors.colormaps.ColorMapRainbow;
import org.jzy3d.maths.Coord3d;
import org.jzy3d.maths.Coordinates;
import org.jzy3d.plot3d.builder.Builder;
import org.jzy3d.plot3d.builder.concrete.OrthonormalTesselator;
import org.jzy3d.plot3d.builder.delaunay.DelaunayCoordinateValidator;
import org.jzy3d.plot3d.builder.delaunay.DelaunayTessellator;
import org.jzy3d.plot3d.builder.delaunay.jdt.Delaunay_Triangulation;
import org.jzy3d.plot3d.primitives.*;
import org.jzy3d.plot3d.rendering.legends.colorbars.ColorbarLegend;
import org.jzy3d.ui.ChartLauncher;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import static org.junit.Assert.*;
import static org.junit.Assert.assertEquals;

public class MarkovErrorModelTest {

    @Test
    public void testThatErrorModelLoadingWorksAlsoForOldModels() {
        // models before refactoring to barna package
        // (BARNA-86)
        try {
            QualityErrorModel model_35 = MarkovErrorModel.loadErrorModel(
                    new File(getClass().getResource("/before_barna_35.error.model.gz").getFile())
            );
            assertNotNull(model_35);
            assertEquals(35, model_35.getQualityModel().getLength());
        } catch (IOException e) {
            e.printStackTrace();
            fail();
        }
        try {
            QualityErrorModel model_76 = MarkovErrorModel.loadErrorModel(
                    new File(getClass().getResource("/before_barna_76.error.model.gz").getFile())
            );
            assertNotNull(model_76);
            assertEquals(76, model_76.getQualityModel().getLength());
        } catch (IOException e) {
            e.printStackTrace();
            fail();
        }
    }
    @Test
    public void testThatDistributedModelsCanBeLoaded() {
        // models before refactoring to barna package
        // (BARNA-86)
        try {
            QualityErrorModel model_35 = MarkovErrorModel.loadErrorModel(
                    new File(getClass().getResource("/35_error.model").getFile())
            );
            assertNotNull(model_35);
        } catch (IOException e) {
            e.printStackTrace();
            fail();
        }
        try {
            QualityErrorModel model_76 = MarkovErrorModel.loadErrorModel(
                    new File(getClass().getResource("/76_error.model").getFile())
            );
            assertNotNull(model_76);
        } catch (IOException e) {
            e.printStackTrace();
            fail();
        }
    }
    @Test
    public void testModelDist() {
        // models before refactoring to barna package
        // (BARNA-86)
        try {
            QualityErrorModel model_76 = MarkovErrorModel.loadErrorModel(
                    //new File(getClass().getResource("/35_error.model").getFile())
                    new File("/tmp/75_errormodel2.xml")
            );
            assertNotNull(model_76);

            System.out.println(model_76.getReadLength());
            int count = 100000;
            int[][] qs = new int[model_76.getReadLength()][count];
            QualityTransitions qm = model_76.getQualityModel();
            Random rndMutator = new Random();
            for(int i=0;i<model_76.getReadLength();i++){
                for(int j=0; j< count;j++){
                    int last = i == 0 ? 0 : qs[i-1][j];
                    qs[i][j] = qm.getQuality(i, last, rndMutator.nextDouble());
                }
            }

            BufferedWriter w = new BufferedWriter(new FileWriter("/tmp/qm.txt"));
            System.out.println(qs.length);
            System.out.println(qs[0].length);
            for (int i = 0; i < qs[0].length; i++) {
                for (int j = 0; j < qs.length; j++) {
                    w.write(""+qs[j][i]);
                    if(j < qs.length-1) w.write("\t");
                }
                w.write("\n");
            }
            w.close();

            long[][][] t = qm.getTransitions();
            w = new BufferedWriter(new FileWriter("/tmp/trans.txt"));
            for (int x = 0; x < t.length; x++) {
                for (int y = 0; y < t[x].length; y++) {
                    for (int z = 0; z < t[x][y].length; z++) {
                        w.write(""+t[x][y][z]);
                        if(z < t[x][y].length-1) w.write("\t");
                    }
                    w.write("\n");
                }
            }
            w.close();


        } catch (IOException e) {
            e.printStackTrace();
            fail();
        }
    }


    @Test
    public void testThatTheErrorModelWriterUsesOnlyTheSimpleClassName(){
        QualityErrorModel model = new QualityErrorModel(Qualities.Technology.Illumina13, 36, new QualityTransitions(10, 36), new CrossTalkModel(4, true));
        XStream streamer = MarkovErrorModel.createXStream();
        String string = streamer.toXML(model);
        assertTrue(string.startsWith("<QualityErrorModel>"));

        // test reloading
        Object loaded_model = streamer.fromXML(string);
        assertNotNull(loaded_model);

    }
    @Test
    public void testThatBARNA106errorModelLoads() throws IOException {
        QualityErrorModel model = MarkovErrorModel.loadErrorModel(
                new File(getClass().getResource("/BARNA-106-errormodel").getFile())
        );
        assertNotNull(model);

    }
}
