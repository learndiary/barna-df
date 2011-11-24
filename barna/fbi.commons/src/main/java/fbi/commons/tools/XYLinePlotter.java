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

import com.lowagie.text.Document;
import com.lowagie.text.Rectangle;
import com.lowagie.text.pdf.DefaultFontMapper;
import com.lowagie.text.pdf.PdfContentByte;
import com.lowagie.text.pdf.PdfTemplate;
import com.lowagie.text.pdf.PdfWriter;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.title.TextTitle;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.RectangleInsets;

import java.awt.*;
import java.awt.geom.Rectangle2D;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.OutputStream;

/**
 * Simple XY Line plotter. The plotter follows the builder pattern and yuo can simply concatenated method calls
 * up to {@link #pdf(java.io.File, int, int)}, which writes the plot to a file.
 *
 * @author Thasso Griebel (Thasso.Griebel@googlemail.com)
 */
public class XYLinePlotter {

    /**
     * The data
     */
    private XYSeriesCollection dataCollection;
    /**
     * The chart
     */
    private JFreeChart jfreechart;
    /**
     * X Axis lower bound
     */
    private double xLowerBound = 0;
    /**
     * X Axis upper bound
     */
    private double xUpperBound = 0;
    /**
     * Y Axis lower bound
     */
    private double yLowerBound = 0;
    /**
     * Y Axis upper bound
     */
    private double yUpperBound = 0;

    /**
     * Create a new plot
     *
     * @param title the title (null permitted)
     * @param xAxis the x axis label
     * @param yAxis the y axis label
     */
    public XYLinePlotter(String title, String xAxis, String yAxis) {
        dataCollection = new XYSeriesCollection();
        jfreechart = ChartFactory.createXYLineChart(title, xAxis, yAxis, dataCollection, PlotOrientation.VERTICAL, true, true, false);
    }

    /**
     * Add data. Note that the {@code xs.length() == ys.length()} must always be true.
     *
     * @param name the name
     * @param xs   the x data
     * @param ys   the y date
     * @return plotter the plotter
     */
    public XYLinePlotter dataset(String name, double[] xs, double[] ys) {
        if (xs.length != ys.length) {
            throw new RuntimeException("xs.lengh != ys.length - You always have to provide x and y data of the same length!");
        }
        XYSeries xySeries = new XYSeries(name, true, false);
        for (int i = 0; i < xs.length; i++) {
            xySeries.add(xs[i], ys[i]);
        }
        dataCollection.addSeries(xySeries);
        return this;
    }

    /**
     * Add data. Note that the array must contain exactly two arrays, the first for the x data, the second for
     * the y data.
     *
     * @param name the dataset name
     * @param data the data
     * @return plotter the plotter
     */
    public XYLinePlotter dataset(String name, double[][] data) {
        return dataset(name, data[0], data[1]);
    }

    /**
     * Add a subtitle
     *
     * @param subTitle the subtitle
     * @param font     the font (null permitted)
     * @return plotter the plotter
     */
    public XYLinePlotter subtitle(String subTitle, Font font) {
        if (font != null) {
            jfreechart.addSubtitle(new TextTitle(subTitle, font));
        } else {
            jfreechart.addSubtitle(new TextTitle(subTitle));
        }
        return this;
    }

    /**
     * Set the X axis range. If {@code min == max}, auto-range is used
     *
     * @param min the lower bound
     * @param max the upper bound
     * @return plotter the plotter
     */
    public XYLinePlotter xBounds(double min, double max) {
        if (min > max) {
            throw new RuntimeException("Lower bound > upper bound");
        }
        this.xLowerBound = min;
        this.xUpperBound = max;
        return this;
    }

    /**
     * Set the Y axis range. If {@code min == max}, auto-range is used
     *
     * @param min the lower bound
     * @param max the upper bound
     * @return plotter the plotter
     */
    public XYLinePlotter yBounds(double min, double max) {
        this.yLowerBound = min;
        this.yUpperBound = max;
        return this;
    }

    /**
     * Create the plot
     *
     * @return plot the plot
     */
    protected XYPlot plot() {
        XYPlot plot = jfreechart.getXYPlot();

        if (xLowerBound != xUpperBound) {
            NumberAxis xaxis = (NumberAxis) plot.getDomainAxis();    // xaxis
            xaxis.setLowerBound(xLowerBound);
            xaxis.setUpperBound(xUpperBound);
        }

        if (yLowerBound != yUpperBound) {
            NumberAxis yaxis = (NumberAxis) plot.getRangeAxis();
            yaxis.setLowerBound(yLowerBound);
            yaxis.setUpperBound(yUpperBound);
        }
        return plot;
    }

    /**
     * Create a PDF document for this plot.
     *
     * @param file   the file
     * @param width  the width
     * @param height the height
     * @throws Exception in case of any errors
     */
    public void pdf(File file, int width, int height) throws Exception {
        plot();
        OutputStream os = new BufferedOutputStream(new FileOutputStream(file));
        Rectangle pageSize = new Rectangle(width, height);
        Document doc = new Document(pageSize, 25, 25, 50, 50);
        PdfWriter writer = PdfWriter.getInstance(doc, os);
        doc.open();
        PdfContentByte contentByte = writer.getDirectContent();
        PdfTemplate temp = contentByte.createTemplate(width, height);
        Graphics2D g2d = temp.createGraphics(width, height, new DefaultFontMapper());
        Rectangle2D r2d = new Rectangle2D.Double(0, 0, width, height);
        jfreechart.setPadding(new RectangleInsets(0, 0, 0, 20));
        jfreechart.setBackgroundPaint(Color.WHITE);
        jfreechart.draw(g2d, r2d);
        g2d.dispose();
        contentByte.addTemplate(temp, 0f, 0f);
        doc.close();
    }
}
