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

package fbi.commons.gui;


import fbi.commons.StringUtils;

import javax.swing.*;
import java.awt.*;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import java.util.Random;

public class SimpleBinPlotterPanel extends JPanel {

    public static Color defaultPanBg = Color.white;

    public static final byte MODE_HISTO = 0, MODE_GEL = 1, MODE_BARPLOT = 2, MODE_LINE = 3;
    public static final byte NORMAL = 0, LOGN = 1, LOG2 = 2, LOG10 = 3;
    protected byte paintMode = MODE_HISTO;
    protected boolean invert = false;

    public static void main(String[] args) {
        int[] a = new int[60000];
        Random rnd = new Random();
        for (int i = 0; i < a.length; i++) {
            a[i] = rnd.nextInt(8000);
            if (a[i] == 0) {
                --i;
            }
        }
        JFrame frame = new JFrame();
        frame.setLayout(new BorderLayout());
        SimpleBinPlotterPanel s = new SimpleBinPlotterPanel();
        frame.add(s, BorderLayout.CENTER);
        frame.setSize(800, 600);
        frame.setVisible(true);
        s.paintOSI(a);
        s.repaint();
    }

    byte log = 0;    // 1= logX, 2= logY
    byte modeX = NORMAL, modeY = NORMAL;
    double minXoffset = 0d;
    int tickNoX = 50, tickNoY = 30;
    Insets axesDim = new Insets(15, 50, 30, 0);
    Dimension ticksDim = new Dimension(5, 5);
    String title, titleX = "values", titleY = "frequency";
    boolean paintAxisX = true, paintAxisY = true;
    double thrUp = -1;
    public BufferedImage offScrImg;
    public static final int[] phageHindIIIEcoRI = new int[]{564, 847, 931, 1375, 1584, 1904, 2027, 3530, 4268, 4973, 5149, 21226},
            phageHindIII = new int[]{564, 2027, 2322, 4361, 6557, 9416, 23130},
            systematicLadder = new int[]{25, 200, 500, 1000, 2000, 4000, 8000};
    protected int[] ladder = systematicLadder;
    protected Color plotBg, plotFg;

    public SimpleBinPlotterPanel(String title) {
        this();
        setTitle(title);
    }

    public SimpleBinPlotterPanel() {
        setInvert(invert);
    }

    /**
     * @param data
     * @param nrBins maximal number, gets smaller if width< nrBins
     */
    protected int osiCalled = 0;
    protected Color[] lineColors = new Color[0];
    protected double[] lastXY = null;

    public void paintOSI(Object data) {

        if (offScrImg == null) {
            offScrImg = new BufferedImage(getWidth() - getInsets().left - getInsets().right,
                    getHeight() - getInsets().top - getInsets().bottom,
                    BufferedImage.TYPE_INT_ARGB);
        }
        Graphics gi = offScrImg.getGraphics();
        gi.setColor(plotBg);
        if ((!(paintMode == MODE_LINE)) || (osiCalled == lineColors.length)) {
            gi.fillRect(0, 0, offScrImg.getWidth(), offScrImg.getHeight());
            if (paintMode == MODE_LINE) {
                osiCalled = 0;
                lastXY = null;
            }
        }

        if (data == null) {
            if (paintMode == MODE_HISTO || paintMode == MODE_BARPLOT || paintMode == MODE_LINE) {
                paintOSIhisto(gi, globalBins, globalMaxXY);
            }
        } else {
            MyArrays a = new MyArrays(data);
            if (paintMode == MODE_HISTO || paintMode == MODE_BARPLOT || paintMode == MODE_LINE) {
                double[] maxXY = new double[2];
                double[] bins = getBins(a, maxXY);
                paintOSIhisto(gi, bins, maxXY);
            }
        }
        ++osiCalled;
    }

    /**
     * @param a
     * @deprecated
     */
    public void paintOSIdownscale(long[] a) {
        double maxY = Double.MIN_VALUE;
        for (int i = 0; i < a.length; i++) {
            if (a[i] > maxY) {
                maxY = a[i];
            }
        }

        if (log >= 1) {
            maxY = Math.log10(maxY);
        }
        double scX = a.length / (double) (getWidth() - (axesDim.left + axesDim.right)), scY = maxY / (getHeight() - (axesDim.top + axesDim.bottom));

        offScrImg = new BufferedImage(getWidth(), getHeight(), BufferedImage.TYPE_BYTE_BINARY);
//		BufferedImage tmpImage= new BufferedImage((int) (addX+ a.length), 
//				(int) (addY+ maxY), BufferedImage.TYPE_BYTE_BINARY);

        double eleH = (offScrImg.getHeight() - axesDim.top - axesDim.bottom) / maxY;
        Graphics gi = offScrImg.getGraphics();
        gi.setColor(getBackground());
        gi.fillRect(0, 0, offScrImg.getWidth(), offScrImg.getHeight());
        for (int i = 14000; i < a.length; i++) {
            float y = 0;
            if (log >= 1) {
                y = (float) (eleH * Math.log10(a[i]));
            } else {
                y = (float) (eleH * a[i]);
            }
            int posX = (int) (axesDim.left + 1 * (i - 14000));
            if (posX > offScrImg.getWidth() - axesDim.right) {
                return;
            }

            int posY = (int) (offScrImg.getHeight() - axesDim.bottom - y);
            if (posY < axesDim.top * scY) {
                gi.setColor(Color.red);
                gi.drawLine(posX - 1, (int) (axesDim.top * scY), posX + 1, (int) (axesDim.top));
                posY = (int) (axesDim.top * scY);
                y = (int) (offScrImg.getHeight() - axesDim.top - axesDim.bottom);
            }
            gi.setColor(plotFg);
            gi.drawRect(posX,
                    posY,
                    (int) 1,
                    (int) y);
        }

//		ImageFilter repl= new ReplicateScaleFilter((int) (tmpImage.getWidth(this)/scX), 
//				(int) (tmpImage.getWidth(this)/scY));
//		ImageProducer prod= new FilteredImageSource(tmpImage.getSource(), repl);
//		Image tmpImage2= createImage(prod);
//		offScrImg= new BufferedImage(getWidth(), getHeight(), BufferedImage.TYPE_4BYTE_ABGR);
//		offScrImg.getGraphics().drawImage(tmpImage2,0,0,this);
    }

    protected void paintOSIhisto(Graphics gi, double[] bins, double[] maxXY) {

        if (bins == null) {
            return;
        }

        int w = offScrImg.getWidth() - axesDim.left - axesDim.right, h = offScrImg.getHeight() - axesDim.bottom;    // h includes upper border
        float binWidth = 1;
        if (paintMode == MODE_LINE) {
            binWidth = w / (float) bins.length;    // (getWidth()- axesDim.left)
        }

        float eleHeight = (float) ((h - axesDim.top) / maxXY[1]);

        // bars
        gi.setColor(plotFg);
        for (int i = 0; maxXY[1] >= 1 && i < bins.length; i++) {    // 100421: changed from maxXY[1]>= 1
            if (paintMode == MODE_LINE) {
                paintOSIValline(gi, bins, i, binWidth, eleHeight);
            } else {
                paintOSIValhisto(gi, bins, i, binWidth, eleHeight);
            }
        }

        // axes and legends
        if (paintAxisY) {
            gi.setColor(plotFg);
//			int tickHY= ((int) (eleHeight* maxY/ tickNoY))* tickNoY;
            int tickHY = h - axesDim.top;        // Y axis
            int tickNoYnow = tickHY / tickNoY;
            gi.drawLine(axesDim.left, h - tickHY, axesDim.left, h);
            for (int i = 0; i < tickNoYnow; i++) {    // add 0, not last (cosmetic)
                double tickPos = ((float) tickHY * i) / (tickNoYnow);
                long val = Math.round(maxXY[1] * tickPos / (eleHeight * maxXY[1]));                    // round to next int
                tickPos = ((double) val * tickHY) / maxXY[1];
                String label = Long.toString(val);
                int tickPosInt = (int) (h - tickPos);
                gi.drawLine(axesDim.left - ticksDim.width, tickPosInt, axesDim.left, tickPosInt);
                Rectangle2D labelDim = gi.getFontMetrics().getStringBounds(label, gi);
                gi.drawString(label,
                        (int) (axesDim.left - labelDim.getWidth() - 4),
                        (int) (tickPosInt + labelDim.getHeight() / 2));
            }
        }
        gi.setColor(plotFg);
        String title = titleY;
        Rectangle2D dimY = gi.getFontMetrics().getStringBounds(title, gi);
        gi.drawString(title, Math.max(0, (int) (axesDim.left - dimY.getWidth())), Math.max((int) dimY.getHeight(), (int) (axesDim.top / 2 - dimY.getHeight())));

        if (paintAxisX) {
            gi.setColor(plotFg);
            gi.drawLine(axesDim.left, h, axesDim.left + w, h);    // X axis
//			int tickWX= ((int) (binWidth* bins.length/ tickNoX))* tickNoX;		
            int tickWX = w;
            int tickNoXnow = tickWX / tickNoX;
            for (int i = 0; i < tickNoXnow; i++) {
                float tickPos = ((float) tickWX * i) / (tickNoXnow);
                int binNr = 0 + ((paintMode == MODE_LINE) ? 1 : 0);
                if (i > 0) {
                    binNr = i * bins.length / tickNoXnow + ((paintMode == MODE_LINE) ? 1 : 0);
                }
                int tickPosInt = (int) (axesDim.left + tickPos);
                gi.drawLine(tickPosInt, h, tickPosInt, h + 5);
                String label = null;
                if (log == 2) {
                    label = StringUtils.fprint((binNr * (maxXY[0] / bins.length)) + minXoffset, 2);
                } else {
                    label = Integer.toString((int) ((binNr * (maxXY[0] / bins.length)) + minXoffset)); // +"-"+ Integer.toString((int) ((binNr+ 1)* (maxXY[0]/ bins.length)))
                }
                if (maxXY[0] == 1) {
                    label = StringUtils.fprint((binNr * (maxXY[0] / bins.length)), 2);
                }
                Rectangle2D labelDim = gi.getFontMetrics().getStringBounds(label, gi);
                gi.drawString(label,
                        (int) (tickPosInt),
                        (int) (h + ticksDim.height + 4 + labelDim.getHeight()));
            }
        }
        gi.setColor(plotFg);
        title = titleX;
        Rectangle2D dimX = gi.getFontMetrics().getStringBounds(title, gi);
        gi.drawString(title, (int) (w + axesDim.left + axesDim.right - dimX.getWidth()), (int) (h + 2 + dimX.getHeight()));

        if (this.title != null) {
            gi.setColor(plotFg);
            title = this.title;
            if (globalCntVals > 0) {
                title += " " + Long.toString(globalCntVals);
            }
            Rectangle2D dim = gi.getFontMetrics().getStringBounds(title, gi);
            gi.drawString(title, (int) (getWidth() - dim.getWidth()), axesDim.top);
        }
    }


    protected void paintComponent(Graphics g) {

        super.paintComponent(g);
        if (offScrImg != null) {
            g.drawImage(offScrImg, getInsets().left, getInsets().top, null);
        } else {
            int w = getWidth(), h = getHeight();
            g.setColor(plotBg);
            g.fillRect(getInsets().left, getInsets().top,
                    w - getInsets().left - getInsets().right,
                    h - getInsets().top - getInsets().bottom);
        }
    }

//		@Override
//		public void setPreferredSize(Dimension d) {
//			super.setPreferredSize(d);
//			System.out.println("pref <- "+d);			
//		}

    @Override
    public Dimension getPreferredSize() {
        //Dimension dim= super.getPreferredSize();	// 10x10
        if (offScrImg != null) {
            return new Dimension(offScrImg.getWidth(), offScrImg.getHeight());
        }
        if (paintMode == MODE_GEL) {
            return new Dimension(getInsets().left + getInsets().right +
                    2 * 30 + 3 * 10, Integer.MAX_VALUE);
        }
        if (paintMode == MODE_HISTO) {
            return new Dimension(getInsets().left + getInsets().right + axesDim.left + 200,
                    getInsets().top + getInsets().bottom + axesDim.top + axesDim.bottom + 200);
        }
        return new Dimension(100, 100);    // something, boxlayout does not like Integer.MAX
    }

    @Override
    public Dimension getMinimumSize() {
        if (offScrImg != null) {
            return new Dimension(offScrImg.getWidth(), offScrImg.getHeight());
        }
        return new Dimension(1, 1);
    }


    protected double[] globalMaxXY;    // [0]= maxX as defined in resetOSI(), [1] maxY determined dynamically
    protected double[] globalBins;
    long globalCntVals;

    public void resetOSI(double maxX) {

/*		if (globalBins!= null&& globalCntVals> 0) {
            offScrImg.getGraphics().setColor(plotBg);
            offScrImg.getGraphics().fillRect(0,0,offScrImg.getWidth(),offScrImg.getHeight());
            if (paintMode== MODE_HISTO) {
                paintOSIhisto(offScrImg.getGraphics(), globalBins, this.globalMaxXY);
            } else if (paintMode== MODE_GEL) {
                paintOSIgel(offScrImg.getGraphics(), globalBins, this.globalMaxXY);
            }
        }
        //repaint();
        paintImmediately(getInsets().left, getInsets().top,
                getWidth()- getInsets().left- getInsets().right,
                getHeight()- getInsets().top- getInsets().bottom);
*/
        this.globalCntVals = 0;
        this.globalMaxXY = new double[]{maxX, Double.MIN_VALUE, 1};
        if (paintMode == MODE_GEL) {
            globalMaxXY[0] = 10000;    // max frag length
        }
/*		if (offScrImg== null) {
			offScrImg= new BufferedImage(getWidth()- getInsets().left- getInsets().right,
					getHeight()-getInsets().top-getInsets().bottom, 
					BufferedImage.TYPE_INT_ARGB);
			offScrImg.getGraphics().setColor(plotBg);
			offScrImg.getGraphics().fillRect(0,0,offScrImg.getWidth(),offScrImg.getHeight());
		} 
*/
        if (globalBins == null) {
            if (paintMode == MODE_GEL) {
                int len = getHeight();
                globalBins = new double[len];
            } else {
                int nr = getWidth() - getInsets().left - getInsets().right - axesDim.left - axesDim.right;
                if (nr < 0) {
                    nr = 0;
                    System.err.println("error in resetOSI()");
                }
                globalBins = new double[nr];
            }
        } else {
            for (int i = 0; i < globalBins.length; i++) {
                globalBins[i] = 0;
            }
        }

    }

    public void addVal(double dd) {
        if (thrUp >= 0 && dd > thrUp) {
            return;
        }
        int pos = (int) (dd * (globalBins.length - 1) / globalMaxXY[0]);    // Math.round() f* slow
        if (pos >= globalBins.length || pos < 0) {
            return;    // skip out of zoom
        }

        ++globalBins[pos];
        ++this.globalCntVals;

        double newMax = 0;
        if (modeY == NORMAL) {
            if (pos != 0 && pos != globalBins.length - 1)    // TODO first currently disregarded for max
            {
                newMax = globalBins[pos];
            }
        } else if (modeY == LOGN) {
            newMax = Math.log(globalBins[pos]);
        } else if (modeY == LOG2) {
            newMax = Math.log(globalBins[pos]) / log2transf;
        } else if (modeY == LOG10) {
            newMax = Math.log10(globalBins[pos]);
        }

        if (newMax > globalMaxXY[1]) {
            globalMaxXY[1] = newMax;
//			if (paintMode== MODE_GEL)
//				paintOSIgel(offScrImg.getGraphics(), globalBins, globalMaxXY);
//			else if (paintMode== MODE_HISTO)
//				paintOSIhisto(offScrImg.getGraphics(), globalBins, globalMaxXY);
        } else {
//			if (paintMode== MODE_GEL) {
//				double minLogX= Math.log10(Math.min(globalMaxXY[2], ladder[0]));
//				double maxLogX= Math.log10(Math.max(globalMaxXY[0], ladder[ladder.length-1]));
//				double m= offScrImg.getHeight()/ (minLogX- maxLogX);
//				double t= offScrImg.getHeight()- (m* minLogX);
//
//				paintOSIValgel(offScrImg.getGraphics(), globalBins, pos, globalMaxXY, m, t);
//			} else if (paintMode== MODE_HISTO) {
//				paintOSIValhisto(offScrImg.getGraphics(), globalBins, pos, 
//						1f, (float) ((offScrImg.getHeight()- axesDim.height)/ globalMaxXY[1]));
//			}
        }
//		repaint();
    }

    private static final int gelBandWidth = 3;
    public static final double log2transf = Math.log(2);

    private void paintOSIValhisto(Graphics gi, double[] bins, int pos,
                                  float w, float h) {
        float y = 0;
        if (log >= 1) {
            if (modeY == LOGN) {
                y = (float) (h * Math.log(bins[pos]));
            } else if (modeY == LOG2) {
                y = (float) (h * Math.log(bins[pos] / log2transf));    // Math.log(x)/Math.log(2)
            } else if (modeY == LOG10) {
                y = (float) (h * Math.log10(bins[pos]));    // Math.log(x)/Math.log(2)
            }
        } else {
            y = (float) (h * bins[pos]);
        }
        int posX = (int) (axesDim.left + w * pos);
        if (posX > offScrImg.getWidth()) {
            return;
        }

        int posY = (int) (offScrImg.getHeight() - axesDim.bottom - y);
        if (posY < axesDim.top) {
            gi.setColor(Color.red);
            gi.drawLine(posX - 1, axesDim.top, posX + 1, axesDim.top);
            posY = axesDim.top;
            y = offScrImg.getHeight() - axesDim.top - axesDim.bottom;
        }
        gi.setColor(plotFg);
        gi.drawRect(posX,
                posY,
                (int) w,
                (int) y);
    }

    protected double[] getBins(MyArrays data, double[] maxXY) {

        int w = offScrImg.getWidth() - axesDim.left - axesDim.right,
                h = offScrImg.getHeight() - axesDim.top - axesDim.bottom;

        // max value of future x coordinates
        double maxX = 0l, minX = Double.MAX_VALUE;
        if (paintMode == MODE_BARPLOT || paintMode == MODE_LINE) {
            maxX = data.length();
            minX = 0;
        } else {
            for (int i = 0; i < data.length(); i++) {
                double dd = data.elementAt(i).doubleValue();
                if (thrUp >= 0 && dd > thrUp) {
                    continue;
                }
                if (dd > maxX) {
                    maxX = dd;
                }
                if (dd < minX) {
                    minX = dd;
                }
            }
        }
        if (log == 2) {
            if (modeY == LOGN) {
                maxX = Math.log(maxX);
                if (minX != 0) {
                    minX = Math.log(minX);
                }
            } else if (modeY == LOG2) {
                maxX = Math.log(maxX) / log2transf;
                if (minX != 0) {
                    minX = Math.log(minX) / log2transf;
                }
            } else if (modeY == LOG10) {
                maxX = Math.log10(maxX);
                if (minX != 0) {
                    minX = Math.log10(minX);
                }
            }
        }
//			if (paintMode== MODE_LINE&& lastXY!= null) {
//				minX= 0; 	
//				maxX= lastXY[0];
//			}

        // bins have to be equal in size
        int dataRef = (paintMode == MODE_GEL) ? h : w;
        int nrBins = 0;
        double binSpan = 0;    // x-interval of a bin
        if (maxX <= dataRef) {
            if (log == 2) {
                nrBins = dataRef;
                binSpan = (dataRef / maxX);
            } else {
                nrBins = (int) Math.ceil(maxX);
                binSpan = 1;
            }
        } else {
            binSpan = (int) (maxX / dataRef);
            nrBins = dataRef;
        }
        double[] bins = new double[nrBins];
        for (int i = 0; i < bins.length; i++) {
            bins[i] = 0;
        }

        // do binning
        double maxY = 0l;
        for (int i = 0; i < data.length(); i++) {
            double dd = data.elementAt(i).doubleValue();
            if (thrUp >= 0 && dd > thrUp) {
                continue;
            }
            int pos = 0;
            if (paintMode == MODE_BARPLOT || paintMode == MODE_LINE) {
                if (log == 2) {
                    if (modeY == LOGN) {
                        pos = (int) (Math.log(1 + i) * binSpan);
                    } else if (modeY == LOG2) {
                        pos = (int) (Math.log(1 + i) * binSpan / log2transf);
                    } else if (modeY == LOG10) {
                        pos = (int) (Math.log10(1 + i) * binSpan);
                    }
                } else {
                    pos = (int) (i / binSpan);
                }
                if (pos >= bins.length) {
                    continue; // TODO
                }
                bins[pos] += dd;
            } else {
                pos = (int) ((dd * (bins.length - 1) / maxX));
                if (pos >= bins.length) {
                    continue;    // skip out of zoom
                }
                ++bins[pos];
            }
            if (bins[pos] > maxY) {
                maxY = bins[pos];
            }
        }
        if (log >= 1) {
            if (modeY == LOGN) {
                maxY = Math.log(maxY);
            } else if (modeY == LOG2) {
                maxY = Math.log(maxY) / log2transf;
            } else if (modeY == LOG10) {
                maxY = Math.log10(maxY);
            }

        }

        maxXY[0] = maxX;
        maxXY[1] = maxY;
        if (maxXY.length > 2) {
            maxXY[2] = minX;
        }

        // common scale for lines
        if (paintMode == MODE_LINE && lastXY != null) {
            maxXY[1] = lastXY[1];
        } else {
            lastXY = maxXY;
        }

        return bins;
    }


    public String getTitleX() {
        return titleX;
    }

    public void setTitleX(String titleX) {
        this.titleX = titleX;
    }

    public String getTitleY() {
        return titleY;
    }

    public void setTitleY(String titleY) {
        this.titleY = titleY;
    }

    public void setThrUp(double thrUp) {
        this.thrUp = thrUp;
    }

    public boolean isInvert() {
        return invert;
    }

    public void setInvert(boolean invert) {
        this.invert = invert;
        if (invert) {
            plotBg = Color.black;
            plotFg = Color.white;
        } else {
            plotBg = defaultPanBg;
            plotFg = Color.black;
        }
        setBackground(plotBg);
    }

    public String getTitle() {
        return title;
    }

    public void setTitle(String title) {
        this.title = title;
    }

    public byte getLog() {
        return log;
    }

    public void setLog(byte log) {
        this.log = log;
        modeY = LOGN;
    }

    public void setLog(byte log, byte modeY) {
        this.log = log;
        this.modeY = modeY;
    }

    public byte getPaintMode() {
        return paintMode;
    }

    public void setPaintMode(byte paintMode) {
        this.paintMode = paintMode;
    }

    private void paintOSIValline(Graphics gi, double[] bins, int pos,
                                 float w, float h) {

        if (pos == 0) {
            return;
        }

        float y = 0, y1 = 0;
        if (log >= 1) {
            if (modeY == LOGN) {
                y = (float) (h * Math.log(bins[pos]));
                y1 = (float) (h * Math.log(bins[pos - 1]));
            } else if (modeY == LOG2) {
                y = (float) (h * Math.log(bins[pos]) / log2transf);
                y1 = (float) (h * Math.log(bins[pos - 1] / log2transf));
            } else if (modeY == LOG10) {
                y = (float) (h * Math.log10(bins[pos]));
                y1 = (float) (h * Math.log10(bins[pos - 1]));
            }
        } else {
            y = (float) (h * bins[pos]);
            y1 = (float) (h * bins[pos - 1]);
        }
        int posX = (int) (axesDim.left + w * pos);
        int posX1 = (int) (axesDim.left + w * (pos - 1));
        if (posX > offScrImg.getWidth()) {
            return;
        }

        int posY = (int) (offScrImg.getHeight() - axesDim.bottom - y);
        int posY1 = (int) (offScrImg.getHeight() - axesDim.bottom - y1);
        if (posY < axesDim.top) {
            gi.setColor(Color.red);
            gi.drawLine(posX - 1, axesDim.top, posX + 1, axesDim.top);
            posY = axesDim.top;
            y = offScrImg.getHeight() - axesDim.top - axesDim.bottom;
        }
        gi.setColor(lineColors[osiCalled]);
        gi.drawLine(posX1, posY1, posX, posY);
    }

    public Color[] getLineColors() {
        return lineColors;
    }

    public void setLineColors(Color[] lineColors) {
        this.lineColors = lineColors;
    }

    public boolean isPaintAxisX() {
        return paintAxisX;
    }

    public void setPaintAxisX(boolean paintAxisX) {
        this.paintAxisX = paintAxisX;
        if (!paintAxisX) {
            axesDim.bottom = 1;
        }
    }

    public boolean isPaintAxisY() {
        return paintAxisY;
    }

    public void setPaintAxisY(boolean paintAxisY) {
        this.paintAxisY = paintAxisY;
        if (!paintAxisY) {
            axesDim.left = 1;
        }
    }

    public double getMinXoffset() {
        return minXoffset;
    }

    public void setMinXoffset(double minXoffset) {
        this.minXoffset = minXoffset;
    }


    protected class MyArrays {
        int[] ia;
        long[] la;
        double[] da;
        float[] fa;
        byte[] ba;
        int len;

        public MyArrays(Object array) {
            if (array instanceof int[]) {
                ia = (int[]) array;
                len = ia.length;
            } else if (array instanceof long[]) {
                la = (long[]) array;
                len = la.length;
            } else if (array instanceof double[]) {
                da = (double[]) array;
                len = da.length;
            } else if (array instanceof float[]) {
                fa = (float[]) array;
                len = fa.length;
            } else if (array instanceof byte[]) {
                ba = (byte[]) array;
                len = ba.length;
            }
        }

        public int length() {
            return len;
        }

        public Number elementAt(int pos) {
            if (ia != null) {
                return ia[pos];
            }
            if (la != null) {
                return la[pos];
            }
            if (da != null) {
                return da[pos];
            }
            if (ba != null) {
                return ba[pos];
            }
            //if (fa!= null)
            return fa[pos];
        }
    }


}
