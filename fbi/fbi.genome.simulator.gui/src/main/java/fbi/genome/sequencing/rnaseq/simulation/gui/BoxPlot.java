package fbi.genome.sequencing.rnaseq.simulation.gui;

import fbi.commons.StringUtils;

import javax.swing.*;
import java.awt.*;
import java.awt.image.BufferedImage;
import java.util.Vector;

public class BoxPlot extends JPanel {

	public static int BOX_WIDTH= 10, BOX_DIST= 5;
	Insets axesDim= new Insets(15, 50, 30, 0);
	
	BufferedImage offScrImg;
	Vector<double[]> v= null;
	Vector<String> vNam= null;
	
	@Override
	public Dimension getMinimumSize() {
		if (offScrImg!= null)
			return new Dimension(offScrImg.getWidth(), offScrImg.getHeight());
		return new Dimension(1,1);
	}
	
	@Override
	public Dimension getPreferredSize() {
		//Dimension dim= super.getPreferredSize();	// 10x10
		if (offScrImg!= null)
			return new Dimension(offScrImg.getWidth(), offScrImg.getHeight());	
		if (v!= null&& v.size()>0)
			return new Dimension(v.size()* (BOX_WIDTH+ BOX_DIST), 200);
		return new Dimension(100,100);	
	}
	
	public void addName(String s) {
		if (vNam== null)
			vNam= new Vector<String>();
		vNam.add(s);
	}
	
	public boolean add(double[] bounds) {
		if (bounds== null|| bounds.length!= 5)
			return false;
		if (v== null)
			v= new Vector<double[]>();
		v.add(bounds);
		return true;
	}
	
	@Override
	protected void paintComponent(Graphics g) {
		super.paintComponent(g);
		if (v== null)
			return;
		int h= getHeight();
		double box_width= Math.min(BOX_WIDTH, getWidth()/ v.size());
		double spacing= Math.min(box_width/ 3, BOX_DIST);
		if (spacing== box_width/ 3)
			box_width= box_width* 2/ 3;
		
		double[] minMax= getMinMax();
		double relMax= minMax[1]- minMax[0];
		for (int i = 0; i < v.size(); i++) {
			int x= axesDim.left+ (i* (BOX_WIDTH+ BOX_DIST));
			double[] a= v.elementAt(i);
			int lastY= -1;
			for (int j = 0; j < a.length; j++) {
				int y= (int) Math.round((a[j]* h)/ relMax);
				g.drawLine(x, y, (int) Math.round(x+box_width), y);
				if (lastY>= 0) {
					if (i== 2|| i== 3) {
						g.drawLine(x, y, x, lastY);
						g.drawLine((int) Math.round(x+ box_width), y, (int) Math.round(x+ box_width), lastY);
					} else
						g.drawLine((int) Math.round((x+ box_width)/2), y, (int) Math.round((x+ box_width)/2), lastY);
				}
			}
		}
		
		// axis
		g.drawLine(axesDim.left, axesDim.bottom, axesDim.left, h);
		for (int i = 0; i < 4; i++) {
			double val= (1d/ (i+1))* relMax;
			String s= StringUtils.fprint(val, 2);
			int y= (int) Math.round(val* h/ relMax);
			int sw= g.getFontMetrics().stringWidth(s);
			g.drawString(s, axesDim.left- sw- 5, y- 5);
			g.drawLine(axesDim.left, y, axesDim.left- 3, y);
		}
	}
	
	double[] getMinMax() {
		double[] minMax= new double[2];
		for (int i = 0; v!= null&& i < v.size(); i++) 
			for (int j = 0; j < v.elementAt(i).length; j++) { 
				minMax[0]= Math.min(minMax[0], v.elementAt(i)[j]);
				minMax[1]= Math.max(minMax[0], v.elementAt(i)[j]);
			}
		return minMax;
	}
}
