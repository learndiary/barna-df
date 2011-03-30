package fbi.genome.sequencing.rnaseq.reconstruction.gui;

import fbi.genome.model.Transcript;
import fbi.genome.model.splicegraph.Edge;
import fbi.genome.model.splicegraph.Graph;
import fbi.genome.model.splicegraph.Node;
import fbi.genome.model.splicegraph.SuperEdge;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.IntervalMarker;
import org.jfree.chart.plot.Marker;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.data.Range;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import javax.swing.*;
import java.awt.*;
import java.util.HashMap;
import java.util.Vector;

public class CoveragePlotter {

	Graph g= null;
	int minMapLen, maxMapLen;
	boolean pairedEnd;
	int seriesCtr= 0;
	double min= Double.MAX_VALUE, max= 0;

	
	public CoveragePlotter(Graph g, boolean pairedEnd, int minMapLen, int maxMapLen) {
		this.g= g;
		this.minMapLen= minMapLen;
		this.maxMapLen= maxMapLen;
		this.pairedEnd= pairedEnd;
	}
	
	public void showPlot(Graph g) {
		Node[] nn= g.getNodesInGenomicOrder();
		XYSeriesCollection allSeries= new XYSeriesCollection();
		HashMap<String, XYSeries> map= new HashMap<String, XYSeries>();
		for (int i = 0; i < nn.length; i++) {
			Vector<Edge> v= nn[i].getOutEdges();
			for (int j = 0; j < v.size(); j++) {
				if (!v.elementAt(j).isExonic())
					continue;
				Edge e= v.elementAt(j);
				if (!pairedEnd)
					addSeries(e, allSeries, map);
				for (int k = 0; e.getSuperEdges()!= null&& k < e.getSuperEdges().size(); k++) {
					SuperEdge se= e.getSuperEdges().elementAt(k);
					if (se.getEdges()[0]!= e&& se.getEdges()[se.getEdges().length- 1]!= e)
						continue;
					if (pairedEnd) {
						if (se.isPend())
							addSeries(se, allSeries, map);
						else {
							for (int m = 0; se.getSuperEdges()!= null&& m < se.getSuperEdges().size(); m++) {
								SuperEdge pe= se.getSuperEdges().elementAt(m);
								assert(pe.isPend());
								if (pe.getEdges()[0]!= se&& pe.getEdges()[pe.getEdges().length- 1]!= se)
									continue;
								addSeries(pe, allSeries, map);
							}
						}
					} else { // not paired-end
						assert(!se.isPend());
						addSeries(se, allSeries, map);
					}
						
				} // end iterate super-edges
			} // end iterate outedges
		} // end iterate nodes
		
		Object[] keys= map.keySet().toArray();
		JFrame myFrame= new JFrame();
		myFrame.getContentPane().setLayout(new BoxLayout(myFrame.getContentPane(), BoxLayout.Y_AXIS));
		System.err.println("min= "+min+", max= "+max);
		Range range= new Range(min, max);
		for (int i = 0; i < keys.length; i++) {
			String s= (String) keys[i];
			if (s.endsWith("_a"))
				continue;			
			XYSeries senseSeries= map.get(s);
			s= s.substring(0, s.length()- 1)+ "a";
			XYSeries asenseSeries= map.get(s);
			XYSeriesCollection bothSeries= new XYSeriesCollection();
			bothSeries.addSeries(senseSeries);
			bothSeries.addSeries(asenseSeries);
			String id= s.substring(0, s.length()- 2);
			JFreeChart chart = ChartFactory.createXYLineChart(id, // Title
                    "genomic position", // x-axis Label
                    "count", // y-axis Label
                    bothSeries, // Dataset
                    PlotOrientation.VERTICAL, // Plot Orientation
                    false, // Show Legend
                    true, // Use tooltips
                    false // Configure chart to generate URLs?
            );
			
			XYPlot plot= chart.getXYPlot();
			for (int x = 0; x < nn.length; x++) {
				Vector<Edge> v= nn[x].getOutEdges();
				for (int j = 0; j < v.size(); j++) {
					if (!v.elementAt(j).isExonic())
						continue;
					Edge e= v.elementAt(j);
					
					addMarker(e, plot, id);
					for (int k = 0; e.getSuperEdges()!= null&& k < e.getSuperEdges().size(); k++) {
						SuperEdge se= e.getSuperEdges().elementAt(k);
						if (se.getEdges()[0]!= e)
							continue;
						if (!se.isPend())
							addMarker(se, plot, id);
/*						if (pairedEnd) {
							if (se.isPend())
								addMarker(se, plot, id);
							else {
								for (int m = 0; se.getSuperEdges()!= null&& m < se.getSuperEdges().size(); m++) {
									SuperEdge pe= se.getSuperEdges().elementAt(m);
									assert(pe.isPend());
									if (pe.getEdges()[0]!= se&& pe.getEdges()[pe.getEdges().length- 1]!= se)
										continue;
									addMarker(pe, plot, id);
								}
							}
						} else { // not paired-end
							assert(!se.isPend());
							addMarker(se, plot, id);
						}
*/						
							
					} // end iterate super-edges
				} // end iterate outedges
			} // end iterate nodes
			
			chart.getXYPlot().getDomainAxis().setAutoRange(false);
			chart.getXYPlot().getDomainAxis().setRange(range, true, false);
			chart.getXYPlot().getDomainAxis().setRange(range);
			ChartPanel panel= new ChartPanel(chart);
			myFrame.getContentPane().add(panel);
		}
		
/*		JFreeChart chart = ChartFactory.createXYLineChart("XY Chart", // Title
					"genomic position", // x-axis Label
					"count", // y-axis Label
					allSeries, // Dataset
					PlotOrientation.VERTICAL, // Plot Orientation
					false, // Show Legend
					true, // Use tooltips
					false // Configure chart to generate URLs?
		);
		XYPlot plot = (XYPlot) chart.getPlot();
		for (int i = 1; i <= seriesCtr; i+=2) {
			Paint paint= plot.getRenderer().getSeriesPaint(i- 1);
			plot.getRenderer().setSeriesPaint(i, 
					paint);
		}
		
		ChartPanel panel= new ChartPanel(chart);
		JFrame myFrame= new JFrame();
		myFrame.getContentPane().add(panel);
*/		
		myFrame.pack();
		myFrame.setVisible(true);
		System.err.println(seriesCtr+ " series");
		
	}
	
	void addSeries_save(Edge e, XYSeriesCollection allSeries, HashMap<String, XYSeries> serMap) {
		int sSense= e.getGpos(true, true, minMapLen, maxMapLen),
			eSense= e.getGpos(true, false, minMapLen, maxMapLen),
			sASense= e.getGpos(false, true, minMapLen, maxMapLen),
			eASense= e.getGpos(false, false, minMapLen, maxMapLen);
		Transcript[] tt= g.decodeTset(e.getTranscripts());
		StringBuilder txString= new StringBuilder(tt.length* 10);
		for (int i = 0; i < tt.length; i++) 
			txString.append(tt[i].getTranscriptID()+",");
		txString.deleteCharAt(txString.length()- 1);
		
		float[] a= e.getCoverageFlat(true, minMapLen, maxMapLen), 
			b= e.getCoverageFlat(false, minMapLen, maxMapLen);
		Vector<SuperEdge> v= e.getSuperEdges();
		HashMap<SuperEdge, String> mapSEID= new HashMap<SuperEdge, String>();
		
		// 
		String seriSname= txString.toString()+ "_s", seriAname= txString.toString()+"_a";
		XYSeries seriS= serMap.get(seriSname), seriA= serMap.get(seriAname);
		if (seriS== null) {
			seriS= new XYSeries(seriSname);
			allSeries.addSeries(seriS);
		}
		if (seriA== null) {
			seriA= new XYSeries(seriAname);
			allSeries.addSeries(seriA);
		}
			
		for (int i = sSense; i <= eASense; ++i) {
			if (i- sSense< eASense) 
				seriS.add(i, a[i- sSense]);
			if (i- sASense>= 0)
				seriA.add(i,b[i- sASense]);
		}
		
	}

	void addMarker(Edge e, XYPlot plot, String id) {
		int sSense= e.getGpos(true, true, minMapLen, maxMapLen),
			eSense= e.getGpos(true, false, minMapLen, maxMapLen),
			sASense= e.getGpos(false, true, minMapLen, maxMapLen),
			eASense= e.getGpos(false, false, minMapLen, maxMapLen);

		Transcript[] tt= g.decodeTset(e.getTranscripts());
		StringBuilder txString= new StringBuilder(tt.length* 10);
//		if (tt.length== 2&& tt[0].getTranscriptID().equals("ENST00000419160")&& tt[1].getTranscriptID().equals("ENST00000425496"))
//			System.currentTimeMillis();
		boolean contained= false;
		for (int i = 0; i < tt.length; i++) { 
			txString.append(tt[i].getTranscriptID()+",");
			if (id.contains(tt[i].getTranscriptID())) {
				contained|= true;
				//break;
			}
		}
		txString.deleteCharAt(txString.length()- 1);
		if (!contained)
			return;
		
//		if (id.equals("ENST00000456623"))
//			System.currentTimeMillis();
		
		Marker m= new IntervalMarker(sSense, eSense);
		if (e instanceof SuperEdge&& ((SuperEdge) e).isPend())
			m.setPaint(Color.blue);
		else
			m.setPaint(Color.green);
		m.setAlpha(0.2f);
		plot.addDomainMarker(m);
		m= new IntervalMarker(sASense, eASense);
		if (e instanceof SuperEdge&& ((SuperEdge) e).isPend())
			m.setPaint(Color.red);
		else
			m.setPaint(Color.yellow);
		m.setAlpha(0.2f);
		plot.addDomainMarker(m);
	}
	
	void addSeries(Edge e, XYSeriesCollection allSeries, HashMap<String, XYSeries> serMap) {
		int sSense= e.getGpos(true, true, minMapLen, maxMapLen),
			eSense= e.getGpos(true, false, minMapLen, maxMapLen),
			sASense= e.getGpos(false, true, minMapLen, maxMapLen),
			eASense= e.getGpos(false, false, minMapLen, maxMapLen);

		min= Math.min(min, sSense);
		min= Math.min(min, sASense);
		max= Math.max(max, eSense);
		max= Math.max(max, eASense);
		
		Transcript[] tt= g.decodeTset(e.getTranscripts());
		StringBuilder txString= new StringBuilder(tt.length* 10);
//		if (tt.length== 2&& tt[0].getTranscriptID().equals("ENST00000419160")&& tt[1].getTranscriptID().equals("ENST00000425496"))
//			System.currentTimeMillis();
		for (int i = 0; i < tt.length; i++) 
			txString.append(tt[i].getTranscriptID()+",");
		txString.deleteCharAt(txString.length()- 1);
		// get series
		String seriSname= txString.toString()+ "_s", seriAname= txString.toString()+"_a";
		XYSeries seriS= serMap.get(seriSname), seriA= serMap.get(seriAname);
		if (seriS== null) {
			seriS= new XYSeries(seriSname);
			allSeries.addSeries(seriS);
			++seriesCtr;
			serMap.put(seriSname, seriS);
		}
		if (seriA== null) {
			seriA= new XYSeries(seriAname);
			allSeries.addSeries(seriA);
			++seriesCtr;
			serMap.put(seriAname, seriA);
		}
		
		float[] a= e.getCoverage(true, minMapLen, maxMapLen), 
			b= e.getCoverage(false, minMapLen, maxMapLen);
		
		// check
		float sums= 0, sumas= 0;
		for (int i = 0; i < a.length; i++) 
			sums+= a[i];
		for (int i = 0; i < b.length; i++) 
			sumas+= b[i];
		if (sums!= sumas)
			System.currentTimeMillis();
		
		for (int i = sSense; i <= eSense; ++i) {
			if (i- sSense< eSense) 
				seriS.add(i, a[i- sSense]);			
		}
		for (int i = sASense; i <= eASense; i++) {
			if (i- sASense< eASense) 
				seriA.add(i, b[i- sASense]);
		}
		
	}
	
}
