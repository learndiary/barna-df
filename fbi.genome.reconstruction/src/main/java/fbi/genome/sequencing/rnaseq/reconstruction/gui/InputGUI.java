package fbi.genome.sequencing.rnaseq.reconstruction.gui;

import fbi.commons.Log;
import fbi.commons.ReadyOrNot;
import fbi.commons.gui.SimpleBinPlotterPanel;
import fbi.commons.thread.StoppableRunnable;
import fbi.genome.io.bed.BEDwrapper;
import fbi.genome.io.gff.GFFReader;
import fbi.genome.model.commons.Distribution;
import fbi.genome.model.constants.Constants;

import javax.swing.*;
import javax.swing.table.DefaultTableModel;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.File;
import java.util.HashMap;
import java.util.Iterator;

public class InputGUI extends JPanel implements ReadyOrNot, StoppableRunnable {

	
	public static void main(String[] args) {
		InputGUI annotGUI= new InputGUI();
		annotGUI.set(null);
		
		JFrame aFrame= new JFrame();
		aFrame.getContentPane().add(annotGUI);
		aFrame.addWindowListener(new WindowAdapter(){
			@Override
			public void windowClosing(WindowEvent e) {
				System.exit(0);
			}
		});
		
		aFrame.setSize(new Dimension(600,400));
		aFrame.setVisible(true);
	}
	
	SimpleBinPlotterPanel lengthDistrPanel, lengthDistrPanel2;
	//AnnotationReader annot;
	
	JTextField refTextField;
	JButton butLoadRef;
	public InputGUI() {
		setLayout(new BorderLayout(10,0));
		add(getPanelRef(), BorderLayout.WEST);
		add(getPanelReads(), BorderLayout.EAST);
	}
	
	static final HashMap<String, Integer> defaultAtrributes= new HashMap<String, Integer>();
	static{
		defaultAtrributes.put("CDS", 43756);
		defaultAtrributes.put("NC", 16312);
	}
	private Object[][] getTableReadLengths(HashMap<Integer, Integer> map) {
		if (map== null)
			return new Object[][] {{36, 98321568}};
		Object[][] obj= new Object[map.size()][];
		Iterator<Integer> iter= map.keySet().iterator();
		int cnt= 0;
		while (iter.hasNext()) {
			Integer x= iter.next();
			obj[cnt]= new Object[2];
			obj[cnt][0]= x;
			obj[cnt][1]= map.get(x);
			++cnt;
		}
		return obj;
	}
	private Object[][] getTableObjects(HashMap<String, Integer> map) {
		Object[][] obj= new Object[map.size()][];
		Iterator<String> iter= map.keySet().iterator();
		for (int i = 0; i < obj.length; i++) {
			obj[i]= new Object[3];
			obj[i][0]= true;
			obj[i][1]= iter.next();
			obj[i][2]= map.get(obj[i][1]);
		}
		return obj;
	}

	JPanel paraPanel, readsPanel;
	JLabel labSize= new JLabel();
	JTextField readsTextField;
	JButton butLoadReads;
	private Component getPanelReads() {
		if (readsPanel == null) {
			readsPanel= new JPanel();
			readsPanel.setLayout(new BorderLayout());
			readsPanel.setBorder(BorderFactory.createTitledBorder("Aligned Reads"));
			
			JPanel panFile= new JPanel();
			panFile.setLayout(new BorderLayout(10,0));
			panFile.add(new JLabel("File"), BorderLayout.WEST);
			readsTextField= new JTextField();
			readsTextField.setEditable(false);
			panFile.add(readsTextField, BorderLayout.CENTER);
			butLoadReads= new JButton("Change");
			butLoadReads.addActionListener(new ActionListener(){
				public void actionPerformed(ActionEvent e) {
				    JFileChooser chooser = FluxCapacitorGUI.singleton.getCommonChooser(); // new JFileChooser();
				    int returnVal = chooser.showOpenDialog(InputGUI.this);
				    if(returnVal == JFileChooser.APPROVE_OPTION) {
				    	File f= chooser.getSelectedFile();
						readsTextField.setText(f.getName());
						FluxCapacitorGUI.singleton.capacitor.fileBED= f;
						new Thread(new Runnable() {
							public void run() {
						    	if (!FluxCapacitorGUI.singleton.capacitor.fileInitBED()) {
                                    Log.progressFinish();
							    	return;
						    	}
                                Log.progressFinish();
						    	BEDwrapper reader= FluxCapacitorGUI.singleton.capacitor.getBedReader();
						    	tfMapAll.setText(Integer.toString(reader.getCountAll()));
						    	tfMapEntire.setText(Integer.toString(reader.getCountEntire()));
						    	tfMapSplit.setText(Integer.toString(reader.getCountSplit()));
						    	
						    	tfReadAll.setText(Integer.toString(reader.getCountReads()));
						    	tfReadMap.setText(Integer.toString(reader.getCountReads()));
						    	tfRfac.setText(Float.toString(reader.getCountAll() / reader.getCountReads()));
						    	
								FluxCapacitorGUI.singleton.uiUpdate(null);
							}
						}).start();
				    }				
				}
			});
			butLoadReads.setEnabled(true);
			panFile.add(butLoadReads, BorderLayout.EAST);
			readsPanel.add(panFile, BorderLayout.NORTH);

			JPanel panAtr= new JPanel(new BorderLayout());
			panAtr.setBorder(BorderFactory.createTitledBorder(BorderFactory.createEmptyBorder(), "File Size"));
			DefaultTableModel dtm= new DefaultTableModel(getTableReadLengths(null), 
					new Object[] {"Length", "Reads"});
			JTable tab= new JTable(dtm);
			JScrollPane attrScroller= new JScrollPane(tab);
			attrScroller.setPreferredSize(new Dimension(200,10));
			panAtr.add(attrScroller, BorderLayout.CENTER);
			readsPanel.add(panAtr, BorderLayout.WEST);

			
			JPanel pan= new JPanel();
			pan.setLayout(new BoxLayout(pan, BoxLayout.Y_AXIS));
			JPanel panFrame= new JPanel();
			panFrame.setBorder(BorderFactory.createTitledBorder(BorderFactory.createEmptyBorder(), "Insert Sizes (genomic)"));
			lengthDistrPanel2= new SimpleBinPlotterPanel();
			lengthDistrPanel2.setTitleX("length [nt]");
			panFrame.add(lengthDistrPanel2);
			pan.add(panFrame);
			pan.add(getPanelReadAttr());
			readsPanel.add(pan, BorderLayout.CENTER);
			
		}
		return readsPanel;
	}
	
	JPanel statsPanel;
	JTextField tfFSize, tfTx, tfLoc, tfLenAvg, tfLenMed, tfASAvg, tfASMed, tfAS1Q, tfAS3Q, tfASSTD, tfLenMin, tfLenMax, tfLen1Q, tfLen3Q, tfLenSTD;
	JCheckBox cbLoadCoding, cbLoadNoncoding;
	private Component getPanelRefAttr() {
		if (statsPanel == null) {
			statsPanel = new JPanel();
			statsPanel.setBorder(BorderFactory.createTitledBorder(BorderFactory.createEmptyBorder(), "Attributes"));
			statsPanel.setLayout(new GridBagLayout());
			GridBagConstraints c= new GridBagConstraints();
			c.anchor= GridBagConstraints.FIRST_LINE_START;
			c.fill = GridBagConstraints.NONE;
			c.insets= new Insets(1,10,1,10);
			c.gridwidth= 1; c.gridheight= 1;

			c.weighty= 0.01;
			c.gridx = 0; 
			c.gridy = 0;
//			statsPanel.add(new JLabel("File Size"), c);
			
			c.gridy = 0;
			statsPanel.add(new JLabel("Spliceforms"), c);
			c.gridy = 1;
			statsPanel.add(new JLabel("Loci"), c);
			c.gridy = 2;
			statsPanel.add(new JLabel("Spliceforms per Locus"), c);
			c.gridy = 3;
			c.insets= new Insets(1,20,1,10);
			statsPanel.add(new JLabel("Average"), c);
			c.gridy = 4;
			c.insets= new Insets(1,20,1,10);
			statsPanel.add(new JLabel("Std Deviation"), c);
			c.gridy = 5;
			c.insets= new Insets(1,20,1,10);
			statsPanel.add(new JLabel("1st Quartile"), c);
			c.gridy = 6;
			c.insets= new Insets(1,20,1,10);
			statsPanel.add(new JLabel("Median"), c);
			c.gridy = 7;
			c.insets= new Insets(1,20,1,10);
			statsPanel.add(new JLabel("3rd Quartile"), c);
			/*c.gridy = 8;
			c.insets= new Insets(1,10,1,10);
			statsPanel.add(new JLabel("Spliceform Length"), c);
			c.gridy = 9;
			c.insets= new Insets(1,20,1,10);
			statsPanel.add(new JLabel("Minimum"), c);
			c.gridy = 10;
			c.insets= new Insets(1,20,1,10);
			statsPanel.add(new JLabel("Average"), c);
			c.gridy = 11;
			c.insets= new Insets(1,20,1,10);
			statsPanel.add(new JLabel("Std Deviation"), c);
			c.gridy = 12;
			c.insets= new Insets(1,20,1,10);
			statsPanel.add(new JLabel("Maximum"), c);
			c.gridy = 13;
			c.insets= new Insets(1,20,1,10);
			statsPanel.add(new JLabel("1st Quartile"), c);
			c.gridy = 14;
			c.insets= new Insets(1,20,1,10);
			statsPanel.add(new JLabel("Median"), c);
			c.gridy = 15;
			c.insets= new Insets(1,20,1,10);
			statsPanel.add(new JLabel("3rd Quartile"), c);
*/			
			int tfWidth= 8;
			c.gridx = 2; 
			c.insets= new Insets(1,10,1,10);
			tfFSize= new JTextField(tfWidth);
			tfFSize.setHorizontalAlignment(JTextField.RIGHT);
			tfFSize.setEditable(false);
			c.gridy = 0;
//			statsPanel.add(tfFSize, c);
			tfTx= new JTextField(tfWidth);
			tfTx.setHorizontalAlignment(JTextField.RIGHT);
			tfTx.setEditable(false);
//			c.gridy = 1;
			statsPanel.add(tfTx, c);
			tfLoc= new JTextField(tfWidth);
			tfLoc.setHorizontalAlignment(JTextField.RIGHT);
			tfLoc.setEditable(false);
			c.gridy = 1;
			statsPanel.add(tfLoc, c);
			
			tfASAvg= new JTextField(tfWidth);
			tfASAvg.setHorizontalAlignment(JTextField.RIGHT);
			tfASAvg.setEditable(false);
			c.gridy = 3;
			statsPanel.add(tfASAvg, c);
			tfASSTD= new JTextField(tfWidth);
			tfASSTD.setHorizontalAlignment(JTextField.RIGHT);
			tfASSTD.setEditable(false);
			c.gridy = 4;
			statsPanel.add(tfASSTD, c);
			tfAS1Q= new JTextField(tfWidth);
			tfAS1Q.setHorizontalAlignment(JTextField.RIGHT);
			tfAS1Q.setEditable(false);
			c.gridy = 5;
			statsPanel.add(tfAS1Q, c);
			tfASMed= new JTextField(tfWidth);
			tfASMed.setHorizontalAlignment(JTextField.RIGHT);
			tfASMed.setEditable(false);
			c.gridy = 6;
			statsPanel.add(tfASMed, c);
			tfAS3Q= new JTextField(tfWidth);
			tfAS3Q.setHorizontalAlignment(JTextField.RIGHT);
			tfAS3Q.setEditable(false);
			c.gridy = 7;
			statsPanel.add(tfAS3Q, c);
			
			tfLenMin= new JTextField(tfWidth);
			tfLenMin.setHorizontalAlignment(JTextField.RIGHT);
			tfLenMin.setEditable(false);
/*			c.gridy = 9;
			statsPanel.add(tfLenMin, c);
			c.gridy = 10;
			tfLenAvg= new JTextField(tfWidth);
			tfLenAvg.setHorizontalAlignment(JTextField.RIGHT);
			tfLenAvg.setEditable(false);
			statsPanel.add(tfLenAvg, c);
			tfLenSTD= new JTextField(tfWidth);
			tfLenSTD.setHorizontalAlignment(JTextField.RIGHT);
			tfLenSTD.setEditable(false);
			c.gridy = 11;
			statsPanel.add(tfLenSTD, c);
			tfLenMax= new JTextField(tfWidth);
			tfLenMax.setHorizontalAlignment(JTextField.RIGHT);
			tfLenMax.setEditable(false);
			c.gridy = 12;
			statsPanel.add(tfLenMax, c);
			tfLen1Q= new JTextField(tfWidth);
			tfLen1Q.setHorizontalAlignment(JTextField.RIGHT);
			tfLen1Q.setEditable(false);
			c.gridy = 13;
			statsPanel.add(tfLen1Q, c);
			tfLenMed= new JTextField(tfWidth);
			tfLenMed.setHorizontalAlignment(JTextField.RIGHT);
			tfLenMed.setEditable(false);
			c.gridy = 14;
			statsPanel.add(tfLenMed, c);
			tfLen3Q= new JTextField(tfWidth);
			tfLen3Q.setHorizontalAlignment(JTextField.RIGHT);
			tfLen3Q.setEditable(false);
			c.gridy = 15;
			statsPanel.add(tfLen3Q, c);
			
			c.insets= new Insets(1,10,1,10);
			c.fill = GridBagConstraints.VERTICAL;
			c.gridx = 0; 
			c.gridy = 16;
			c.weighty= 1;
			c.gridwidth= 3;
			statsPanel.add(new JPanel(), c);
*/
		}

		return statsPanel;
	}

	public boolean isReady() {
		//if (settings== null|| settings.getProfiler()== null) 
//		if (settings== null|| settings.getParFile()== null|| (!settings.getParFile().exists())
//				|| (!settings.getParFile().canRead()))
//			return false;
		return true;	// settings.getProfiler().isReady()
	}
	
	public void set(Object o) {
		
	}

	public boolean setStop() {
//		if (settings.getProfiler()== null)
//			return false;
//		return settings.getProfiler().setStop(true);
		return true;
	}

	public void run() {		
				
//		settings.getProfiler().setStop(false);
//		if (loadStats) {
//			loadStats();
//			setLoadStats(false);
//			return;
//		}
//		
//		if (!settings.getProfiler().readAnnotation())
//			return;
//		
//		if (!settings.getProfiler().isStop()) {
//			settings.getProfiler().writeProfile();
//			settings.save();
//		}
		
		uiUpdateStats();
	}
	
	private void uiUpdatePars() {
//		if (settings!= null&& settings.getRefFile()!= null) {
//			refTextField.setText(settings.getRefFile().getAbsolutePath());
//			//butLoadRef.setEnabled(true);
//			
//			cbLoadCoding.repaint();
////			cbLoadCoding.setEnabled(true);
//			cbLoadCoding.revalidate();
//			cbLoadNoncoding.revalidate();
//		} else
			refTextField.setText(Constants.EMPTYSTRING);
		
		FluxCapacitorGUI.repaintEverywhere(this, false);
	}
	
	private void uiUpdateStats() {
		
		//if (settings== null|| settings.getProfiler()== null|| (!settings.getProfiler().isFinishedReadAnnotation())) {
		if (true) {
			tfFSize.setText(Constants.EMPTYSTRING);
			tfTx.setText(Constants.EMPTYSTRING);
			tfLoc.setText(Constants.EMPTYSTRING);
			tfASAvg.setText(Constants.EMPTYSTRING);
			tfASMed.setText(Constants.EMPTYSTRING);
			tfAS1Q.setText(Constants.EMPTYSTRING);
			tfAS3Q.setText(Constants.EMPTYSTRING);
			tfASSTD.setText(Constants.EMPTYSTRING);
			tfLenMin.setText(Constants.EMPTYSTRING);
			tfLenMax.setText(Constants.EMPTYSTRING);
			tfLenMed.setText(Constants.EMPTYSTRING);
			tfLenSTD.setText(Constants.EMPTYSTRING);
			tfLen1Q.setText(Constants.EMPTYSTRING);
			tfLen3Q.setText(Constants.EMPTYSTRING);
			tfLenAvg.setText(Constants.EMPTYSTRING);
			lengthDistrPanel.offScrImg= null;
			lengthDistrPanel2.offScrImg= null;

		} else {
//			tfFSize.setText(MyFile.humanReadableSize(settings.getRefFile().length()));
//			tfTx.setText(Integer.toString(settings.getProfiler().getLen().length));
//			tfLoc.setText(Integer.toString(settings.getProfiler().getCntLoci()));	//MyFormatter.american(
//			tfASAvg.setText(MyFormatter.fprint(settings.getProfiler().getTxLocAvg(), 2));
//	//		String s;
//	//		int p;
//	//		s= MyFormatter.fprint(profiler.getTxLocMed(), 1);
//	//		p= s.indexOf(".");
//	//		if (p>= 0)
//	//			s= MyFormatter.american(Integer.parseInt(s.substring(0,p)))+ "."+ s.substring(p+1);
//			tfASMed.setText(MyFormatter.fprint(settings.getProfiler().getTxLocMed(), 2));
//			tfAS1Q.setText(MyFormatter.fprint(settings.getProfiler().getTxLoc1Q(), 2));
//			tfAS3Q.setText(MyFormatter.fprint(settings.getProfiler().getTxLoc3Q(), 2));
//			tfASSTD.setText(MyFormatter.fprint(settings.getProfiler().getTxLocSTD(), 2));
//			tfLenMin.setText(MyFormatter.fprint(settings.getProfiler().getLenMin(), 2));
//			tfLenMax.setText(MyFormatter.fprint(settings.getProfiler().getLenMax(), 2));
//			tfLenMed.setText(MyFormatter.fprint(settings.getProfiler().getLenMed(), 2));
//			tfLenSTD.setText(MyFormatter.fprint(settings.getProfiler().getLenSTD(), 2));
//			tfLen1Q.setText(MyFormatter.fprint(settings.getProfiler().getLen1Q(), 2));
//			tfLen3Q.setText(MyFormatter.fprint(settings.getProfiler().getLen3Q(), 2));
//			tfLenAvg.setText(MyFormatter.fprint(settings.getProfiler().getLenAvg(), 2));
//			lengthDistrPanel.paintOSI(settings.getProfiler().getLen());
//			lengthDistrPanel.repaint();
//			lengthDistrPanel2.setThrUp(settings.getProfiler().getLen3Q());
//			lengthDistrPanel2.paintOSI(settings.getProfiler().getLen());
//			lengthDistrPanel2.repaint();
		}
		
		FluxCapacitorGUI.repaintEverywhere(this, false); 
	}
	
	public boolean loadStats() {
//		if (settings.getProfiler().getIds()== null|| settings.getProfiler().getLen()== null) {
//			settings.getProfiler().setStop(false);
//			if (!settings.getProfiler().loadStats())
//				return false;
//		}
//		if (settings.getProfiler().isFinishedReadAnnotation()) {
//			uiUpdateStats();
//		}
		return true;
	}
	
	public boolean isFinished() {		
		return FluxCapacitorGUI.singleton.capacitor.isInputReady();
	}

	public boolean getLoadStats() {
		return loadStats;
	}

	boolean loadStats= false;
	public void setLoadStats(boolean val) {
		loadStats= val;
	}

	public void killResult() {
//		if (settings!= null&& settings.getProfiler()!= null) {
//			settings.setProfiler(null);
//		}
//		if (settings!= null&& settings.getProFile()!= null&& settings.getProFile().exists()) {
//			settings.getProFile().delete();
//		}
//		set(settings);
		lengthDistrPanel.offScrImg= null;
		lengthDistrPanel2.offScrImg= null;
		//uiUpdate();
	}

	public boolean isStop() {		
		return false; //settings== null|| settings.getProfiler()== null|| settings.getProfiler().isStop();
	}

	public boolean setStop(boolean stop) {
//		if (settings.getProfiler()== null)
//			return false;
//		return settings.getProfiler().setStop(stop);
		return true;
	}

	private Component getPanelRef() {
			if (paraPanel == null) {								
				paraPanel= new JPanel();
				paraPanel.setBorder(BorderFactory.createTitledBorder("Reference Annotation"));
				paraPanel.setLayout(new BorderLayout());
				
				JPanel panFile= new JPanel();
				panFile.setLayout(new BorderLayout(10,0));
				panFile.add(new JLabel("File"), BorderLayout.WEST);
				refTextField= new JTextField();
				refTextField.setEditable(false);
				panFile.add(refTextField, BorderLayout.CENTER);
				butLoadRef= new JButton("Change");
				butLoadRef.addActionListener(new ActionListener(){
					public void actionPerformed(ActionEvent e) {
					    JFileChooser chooser = FluxCapacitorGUI.singleton.getCommonChooser(); // new JFileChooser();
					    int returnVal = chooser.showOpenDialog(InputGUI.this);
					    if(returnVal == JFileChooser.APPROVE_OPTION) {
					    	File f= chooser.getSelectedFile();
							refTextField.setText(f.getName());
							refTextField.repaint();
							FluxCapacitorGUI.singleton.capacitor.fileGTF= f;
							new Thread(new Runnable() {
								public void run() {
									if (!FluxCapacitorGUI.singleton.capacitor.fileInitReference()) {
                                        Log.progressFinish();
								    	return;
									}
                                    Log.progressFinish();

							    	GFFReader reader= FluxCapacitorGUI.singleton.capacitor.getGTFreader();
							    	tfLoc.setText(Integer.toString(reader.getReadGenes()));
							    	tfTx.setText(Integer.toString(reader.getReadTranscripts()));
							    	
							    	Distribution dist= new Distribution(reader.getTxPerLocus());
							    	tfASAvg.setText(Float.toString(((int) (dist.getMean()* 10))/ 10f));
							    	tfASSTD.setText(Float.toString(((int) (dist.getStandardDeviation()* 10)/ 10f)));
							    	tfAS1Q.setText(Float.toString((float) dist.get1stQuart())); 
							    	tfASMed.setText(Float.toString((float) dist.getMedian())); 
							    	tfAS3Q.setText(Float.toString((float) dist.get3rdQuart()));
							    	
							    	lengthDistrPanel.offScrImg= null;
							    	lengthDistrPanel.paintOSI(reader.getTxLengths());
									lengthDistrPanel.repaint();
							    	
									FluxCapacitorGUI.singleton.uiUpdate(null);
								}
							}).start();
					    }				
					}
				});
				butLoadRef.setEnabled(true);
				panFile.add(butLoadRef, BorderLayout.EAST);
				paraPanel.add(panFile, BorderLayout.NORTH);

				
				JPanel panAtr= new JPanel(new BorderLayout());
				panAtr.setBorder(BorderFactory.createTitledBorder(BorderFactory.createEmptyBorder(), "File Size"));
				DefaultTableModel dtm= new DefaultTableModel(getTableObjects(defaultAtrributes), 
						new Object[] {"Load", "Attribute", "Spliceforms"});
				JTable tab= new JTable(dtm);
				//for (int i = 0; i < tab.getColumnCount(); i++) {
//				tab.getColumnModel().getColumn(0).setPreferredWidth(15);
//				tab.getColumnModel().getColumn(1).setPreferredWidth(25);
//				tab.getColumnModel().getColumn(2).setPreferredWidth(35);
				//}
				JScrollPane attrScroller= new JScrollPane(tab);
				attrScroller.setPreferredSize(new Dimension(230,10));
				panAtr.add(attrScroller, BorderLayout.CENTER);
				paraPanel.add(panAtr, BorderLayout.WEST);
				
				JPanel pan= new JPanel();
				pan.setLayout(new BoxLayout(pan, BoxLayout.Y_AXIS));
				JPanel panFrame= new JPanel();
				panFrame.setBorder(BorderFactory.createTitledBorder(BorderFactory.createEmptyBorder(), "Spliceform Lengths"));
				lengthDistrPanel= new SimpleBinPlotterPanel();
				lengthDistrPanel.setTitleX("length [nt]");
				panFrame.add(lengthDistrPanel);
				pan.add(panFrame);
				pan.add(getPanelRefAttr());
				paraPanel.add(pan, BorderLayout.CENTER);

			}
			return paraPanel;
		}

	JPanel readStatPanel= null;
	JTextField tfMapAll, tfMapEntire, tfMapSplit, tfReadAll, tfReadMap, tfRfac;
	private Component getPanelReadAttr() {
			if (readStatPanel == null) {
				readStatPanel = new JPanel();
				readStatPanel.setBorder(BorderFactory.createTitledBorder(BorderFactory.createEmptyBorder(), "Attributes"));
				readStatPanel.setLayout(new GridBagLayout());
				GridBagConstraints c= new GridBagConstraints();
				c.anchor= GridBagConstraints.FIRST_LINE_START;
				c.fill = GridBagConstraints.NONE;
				c.insets= new Insets(1,10,1,10);
				c.gridwidth= 1; c.gridheight= 1;
	
				c.weighty= 0.01;
				c.gridx = 0; 
				
				c.gridy = 0;
				readStatPanel.add(new JLabel("Mappings"), c);
				c.insets= new Insets(1,20,1,10);
				c.gridy = 1;
				readStatPanel.add(new JLabel("all"), c);
				c.gridy = 2;
				readStatPanel.add(new JLabel("entire"), c);
				c.gridy = 3;
				readStatPanel.add(new JLabel("split"), c);
				c.gridy = 4;
				c.insets= new Insets(1,10,1,10);
				readStatPanel.add(new JLabel("Reads"), c);
				c.gridy = 5;
				c.insets= new Insets(1,20,1,10);
				readStatPanel.add(new JLabel("all"), c);
				c.gridy = 6;
				c.insets= new Insets(1,20,1,10);
				readStatPanel.add(new JLabel("mapped"), c);
				c.gridy = 7;
				c.insets= new Insets(1,20,1,10);
				readStatPanel.add(new JLabel("R-factor"), c);

				int tfWidth= 8;
				c.gridx = 1; 
				c.gridy = 1;
				c.insets= new Insets(1,10,1,10);
				tfMapAll= new JTextField(tfWidth);
				tfMapAll.setHorizontalAlignment(JTextField.RIGHT);
				tfMapAll.setEditable(false);
				readStatPanel.add(tfMapAll, c);
				c.gridy = 2;
				tfMapEntire= new JTextField(tfWidth);
				tfMapEntire.setHorizontalAlignment(JTextField.RIGHT);
				tfMapEntire.setEditable(false);
				readStatPanel.add(tfMapEntire, c);
				c.gridy = 3;
				tfMapSplit= new JTextField(tfWidth);
				tfMapSplit.setHorizontalAlignment(JTextField.RIGHT);
				tfMapSplit.setEditable(false);
				readStatPanel.add(tfMapSplit, c);
				
				c.gridy = 5;
				tfReadAll= new JTextField(tfWidth);
				tfReadAll.setHorizontalAlignment(JTextField.RIGHT);
				tfReadAll.setEditable(false);
				readStatPanel.add(tfReadAll, c);
				c.gridy = 6;
				tfReadMap= new JTextField(tfWidth);
				tfReadMap.setHorizontalAlignment(JTextField.RIGHT);
				tfReadMap.setEditable(false);
				readStatPanel.add(tfReadMap, c);
				c.gridy = 7;
				tfRfac= new JTextField(tfWidth);
				tfRfac.setHorizontalAlignment(JTextField.RIGHT);
				tfRfac.setEditable(false);
				readStatPanel.add(tfRfac, c);

			}
	
			return readStatPanel;
		}

}
