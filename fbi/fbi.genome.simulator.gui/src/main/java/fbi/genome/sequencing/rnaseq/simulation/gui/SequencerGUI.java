package fbi.genome.sequencing.rnaseq.simulation.gui;

import fbi.commons.ReadyOrNot;
import fbi.commons.StringUtils;
import fbi.commons.gui.MyProgressBar;
import fbi.commons.gui.SimpleBinPlotterPanel;
import fbi.commons.thread.StoppableRunnable;
import fbi.genome.sequencing.rnaseq.simulation.FluxSimulatorSettings;
import fbi.genome.sequencing.rnaseq.simulation.Sequencer;
import fbi.genome.sequencing.rnaseq.simulation.error.CrossTalkTable;
import fbi.genome.sequencing.rnaseq.simulation.error.ErrorModel;
import fbi.genome.sequencing.rnaseq.simulation.error.ModelPool;
import fbi.genome.sequencing.rnaseq.simulation.error.PositionErrorModel;

import javax.swing.*;
import javax.swing.border.BevelBorder;
import javax.swing.filechooser.FileFilter;
import javax.swing.table.DefaultTableModel;
import javax.swing.table.TableCellRenderer;
import java.awt.*;
import java.awt.event.*;
import java.io.File;
import java.util.Arrays;
import java.util.Random;
import java.util.Vector;

public class SequencerGUI extends JPanel implements ReadyOrNot, StoppableRunnable  {
	
	static final String labelTextReads= "Reads", labelTextPairs= "Read Pairs", labelTextPairsOnly= "Pairs", labelTextNr= "Number of ";
	
	static class XtableCellRenderer implements TableCellRenderer {

		static final Color[] BASE_COLS= new Color[] {
			Color.green.darker(), Color.blue.darker(), Color.yellow.darker(), Color.black, Color.red.darker()
		};

		class PercentageLabel extends JPanel {
			double[] perc;
			boolean sel;
			public PercentageLabel(double[] a, boolean selected) {
				perc= a;
				sel= selected;
				assert(a.length== BASE_COLS.length);
			}
			@Override
			protected void paintComponent(Graphics g) {
				
				super.paintComponent(g);
				
				int w= getWidth()- getInsets().left- getInsets().right;
				for (int i = perc.length- 1; i>= 0; --i) {
					if (sel)
						g.setColor(BASE_COLS[i].brighter());
					else
						g.setColor(BASE_COLS[i]);
					g.fillRect(0, 0, 
							(int) (perc[i]* w), 
							getHeight());
				}
			}
	        // The following methods override the defaults for performance reasons
	        public void validate() {}
	        public void revalidate() {}
	        protected void firePropertyChange(String propertyName, Object oldValue, Object newValue) {}
	        public void firePropertyChange(String propertyName, boolean oldValue, boolean newValue) {}

		}
		
		public Component getTableCellRendererComponent(JTable table,
				Object value, boolean isSelected, boolean hasFocus, int row,
				int column) {
			
			
			if (row== -1) {	// && column> 0
				String s= value.toString();
				JLabel lab= new JLabel(s, JLabel.CENTER);
				Font f = lab.getFont();
				lab.setFont(f.deriveFont(f.getStyle() | Font.BOLD));
				int p= Arrays.binarySearch(CrossTalkTable.SYMBOLS, s.charAt(0));
				if (p>= 0)
					lab.setBackground(BASE_COLS[p]);
				lab.setForeground(Color.white);
				lab.setOpaque(true);
				lab.setBorder(BorderFactory.createBevelBorder(BevelBorder.RAISED));
				return lab;
			}
			
			if (value instanceof double[]) {
				PercentageLabel lab= new PercentageLabel((double[]) value, isSelected);
				lab.setBorder(BorderFactory.createBevelBorder(BevelBorder.LOWERED));
				return lab;
			} else {
				
				JLabel lab= null;
				if (value!= null) {
					String s= value.toString();
					lab= new JLabel(s);
				} else
					lab= new JLabel();
				
				if (isSelected) {
					lab.setBackground(FluxSimulatorGUI.COL_BG_YEL);
					lab.setForeground(Color.black);
				} else {
					lab.setBackground(FluxSimulatorGUI.COL_BG_GREEN);
					lab.setForeground(Color.white);
				}
				lab.setOpaque(true);
				lab.setBorder(BorderFactory.createBevelBorder(BevelBorder.LOWERED));
				return lab;
			}
			
			//return new JPanel();
		}
		
	}
	
	
	JTable statsTable;
	FluxSimulatorSettings settings;
	JTextField tfIn, tfOut;
	JButton butChangeIn, butChangeOut;
	SimpleBinPlotterPanel panRankSL, panRankSM, panRankSH,
		panRankML, panRankMM, panRankMH, 
		panRankLL, panRankLM, panRankLH;

	public SequencerGUI() {
		
		setLayout(new BorderLayout());
		
		JPanel panInPan= new JPanel();
		panInPan.setBorder(BorderFactory.createTitledBorder("Input / Output")); // or on lineborder black
		panInPan.setLayout(new GridBagLayout());
		GridBagConstraints c= new GridBagConstraints();
		c.anchor= GridBagConstraints.FIRST_LINE_START;
		c.fill = GridBagConstraints.NONE;
		c.insets= new Insets(1,10,1,10);
		c.gridwidth= 1; c.gridheight= 1;			
		c.weightx= 0; c.weighty= 0;
		
		c.gridx = 0; 
		c.gridy= 0;
		panInPan.add(new JLabel("Intput"), c);
		c.gridy= 1;
		panInPan.add(new JLabel("Output"), c);
		c.gridx = 1; 
		c.gridy= 0;
		c.weightx= 1;
		c.fill = GridBagConstraints.HORIZONTAL;
		tfIn= new JTextField();
		tfIn.setEditable(false);
		panInPan.add(tfIn, c);
		c.gridy= 1;
		tfOut= new JTextField();
		tfOut.setEditable(false);
		panInPan.add(tfOut, c);
		c.weightx= 0;
		c.fill = GridBagConstraints.NONE;
		c.gridx = 2; 
		c.gridy= 0;
		butChangeIn= new JButton("Change");
		butChangeIn.setEnabled(false);
		panInPan.add(butChangeIn, c);
		c.gridy= 1;
		butChangeOut= new JButton("Change");
		butChangeOut.setEnabled(false);
		panInPan.add(butChangeOut, c);
		add(panInPan, BorderLayout.NORTH);		
		
		JPanel panSide= new JPanel();
		panSide.setLayout(new BoxLayout(panSide, BoxLayout.Y_AXIS));
		panSide.add(getPanelParameters());
		panSide.add(getPanelStatistics());
		//panSide.add(new Box.Filler(new Dimension(1,1), new Dimension(1,600), new Dimension(1,600)));
		add(panSide, BorderLayout.WEST);

		add(getPanViz(), BorderLayout.CENTER);
		
	}
	
	private void setXtalkColumnHeaders(boolean hasQual, TableCellRenderer renderer) {
		if (hasQual)
			xtalkJT.getColumnModel().getColumn(0).setHeaderValue("Qual");
		else
			xtalkJT.getColumnModel().getColumn(0).setHeaderValue("Real");

		for (int i = 0; i < xtalkJT.getColumnCount()- 1; i++) {				
			xtalkJT.getColumnModel().getColumn(i+1).setHeaderRenderer(renderer);

			if (i< CrossTalkTable.SYMBOLS.length- 2)
				xtalkJT.getColumnModel().getColumn(i+1).setHeaderValue(CrossTalkTable.SYMBOLS[i]);
			else
				xtalkJT.getColumnModel().getColumn(i+1).setHeaderValue(CrossTalkTable.SYMBOLS[i+1]);
			for (int j = 0; hasQual&& j < xtalkJT.getRowCount(); j++) 
				xtalkJT.setValueAt(Integer.toString(j+ ModelPool.qualLevels[0]), j, 0);
			if (!hasQual)
				xtalkJT.setValueAt("Called", 0, 0);
		}

	}
	
	
	JPanel panViz;
	JTable xtalkJT;
	SimpleBinPlotterPanel panPEMLoc, panPEMQual;
	JLabel labCases;
	JTable modJT;	
	TableCellRenderer renderer; 
	private JPanel getPanViz() {
		if (panViz == null) {
			panViz= new JPanel();
			panViz.setLayout(new GridLayout(2,1));
			panViz.setOpaque(true);
			panViz.setForeground(Color.cyan);
			
			JPanel panErrPos= new JPanel();
			panViz.add(panErrPos);
			panErrPos.setBorder(BorderFactory.createTitledBorder("Position error model"));
			panErrPos.setLayout(new GridLayout(1,3,10,10));
			
			JPanel panTBtot= new JPanel();
			panTBtot.setLayout(new BoxLayout(panTBtot, BoxLayout.Y_AXIS));			
			panErrPos.add(panTBtot);
			
			JLabel lab= new JLabel("Crosstalk table");
//			lab.setOpaque(false);
			panTBtot.add(lab);
			
			renderer= new XtableCellRenderer();
			xtalkJT = new JTable(new DefaultTableModel(
					1,
					CrossTalkTable.SYMBOLS.length){	// 2nd col, wo N
				@Override
				public boolean isCellEditable(int row, int column) {					
					return false;
				}
			}){
				@Override
				public TableCellRenderer getCellRenderer(int row, int column) {
//					if (column> 0)
						return renderer;
//					return super.getCellRenderer(row, column);
				}
			};

			setXtalkColumnHeaders(false, renderer);
			
			JScrollPane scroll= new JScrollPane(xtalkJT);
			scroll.getViewport().setBackground(FluxSimulatorGUI.COL_BG_YEL);
			//scroll.getViewport().setOpaque(false);
			panTBtot.add(scroll);
			
			labCases= new JLabel("Model length distribution");
//			labCases.setOpaque(false);
			panTBtot.add(labCases);		
			modJT= new JTable();
			int defLines= 7;
			if (settings!= null)
				defLines= settings.getReadLength();
			modJT.setModel(new DefaultTableModel(defLines,2){
				@Override
				public boolean isCellEditable(int row, int column) {
					return false;
				}
			});
			modJT.getColumnModel().getColumn(0).setHeaderValue("length");
			modJT.getColumnModel().getColumn(1).setHeaderValue("occurrences");
			for (int i = 0; i < 7; i++) {
				modJT.setValueAt(Integer.toString(i+1), i, 0);
				modJT.setValueAt(Integer.toString(0), i, 1);
			}
			JScrollPane scroller= new JScrollPane(modJT);
			scroller.getViewport().setBackground(FluxSimulatorGUI.COL_BG_YEL);
			panTBtot.add(scroller);
			
			panPEMLoc= new SimpleBinPlotterPanel("positional distribution");
			panPEMLoc.setPaintMode(SimpleBinPlotterPanel.MODE_LINE);
			//panPEMLoc.setPaintAxisY(false);
			panPEMLoc.setTitleX("position");
			panPEMLoc.setTitleY("models");
			panPEMLoc.setLineColors(new Color[]{
				Color.black.darker(), Color.blue.darker(), Color.red.darker() 
			});
			panErrPos.add(panPEMLoc);
			
			//panBox= new BoxPlot();
			panPEMQual= new SimpleBinPlotterPanel("quality distribution");
			panPEMQual.setPaintMode(SimpleBinPlotterPanel.MODE_LINE);
			panPEMQual.setTitleX("position");
			panPEMQual.setTitleY("quality");
			panPEMQual.setLineColors(new Color[]{
					Color.black.darker(), Color.blue.darker(), 
					Color.green.darker(), Color.orange.darker(), Color.red.darker()
			});
			panErrPos.add(panPEMQual);
		}

		return panViz;
	}
	
	JPanel panPars;
	JCheckBox cbErrPos, cbFastQ; 
	JTextField tfQThres;
	JLabel labErrFile;
	JButton butErrFile, butErrStop;
	StoppableRunnable parseRunner;
	private JLabel labNrReads, labTotReads;
	private JLabel getLabNrReads() {
		if (labNrReads == null) {
			labNrReads = new JLabel();
		}
		if (settings!= null&& settings.isPairedEnd())
			labNrReads.setText(labelTextNr+ labelTextPairsOnly);
		else
			labNrReads.setText(labelTextNr+ labelTextReads);
		
		return labNrReads;
	}
	
	private JLabel getLabTotReads() {
		if (labTotReads == null) {
			labTotReads = new JLabel();
		}
		if (settings!= null&& settings.isPairedEnd())
			labTotReads.setText(labelTextPairs);
		else
			labTotReads.setText(labelTextReads);

		return labTotReads;
	}
	
	private JTextField tfErrFname;
	private Component getPanelParameters() {
		if (panPars == null) {
			panPars = new JPanel();
			panPars.setBorder(BorderFactory.createTitledBorder("Parameters"));
			panPars.setLayout(new GridBagLayout());
			int tfWidth= 7;
			GridBagConstraints c= new GridBagConstraints();
			c.anchor= GridBagConstraints.FIRST_LINE_START;
			c.fill = GridBagConstraints.HORIZONTAL;
			c.gridwidth= 1; c.gridheight= 1;			

			c.weightx= 1;
			c.weighty= 1;

			c.gridx = 0; 
			c.gridy = 0;
			panPars.add(getLabNrReads(), c);
			c.gridy = 1; 
			panPars.add(new JLabel("Readlength"), c);
			
			c.gridy = 2;
			c.gridx = 0; 
			panPars.add(new JLabel("Output options"), c);
			c.gridx = 1; c.gridwidth= 2;
			checkPairedEnd= new JCheckBox("Enable Paired End");
			checkPairedEnd.setOpaque(false);
			checkPairedEnd.setEnabled(false);
			checkPairedEnd.addActionListener(new ActionListener(){
				public void actionPerformed(ActionEvent e) {
					settings.setPairedEnd(checkPairedEnd.isSelected());
					getLabNrReads().repaint();
					getLabTotReads().repaint();
				}
			});
			panPars.add(checkPairedEnd, c);
			
			c.gridy = 3;
			cbFastQ= new JCheckBox("Enable FastA/FastQ Output");
			cbFastQ.setOpaque(false);
			cbFastQ.addActionListener(new ActionListener(){
				public void actionPerformed(ActionEvent e) {
					if (settings!= null) {
						if (cbFastQ.isSelected()&& settings.getGenDir()== null) {
//							getChooseGenomeDirDia().setSize(dia.getPreferredSize());
//							getChooseGenomeDirDia().setLocation(200,200);
//							System.out.println(getChooseGenomeDirDia().getBounds());
//							getChooseGenomeDirDia().setVisible(true);
							cbFastQ.setSelected(false);
							JOptionPane.showMessageDialog(
									FluxSimulatorGUI.singleton, 
									"Please choose Genome folder first!\n" +
									"Menu: Project > Genome",
									"Choose a Genome folder",
									JOptionPane.WARNING_MESSAGE, 
									FluxSimulatorGUI.singleton.favicon
							);
						} else
							settings.setFastQ(cbFastQ.isSelected());
					}
				}
			});
			cbFastQ.setEnabled(false);
			panPars.add(cbFastQ, c);
			
			c.gridx = 1; 
			c.gridy = 0; //c.gridwidth= 2;
			tfNrReads= new JTextField(6);
			tfNrReads.setEditable(false);
			tfNrReads.setEnabled(false);
			tfNrReads.addActionListener(new ActionListener(){
				public void actionPerformed(ActionEvent e) {
					if (settings!= null)
						settings.setReadNr(FluxSimulatorGUI.parseLong(tfNrReads, settings.getReadNr()));
				}
			});
			tfNrReads.addFocusListener(new FocusAdapter(){
				public void focusLost(FocusEvent e) {
					if (settings!= null)
						settings.setReadNr(FluxSimulatorGUI.parseLong(tfNrReads, settings.getReadNr()));
				}
			});
			panPars.add(tfNrReads, c);
			
			c.gridy = 1; c.gridwidth= 1;
			tfReadLen= new JTextField(4);
			tfReadLen.setEditable(false);
			tfReadLen.setEnabled(false);
			tfReadLen.addActionListener(new ActionListener(){
				public void actionPerformed(ActionEvent e) {
					int x= FluxSimulatorGUI.parseInt(tfReadLen, settings.getReadLength());
					if (cbErrPos.isSelected()&& sequencer.getBabes()!= null&& sequencer.getBabes().getReadLength()!= x) {
						tfReadLen.setText(Integer.toString(sequencer.getBabes().getReadLength()));
						JOptionPane.showMessageDialog(
								FluxSimulatorGUI.singleton, 
								"The entered readlength does not correspond to the length of the error model!\n" +
								"Use another error model or deactivate the model to change the read-length.",
								"Incompatible Read Length",
								JOptionPane.WARNING_MESSAGE, 
								FluxSimulatorGUI.singleton.favicon
						);
						return;
					}
					
					if (settings!= null) 
						settings.setReadLength(x);
					
				}
			});
			tfReadLen.addFocusListener(new FocusAdapter(){
				public void focusLost(FocusEvent e) {
					int x= FluxSimulatorGUI.parseInt(tfReadLen, settings.getReadLength());
					if (cbErrPos.isSelected()&& sequencer.getBabes()!= null&& sequencer.getBabes().getReadLength()!= x) {
						tfReadLen.setText(Integer.toString(sequencer.getBabes().getReadLength()));
						JOptionPane.showMessageDialog(
								FluxSimulatorGUI.singleton, 
								"The entered readlength does not correspond to the length of the error model!\n" +
								"Use another error model or deactivate the model to change the read-length.",
								"Incompatible Read Length",
								JOptionPane.WARNING_MESSAGE, 
								FluxSimulatorGUI.singleton.favicon
						);
						return;
					}
					
					if (settings!= null)
						settings.setReadLength(x);
				}
			});
			panPars.add(tfReadLen, c);
			
			// error models
			c.gridy = 4;
			c.gridx= 0; c.gridwidth= 3;
			panPars.add(new JSeparator(JSeparator.HORIZONTAL), c);
			
			c.gridy++; 
			cbErrPos= new JCheckBox("Enable Positional Errors");
			cbErrPos.setOpaque(false);
			cbErrPos.setEnabled(false);
			cbErrPos.addActionListener(new ActionListener(){
				public void actionPerformed(ActionEvent e) {
					uiUpdatePars();
				}
			});
			panPars.add(cbErrPos, c);
			
			c.gridy++; c.gridwidth= 1;
			panPars.add(new JLabel("Quality Threshold"), c);
			c.gridx= 1; 
			c.fill= GridBagConstraints.HORIZONTAL;
			tfQThres= new JTextField(4);
			tfQThres.setEditable(false);
			panPars.add(tfQThres, c);
			tfQThres.addActionListener(new ActionListener(){
				public void actionPerformed(ActionEvent e) {
					if (settings!= null&& sequencer.getBabes()!= null&& sequencer.getBabes().hasQualities()) {						
						//settings.setQthold((float) FluxSimulatorGUI.parseDouble(tfQThres, settings.getQthold()));
					}
				}
			});
			tfQThres.addFocusListener(new FocusAdapter(){
				public void focusLost(FocusEvent e) {
					if (settings!= null&& sequencer.getBabes()!= null&& sequencer.getBabes().hasQualities())
						//settings.setQthold((float) FluxSimulatorGUI.parseDouble(tfQThres, settings.getQthold()));
				}
			});
			
			c.gridy++; c.gridx= 0;
			c.fill= GridBagConstraints.NONE;
			c.anchor= GridBagConstraints.WEST;
			labErrFile= new JLabel("Error File");
			labErrFile.setIcon(FluxSimulatorGUI.imgGrey);
			panPars.add(labErrFile, c);

			c.gridx= 1; c.gridwidth= 1;
			//c.weightx= 1; c.anchor= GridBagConstraints.EAST;
			//c.fill= GridBagConstraints.NONE;
			butErrFile= new JButton("Load");			
			panPars.add(butErrFile, c);
			butErrFile.setEnabled(false);
			tfErrFname= new JTextField(25);
			butErrFile.addActionListener(new ActionListener(){
				public void actionPerformed(ActionEvent e) {
					JFileChooser chooser= getFileChooser();
					butErrStop.setEnabled(true);
					int rval= chooser.showOpenDialog(getParent());
					if (rval== JFileChooser.CANCEL_OPTION)
						return;
					final File f= chooser.getSelectedFile(); // dia.getFile()
					if (f!= null&& f.exists()) {
						//getProgressFrame().setVisible(true);	// deadlocks when modal
						//getProgressFrame().setModal(true);
						if (f.getName().toLowerCase().endsWith(".map")) 
							parseRunner= new StoppableRunnable() {
								public void run() {
									ModelPool.parsingStopped= false;
									ModelPool babe= null;
									try {
										babe= ModelPool.readlFromGEM(f,
												-1,//settings.getQthold(),
												FluxSimulatorGUI.singleton.labStatus);
									} catch (Exception e2) {
										JOptionPane.showMessageDialog(
											FluxSimulatorGUI.singleton, 
											e2.getMessage(),
											"Incompatible alignment file",
											JOptionPane.ERROR_MESSAGE);
										return;
									}
									if (babe!= null)
										SequencerGUI.this.sequencer.setBabes(babe);
									
									parseRunner= null;
									butErrStop.setEnabled(false);
									butErrFile.setEnabled(true);
									FluxSimulatorGUI.singleton.butRun.setEnabled(true);
									
									if (ModelPool.parsingStopped|| 
											SequencerGUI.this.sequencer.getBabes()== null) {
										settings.setErrFile(null);
										return;
									}
									if (SequencerGUI.this.sequencer.getBabes().getReadLength()!= settings.getReadLength()) {
										int res= JOptionPane.showOptionDialog(
												FluxSimulatorGUI.singleton, 
												"Sequenced read length "+settings.getReadLength()+" is incompatible\n" +
														"with error model read length "+SequencerGUI.this.sequencer.getBabes().getReadLength()+"!",
												"Incompatible error model",
												JOptionPane.YES_NO_OPTION,
												JOptionPane.QUESTION_MESSAGE,
												null,
												new Object[] {"Adapt readlength", "Send error model to hell"},
												"Adapt read length");
										if (res== JOptionPane.NO_OPTION) {
											settings.setErrFile(null);
											return;
										} else {
											settings.setReadLength(babe.getReadLength());
											tfReadLen.setText(Integer.toString(babe.getReadLength()));
										}
									}
									String fName= settings.getSeqFile().getAbsolutePath();
									fName= fName.substring(0, fName.lastIndexOf(File.separator));
									fName+= File.separator+ f.getName()+ settings.DEF_SFX_ERR;
									File f1= new File(fName);
									sequencer.getBabes().write(f1);
									settings.setErrFile(f1);
									settings.save();
									plot();
									
									if (!ModelPool.parsingStopped) {
										settings.setErrFile(f);
										uiUpdatePars();
									}
									FluxSimulatorGUI.singleton.labStatus.clear();
								}
	
								public boolean isStop() {									
									return ModelPool.parsingStopped;
								}
	
								public boolean setStop() {
									ModelPool.parsingStopped= true;
									return true;
								}
								
								public boolean setStop(boolean stop) {
									if (stop)
										return setStop();
									ModelPool.parsingStopped= stop;
									return true;
								}
							};
					
						else if (f.getName().toLowerCase().endsWith(FluxSimulatorSettings.DEF_SFX_ERR))
							parseRunner= new StoppableRunnable() {
								public void run() {
									ModelPool.parsingStopped= false;
									ModelPool babe= null;
									try {
										babe= ModelPool.read(f, settings);
									} catch (Exception e2) {
										JOptionPane.showMessageDialog(
											FluxSimulatorGUI.singleton, 
											e2.getMessage(),
											"Incompatible alignment file",
											JOptionPane.ERROR_MESSAGE);
										settings.setErrFile(null);
										tfErrFname.setText(FluxSimulatorGUI.emptyString);
										return;
									}
									if (babe!= null)
										SequencerGUI.this.sequencer.setBabes(babe);

									parseRunner= null;
									butErrStop.setEnabled(false);
									butErrFile.setEnabled(true);
									FluxSimulatorGUI.singleton.butRun.setEnabled(true);
									FluxSimulatorGUI.singleton.uiUpdate("done");
									
									if (ModelPool.parsingStopped|| SequencerGUI.this.sequencer.getBabes()== null) {
										settings.setErrFile(null);
										return;
									}
									if (SequencerGUI.this.sequencer.getBabes().getReadLength()!= settings.getReadLength()) {
										int res= JOptionPane.showOptionDialog(
												FluxSimulatorGUI.singleton, 
												"Sequenced read length "+settings.getReadLength()+" is incompatible\n" +
														"with error model read length "+SequencerGUI.this.sequencer.getBabes().getReadLength()+"!",
												"Incompatible error model",
												JOptionPane.YES_NO_OPTION,
												JOptionPane.QUESTION_MESSAGE,
												null,
												new Object[] {"Adapt read length", "Send error model to hell"},
												"Adapt read length");
										if (res== JOptionPane.NO_OPTION) {
											settings.setErrFile(null);
											return;
										} else {
											settings.setReadLength(babe.getReadLength());
											tfReadLen.setText(Integer.toString(babe.getReadLength()));
										}									}
									
									File f1= new File(f.getAbsolutePath());
									settings.setErrFile(f1);
									settings.save();
									plot();
									
									if (!ModelPool.parsingStopped) {
										tfErrFname.setText(f.getName());
										tfErrFname.setEnabled(true);
										uiUpdatePars();
									}
								}
	
								public boolean isStop() {									
									return ModelPool.parsingStopped;
								}
	
								public boolean setStop(boolean stop) {
									if (stop) 
										return setStop();
									ModelPool.parsingStopped= stop;
									return true;
								}
								
								public boolean setStop() {
									ModelPool.parsingStopped= true;
									return true;
								}
							};
				

						Thread t= new Thread(parseRunner);
						butErrFile.setEnabled(false);
						butErrStop.setEnabled(true);
						FluxSimulatorGUI.singleton.butRun.setEnabled(false);
						t.start();
					} else {
						settings.setErrFile(null);
						tfErrFname.setText("");
					}
				}
			});
			c.gridx= 2;
			butErrStop= new JButton("Stop");			
			panPars.add(butErrStop, c);
			butErrStop.setEnabled(false);
			butErrStop.addActionListener(new ActionListener(){
				public void actionPerformed(ActionEvent e) {
					if (parseRunner!= null)
						parseRunner.setStop();
				}
			});

			
			c.gridy++; c.gridx= 0; c.gridwidth= 3;
			c.gridwidth= GridBagConstraints.NONE;
			c.weightx= 0.1;
			//tfErrFname.setEnabled(false);
			tfErrFname.setEditable(false);
			panPars.add(tfErrFname, c);

		}			

		return panPars;
	}
	
	JDialog prgFrame;
	MyProgressBar prgBar;
	JDialog getProgressFrame() {
		if (prgFrame == null) {
			//prgFrame = new JInternalFrame("",false,false,false,false);
			//prgFrame= new JFrame();
			//prgFrame.setUndecorated(true);
			prgFrame= new JDialog(
					FluxSimulatorGUI.singleton,// SwingUtilities.getWindowAncestor(this),
					"Progress"); 
			//prgFrame.setUndecorated(true);
			prgFrame.setResizable(false);
			prgFrame.setDefaultCloseOperation(
				    JDialog.DO_NOTHING_ON_CLOSE);


			prgBar= new MyProgressBar();
			prgBar.setMaximum(10);
			
			//prgFrame.getContentPane().setLayout(new BorderLayout());
			prgFrame.getContentPane().add(prgBar);
			prgBar.setSize(prgBar.getPreferredSize().width, prgBar.getPreferredSize().height);
			prgFrame.pack();
			Dimension dim= SwingUtilities.getWindowAncestor(this).getSize();
			prgFrame.setLocation(
					dim.width/ 2- prgFrame.getWidth()/ 2, 
					dim.height/ 2- prgFrame.getHeight()/ 2
			);
			//prgFrame.setLocation(500, 500);
			prgFrame.setSize(400, 80);
		}

		return prgFrame;

	}
	
	private JFileChooser chooser;
	private JFileChooser getFileChooser() {
		if (chooser == null) {
			chooser= new JFileChooser();	// not in the EDT !!!
			if (FluxSimulatorGUI.singleton.parFile!= null)
				chooser.setCurrentDirectory(FluxSimulatorGUI.singleton.parFile.getParentFile());
			chooser.setMultiSelectionEnabled(false);
			chooser.setPreferredSize(new Dimension(chooser.getPreferredSize().width, 600));
			chooser.setDialogTitle("Open Error File");	//dia.setTitle()
			chooser.setFileSelectionMode(JFileChooser.FILES_ONLY);
			chooser.addChoosableFileFilter(new FileFilter() {
				public boolean accept(File f) {
			        return f.isDirectory()|| f.getName().toLowerCase().endsWith(FluxSimulatorSettings.DEF_SFX_ERR)
			         || f.getName().toLowerCase().endsWith(FluxSimulatorSettings.DEF_SFX_GEM);
			    }
			    
			    public String getDescription() {
			        return "Error profile (*.err), GEM alignment (*.map)";
			    }
			});
			try {
				chooser.setCurrentDirectory(settings.getParFile().getParentFile());
			} catch (Exception ex) {
				; // :)
			}
//			chooser.addChoosableFileFilter(new FileFilter() {
//				public boolean accept(File f) {
//			        return f.isDirectory()|| f.getName().toLowerCase().endsWith(".map");
//			    }
//			    
//			    public String getDescription() {
//			        return "GEM alignment (*.map)";
//			    }
//			});
		}

		return chooser;
	}

	JPanel panStats;
	JCheckBox checkPairedEnd;
	JTextField tfNrReads, tfReadLen;
	JTextField tfLociSeq, tfLociSeqPerc, tfTrptSeq, tfTrptSeqPerc, tfReadSeq, tfReadSeqPerc;
	private Component getPanelStatistics() {
		if (panStats == null) {
			panStats = new JPanel();
			panStats.setBorder(BorderFactory.createTitledBorder("Characteristics"));
			panStats.setLayout(new GridBagLayout());
			int tfWidthNr= 5;
			GridBagConstraints c= new GridBagConstraints();
			c.anchor= GridBagConstraints.FIRST_LINE_START;
			c.fill = GridBagConstraints.NONE;
			c.weightx= 1;
			c.weighty= 1;
			c.insets= new Insets(1,10,1,10);
			c.gridwidth= 1; c.gridheight= 1;			
			
			c.gridx = 0; 
			c.gridy= 0;
			panStats.add(new JLabel("Sequenced"), c);
			c.gridy = 1;
			c.insets= new Insets(1,20,1,10);
			panStats.add(getLabTotReads(), c);
			c.gridy = 2;
			panStats.add(new JLabel("Loci"), c);
			c.gridy = 3;
			panStats.add(new JLabel("Transcripts"), c);

			c.gridy = 1;
			c.gridx= 1;
			tfReadSeq= new JTextField(tfWidthNr);
			tfReadSeq.setHorizontalAlignment(JTextField.RIGHT);
			tfReadSeq.setEnabled(false);
			tfReadSeq.setEditable(false);
			panStats.add(tfReadSeq, c);
			c.gridx= 2;
			tfReadSeqPerc= new JTextField(tfWidthNr);
			tfReadSeqPerc.setHorizontalAlignment(JTextField.RIGHT);
			tfReadSeqPerc.setEnabled(false);
			tfReadSeqPerc.setEditable(false);
			panStats.add(tfReadSeqPerc, c);
			
			c.gridy = 2;
			c.gridx= 1;
			tfLociSeq= new JTextField(tfWidthNr);
			tfLociSeq.setHorizontalAlignment(JTextField.RIGHT);
			tfLociSeq.setEditable(false);
			tfLociSeq.setEnabled(false);
			panStats.add(tfLociSeq, c);
			c.gridx= 2;
			tfLociSeqPerc= new JTextField(tfWidthNr);
			tfLociSeqPerc.setHorizontalAlignment(JTextField.RIGHT);
			tfLociSeqPerc.setEditable(false);
			tfLociSeqPerc.setEnabled(false);
			panStats.add(tfLociSeqPerc, c);
			
			c.gridy = 3;
			c.gridx= 1;
			tfTrptSeq= new JTextField(tfWidthNr);
			tfTrptSeq.setHorizontalAlignment(JTextField.RIGHT);
			tfTrptSeq.setEnabled(false);
			tfTrptSeq.setEditable(false);
			panStats.add(tfTrptSeq, c);
			c.gridx= 2;
			tfTrptSeqPerc= new JTextField(tfWidthNr);
			tfTrptSeqPerc.setHorizontalAlignment(JTextField.RIGHT);
			tfTrptSeqPerc.setEnabled(false);
			tfTrptSeqPerc.setEditable(false);
			panStats.add(tfTrptSeqPerc, c);
		}

		return panStats;
	}

	public static void main(String[] args) {
		SequencerGUI gui= new SequencerGUI();
		
		
		JFrame aFrame= new JFrame();
		aFrame.getContentPane().add(gui);
		aFrame.addWindowListener(new WindowAdapter(){
			@Override
			public void windowClosing(WindowEvent e) {
				System.exit(0);
			}
		});
		
		aFrame.setSize(new Dimension(800,600));
		aFrame.setVisible(true);
		
		int[] rnd= new int[30000];
		Random rndm= new Random();
		for (int i = 0; i < rnd.length; i++) {
			rnd[i]= rndm.nextInt(20000);
			if (rnd[i]== 0)
				--i;
		}
//		gui.panGel.paintOSI(rnd);
//		gui.panGel.repaint();
		
	}

	public boolean isReady() {
		if (sequencer== null)
			return false;
		return sequencer.isReady() == null;
	}

	Sequencer sequencer;
	public void set(Object o) {
		if (o== null|| !(o instanceof FluxSimulatorSettings))
			return;
		FluxSimulatorSettings settings= (FluxSimulatorSettings) o;
		this.settings= settings;
		if (settings== null) {
			sequencer= null;
			uiUpdateStats();
		} else {
			// save old error model, why?
//			ModelPool babe= null;
//			if (sequencer!= null)
//				babe= sequencer.getBabes();
			sequencer= new Sequencer(settings);
			sequencer.loadErrors();
			if (sequencer.getBabes()!= null)
				cbErrPos.setSelected(true);
//			sequencer.setBabes(babe);
		}
		uiUpdatePars();
	}
	
	private void uiUpdateStats() {
		if (sequencer== null) {
			tfLociSeq.setText(FluxSimulatorGUI.emptyString);
			tfLociSeq.setEnabled(false);
			tfTrptSeq.setText(FluxSimulatorGUI.emptyString);
			tfTrptSeq.setEnabled(false);
			tfReadSeq.setText(FluxSimulatorGUI.emptyString);
			tfReadSeq.setEnabled(false);
			tfLociSeqPerc.setText(FluxSimulatorGUI.emptyString);
			tfLociSeqPerc.setEnabled(false);
			tfTrptSeqPerc.setText(FluxSimulatorGUI.emptyString);
			tfTrptSeqPerc.setEnabled(false);
			tfReadSeqPerc.setText(FluxSimulatorGUI.emptyString);
			tfReadSeqPerc.setEnabled(false);
			cbFastQ.setEnabled(false);
			cbErrPos.setEnabled(false);
			
		} else {
			if (sequencer.getCntLoci()!= 0) {
				tfLociSeq.setText(Integer.toString(sequencer.getCntLoci()));
				tfLociSeq.setEnabled(true);
				tfLociSeqPerc.setText(StringUtils.fprint(sequencer.getCntLoci() * 100d / settings.getProfiler().getCntLoci(), 2)+"%");
				tfLociSeqPerc.setEnabled(true);
			}

			if (sequencer!= null&& settings.getProfiler()!= null&& settings.getProfiler().getIds()!= null&&
					sequencer.getCntTrpts()!= 0) {
				tfTrptSeq.setText(Integer.toString(sequencer.getCntTrpts()));
				tfTrptSeq.setEnabled(true);
				tfTrptSeqPerc.setText(StringUtils.fprint(sequencer.getCntTrpts() * 100d / settings.getProfiler().getIds().length, 2)+"%");
				tfTrptSeqPerc.setEnabled(true);
			}
			if (sequencer.getTotalReads()!= 0) {
				tfReadSeq.setText(Long.toString(sequencer.getTotalReads()));
				tfReadSeq.setEnabled(true);
				tfReadSeqPerc.setText(StringUtils.fprint(sequencer.getTotalReads() * 100d / (settings.getReadNr() * (settings.isPairedEnd() ? 2 : 1)), 2)+"%");
				tfReadSeqPerc.setEnabled(true);
			}
		}
		FluxSimulatorGUI.repaintEverywhere(this, false);
	}
	
	private void uiUpdatePars() {
		
		if (settings== null) {
			tfIn.setText(FluxSimulatorGUI.emptyString);
			tfOut.setText(FluxSimulatorGUI.emptyString);

			tfNrReads.setText(FluxSimulatorGUI.emptyString);
			tfNrReads.setEditable(false);
			//tfNrReads.setEnabled(false);
			tfReadLen.setText(FluxSimulatorGUI.emptyString);
			tfReadLen.setEditable(false);
			//tfReadLen.setEnabled(false);
			checkPairedEnd.setSelected(false);
			checkPairedEnd.setEnabled(false);
			
			labErrFile.setIcon(FluxSimulatorGUI.imgGrey);
			butErrFile.setEnabled(false);
			cbFastQ.setEnabled(false);
			cbErrPos.setEnabled(false);
			
		} else {
			if (settings.getFrgFile()!= null)
				tfIn.setText(settings.getFrgFile().getAbsolutePath());
			if (settings.getSeqFile()!= null)
				tfOut.setText(settings.getSeqFile().getAbsolutePath());

			tfNrReads.setText(Long.toString(settings.getReadNr()));
			tfNrReads.setEditable(true);
			tfNrReads.setEnabled(true);
			tfReadLen.setText(Integer.toString(settings.getReadLength()));
			tfReadLen.setEditable(true);
			tfReadLen.setEnabled(true);
			checkPairedEnd.setSelected(settings.isPairedEnd());
			checkPairedEnd.setEnabled(true);
			cbErrPos.setEnabled(true);

			cbFastQ.setSelected(settings.isFastQ());
			if (settings.getGenDir()!= null) {
				cbFastQ.setEnabled(true);
				//cbFastQ.setSelected(true); // only init
				
				cbErrPos.setEnabled(true);
				if (cbErrPos.isSelected())
					butErrFile.setEnabled(true);
				else
					butErrFile.setEnabled(false);
			} else {
				cbFastQ.setEnabled(false);
				cbFastQ.setSelected(false);
				cbErrPos.setEnabled(false);
				cbErrPos.setSelected(false);
				butErrFile.setEnabled(false);
			}
			
			if (settings.getErrFile()!= null)	// these are 2 things
				tfErrFname.setText(settings.getErrFile().getName());
			else
				tfErrFname.setText(FluxSimulatorGUI.emptyString);
			tfErrFname.revalidate();
			tfErrFname.repaint();
			if (sequencer.getBabes()!= null) {
				plot();
				if (settings.getGenDir()!= null) {
					//cbErrPos.setSelected(true);	// only init
					cbErrPos.setEnabled(true);
				}
				labErrFile.setIcon(FluxSimulatorGUI.imgGreen);
				if (sequencer.getBabes().hasQualities())
					tfQThres.setEditable(true);
				labErrFile.revalidate();
				labErrFile.repaint();
			} else {
				if (settings.getErrFile()== null)
					labErrFile.setIcon(FluxSimulatorGUI.imgGrey);
				else 
					labErrFile.setIcon(FluxSimulatorGUI.imgRed);
			}

		}
		FluxSimulatorGUI.repaintEverywhere(this, false);
	}

	public boolean setStop() {
		if (sequencer== null)
			return false;		
		return sequencer.setStop(true);
	}

	private void loadTextFields() {
		settings.setReadNr(FluxSimulatorGUI.parseLong(tfNrReads,settings.getReadNr()));
		settings.setReadLength(FluxSimulatorGUI.parseInt(tfReadLen,settings.getReadLength()));
	}
	
	BoxPlot panBox;
	public void plot() {
		
		if (sequencer.getBabes()== null)
			return;
		
		// x-talk
		// TODO colors
		
		boolean qOn= sequencer!= null&& sequencer.getBabes()!= null&& sequencer.getBabes().hasQualities();
		int rowCount= qOn?
				sequencer.getBabes().qualLevels[1]- sequencer.getBabes().qualLevels[0]+ 1: 1;
		xtalkJT.setModel(new DefaultTableModel(
				rowCount, 
				CrossTalkTable.SYMBOLS.length){
			@Override
			public boolean isCellEditable(int row, int column) {
				return false;
			}
		});
		setXtalkColumnHeaders(qOn, renderer);
		if (qOn) {
			double[][][] a= (double[][][]) sequencer.getBabes().getXTalkTable().getTable();
			for (int i = 0; i < (a.length- 1); i++) { 
				int sub= 0;	// col change
				if (CrossTalkTable.SYMBOLS[i]== 'N') {
					sub= 1;
					//continue;
				}
				for (int j = 0; qOn&& j < a[i].length; j++)					
					xtalkJT.setValueAt(a[i+sub][j], j, i+1);	// i+1-sub
			}
		} else {
			double[][] a= (double[][]) sequencer.getBabes().getXTalkTable().getTable();
			for (int i = 0; i < (a.length- 1); i++) { 
				int sub= 0;	// col change
				if (CrossTalkTable.SYMBOLS[i]== 'N') {
					sub= 1;
					//continue;
				}
				for (int j = 0; j < a[i].length; j++) 
				xtalkJT.setValueAt(a[i+sub], 0, i+1);	// i+1-sub
			}
		}

		// lengths
		int[] len= sequencer.getBabes().getLengthDistr();
//		if (len.length!= modJT.getRowCount()- 1) {
//			for (int i = 0; i < Math.abs(len.length- (modJT.getRowCount()-1)); i++) {
//				if (len.length> (modJT.getRowCount()- 1))
//					((DefaultTableModel) modJT.getModel()).addRow(new String[] {Integer.toString(modJT.getRowCount()+1), ""});
//				else 
//					((DefaultTableModel) modJT.getModel()).removeRow(modJT.getRowCount()- 1);
//			}
//		}
		modJT.setModel(new DefaultTableModel(len.length,2){
			@Override
			public boolean isCellEditable(int row, int column) {
				return false;
			}
		});
		modJT.getColumnModel().getColumn(0).setHeaderValue("length");
		modJT.getColumnModel().getColumn(1).setHeaderValue("occurrences");
		for (int i = 0; i < modJT.getRowCount(); i++) {
			modJT.setValueAt(Integer.toString(i+1), i, 0);
			modJT.setValueAt(Integer.toString(len[i]), i, 1);
		}
//		for (int i = 0; i < (len.length- modJT.getRowCount()); i++) 
//			((DefaultTableModel) modJT.getModel()).addRow(
//					new Object[]{
//							Integer.toString(modJT.getRowCount()+i),
//							Integer.toString(len[modJT.getRowCount()-1+i])});
		
		// starts, ends, overlaps
		long[] ovl= sequencer.getBabes().getPEMoverlap(settings.getReadLength());
		panPEMLoc.paintOSI(ovl); 
		panPEMLoc.paintOSI(sequencer.getBabes().getPEMends(settings.getReadLength()));
		panPEMLoc.paintOSI(sequencer.getBabes().getPEMstarts(settings.getReadLength()));
		

		if (sequencer.getBabes().hasQualities()) {
			Vector<ErrorModel> v= sequencer.getBabes().getModelV();
			System.err.println();
			long[][] allQ= new long[ovl.length][];
			for (int i = 0; i < allQ.length; i++) 
				allQ[i]= new long[((PositionErrorModel) v.elementAt(0)).getQuals()[0].length];	// NOT .getExtension()];
			for (int i = 0; i < v.size(); i++) {
				if (!(v.elementAt(i) instanceof PositionErrorModel))
					continue;
				double[][] q= ((PositionErrorModel) v.elementAt(i)).getQuals();
				double f= ((PositionErrorModel) v.elementAt(i)).getBaseProbability()* 
					sequencer.getBabes().getTotalCases();
				int offs= ((PositionErrorModel) v.elementAt(i)).getStart();
				for (int j = 0; q!= null&& j < q.length; j++) 
					for (int h = 0; h < q[j].length; h++) 
						allQ[offs+ j][h]+= f* (q[j][h]-((h>0)?q[j][h-1]:0));
			} 
			int[][] vals= new int[5][];
			for (int i = 0; i < vals.length; i++) {
				vals[i]= new int[allQ.length];
				for (int j = 0; j < vals[i].length; j++) 
					vals[i][j]= Integer.MIN_VALUE;
			}
			for (int i = 0; i < allQ.length; i++) { // pos
				double sum;
				int j;
				for (j = 0, sum= 0; j < allQ[i].length; sum+= allQ[i][j++]); 
				double med= sum/ 2d;
				double first= sum/ 4d;
				double third= 3* first;			
				for (j = 0, sum= 0; j < allQ[i].length; j++) {				
					if (allQ[i][j]== 0)
						continue;
					sum+= allQ[i][j];
					int val= ModelPool.qualLevels[0]+ j;
					if (vals[0][i]== Integer.MIN_VALUE) 
						vals[0][i]= val;
					if (vals[1][i]== Integer.MIN_VALUE&& sum>= first)
						vals[1][i]= val;
					if (vals[2][i]== Integer.MIN_VALUE&& sum>= med)
						vals[2][i]= val;
					if (vals[3][i]== Integer.MIN_VALUE&& sum>= third)
						vals[3][i]= val;
					vals[4][i]= val;
				}
	//			double[] pts= new double[5];
	//			for (int j = 0; j < pts.length; j++) 
	//				panBox.add(pts);
			}
			for (int i = vals.length-1; i >=0; --i) 
				panPEMQual.paintOSI(vals[i]);
		
		} else {
			panPEMQual.paintOSI(null);
		}
			
		
		SwingUtilities.invokeLater(
			new Runnable() {
				public void run() {
					getPanViz().repaint();
				}
		});
		
	}
	
	public void run() {
		
		if (loadStats) {
			loadStats();
			setLoadStats(false);
			return;
		} else {
			loadTextFields();
			sequencer.setStop(false);
			try {
				sequencer.run();
			} catch (Exception e) {
				FluxSimulatorGUI.dialogError(e);
			}
		}
		
		if (sequencer.isStop()) {
			if (settings.getSeqFile().exists())
				settings.getSeqFile().delete();
		} else {
			settings.save();
			uiUpdateStats();
		}
	}

	public boolean loadStats() {
		if (sequencer== null)
			return false;
		if (!sequencer.loadStats())
			return false;
		uiUpdatePars();
		uiUpdateStats();
		if (sequencer.loadErrors())	// sequencer.getBabes()!= null
			plot();
		return true;
	} 
	
	public boolean isFinished() {
		if (sequencer== null)
			return false;
		return sequencer.isFinished();
	}

	public boolean getLoadStats() {
		return loadStats;
	}

	boolean loadStats= false;
	public void setLoadStats(boolean val) {
		loadStats= val;
	}

	public void killResult() {
		if (settings!= null&& settings.getSeqFile()!= null&& settings.getSeqFile().exists())
			settings.getSeqFile().delete();
//		ModelPool babes= sequencer.getBabes();
		if (sequencer!= null) {
			sequencer= null;
			uiUpdateStats();
		}
		sequencer= new Sequencer(settings);
		sequencer.loadErrors();
//		sequencer.setBabes(babes);
	}

	public boolean isStop() {		
		return sequencer== null|| sequencer.isStop();
	}

	public boolean setStop(boolean stop) {
		if (stop)
			return setStop();
		else 
			return sequencer.setStop(stop);
	}
}
