package fbi.genome.sequencing.rnaseq.simulation.gui;

import fbi.commons.ByteArrayCharSequence;
import fbi.commons.ReadyOrNot;
import fbi.commons.gui.SimpleBinPlotterPanel;
import fbi.commons.thread.StoppableRunnable;
import fbi.genome.sequencing.rnaseq.simulation.FluxSimulatorSettings;
import fbi.genome.sequencing.rnaseq.simulation.Fragmenter;
import fbi.genome.sequencing.rnaseq.simulation.ProgressablePlotter;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
import java.util.Hashtable;
import java.util.Random;

public class FragmenterGUI extends JPanel implements ReadyOrNot, StoppableRunnable  {
	public static final int FRAG_LEN_INTEREST_MAX= 10000;
	public class Plotter9 implements ProgressablePlotter {

		// bad idea, hashCode incompatibilities String<>BACharSequence JDK dependent
		Hashtable<ByteArrayCharSequence,int[]> hash;
		public Plotter9() {
			hash= new Hashtable<ByteArrayCharSequence, int[]>();
		}
		
		public void addBase(ByteArrayCharSequence ID, int length, int mol) {
			hash.put(ID, new int[]{length, mol});
		}

		long moInit= -1;
		public void setMolTot(final long value) {
			if (moInit< 0)
				moInit= value;
			try {
				SwingUtilities.invokeLater(new Runnable(){
					public void run() {
						if (FragmenterGUI.this.tabPlots.getTabCount()== 0) {
							tfMolRNA.setText(Long.toString(value));
							tfMolRNA.paintImmediately(tfMolRNA.getBounds());
						} else {
							String s= FragmenterGUI.this.tabPlots.getTitleAt(FragmenterGUI.this.tabPlots.getTabCount()-1);
							if (s.startsWith(Fragmenter.MODE_RT_MESSAGE)) {
								tfMolRT.setText(Long.toString(value));
								tfMolRT.paintImmediately(tfMolRT.getBounds());
							} else if (s.startsWith(Fragmenter.MODE_NEBU_MESSAGE)|| s.startsWith(Fragmenter.MODE_FRAG_MESSAGE)) {
								tfMolFrag.setText(Long.toString(value));
								tfMolFrag.paintImmediately(tfMolFrag.getBounds());
							} else if (s.startsWith(Fragmenter.MODE_FILT_MESSAGE)) {
								tfMolFilt.setText(Long.toString(value));
								tfMolFilt.paintImmediately(tfMolFilt.getBounds());
							}
						}
					}
				});
			} catch (Exception e) {
				System.out.println(e);; // :)
			}
		}

		public void reset(final String message) {
			try {
				SwingUtilities.invokeAndWait(new Runnable(){
					public void run() {
						FragmenterGUI.this.addTab(message);
						FragmenterGUI.this.getTabPlots().setSelectedComponent(
								FragmenterGUI.this.getTabPlots().getComponent(
										FragmenterGUI.this.getTabPlots().getTabCount()-1));
						FragmenterGUI.this.getTabPlots().paintImmediately(
								FragmenterGUI.this.getTabPlots().getBounds());	// need size for resetOSI()
						FragmenterGUI.this.getTabPlots().revalidate();
						
						FragmenterGUI.this.panRankSL.resetOSI(1);
						//FragmenterGUI.this.panRankSM.resetOSI();
						FragmenterGUI.this.panRankSH.resetOSI(1);
						FragmenterGUI.this.panRankML.resetOSI(1);
						//FragmenterGUI.this.panRankMM.resetOSI();
						FragmenterGUI.this.panRankMH.resetOSI(1);
						FragmenterGUI.this.panRankLL.resetOSI(1);
						//FragmenterGUI.this.panRankLM.resetOSI();
						FragmenterGUI.this.panRankLH.resetOSI(1);
						FragmenterGUI.this.panGel.resetOSI(10000);
						
						FragmenterGUI.this.panFragLog.resetOSI(Math.log10(FRAG_LEN_INTEREST_MAX));
						FragmenterGUI.this.panFragLen.resetOSI(FRAG_LEN_INTEREST_MAX);
						if (settings!= null)
							FragmenterGUI.this.panFragSel.resetOSI(fragmenter.getDeltaFiltMax()- fragmenter.getDeltaFiltMin()+ 1);
						else
							FragmenterGUI.this.panFragSel.resetOSI(1);
					}
				});
			} catch (Exception e) {
				e.printStackTrace();	// :)
			}
			
		}
		
		public void paint() {
			try {
				SwingUtilities.invokeAndWait(new Runnable(){
					public void run() {
						FragmenterGUI.this.panRankSL.paintOSI(null);
						//FragmenterGUI.this.panRankSM.paintOSI(null);
						FragmenterGUI.this.panRankSH.paintOSI(null);
						FragmenterGUI.this.panRankML.paintOSI(null);
						//FragmenterGUI.this.panRankMM.paintOSI(null);
						FragmenterGUI.this.panRankMH.paintOSI(null);
						FragmenterGUI.this.panRankLL.paintOSI(null);
						//FragmenterGUI.this.panRankLM.paintOSI(null);
						FragmenterGUI.this.panRankLH.paintOSI(null);
						FragmenterGUI.this.panGel.paintOSI(null);
						
						FragmenterGUI.this.panFragLog.paintOSI(null);
						FragmenterGUI.this.panFragLen.paintOSI(null);
						FragmenterGUI.this.panFragSel.paintOSI(null);
						
						
						FragmenterGUI.this.panRankSL.repaint();
						//FragmenterGUI.this.panRankSM.repaint();
						FragmenterGUI.this.panRankSH.repaint();
						FragmenterGUI.this.panRankML.repaint();
						//FragmenterGUI.this.panRankMM.repaint();
						FragmenterGUI.this.panRankMH.repaint();
						FragmenterGUI.this.panRankLL.repaint();
						//FragmenterGUI.this.panRankLM.repaint();
						FragmenterGUI.this.panRankLH.repaint();
						FragmenterGUI.this.panGel.repaint();
						
						FragmenterGUI.this.panFragLog.repaint();
						FragmenterGUI.this.panFragLen.repaint();
						FragmenterGUI.this.panFragSel.repaint();

					}
				});
			} catch (Exception e) {
				; // :)
			}
		}

		private int[] tmp;
		public boolean plot(int start, int end, int segLen, ByteArrayCharSequence ID) {
			//String s= ID.toString();
			tmp= hash.get(ID);
			double molPerCell= tmp[1]/ ((moInit/ FluxSimulatorSettings.AVG_MOL_CELL)+ 1);	// +1 avoids /0
			double fracStart= ((double) start)/ tmp[0],
				fracEnd= ((double) end)/ tmp[0];
			if (tmp[0]< FluxSimulatorSettings.SF_SHORT) {
				if (molPerCell>= FluxSimulatorSettings.SF_CELL_HI) {
					if (fracStart>= 0&& fracStart<= 1)
						FragmenterGUI.this.panRankSH.addVal(fracStart);
					if (fracEnd>= 0&& fracEnd<= 1)
						FragmenterGUI.this.panRankSH.addVal(fracEnd);
/*				} else if (molPerCell>= FluxSimulatorSettings.SF_CELL_MED) {
					if (fracStart>= 0&& fracStart<= 1)
						FragmenterGUI.this.panRankSM.addVal(fracStart);
					if (fracEnd>= 0&& fracEnd<= 1)
						FragmenterGUI.this.panRankSM.addVal(fracEnd);
*/				} else {
					if (fracStart>= 0&& fracStart<= 1)
						FragmenterGUI.this.panRankSL.addVal(fracStart);
					if (fracEnd>= 0&& fracEnd<= 1)
						FragmenterGUI.this.panRankSL.addVal(fracEnd);
				}
			} else if (tmp[0]< FluxSimulatorSettings.SF_MEDIUM) {
				if (molPerCell>= FluxSimulatorSettings.SF_CELL_HI) {
					if (fracStart>= 0&& fracStart<= 1)
						FragmenterGUI.this.panRankMH.addVal(fracStart);
					if (fracEnd>= 0&& fracEnd<= 1)
						FragmenterGUI.this.panRankMH.addVal(fracEnd);
/*				} else if (molPerCell>= FluxSimulatorSettings.SF_CELL_MED) {
					if (fracStart>= 0&& fracStart<= 1)
						FragmenterGUI.this.panRankMM.addVal(fracStart);
					if (fracEnd>= 0&& fracEnd<= 1)
						FragmenterGUI.this.panRankMM.addVal(fracEnd);
*/				} else {
					if (fracStart>= 0&& fracStart<= 1)
						FragmenterGUI.this.panRankML.addVal(fracStart);
					if (fracEnd>= 0&& fracEnd<= 1)
						FragmenterGUI.this.panRankML.addVal(fracEnd);
				}
			} else {
				if (molPerCell>= FluxSimulatorSettings.SF_CELL_HI) {
					if (fracStart>= 0&& fracStart<= 1)
						FragmenterGUI.this.panRankLH.addVal(fracStart);
					if (fracEnd>= 0&& fracEnd<= 1)
						FragmenterGUI.this.panRankLH.addVal(fracEnd);
/*				} else if (molPerCell>= FluxSimulatorSettings.SF_CELL_MED) {
					if (fracStart>= 0&& fracStart<= 1)
						FragmenterGUI.this.panRankLM.addVal(fracStart);
					if (fracEnd>= 0&& fracEnd<= 1)
						FragmenterGUI.this.panRankLM.addVal(fracEnd);
*/				} else {
					if (fracStart>= 0&& fracStart<= 1)
						FragmenterGUI.this.panRankLL.addVal(fracStart);
					if (fracEnd>= 0&& fracEnd<= 1)
						FragmenterGUI.this.panRankLL.addVal(fracEnd);
				}
			}
			
			int len= end- start+ 1;
			FragmenterGUI.this.panGel.addVal(len);
			if (segLen< 0)
				segLen= Fragmenter.getSegregatedLength(len, true); 
			if (len<= FRAG_LEN_INTEREST_MAX) {
				FragmenterGUI.this.panFragLen.addVal(len);
				FragmenterGUI.this.panFragLog.addVal(Math.log10(segLen));
			}
			if (settings!= null&& settings.isFilter() 
					&& segLen>= fragmenter.getDeltaFiltMin()
					&& segLen<= fragmenter.getDeltaFiltMax()) { 				
				FragmenterGUI.this.panFragSel.addVal(len- fragmenter.getDeltaFiltMin());
				return true;
			}
			// len= SimpleGelBinPlotterPanel.getEffectiveLength(len);	// f* slow
			return false;
		}
		
	}
	
	boolean loadStats= false;
	JTable statsTable;
	FluxSimulatorSettings settings;
	JTextField tfIn, tfOut;
	JButton butChangeIn, butChangeOut;
	SimpleGelBinPlotterPanel panGel;
	SimpleBinPlotterPanel panRankSL, panFragLog, panRankSH,
		panRankML, panFragLen, panRankMH, 
		panRankLL, panFragSel, panRankLH;
	
	public FragmenterGUI() {
		
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
		
		JPanel panSide= new JPanel();
		panSide.setLayout(new BoxLayout(panSide, BoxLayout.Y_AXIS));
		panSide.add(getPanelParameters());
		panSide.add(getPanelStatistics());
		//panSide.add(new Box.Filler(new Dimension(1,1), new Dimension(1,600), new Dimension(1,600)));
		add(panSide, BorderLayout.WEST);

		JPanel panViz= new JPanel();
		panViz.setLayout(new BorderLayout());
		//panViz.setBorder(BorderFactory.createTitledBorder("Fragment Distribution"));
//		panGel= new SimpleBinPlotterPanel();
//		panGel.paintMode= SimpleBinPlotterPanel.MODE_GEL;
//		panGel.setInvert(true);
//		panGel.setBorder(BorderFactory.createTitledBorder("N/S blot"));
//		panViz.add(panGel, BorderLayout.EAST);
		panViz.add(getTabPlots(), BorderLayout.CENTER);
		JPanel panHelp= new JPanel(new BorderLayout());
		add(panHelp, BorderLayout.CENTER);
		panHelp.add(panInPan, BorderLayout.NORTH);		
		panHelp.add(panViz, BorderLayout.CENTER);
		
	}
	
	JTabbedPane tabPlots;
	public JTabbedPane getTabPlots() {
		if (tabPlots == null) {
			tabPlots = new JTabbedPane(JTabbedPane.BOTTOM);
			// didnt solve repaint pb when tab opens
//			tabPlots.addChangeListener( new ChangeListener(){
//				public void stateChanged(ChangeEvent e) {
//					tabPlots.repaint();
//				}
//			});

		}

		return tabPlots;
	}
	
	public void addTab(String message) {
		
		if (getTabPlots().getTabCount()> 0)
			try {	// eye
				Thread.sleep(1000);
			} catch (Exception e) {
				; // :)
			}
		
		/*{
			@Override
			protected void paintComponent(Graphics g) {
				super.paintComponent(g);
				int w= getWidth()- getInsets().left- getInsets().right, 
					h= getHeight()- getInsets().top- getInsets().bottom;
				Dimension d= new Dimension((w- 2*10)/3,(h- 2*10)/3);
				panRankSL.setPreferredSize(d);
				panRankSM.setPreferredSize(d);
				panRankSH.setPreferredSize(d);
				panRankML.setPreferredSize(d);
				panRankMM.setPreferredSize(d);
				panRankMH.setPreferredSize(d);
				panRankLL.setPreferredSize(d);
				panRankLM.setPreferredSize(d);
				panRankLH.setPreferredSize(d);
				panGel.setPreferredSize(new Dimension(panGel.getPreferredSize().width, h));
			}
		};*/
		JPanel panPlots= new JPanel(); 
		panPlots.setBorder(BorderFactory.createTitledBorder("Breakpoint Distribution"));
		panPlots.setLayout(new GridLayout(2,3,10,10));
		//Dimension d= new Dimension(203, 156);
		panRankSL= new SimpleBinPlotterPanel("short-low");
		//panRankSL.setPreferredSize(d);
		panPlots.add(panRankSL);
		panRankML= new SimpleBinPlotterPanel("medium-low");
		//panRankML.setPreferredSize(d);
		panPlots.add(panRankML);
		panRankLL= new SimpleBinPlotterPanel("long-low");
		//panRankLL.setPreferredSize(d);
		panPlots.add(panRankLL);
		panRankSH= new SimpleBinPlotterPanel("short-high");
		//panRankSH.setPreferredSize(d);
		panPlots.add(panRankSH);
		panRankMH= new SimpleBinPlotterPanel("medium-high");
		//panRankMH.setPreferredSize(d);
		panPlots.add(panRankMH);
		panRankLH= new SimpleBinPlotterPanel("long-high");
		//panRankLH.setPreferredSize(d);
		panPlots.add(panRankLH);
		
		JPanel panPlots2= new JPanel(); 
		panPlots2.setBorder(BorderFactory.createTitledBorder("Fragment Lengths"));
		panPlots2.setLayout(new GridLayout(1,3,10,10));
		panFragLen= new SimpleBinPlotterPanel("real");
		//panRankMM.setPreferredSize(d);
		panPlots2.add(panFragLen);
		panFragLog= new SimpleBinPlotterPanel("gel");
		panFragLog.setLog((byte) 2, SimpleBinPlotterPanel.LOG2);
		//panRankSM.setPreferredSize(d);
		panPlots2.add(panFragLog);
		panFragSel= new SimpleBinPlotterPanel("selected");
		//panRankLM.setPreferredSize(d);
		if (settings!= null)
			panFragSel.setMinXoffset(fragmenter.getDeltaFiltMin());	
		panPlots2.add(panFragSel);

		JSplitPane vpane = new JSplitPane(JSplitPane.VERTICAL_SPLIT, panPlots, panPlots2);
		//vpane.setDividerLocation(0.666d);	// no fx, not painted
		vpane.setDividerLocation(350);
		vpane.setDividerSize(0);
		
		JPanel panPlotsNgel= new JPanel();
		panPlotsNgel.setLayout(new BorderLayout());
		panGel= new SimpleGelBinPlotterPanel();
		//panGel.setPreferredSize(new Dimension(80,489));
		panGel.setPaintMode(SimpleBinPlotterPanel.MODE_GEL);
		panGel.setInvert(true);
		JPanel gelWrapper= new JPanel();
		gelWrapper.setLayout(new BorderLayout());
		gelWrapper.setBorder(BorderFactory.createEmptyBorder(0,10,0,0));
		gelWrapper.add(panGel, BorderLayout.CENTER);
		panPlotsNgel.add(vpane, BorderLayout.CENTER);
		panPlotsNgel.add(gelWrapper, BorderLayout.EAST);
		
		//getTabPlots().setEnabled(false);
		getTabPlots().addTab(message, panPlotsNgel);
		// didnt solve the repaint pb
//		if (getTabPlots().getTabCount()>1&& getTabPlots().getSelectedIndex()== getTabPlots().getTabCount()-1)
//			getTabPlots().setSelectedIndex(getTabPlots().getTabCount()- 2);
//		getTabPlots().invalidate();
//		getTabPlots().revalidate();
//		getTabPlots().paintImmediately(getTabPlots().getBounds());
//		getTabPlots().repaint();
	}
	
	JPanel panPars;
	JTextField tfTSSmean, tfPolyAshape, tfPolyAscale;
	JTextField tfRTmax, tfRTmin, tfNebLamda, tfNebSigma, tfNebThres;
	JCheckBox checkNebFirst, cbFragmentation, cbFiltering;
	JRadioButton butFragPhys, butFragChem, butRTrnd, butRTpdT;
	static double parseDouble(JTextField tf, double defVal) {
		try {
			return Double.parseDouble(tf.getText());
		} catch (NumberFormatException xxx) {
			tf.setText(Double.toString(defVal));
		}
		return defVal;
	}
	
	static int parseInt(JTextField tf, int defVal) {
		try {
			return Integer.parseInt(tf.getText());
		} catch (NumberFormatException xxx) {
			tf.setText(Integer.toString(defVal));
		}
		return defVal;
	}
	
	private Component getPanelParameters() {
		if (panPars == null) {
			panPars = new JPanel();
			panPars.setBorder(BorderFactory.createTitledBorder("Parameters"));
			panPars.setLayout(new GridBagLayout());
			int tfWidth= 7;
			GridBagConstraints c= new GridBagConstraints();
			c.anchor= GridBagConstraints.FIRST_LINE_START;
			c.fill = GridBagConstraints.NONE;
			Insets ins= c.insets;
			c.gridwidth= 1; c.gridheight= 1;			

			c.weightx= 1;
			c.weighty= 1;
			c.gridwidth = GridBagConstraints.RELATIVE;

			c.gridx = 0; c.gridy = 0;
			c.gridwidth= 2;
			c.insets= new Insets(1,10,1,10);
			panPars.add(new JLabel("Transcript Variation"), c);
			c.gridy++;
			c.gridwidth= 1;
			c.insets= new Insets(1,30,1,10);
			panPars.add(new JLabel("TSS mean"), c);
			c.gridx= 1;
			c.insets= new Insets(1,10,1,10);
			tfTSSmean= new JTextField(tfWidth);
			tfTSSmean.setHorizontalAlignment(JTextField.RIGHT);
			tfTSSmean.setEditable(false);
			tfTSSmean.addActionListener(new ActionListener(){
				public void actionPerformed(ActionEvent e) {
					if (settings!= null) {
						if (tfTSSmean.getText().equals(""))
							settings.setTssVar(false);
						else {
							settings.setTssMean(parseDouble(tfTSSmean, settings.getTssMean()));
							if (Double.isNaN(settings.getTssMean()))
								tfTSSmean.setText("");
						}
					}
				}
			});
			tfTSSmean.addFocusListener(new FocusAdapter(){
				public void focusLost(FocusEvent e) {
					if (settings!= null) {
						if (tfTSSmean.getText().equals(""))
							settings.setTssVar(false);
						else {
							settings.setTssMean(parseDouble(tfTSSmean, settings.getTssMean()));
							if (Double.isNaN(settings.getTssMean()))
								tfTSSmean.setText("");
						}
					}
				}
			});
			panPars.add(tfTSSmean, c);

			c.gridy++; c.gridx= 0;
			c.insets= new Insets(1,30,1,10);
			panPars.add(new JLabel("polyA shape"), c);
			c.gridx= 1;
			c.insets= new Insets(1,10,1,10);
			tfPolyAshape= new JTextField(tfWidth);
			tfPolyAshape.setHorizontalAlignment(JTextField.RIGHT);
			tfPolyAshape.setEditable(false);
			tfPolyAshape.addActionListener(new ActionListener(){
				public void actionPerformed(ActionEvent e) {
					if (settings!= null) {
						if (tfPolyAshape.getText().trim().equals("")) {
							settings.setPolyAVar(false);
							tfPolyAscale.setText("");
						} else {
							settings.setPolyAshape(parseDouble(tfPolyAshape, settings.getPolyAshape()));
							if (Double.isNaN(settings.getPolyAshape()))
								tfPolyAshape.setText("");
						}
					}
				}
			});
			tfPolyAshape.addFocusListener(new FocusAdapter(){
				public void focusLost(FocusEvent e) {
					if (settings!= null) {
						if (tfPolyAshape.getText().trim().equals("")) {
							settings.setPolyAVar(false);
						} else {
							settings.setPolyAshape(parseDouble(tfPolyAshape, settings.getPolyAshape()));
							if (Double.isNaN(settings.getPolyAshape()))
								tfPolyAshape.setText("");
						}
					}
				}
			});
			panPars.add(tfPolyAshape, c);
			
			c.gridy++; c.gridx= 0;
			c.insets= new Insets(1,30,1,10);
			panPars.add(new JLabel("polyA scale"), c);
			c.gridx= 1;
			c.insets= new Insets(1,10,1,10);
			tfPolyAscale= new JTextField(tfWidth);
			tfPolyAscale.setHorizontalAlignment(JTextField.RIGHT);
			tfPolyAscale.setEditable(false);
			tfPolyAscale.addActionListener(new ActionListener(){
				public void actionPerformed(ActionEvent e) {
					if (settings!= null) {
						if (tfPolyAscale.getText().trim().equals("")) 
							settings.setPolyAVar(false);
						else {
							settings.setPolyAscale(parseDouble(tfPolyAscale, settings.getPolyAscale()));
							if (Double.isNaN(settings.getPolyAscale()))
								tfPolyAscale.setText("");
						}
					}
				}
			});
			tfPolyAscale.addFocusListener(new FocusAdapter(){
				public void focusLost(FocusEvent e) {
					if (settings!= null) {
						if (tfPolyAscale.getText().trim().equals("")) 
							settings.setPolyAVar(false);
						else {
							settings.setPolyAscale(parseDouble(tfPolyAscale, settings.getPolyAscale()));
							if (Double.isNaN(settings.getPolyAscale()))
								tfPolyAscale.setText("");
						}
					}
				}
			});
			panPars.add(tfPolyAscale, c);
			
			c.gridy++; c.gridx= 0;
			c.gridwidth= 2;
			c.fill= GridBagConstraints.HORIZONTAL;
			panPars.add(new JSeparator(JSeparator.HORIZONTAL), c);
			
			c.gridy++; c.gridx= 0;
			c.gridwidth= 2;
			c.insets= new Insets(1,10,1,10);
			ButtonGroup groupRT= new ButtonGroup(), groupFrag= new ButtonGroup();
			panPars.add(new JLabel("Reverse Transcription"), c);
			c.gridy++ ;
			c.gridwidth= 1;
			//panPars.add(new JLabel("Falloff Length"), c);
			
			c.gridy++; c.gridx= 0; 
			c.gridwidth= 1;
			c.insets= new Insets(1,30,1,10);
			panPars.add(new JLabel("RT Min Length"), c);			
			c.gridx = 1; 
			//c.fill= GridBagConstraints.HORIZONTAL;
			tfRTmin= new JTextField(tfWidth);
			tfRTmin.setHorizontalAlignment(JTextField.RIGHT);
			tfRTmin.setEditable(false);
			tfRTmin.addActionListener(new ActionListener(){
				public void actionPerformed(ActionEvent e) {
					if (settings!= null)
						settings.setRTminLen(parseInt(tfRTmin, settings.getRTminLen()));
				} 
			});
			tfRTmin.addFocusListener(new FocusAdapter(){
				public void focusLost(FocusEvent e) {
					if (settings!= null)
						settings.setRTminLen(parseInt(tfRTmin, settings.getRTminLen()));
				}
			});
			panPars.add(tfRTmin, c);

			c.gridy++; c.gridx= 0;
			panPars.add(new JLabel("RT Max Length"), c);			
			c.gridx= 1;
			tfRTmax= new JTextField(tfWidth);
			tfRTmax.setHorizontalAlignment(JTextField.RIGHT);
			tfRTmax.setEditable(false);
			tfRTmax.addActionListener(new ActionListener(){
				public void actionPerformed(ActionEvent e) {
					if (settings!= null)
						settings.setRTmaxLen(parseInt(tfRTmax, settings.getRTmaxLen()));
				}
			});
			tfRTmax.addFocusListener(new FocusAdapter(){
				public void focusLost(FocusEvent e) {
					if (settings!= null)
						settings.setRTmaxLen(parseInt(tfRTmax, settings.getRTmaxLen()));
				}
			});
			panPars.add(tfRTmax, c);

			c.gridy++; c.gridx= 0;
			c.gridwidth= 2;
			c.insets= new Insets(1,10,1,10);
			butRTrnd= new JRadioButton("Random Priming");
			butRTrnd.setEnabled(false);
			butRTrnd.setOpaque(false);
			butRTrnd.addActionListener(new ActionListener(){
				public void actionPerformed(ActionEvent e) {
					if (settings!= null)
						settings.setRtMode(FluxSimulatorSettings.PAR_RT_MODE_RANDOM);
				}
			});
			groupRT.add(butRTrnd);
			panPars.add(butRTrnd, c);

			c.gridy++; c.gridx= 0; 
			c.gridwidth= 1;
			c.insets= new Insets(1,10,1,10);
			butRTpdT= new JRadioButton("poly-dT");
			butRTpdT.setEnabled(false);
			butRTpdT.setOpaque(false);
			butRTpdT.addActionListener(new ActionListener(){
				public void actionPerformed(ActionEvent e) {
					if (settings!= null)
						settings.setRtMode(FluxSimulatorSettings.PAR_RT_MODE_POLY_DT);
				}
			});
			groupRT.add(butRTpdT);
			panPars.add(butRTpdT, c);


			c.gridy++;
			c.gridx = 0; 
			c.gridwidth= 2;
//			c.weightx= 0.1;
//			c.weighty= 0.1;
			c.anchor= GridBagConstraints.CENTER;
			c.fill= GridBagConstraints.HORIZONTAL;
			JSeparator sep= new JSeparator(JSeparator.HORIZONTAL);
//			sep.setPreferredSize(new Dimension(20,5));
//			sep.setMaximumSize(new Dimension(20,5));
			panPars.add(sep, c);
			c.gridy++;
			//c.weightx= 1;
			//c.weighty= 1;
			c.fill= GridBagConstraints.NONE;
			c.anchor= GridBagConstraints.WEST;
			cbFragmentation= new JCheckBox("Fragmentation");
			cbFragmentation.setOpaque(false);
			cbFragmentation.addActionListener(new ActionListener(){
				public void actionPerformed(ActionEvent e) {
					settings.setFragment(cbFragmentation.isSelected());
					if (settings.isFragment()) {
						tfNebLamda.setEnabled(true);
						tfNebLamda.setEditable(true);
						tfNebLamda.setText(Integer.toString((int) settings.getFragNBlambda()));
						butFragChem.setEnabled(true);
						butFragChem.setSelected(settings.getFragMode().equals(FluxSimulatorSettings.PAR_FRAG_METHOD_UR));
						butFragPhys.setEnabled(true);
						butFragPhys.setSelected(settings.getFragMode().equals(FluxSimulatorSettings.PAR_FRAG_METHOD_NB));
						cbFragB4RT.setEnabled(true);						
						cbFragB4RT.setSelected(settings.isFragB4RT());
					} else {
						tfNebLamda.setEnabled(false);
						tfNebLamda.setEditable(false);
						butFragChem.setEnabled(false);
						butFragPhys.setEnabled(false);
					}
				}
			});
			cbFragmentation.setEnabled(false);
			panPars.add(cbFragmentation, c);
			
			c.gridy++; 
			c.gridwidth= 1;
			c.insets= new Insets(1,30,1,10);
			panPars.add(new JLabel("Lambda"), c);
			c.gridx = 1;
			c.gridwidth= 1;
			tfNebLamda= new JTextField(tfWidth);
			tfNebLamda.setEditable(false);
			tfNebLamda.setHorizontalAlignment(JTextField.RIGHT);
			tfNebLamda.addActionListener(new ActionListener(){
				public void actionPerformed(ActionEvent e) {
					if (settings!= null)
						settings.setLambda(parseInt(tfNebLamda, (int) settings.getFragNBlambda()));
				}
			});
			tfNebLamda.addFocusListener(new FocusAdapter(){
				public void focusLost(FocusEvent e) {
					if (settings!= null)
						settings.setLambda(parseInt(tfNebLamda, (int) settings.getFragNBlambda()));
				}
			});
			panPars.add(tfNebLamda, c);
			/*c.gridy = 8; 
			panPars.add(new JLabel("Sigma"), c);
			c.gridy = 9; 
			panPars.add(new JLabel("Threshold"), c);*/
			
			c.gridy++; c.gridx= 0; 
			c.gridwidth= 2;
			c.insets= new Insets(1,10,1,10);
			butFragPhys= new JRadioButton("Physical Shearing");
			butFragPhys.setEnabled(false);
			butFragPhys.setOpaque(false);
			butFragPhys.addActionListener(new ActionListener(){
				public void actionPerformed(ActionEvent e) {
					if (settings!= null) {
						settings.setFragMode(FluxSimulatorSettings.PAR_FRAG_METHOD_NB);
						settings.setLambda(FluxSimulatorSettings.DEF_FRAG_LAMBDA);
						tfNebLamda.setText(Integer.toString((int) settings.getFragNBlambda()));
					}
				}
			});
			groupFrag.add(butFragPhys);
			panPars.add(butFragPhys, c);
			

			
			c.gridy++; 
			c.gridx= 0;
			c.gridwidth= 2;
			c.insets= new Insets(1,10,1,10);
			butFragChem= new JRadioButton("Chemical Cleavage/AFA");
			butFragChem.setEnabled(false);
			butFragChem.setOpaque(false);
			butFragChem.addActionListener(new ActionListener(){
				public void actionPerformed(ActionEvent e) {
					if (settings!= null) {
						settings.setFragMode(FluxSimulatorSettings.PAR_FRAG_METHOD_UR);
						settings.setLambda(FluxSimulatorSettings.DEF_FRG_LAMBDA);
						tfNebLamda.setText(Integer.toString((int) settings.getFragNBlambda()));
					}
				}
			});
			groupFrag.add(butFragChem);
			panPars.add(butFragChem, c);			
			c.gridy++; 
			c.gridwidth= 2;
			//c.fill= GridBagConstraints.HORIZONTAL;
			cbFragB4RT= new JCheckBox("Fragmentation before RT");
			cbFragB4RT.addActionListener(new ActionListener(){
				public void actionPerformed(ActionEvent e) {
					if (settings!= null) {
						settings.setFragB4RT(cbFragB4RT.isSelected());
						if (cbFragB4RT.isSelected()) {
							if (butRTpdT.isSelected()) {
								butRTrnd.setSelected(true);
								settings.setRtMode(FluxSimulatorSettings.PAR_RT_MODE_RANDOM);							
							}
							butRTpdT.setEnabled(false);
						} else {
							butRTpdT.setEnabled(true);
						}
					}
				}
			});
			cbFragB4RT.setEnabled(false);
			cbFragB4RT.setOpaque(false);
			panPars.add(cbFragB4RT, c);
			
			c.gridy++; c.gridx= 0;
			c.gridwidth= 2;
			c.fill= GridBagConstraints.HORIZONTAL;
			panPars.add(new JSeparator(JSeparator.HORIZONTAL), c);
			
			
			c.gridy++;
			c.fill= GridBagConstraints.NONE;
			c.insets= new Insets(1,10,1,10);
			cbFiltering= new JCheckBox(" cutoff");
			cbFiltering.setOpaque(false);
			cbFiltering.addActionListener(new ActionListener(){
				public void actionPerformed(ActionEvent e) {
					settings.setFilter(cbFiltering.isSelected());
					if (settings.isFilter()) {
						tfFiltMin.setEnabled(true);
						tfFiltMin.setEditable(true);
						tfFiltMin.setText(Integer.toString(settings.getFiltMin()));
						tfFiltMax.setEnabled(true);
						tfFiltMax.setEditable(true);
						tfFiltMax.setText(Integer.toString(settings.getFiltMax()));
					} else {
						tfFiltMin.setEnabled(false);
						tfFiltMin.setEditable(false);
						tfFiltMin.setText(Integer.toString(settings.getFiltMin()));
						tfFiltMax.setEnabled(false);
						tfFiltMax.setEditable(false);
						tfFiltMax.setText(Integer.toString(settings.getFiltMax()));
					}
				}
			});
			cbFiltering.setEnabled(false);
			panPars.add(cbFiltering, c);
			
			c.gridy++; c.gridx= 0;
			c.gridwidth= 1;
			c.insets= new Insets(1,20,1,10);
			panPars.add(new JLabel("min"), c);
			c.gridx= 1;
			tfFiltMin= new JTextField(tfWidth);
			tfFiltMin.setEditable(false);
			tfFiltMin.setHorizontalAlignment(JTextField.RIGHT);
			tfFiltMin.addActionListener(new ActionListener(){
				public void actionPerformed(ActionEvent e) {
					if (settings!= null)
						settings.setFiltMin(parseInt(tfFiltMin, settings.getFiltMin()));
				}
			});
			tfFiltMin.addFocusListener(new FocusAdapter(){
				public void focusLost(FocusEvent e) {
					if (settings!= null)
						settings.setFiltMin(parseInt(tfFiltMin, settings.getFiltMin()));
				}
			});
			panPars.add(tfFiltMin, c);
			
			c.gridy++; c.gridx= 0; 
			c.gridwidth= 1;
			panPars.add(new JLabel("max"), c);			
			c.gridx = 1;
			tfFiltMax= new JTextField(tfWidth);
			panPars.add(tfFiltMax, c);
			tfFiltMax.setEditable(false);
			tfFiltMax.setHorizontalAlignment(JTextField.RIGHT);
			tfFiltMax.addActionListener(new ActionListener(){
				public void actionPerformed(ActionEvent e) {
					if (settings!= null)
						settings.setFiltMax(parseInt(tfFiltMax, settings.getFiltMax()));
				}
			});
			tfFiltMax.addFocusListener(new FocusAdapter(){
				public void focusLost(FocusEvent e) {
					if (settings!= null)
						settings.setFiltMax(parseInt(tfFiltMax, settings.getFiltMax()));
				}
			});			

						
			
			/*c.gridy = 8;
			tfNebSigma= new JTextField(tfWidth);
			tfNebSigma.setEditable(false);
			tfNebSigma.addActionListener(new ActionListener(){
				public void actionPerformed(ActionEvent e) {
					if (settings!= null)
						settings.setSigma(parseDouble(tfNebSigma, settings.getSigma()));
				}
			});
			tfNebSigma.addFocusListener(new FocusAdapter(){
				public void focusLost(FocusEvent e) {
					if (settings!= null)
						settings.setSigma(parseDouble(tfNebSigma, settings.getSigma()));
				}
			});
			panPars.add(tfNebSigma, c);
			
			c.gridy = 9;
			tfNebThres= new JTextField(tfWidth);
			tfNebThres.setEditable(false);
			tfNebThres.addActionListener(new ActionListener(){
				public void actionPerformed(ActionEvent e) {
					if (settings!= null)
						settings.setThold(parseDouble(tfNebThres, settings.getThold()));
				}
			});
			tfNebThres.addFocusListener(new FocusAdapter(){
				public void focusLost(FocusEvent e) {
					if (settings!= null)
						settings.setThold(parseDouble(tfNebThres, settings.getThold()));
				}
			});
			panPars.add(tfNebThres, c);*/

			
		}

		return panPars;
	}

	JPanel panStats;
	JCheckBox cbFragB4RT;
	JTextField tfFiltMin, tfFiltMax;
	JTextField tfMolRNA, tfMolRT, tfMolFrag, tfMolFilt;
	private Component getPanelStatistics() {
		if (panStats == null) {
			panStats = new JPanel();
			panStats.setBorder(BorderFactory.createTitledBorder("Characteristics"));
			panStats.setLayout(new GridBagLayout());
			int tfWidthNr= 8;
			GridBagConstraints c= new GridBagConstraints();
			c.anchor= GridBagConstraints.FIRST_LINE_START;
			c.fill = GridBagConstraints.NONE;
			c.weightx= 1;
			c.weighty= 1;
			c.insets= new Insets(1,10,1,10);
			c.gridwidth= 1; c.gridheight= 1;			
			
			c.gridx = 0; 
			c.gridy= 0;
			panStats.add(new JLabel("Molecules"), c);
			c.gridy = 1;
			c.insets= new Insets(1,20,1,10);
			panStats.add(new JLabel("Initial RNA"), c);
			c.gridx = 1; 
			//c.gridy = 2;
			tfMolRNA= new JTextField(tfWidthNr);
			tfMolRNA.setHorizontalAlignment(JTextField.RIGHT);
			tfMolRNA.setEditable(false);
			panStats.add(tfMolRNA, c);
			c.gridx = 0; 
			c.gridy= 3;
			panStats.add(new JLabel(Fragmenter.MODE_RT_MESSAGE), c);
			c.gridx = 1; 
			//c.gridy= 4;
			tfMolRT= new JTextField(tfWidthNr);
			tfMolRT.setHorizontalAlignment(JTextField.RIGHT);
			tfMolRT.setEditable(false);
			panStats.add(tfMolRT, c);
			c.gridx = 0; 
			c.gridy = 5;
			panStats.add(new JLabel(Fragmenter.MODE_NEBU_MESSAGE), c);
			c.gridx = 1; 
			//c.gridy = 6;
			tfMolFrag= new JTextField(tfWidthNr);
			tfMolFrag.setHorizontalAlignment(JTextField.RIGHT);
			tfMolFrag.setEditable(false);
			panStats.add(tfMolFrag, c);
			c.gridx = 0; 
			c.gridy = 7;
			panStats.add(new JLabel(Fragmenter.MODE_FILT_MESSAGE), c);
			c.gridx = 1; 
			//c.gridy = 8;
			tfMolFilt= new JTextField(tfWidthNr);
			tfMolFilt.setHorizontalAlignment(JTextField.RIGHT);
			tfMolFilt.setEditable(false);
			panStats.add(tfMolFilt, c);

		}

		return panStats;
	}

	public static void main(String[] args) {
		FragmenterGUI gui= new FragmenterGUI();
		
		
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
		if (fragmenter== null)
			return false;
		return fragmenter.isReady() == null;
	}

	Fragmenter fragmenter;
	public void set(Object o) {
		if (o== null|| !(o instanceof FluxSimulatorSettings))
			return;
		FluxSimulatorSettings settings= (FluxSimulatorSettings) o;

		this.settings= settings;
		if (settings== null) {
			fragmenter= null;
			uiUpdateStats();
		} else
			fragmenter= new Fragmenter(settings);
		uiUpdatePars();
	}

	private void uiUpdateStats() {
		if (fragmenter== null) {
			if (tabPlots!= null&& tabPlots.getTabCount()> 0)
				tabPlots.removeAll();
			tfMolRNA.setText(FluxSimulatorGUI.emptyString);
			tfMolFrag.setText(FluxSimulatorGUI.emptyString);
			tfMolRT.setText(FluxSimulatorGUI.emptyString);
			tfMolFilt.setText(FluxSimulatorGUI.emptyString);
		} else {
			// plots draw themselves dynamically
			// also mol counts are set dynamically
		}
		FluxSimulatorGUI.repaintEverywhere(this, false);	
	}

	private void uiUpdatePars() {		
		if (settings== null) {
			tfIn.setText(FluxSimulatorGUI.emptyString);
			tfOut.setText(FluxSimulatorGUI.emptyString);
			tfRTmax.setText(FluxSimulatorGUI.emptyString);
			tfRTmax.setEditable(false);
			tfRTmin.setText(FluxSimulatorGUI.emptyString);
			tfRTmin.setEditable(false);
			tfNebLamda.setText(FluxSimulatorGUI.emptyString);
			tfNebLamda.setEditable(false);
			//tfNebSigma.setText(FluxSimulatorGUI.emptyString);
			//tfNebSigma.setEditable(false);
			//tfNebThres.setText(FluxSimulatorGUI.emptyString);
			//tfNebThres.setEditable(false);
			tfFiltMin.setText(FluxSimulatorGUI.emptyString);
			tfFiltMin.setEditable(false);
			tfFiltMax.setText(FluxSimulatorGUI.emptyString);
			tfFiltMax.setEditable(false);
			tfPolyAscale.setEditable(false);
			tfPolyAshape.setEditable(false);
			tfTSSmean.setEditable(false);
			cbFragB4RT.setSelected(false);
			cbFragB4RT.setEnabled(false);
			butFragChem.setEnabled(false);
			butFragChem.setSelected(false);
			butFragPhys.setEnabled(false);
			butFragPhys.setSelected(false);
			butRTpdT.setEnabled(false);
			butRTpdT.setSelected(false);
			butRTrnd.setEnabled(false);	
			butRTrnd.setSelected(false);
			cbFragmentation.setEnabled(false);
			cbFiltering.setEnabled(false);
		} else {
			if (settings.getProFile()!= null)
				tfIn.setText(settings.getProFile().getAbsolutePath());
			if (settings.getFrgFile()!= null)
				tfOut.setText(settings.getFrgFile().getAbsolutePath());
			tfRTmax.setText(Integer.toString(settings.getRTmaxLen()));
			tfRTmax.setEditable(true);
			tfRTmin.setText(Integer.toString(settings.getRTminLen()));
			tfRTmin.setEditable(true);
			tfNebLamda.setText(Integer.toString((int) settings.getFragNBlambda()));
			tfNebLamda.setEditable(false);	// deactivated
//			tfNebSigma.setText(Double.toString(settings.getSigma()));
//			tfNebSigma.setEditable(true);
//			tfNebThres.setText(Double.toString(settings.getThold()));
//			tfNebThres.setEditable(true);
			butFragChem.setEnabled(true);
			butFragPhys.setEnabled(true);
			butRTpdT.setEnabled(true);
			butRTrnd.setEnabled(true);
			if (settings.getRtMode().equals(FluxSimulatorSettings.PAR_RT_MODE_RANDOM))
				butRTrnd.setSelected(true);
			else if (settings.getRtMode().equals(FluxSimulatorSettings.PAR_RT_MODE_POLY_DT))
				butRTpdT.setSelected(true);
			
			cbFragmentation.setEnabled(true);
			cbFragmentation.setSelected(settings.isFragment());
			if (settings.isFragment()) {
				butFragPhys.setEnabled(true);
				butFragChem.setEnabled(true);
				if (settings.getFragMode().equals(FluxSimulatorSettings.PAR_FRAG_METHOD_NB))
					butFragPhys.setSelected(true);
				else if (settings.getFragMode().equals(FluxSimulatorSettings.PAR_FRAG_METHOD_UR))
					butFragChem.setSelected(true);
				cbFragB4RT.setSelected(settings.isFragB4RT());
				cbFragB4RT.setEnabled(true);
			} else {
				butFragPhys.setEnabled(false);
				butFragChem.setEnabled(false);
				butFragPhys.setSelected(false);
				butFragChem.setSelected(false);
				cbFragB4RT.setSelected(false);
				cbFragB4RT.setEnabled(false);
			}
			
			cbFiltering.setEnabled(true);
			cbFiltering.setSelected(settings.isFilter());
			if (settings.isFilter()) {
				tfFiltMin.setText(Integer.toString(settings.getFiltMin()));
				tfFiltMin.setEditable(true);
				tfFiltMax.setText(Integer.toString(settings.getFiltMax()));
				tfFiltMax.setEditable(true);
			} else {
				tfFiltMin.setEditable(false);
				tfFiltMax.setEditable(false);
			}
			
			if (settings.isPolyAVar()) {
				tfPolyAscale.setText(Double.toString(settings.getPolyAscale()));
				tfPolyAscale.setEditable(true);
				tfPolyAshape.setText(Double.toString(settings.getPolyAshape()));
				tfPolyAshape.setEditable(true);
			} else {
				tfPolyAscale.setText("");
				tfPolyAscale.setEditable(true);
				tfPolyAshape.setText("");
				tfPolyAshape.setEditable(true);
			}
			if (settings.isTssVar()) {
				tfTSSmean.setText(Double.toString(settings.getTssMean()));
				tfTSSmean.setEditable(true);
			} else {
				tfTSSmean.setText("");
				tfTSSmean.setEditable(true);
			}

		}
		
		FluxSimulatorGUI.repaintEverywhere(this, false);
	} 

	public boolean setStop() {
		if (fragmenter== null)
			return false;
		return fragmenter.setStop();
	}

	private void loadTextFields() {
		settings.setRTminLen(FluxSimulatorGUI.parseInt(tfRTmin,settings.getRTminLen()));
		settings.setRTmaxLen(FluxSimulatorGUI.parseInt(tfRTmax,settings.getRTmaxLen()));
		settings.setLambda(FluxSimulatorGUI.parseDouble(tfNebLamda,settings.getFragNBlambda()));
		settings.setFiltMin(FluxSimulatorGUI.parseInt(tfFiltMin,settings.getFiltMin()));
		settings.setFiltMax(FluxSimulatorGUI.parseInt(tfFiltMax,settings.getFiltMax()));
	}
	
	public void run() {		
		fragmenter.setStop(false);
		if (tabPlots.getTabCount()> 0) {
			fragmenter= null;
			uiUpdateStats();
			fragmenter= new Fragmenter(settings);
		}

		loadTextFields();
		if (settings.getProfiler()== null|| settings.getProfiler().getIds()== null)
			return;
		
		if (loadStats) {
			loadStats();
			setLoadStats(false);
		} else {
			fragmenter.setPlotter(new Plotter9());
			fragmenter.run();
		}
		
		if (fragmenter.isStop()) {
			if (settings.getFrgFile().exists())
				settings.getFrgFile().delete();
			settings.save();
			uiUpdateStats();
		} else
			settings.setReadNr(Math.min(settings.getReadNr(),fragmenter.getMoleculeNr()));
	}

	public boolean loadStats() {
		
		if (fragmenter== null|| (fragmenter.isReady() != null )|| (!settings.getFrgFile().exists()))
			return false;
		
		Plotter9 plotter= new Plotter9();
		plotter.reset(Fragmenter.MODE_FILT_MESSAGE);
		
		for (int i = 0; (!isStop())&& i < settings.getProfiler().getIds().length; i++) {
			plotter.addBase(settings.getProfiler().getCombinedID(i),
					settings.getProfiler().getLen()[i],
					(int) settings.getProfiler().getMolecules()[i]);
		}
		fragmenter.setPlotter(plotter);		
		
		if (!fragmenter.loadStats()) {
			getTabPlots().removeAll();
			return false;
		} 
		uiUpdateStats();
		return true;
	}
	
	public boolean isFinished() {
		if (fragmenter== null)
			return false;
		return fragmenter.isFinished();
	}

	public boolean getLoadStats() {
		return loadStats;
	}

	public void setLoadStats(boolean val) {
		loadStats= val;		
	}

	public void killResult() {
		if (settings!= null&& settings.getFrgFile()!= null&& settings.getFrgFile().exists())
			settings.getFrgFile().delete();
		if (fragmenter!= null) {			
			fragmenter= null;
			uiUpdateStats();
		}
		tabPlots.removeAll();
		this.fragmenter= new Fragmenter(settings);
	}

	public boolean isStop() {
		return fragmenter== null|| fragmenter.isStop();
	}

	public boolean setStop(boolean stop) {
		if (stop)
			return setStop();
		if (fragmenter== null)
			return false;
		return fragmenter.setStop(stop);
	}

}
