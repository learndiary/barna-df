package fbi.genome.sequencing.rnaseq.reconstruction.gui;

import fbi.genome.model.constants.Constants;
import fbi.genome.sequencing.rnaseq.reconstruction.FluxCapacitor;
import fbi.genome.sequencing.rnaseq.reconstruction.TProfile;
import fbi.genome.sequencing.rnaseq.reconstruction.TProfileFunction;


import java.awt.BorderLayout;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.GridLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.FocusAdapter;
import java.awt.event.FocusEvent;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.File;
import java.util.Hashtable;
import java.util.Random;

import javax.swing.BorderFactory;
import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JSeparator;
import javax.swing.JTable;
import javax.swing.JTextField;
import javax.swing.SwingUtilities;

import commons.ReadyOrNot;
import commons.gui.SimpleBinPlotterPanel;
import commons.thread.StoppableRunnable;
import fbi.genome.sequencing.rnaseq.simulation.Fragmenter;

public class ProfilerGUI extends JPanel implements ReadyOrNot, StoppableRunnable  {
	public class Plotter9 {

		Hashtable<String,int[]> hash;
		public Plotter9() {
			hash= new Hashtable<String, int[]>();
		}
		
		public void addBase(String ID, int length, int mol) {
			hash.put(ID, new int[]{length, mol});
		}

		long moInit= -1;
		public void setMolTot(final long value) {
			if (moInit< 0)
				moInit= value;
			try {
				SwingUtilities.invokeLater(new Runnable(){
					public void run() {
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
						ProfilerGUI.this.getTabPlots().paintImmediately(
								ProfilerGUI.this.getTabPlots().getBounds());	// need size for resetOSI()
						ProfilerGUI.this.getTabPlots().revalidate();
						
						ProfilerGUI.this.panRankSL.resetOSI(1);
						ProfilerGUI.this.panRankSH.resetOSI(1);
						ProfilerGUI.this.panRankML.resetOSI(1);
						ProfilerGUI.this.panRankMH.resetOSI(1);
						ProfilerGUI.this.panRankLL.resetOSI(1);
						ProfilerGUI.this.panRankLH.resetOSI(1);
						ProfilerGUI.this.panGel.resetOSI(10000);
						ProfilerGUI.this.panRankSM.resetOSI(1);
						ProfilerGUI.this.panRankMM.resetOSI(1);
						ProfilerGUI.this.panRankLM.resetOSI(1);
					}
				});
			} catch (Exception e) {
				;	// :)
			}
			
		}
		
		public void paint() {
			try {
				SwingUtilities.invokeAndWait(new Runnable(){
					public void run() {
						ProfilerGUI.this.panRankSL.paintOSI(null);
						ProfilerGUI.this.panRankSM.paintOSI(null);
						ProfilerGUI.this.panRankSH.paintOSI(null);
						ProfilerGUI.this.panRankML.paintOSI(null);
						ProfilerGUI.this.panRankMM.paintOSI(null);
						ProfilerGUI.this.panRankMH.paintOSI(null);
						ProfilerGUI.this.panRankLL.paintOSI(null);
						ProfilerGUI.this.panRankLM.paintOSI(null);
						ProfilerGUI.this.panRankLH.paintOSI(null);
						ProfilerGUI.this.panGel.paintOSI(null);
						
						ProfilerGUI.this.panRankSL.repaint();
						ProfilerGUI.this.panRankSM.repaint();
						ProfilerGUI.this.panRankSH.repaint();
						ProfilerGUI.this.panRankML.repaint();
						ProfilerGUI.this.panRankMM.repaint();
						ProfilerGUI.this.panRankMH.repaint();
						ProfilerGUI.this.panRankLL.repaint();
						ProfilerGUI.this.panRankLM.repaint();
						ProfilerGUI.this.panRankLH.repaint();
						ProfilerGUI.this.panGel.repaint();

					}
				});
			} catch (Exception e) {
				; // :)
			}
		}

		private int[] tmp;
		public void plot(int pos, int len, float rpkm) {
			double fracStart= ((double) pos)/ len;
			
			if (len< TProfileFunction.LEN_LO) {
				if (rpkm>= TProfileFunction.EXP_UP) {
					if (fracStart>= 0&& fracStart<= 1)
						ProfilerGUI.this.panRankSH.addVal(fracStart);
				} else if (rpkm>= TProfileFunction.EXP_LO) {
					if (fracStart>= 0&& fracStart<= 1)
						ProfilerGUI.this.panRankSM.addVal(fracStart);
				} else {
					if (fracStart>= 0&& fracStart<= 1)
						ProfilerGUI.this.panRankSL.addVal(fracStart);
				}
			} else if (len< TProfileFunction.LEN_UP) {
				if (rpkm>= TProfileFunction.EXP_UP) {
					if (fracStart>= 0&& fracStart<= 1)
						ProfilerGUI.this.panRankMH.addVal(fracStart);
				} else if (rpkm>= TProfileFunction.EXP_LO) {
					if (fracStart>= 0&& fracStart<= 1)
						ProfilerGUI.this.panRankMM.addVal(fracStart);
				} else {
					if (fracStart>= 0&& fracStart<= 1)
						ProfilerGUI.this.panRankML.addVal(fracStart);
				}
			} else {
				if (rpkm>= TProfileFunction.EXP_UP) {
					if (fracStart>= 0&& fracStart<= 1)
						ProfilerGUI.this.panRankLH.addVal(fracStart);
				} else if (rpkm>= TProfileFunction.EXP_LO) {
					if (fracStart>= 0&& fracStart<= 1)
						ProfilerGUI.this.panRankLM.addVal(fracStart);
				} else {
					if (fracStart>= 0&& fracStart<= 1)
						ProfilerGUI.this.panRankLL.addVal(fracStart);
				}
			}
			
		}
		
	}
	
	boolean loadStats= false;
	JTable statsTable;
	JTextField tfIn, tfOut;
	JButton butChangeIn;
	SimpleBinPlotterPanel panGel;
	SimpleBinPlotterPanel panRankSL, panRankSM, panRankSH,
		panRankML, panRankMM, panRankMH, 
		panRankLL, panRankLM, panRankLH;
	public ProfilerGUI() {
		
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
		panInPan.add(new JLabel("File"), c);
		c.gridx = 1; 
		c.weightx= 1;
		c.fill = GridBagConstraints.HORIZONTAL;
		tfIn= new JTextField();
		tfIn.setEditable(false);
		panInPan.add(tfIn, c);
		c.weightx= 0;
		c.fill = GridBagConstraints.NONE;
		c.gridx = 2; 
		butChangeIn= new JButton("Change");
		butChangeIn.setEnabled(false);
		panInPan.add(butChangeIn, c);
		
		JPanel panSide= new JPanel();
		panSide.setLayout(new BoxLayout(panSide, BoxLayout.Y_AXIS));
		panSide.add(getPanelAttributes());
		panSide.add(getPanelISize());
		//panSide.add(new Box.Filler(new Dimension(1,1), new Dimension(1,600), new Dimension(1,600)));
		add(panSide, BorderLayout.WEST);

		JPanel panViz= getTabPanel();
		JPanel panHelp= new JPanel(new BorderLayout());
		add(panHelp, BorderLayout.CENTER);
		panHelp.add(panInPan, BorderLayout.NORTH);		
		panHelp.add(panViz, BorderLayout.CENTER);
		
	}
	
	JPanel tabPlots;
	public JPanel getTabPlots() {
		if (tabPlots == null) {
			tabPlots = new JPanel();
			// didnt solve repaint pb when tab opens
//			tabPlots.addChangeListener( new ChangeListener(){
//				public void stateChanged(ChangeEvent e) {
//					tabPlots.repaint();
//				}
//			});

		}

		return tabPlots;
	}
	
	JPanel tabPanel;
	private JPanel getTabPanel() {
		
		if (tabPanel == null) {
			tabPanel= new JPanel(); 
			tabPanel.setBorder(BorderFactory.createTitledBorder("Read Distribution"));
			tabPanel.setLayout(new GridLayout(3,3,10,10));
			Dimension d= new Dimension(203, 156);
			panRankSL= new SimpleBinPlotterPanel("short-low");
			//panRankSL.setPreferredSize(d);
			tabPanel.add(panRankSL);
			panRankSM= new SimpleBinPlotterPanel("short-medium");
			//panRankSM.setPreferredSize(d);
			tabPanel.add(panRankSM);
			panRankSH= new SimpleBinPlotterPanel("short-high");
			//panRankSH.setPreferredSize(d);
			tabPanel.add(panRankSH);
			panRankML= new SimpleBinPlotterPanel("medium-low");
			//panRankML.setPreferredSize(d);
			tabPanel.add(panRankML);
			panRankMM= new SimpleBinPlotterPanel("medium-medium");
			//panRankMM.setPreferredSize(d);
			tabPanel.add(panRankMM);
			panRankMH= new SimpleBinPlotterPanel("medium-high");
			//panRankMH.setPreferredSize(d);
			tabPanel.add(panRankMH);
			panRankLL= new SimpleBinPlotterPanel("long-low");
			//panRankLL.setPreferredSize(d);
			tabPanel.add(panRankLL);
			panRankLM= new SimpleBinPlotterPanel("long-medium");
			//panRankLM.setPreferredSize(d);
			tabPanel.add(panRankLM);
			panRankLH= new SimpleBinPlotterPanel("long-high");
			//panRankLH.setPreferredSize(d);
			tabPanel.add(panRankLH);

		}

		return tabPanel;
	}
	
	JPanel panAttrib;
	JTextField tfBinLenLo, tfBinLenUp, tfBinExpUp, tfBinExpLo;
	JTextField tfNrLoci, tfNrProfiles, tfNrMappings, tfNrMapAnnotation;
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
	
	private Component getPanelAttributes() {
		if (panAttrib == null) {
			panAttrib = new JPanel();
			panAttrib.setBorder(BorderFactory.createTitledBorder("Attributes"));
			panAttrib.setLayout(new GridBagLayout());
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
			panAttrib.add(new JLabel("Profile"), c);
			c.gridy++;
			c.gridwidth= 1;
			c.insets= new Insets(1,30,1,10);
			panAttrib.add(new JLabel("Loci"), c);
			c.gridx= 1;
			c.insets= new Insets(1,10,1,10);
			tfNrLoci= new JTextField(tfWidth);
			tfNrLoci.setHorizontalAlignment(JTextField.RIGHT);
			tfNrLoci.setEditable(false);
			panAttrib.add(tfNrLoci, c);
			c.gridy++;
			c.gridx= 0;
			c.insets= new Insets(1,30,1,10);
			panAttrib.add(new JLabel("Profiles"), c);
			c.gridx= 1;
			c.insets= new Insets(1,10,1,10);
			tfNrProfiles= new JTextField(tfWidth);
			tfNrProfiles.setHorizontalAlignment(JTextField.RIGHT);
			tfNrProfiles.setEditable(false);
			panAttrib.add(tfNrProfiles, c);
			c.gridy++;
			c.gridx= 0;
			c.insets= new Insets(1,30,1,10);
			panAttrib.add(new JLabel("Mappings"), c);
			c.gridx= 1;
			c.insets= new Insets(1,10,1,10);
			tfNrMappings= new JTextField(tfWidth);
			tfNrMappings.setHorizontalAlignment(JTextField.RIGHT);
			tfNrMappings.setEditable(false);
			panAttrib.add(tfNrMappings, c);
			c.gridy++;
			c.gridx= 0;
			c.insets= new Insets(1,30,1,10);
			panAttrib.add(new JLabel("Annotated"), c);
			c.gridx= 1;
			c.insets= new Insets(1,10,1,10);
			tfNrMapAnnotation= new JTextField(tfWidth);
			tfNrMapAnnotation.setHorizontalAlignment(JTextField.RIGHT);
			tfNrMapAnnotation.setEditable(false);
			panAttrib.add(tfNrMapAnnotation, c);
			
			c.gridy++; c.gridx= 0;
			c.gridwidth= 2;
			c.fill= GridBagConstraints.HORIZONTAL;
			panAttrib.add(new JSeparator(JSeparator.HORIZONTAL), c);
			

			c.fill = GridBagConstraints.NONE;
			c.gridy++;
			c.gridx= 0;
			c.insets= new Insets(1,10,1,10);
			c.gridwidth= 2;
			panAttrib.add(new JLabel("Length Bins"), c);
			c.gridy++;
			c.gridwidth= 1;
			c.insets= new Insets(1,30,1,10);
			panAttrib.add(new JLabel("lower"), c);
			c.gridx= 1;
			c.insets= new Insets(1,10,1,10);
			tfBinLenLo= new JTextField(tfWidth);
			tfBinLenLo.setHorizontalAlignment(JTextField.RIGHT);
			tfBinLenLo.setEditable(false);
			tfBinLenLo.addActionListener(new ActionListener(){
				public void actionPerformed(ActionEvent e) {
				}
			});
			tfBinLenLo.addFocusListener(new FocusAdapter(){
				public void focusLost(FocusEvent e) {
				}
			});
			panAttrib.add(tfBinLenLo, c);

			c.gridy++; c.gridx= 0;
			c.insets= new Insets(1,30,1,10);
			panAttrib.add(new JLabel("upper"), c);
			c.gridx= 1;
			c.insets= new Insets(1,10,1,10);
			tfBinLenUp= new JTextField(tfWidth);
			tfBinLenUp.setHorizontalAlignment(JTextField.RIGHT);
			tfBinLenUp.setEditable(false);
			tfBinLenUp.addActionListener(new ActionListener(){
				public void actionPerformed(ActionEvent e) {
				}
			});
			tfBinLenUp.addFocusListener(new FocusAdapter(){
				public void focusLost(FocusEvent e) {
				}
			});
			panAttrib.add(tfBinLenUp, c);
			
			
			c.gridy++; c.gridx= 0;
			c.gridwidth= 2;
			c.fill= GridBagConstraints.HORIZONTAL;
			panAttrib.add(new JSeparator(JSeparator.HORIZONTAL), c);
			
			c.fill = GridBagConstraints.NONE;
			c.gridy++; c.gridx= 0;
			c.gridwidth= 2;
			c.insets= new Insets(1,10,1,10);
			panAttrib.add(new JLabel("Expression bins"), c);
			c.gridy++ ;
			c.gridwidth= 1;
			//panPars.add(new JLabel("Falloff Length"), c);
			
			c.gridy++; c.gridx= 0; 
			c.gridwidth= 1;
			c.insets= new Insets(1,30,1,10);
			panAttrib.add(new JLabel("lower"), c);			
			c.gridx = 1; 
			c.insets= new Insets(1,10,1,10);
			tfBinExpLo= new JTextField(tfWidth);
			tfBinExpLo.setHorizontalAlignment(JTextField.RIGHT);
			tfBinExpLo.setEditable(false);
			tfBinExpLo.addActionListener(new ActionListener(){
				public void actionPerformed(ActionEvent e) {
				} 
			});
			tfBinExpLo.addFocusListener(new FocusAdapter(){
				public void focusLost(FocusEvent e) {
				}
			});
			panAttrib.add(tfBinExpLo, c);

			c.gridy++; c.gridx= 0;
			c.insets= new Insets(1,30,1,10);
			panAttrib.add(new JLabel("upper"), c);			
			c.gridx= 1;
			c.insets= new Insets(1,10,1,10);
			tfBinExpUp= new JTextField(tfWidth);
			tfBinExpUp.setHorizontalAlignment(JTextField.RIGHT);
			tfBinExpUp.setEditable(false);
			panAttrib.add(tfBinExpUp, c);
		}

		return panAttrib;
	}

	JPanel panISize;
	SimpleBinPlotterPanel lengthDistrPanel;
	private Component getPanelISize() {
		if (panISize == null) {
			panISize = new JPanel();
			panISize.setBorder(BorderFactory.createTitledBorder("Insert Sizes (transcriptomic)"));
			
			lengthDistrPanel= new SimpleBinPlotterPanel();
			lengthDistrPanel.setTitleX("length [nt]");
			panISize.add(lengthDistrPanel);
			
		}

		return panISize;
	}

	public static void main(String[] args) {
		ProfilerGUI gui= new ProfilerGUI();
		
		
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
		return FluxCapacitorGUI.singleton.capacitor.isInputReady();
	}

	private void uiUpdate() {
		File f= FluxCapacitorGUI.singleton.capacitor.getFileProfile();
		if (f== null) {
			tfIn.setText("");
			tfIn.setToolTipText(null);
		} else {
			tfIn.setText(f.getName());
			tfIn.setToolTipText(f.getPath());
		}
		
		FluxCapacitorGUI.repaintEverywhere(this, false);	
	}

	public boolean setStop() {
		return false;
	}

	private void loadTextFields() {
		
	}
	
	public void run() {
		
		FluxCapacitor capacitor= FluxCapacitorGUI.singleton.capacitor;
		try {
			capacitor.explore(FluxCapacitor.MODE_LEARN);
		} catch (Throwable e) {
			FluxCapacitorGUI.singleton.dialogError(e);
			if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
				System.err.println("[FATAL] Error occured during scanning\n\t"+ e.getMessage());
		}

		tfBinLenLo.setText(Integer.toString(TProfileFunction.LEN_LO));
		tfBinLenUp.setText(Integer.toString(TProfileFunction.LEN_UP));
		tfBinExpLo.setText(Integer.toString(TProfileFunction.EXP_LO));
		tfBinExpUp.setText(Integer.toString(TProfileFunction.EXP_UP));
		
		TProfileFunction func= FluxCapacitorGUI.singleton.capacitor.getFunc();
		tfNrLoci.setText(Integer.toString(capacitor.getNrSingleTranscriptLoci()));
		tfNrProfiles.setText(Integer.toString(func.getNrProfiles()));
		tfNrMappings.setText(Integer.toString(capacitor.getNrReadsSingleLoci()));
		tfNrMapAnnotation.setText(Integer.toString(capacitor.getNrReadsSingleLociMapped()));
		
		
		TProfile[] profis= func.getProfis();
		Plotter9 plotter= new Plotter9();
		plotter.reset(Fragmenter.MODE_FILT_MESSAGE);
		for (int i = 0; i < profis.length; i++) {
			for (int j = 0; j < profis[i].getLength(); j++) {
				int[] a= new int[] {j, j+1};
				int x= (int) profis[i].getArea(a, capacitor.getReadLength(), capacitor.getInsertMinMax(), (byte) 0);
				for (int k = 0; k < x; k++) {
					plotter.plot(j, profis[i].getLength(), 
							FluxCapacitorGUI.singleton.capacitor.calcRPKM(profis[i].getReads(), 
									profis[i].getLength()));
				}
			}
		}
		plotter.paint();
	}

	public boolean loadStats() {
		
		if (true)
			return false;
		
		Plotter9 plotter= new Plotter9();
		plotter.reset(Fragmenter.MODE_FILT_MESSAGE);
		
		uiUpdate();
		return true;
	}
	
	public boolean isFinished() {		
		return (FluxCapacitorGUI.singleton.capacitor.getFunc()!= null
				&& FluxCapacitorGUI.singleton.capacitor.getFunc().getProfis()!= null
				&& FluxCapacitorGUI.singleton.capacitor.getFunc().getProfis().length> 0);
	}

	public boolean getLoadStats() {
		return loadStats;
	}

	public void setLoadStats(boolean val) {
		loadStats= val;		
	}

	public void killResult() {
		if (true)
			FluxCapacitorGUI.singleton.capacitor.getFileProfile().delete();
		tabPlots.removeAll();
	}

	public boolean isStop() {
		return false;
	}

	public boolean setStop(boolean stop) {
		return false;
	}

	public void set(Object o) {
		uiUpdate();
	}
}
