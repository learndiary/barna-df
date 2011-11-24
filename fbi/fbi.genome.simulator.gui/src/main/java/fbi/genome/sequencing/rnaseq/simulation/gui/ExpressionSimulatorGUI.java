package fbi.genome.sequencing.rnaseq.simulation.gui;

import fbi.commons.ReadyOrNot;
import fbi.commons.StringUtils;
import fbi.commons.gui.SimpleBinPlotterPanel;
import fbi.commons.thread.StoppableRunnable;
import fbi.genome.sequencing.rnaseq.simulation.FluxSimulatorSettings;
import fbi.genome.sequencing.rnaseq.simulation.Profiler;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
import java.io.File;
import java.util.*;

public class ExpressionSimulatorGUI extends JPanel implements ReadyOrNot, StoppableRunnable  {
	
	JTable statsTable;
	FluxSimulatorSettings settings;
	JTextField tfIO;
	JButton butChangeIO;
	SimpleGelBinPlotterPanel panGel;
	SimpleBinPlotterPanel panRank1, panRank2;
	public ExpressionSimulatorGUI() {
		
		setLayout(new BorderLayout());
		
		JPanel panInPan= new JPanel();
		panInPan.setBorder(BorderFactory.createTitledBorder("Input / Output")); // or on lineborder black
		panInPan.setLayout(new BorderLayout(10,0));
		panInPan.add(new JLabel("File"), BorderLayout.WEST);
		tfIO= new JTextField();
		tfIO.setEditable(false);
		panInPan.add(tfIO, BorderLayout.CENTER);
		butChangeIO= new JButton("Change");
		butChangeIO.addActionListener(new ActionListener(){
			public void actionPerformed(ActionEvent e) {
			    JFileChooser chooser = new JFileChooser();
			    int returnVal = chooser.showOpenDialog(ExpressionSimulatorGUI.this);
			    if(returnVal == JFileChooser.APPROVE_OPTION) {
			    	String absPath= chooser.getSelectedFile().getAbsolutePath();
					File f= new File(absPath);
					ExpressionSimulatorGUI.this.settings.setRefFile(f);
					tfIO.setText(f.getName());
					tfIO.repaint();
			    }				
			}
		});
		butChangeIO.setEnabled(false);
		panInPan.add(butChangeIO, BorderLayout.EAST);
		add(panInPan, BorderLayout.NORTH);		
		
		JPanel panSide= new JPanel();
		panSide.setLayout(new BoxLayout(panSide, BoxLayout.Y_AXIS));
		panSide.add(getPanelParameters());
		panSide.add(getPanelStatistics());
		add(panSide, BorderLayout.WEST);

		JPanel panViz= new JPanel();
		panViz.setLayout(new BorderLayout());
		panGel= new SimpleGelBinPlotterPanel();
		panGel.setPaintMode(SimpleBinPlotterPanel.MODE_GEL);
		panGel.setInvert(true);
		panGel.setBorder(BorderFactory.createTitledBorder("Northern blot"));
		panViz.add(panGel, BorderLayout.EAST);
		JPanel panPlots= new JPanel();
		panPlots.setBorder(BorderFactory.createTitledBorder("Expression Distribution"));
		panPlots.setLayout(new BoxLayout(panPlots, BoxLayout.Y_AXIS));
		panRank1= new SimpleBinPlotterPanel();
		panRank1.setLog((byte) 1);
		panRank1.setTitleY("ln(reads)");
		panRank1.setTitleX("rank");
		panRank1.setPaintMode(SimpleBinPlotterPanel.MODE_BARPLOT);
		panPlots.add(panRank1);
		panPlots.add(Box.createRigidArea(new Dimension(Integer.MAX_VALUE, 10)));	
		panRank2= new SimpleBinPlotterPanel();
		panRank2.setPaintMode(SimpleBinPlotterPanel.MODE_BARPLOT);
		panRank2.setLog((byte) 2);
		panRank2.setTitleX("ln(reads)");
		panRank2.setTitleY("ln(1+frequency)");
		panPlots.add(panRank2);
		panViz.add(panPlots, BorderLayout.CENTER);
		add(panViz, BorderLayout.CENTER);
		
	}
	
	JPanel panPars;
	JTextField tfNrCells, tfNrMol, tfParX, tfParK, tfParExp;
	private Component getPanelParameters() {
		if (panPars == null) {
			panPars = new JPanel();
			panPars.setBorder(BorderFactory.createTitledBorder("Parameters"));
			panPars.setLayout(new GridBagLayout());
			int tfWidth= 14;
			GridBagConstraints c= new GridBagConstraints();
			c.anchor= GridBagConstraints.FIRST_LINE_START;
			c.fill = GridBagConstraints.NONE;
			c.insets= new Insets(1,10,1,10);
			c.gridwidth= 1; c.gridheight= 1;			
			
			c.weighty= 0.01;
			c.gridx = 0; 
			c.gridy = 0; 
			panPars.add(new JLabel("Number of cells"), c);
			c.gridy = 1; 
			panPars.add(new JLabel("Number of molecules"), c);
			c.gridy = 2; 
			panPars.add(new JLabel("Pareto distribution"), c);
			c.gridy = 3; 
			c.insets= new Insets(1,20,1,10);
			panPars.add(new JLabel("k"), c);
			c.gridy = 4; 
			panPars.add(new JLabel("x_m"), c);
			c.insets= new Insets(1,10,1,10);
			c.gridy = 5; 
			panPars.add(new JLabel("Exponential decay"), c);


			c.gridx = 1; 
			c.gridy = 0;
			tfNrCells= new JTextField(tfWidth);
			tfNrCells.setEditable(false);
			tfNrCells.addKeyListener(new KeyAdapter(){
				@Override
				public void keyTyped(KeyEvent e) {
					try {
						long nr= Long.parseLong(tfNrCells.getText());
						settings.setNbCells(nr);
						nr*= FluxSimulatorSettings.AVG_MOL_CELL;
						tfNrMol.setText(Long.toString(nr));
						settings.setNbMolecules(nr);
					} catch (NumberFormatException xxx) {
						tfNrCells.setText(Long.toString(settings.getNbCells()));
					}
				}
			});
			
			panPars.add(tfNrCells, c);
			c.gridy = 1;
			tfNrMol= new JTextField(tfWidth);
			tfNrMol.setEditable(false);
			tfNrMol.addKeyListener(new KeyAdapter(){
				@Override
				public void keyTyped(KeyEvent e) {
					try {						
						String s= tfNrMol.getText();
						if (Character.isDigit(e.getKeyChar())) 
							s+= e.getKeyChar();
						else {
							tfNrMol.setText(s);
							return;
						}
						long nr= Long.parseLong(s);
						settings.setNbMolecules(nr);
						nr/= FluxSimulatorSettings.AVG_MOL_CELL;
						tfNrCells.setText(Long.toString(nr));
						settings.setNbCells(nr);
					} catch (NumberFormatException xxx) {
						tfNrMol.setText(Long.toString(settings.getNbMolecules()));
					}
				}
			});
			panPars.add(tfNrMol, c);
			c.gridy = 3;
			tfParX= new JTextField(tfWidth);
			tfParX.setEditable(false);
			tfParX.addKeyListener(new KeyAdapter(){
				@Override
				public void keyTyped(KeyEvent e) {
					try {
						double x= Double.parseDouble(tfParX.getText());
						settings.setExpDistrP1(x);
					} catch (NumberFormatException xxx) {
						tfParX.setText(Double.toString(settings.getExpDistrP1()));
					}
				}
			});
			panPars.add(tfParX, c);
			c.gridy = 4;
			tfParK= new JTextField(tfWidth);
			tfParK.setEditable(false);
			tfParK.addKeyListener(new KeyAdapter(){
				@Override
				public void keyTyped(KeyEvent e) {
					try {
						double x= Double.parseDouble(tfParK.getText());
						settings.setExpDistrP2(x);
					} catch (NumberFormatException xxx) {
						tfParK.setText(Double.toString(settings.getExpDistrP2()));
					}
				}
			});
			panPars.add(tfParK, c);
			c.gridy = 5;
			tfParExp= new JTextField(tfWidth);
			tfParExp.setEditable(false);
			tfParExp.addKeyListener(new KeyAdapter(){
				@Override
				public void keyTyped(KeyEvent e) {
					try {
						double x= Double.parseDouble(tfParExp.getText());
						settings.setDecDistrP1(x);
					} catch (NumberFormatException xxx) {
						tfParExp.setText(Double.toString(settings.getDecDistrP1()));
					}
				}
			});
			panPars.add(tfParExp, c);

			// fill
//			c.gridx=2;
//			c.weightx= 1;
//			for (int i = 0; i < 6; i++) {
//				c.gridy=i;
//				panPars.add(new JPanel(), c);
//			}
//			c.gridx=0;
//			c.gridy= 6;
//			c.gridwidth=3;
//			c.weighty= 1;
//			panPars.add(new JPanel(), c);
		}

		return panPars;
	}

	JPanel panStats;
	JTextField tfTxTot, tfRNAtot;
	JTextField tfTxHi, tfTxMed, tfTxLo, tfTxNo, tfRNAhi, tfRNAmed, tfRNAlo;
	JTextField tfTxHiFrac, tfTxMedFrac, tfTxLoFrac, tfTxNoFrac, tfRNAhiFrac, tfRNAmedFrac, tfRNAloFrac;
	
	private Component getPanelStatistics() {
		if (panStats == null) {
			panStats = new JPanel();
			panStats.setBorder(BorderFactory.createTitledBorder("Characteristics"));
			panStats.setLayout(new GridBagLayout());
			int tfWidthNr= 10, tfWidthPercNr= 5;
			GridBagConstraints c= new GridBagConstraints();
			c.anchor= GridBagConstraints.FIRST_LINE_START;
			c.fill = GridBagConstraints.NONE;
			c.insets= new Insets(1,10,1,10);
			c.gridwidth= 1; c.gridheight= 1;			
			
			c.weightx= 0.01;
			c.weighty= 0.01;
			c.gridx = 0; 
			c.gridy= 0;
			panStats.add(new JLabel("Spliceforms"), c);
			c.gridy = 1;
			c.insets= new Insets(1,20,1,10);
			panStats.add(new JLabel("high"), c);
			c.gridy = 2;
			panStats.add(new JLabel("medium"), c);
			c.gridy = 3;
			panStats.add(new JLabel("low"), c);
			c.gridy = 4;
			panStats.add(new JLabel("not"), c);
			c.insets= new Insets(1,10,1,10);
			c.gridy = 5;
			panStats.add(new JLabel("RNA mass"), c);
			c.gridy = 6;
			c.insets= new Insets(1,20,1,10);
			panStats.add(new JLabel("high"), c);
			c.gridy = 7;
			panStats.add(new JLabel("medium"), c);
			c.gridy = 8;
			panStats.add(new JLabel("low"), c);
			c.insets= new Insets(1,10,1,10);
			
			c.insets= new Insets(1,5,1,5);
			c.gridx = 1; 
			c.gridy = 0;
			tfTxTot= new JTextField(tfWidthNr);
			tfTxTot.setHorizontalAlignment(JTextField.RIGHT);
			tfTxTot.setEditable(false);
			panStats.add(tfTxTot, c);
			c.gridy = 1;
			tfTxHi= new JTextField(tfWidthNr);
			tfTxHi.setHorizontalAlignment(JTextField.RIGHT);
			tfTxHi.setEditable(false);
			panStats.add(tfTxHi, c);
			c.gridy = 2;
			tfTxMed= new JTextField(tfWidthNr);
			tfTxMed.setHorizontalAlignment(JTextField.RIGHT);
			tfTxMed.setEditable(false);
			panStats.add(tfTxMed, c);
			c.gridy = 3;
			tfTxLo= new JTextField(tfWidthNr);
			tfTxLo.setHorizontalAlignment(JTextField.RIGHT);
			tfTxLo.setEditable(false);
			panStats.add(tfTxLo, c);
			c.gridy = 4;
			tfTxNo= new JTextField(tfWidthNr);
			tfTxNo.setHorizontalAlignment(JTextField.RIGHT);
			tfTxNo.setEditable(false);
			panStats.add(tfTxNo, c);			
			c.gridy = 5;
			tfRNAtot= new JTextField(tfWidthNr);
			tfRNAtot.setHorizontalAlignment(JTextField.RIGHT);
			tfRNAtot.setEditable(false);
			panStats.add(tfRNAtot, c);
			c.gridy = 6;
			tfRNAhi= new JTextField(tfWidthNr);
			tfRNAhi.setHorizontalAlignment(JTextField.RIGHT);
			tfRNAhi.setEditable(false);
			panStats.add(tfRNAhi, c);
			c.gridy = 7;
			tfRNAmed= new JTextField(tfWidthNr);
			tfRNAmed.setHorizontalAlignment(JTextField.RIGHT);
			tfRNAmed.setEditable(false);
			panStats.add(tfRNAmed, c);
			c.gridy = 8;
			tfRNAlo= new JTextField(tfWidthNr);
			tfRNAlo.setHorizontalAlignment(JTextField.RIGHT);
			tfRNAlo.setEditable(false);
			panStats.add(tfRNAlo, c);
			
			
			c.gridx = 2; 
			c.gridy = 1;
			tfTxHiFrac= new JTextField(tfWidthPercNr);
			tfTxHiFrac.setHorizontalAlignment(JTextField.RIGHT);
			tfTxHiFrac.setEditable(false);
			panStats.add(tfTxHiFrac, c);
			c.gridy = 2;
			tfTxMedFrac= new JTextField(tfWidthPercNr);
			tfTxMedFrac.setHorizontalAlignment(JTextField.RIGHT);
			tfTxMedFrac.setEditable(false);
			panStats.add(tfTxMedFrac, c);
			c.gridy = 3;
			tfTxLoFrac= new JTextField(tfWidthPercNr);
			tfTxLoFrac.setHorizontalAlignment(JTextField.RIGHT);
			tfTxLoFrac.setEditable(false);
			panStats.add(tfTxLoFrac, c);
			c.gridy = 4;
			tfTxNoFrac= new JTextField(tfWidthPercNr);
			tfTxNoFrac.setHorizontalAlignment(JTextField.RIGHT);
			tfTxNoFrac.setEditable(false);
			panStats.add(tfTxNoFrac, c);			
			c.gridy = 6;
			tfRNAhiFrac= new JTextField(tfWidthPercNr);
			tfRNAhiFrac.setHorizontalAlignment(JTextField.RIGHT);
			tfRNAhiFrac.setEditable(false);
			panStats.add(tfRNAhiFrac, c);
			c.gridy = 7;
			tfRNAmedFrac= new JTextField(tfWidthPercNr);
			tfRNAmedFrac.setHorizontalAlignment(JTextField.RIGHT);
			tfRNAmedFrac.setEditable(false);
			panStats.add(tfRNAmedFrac, c);
			c.gridy = 8;
			tfRNAloFrac= new JTextField(tfWidthPercNr);
			tfRNAloFrac.setHorizontalAlignment(JTextField.RIGHT);
			tfRNAloFrac.setEditable(false);
			panStats.add(tfRNAloFrac, c);			
			
//			c.gridx=3;
//			c.weightx= 1;
//			for (int i = 0; i < 9; i++) {
//				c.gridy=i;
//				panPars.add(new JPanel(), c);
//			}
//			c.gridx=0;
//			c.gridy= 9;
//			c.gridwidth=4;
//			c.weighty= 1;
//			panStats.add(new JPanel(), c);

		}

		return panStats;
	}

	public static void main(String[] args) {
		ExpressionSimulatorGUI gui= new ExpressionSimulatorGUI();
		
		
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
		
		if (settings!= null&& settings.getProfiler()!= null&& settings.getProfiler().getIds()!= null&& settings.getProfiler().getIds().length> 0)
			return true;
		return false;
	}
	
	public boolean loadStats() {
		if (settings.getProfiler().getMolecules()== null) {
			settings.getProfiler().setStop(false);
			if (!settings.getProfiler().loadStats())
				return false;
		}
		if (settings.getProfiler().isFinishedExpression()) {
			uiUpdateStats();
		}
		return true;
	}

	public void set(Object o) {
		if (o== null|| !(o instanceof FluxSimulatorSettings))
			return;
		FluxSimulatorSettings settings= (FluxSimulatorSettings) o;
		this.settings= settings;
		if (settings== null) 
			uiUpdateStats();
		else if (settings.getProfiler()== null)
			settings.setProfiler(new Profiler(settings));
		uiUpdatePars();
	}

	public boolean setStop() {
		if (settings.getProfiler()== null)
			return false;
		return settings.getProfiler().setStop(true);		
	}

	private void loadTextFields() {
		settings.setNbMolecules(FluxSimulatorGUI.parseLong(tfNrMol, settings.getNbMolecules()));
		settings.setExpDistrP1(FluxSimulatorGUI.parseDouble(tfParX, settings.getExpDistrP1()));
		settings.setExpDistrP2(FluxSimulatorGUI.parseDouble(tfParK, settings.getExpDistrP2()));
		settings.setDecDistrP1(FluxSimulatorGUI.parseDouble(tfParExp, settings.getDecDistrP1()));
	}
	
	public void run() {
		settings.getProfiler().setStop(false);
		if (loadStats) {
			loadStats();
			setLoadStats(false);
			return;
		}
		
		loadTextFields();
		settings.getProfiler().profile();
		
		// stats
		if (!settings.getProfiler().isStop()) {
			uiUpdateStats();
			settings.save();
		}
		
		System.gc();
	}

	
	private void uiUpdatePars() {
		if (settings== null) {
			tfIO.setText(FluxSimulatorGUI.emptyString);
			butChangeIO.setEnabled(false);
			tfNrCells.setText(FluxSimulatorGUI.emptyString);
			tfNrCells.setEditable(false);
			tfNrMol.setText(FluxSimulatorGUI.emptyString);
			tfNrMol.setEditable(false);
			tfParK.setText(FluxSimulatorGUI.emptyString);
			tfParK.setEditable(false);
			tfParX.setText(FluxSimulatorGUI.emptyString);
			tfParX.setEditable(false);
			tfParExp.setText(FluxSimulatorGUI.emptyString);
			tfParExp.setEditable(false);
		} else {
			if (settings.getProFile()!= null)
				tfIO.setText(settings.getProFile().getAbsolutePath());
			//butChangeIO.setEnabled(true);
			tfNrCells.setText(Long.toString(settings.getNbCells()));
			tfNrCells.setEditable(true);
			tfNrMol.setText(Long.toString(settings.getNbMolecules()));
			tfNrMol.setEditable(true);
			tfParK.setText(Double.toString(settings.getExpDistrP2()));
			tfParK.setEditable(true);
			tfParX.setText(Double.toString(settings.getExpDistrP1()));
			tfParX.setEditable(true);
			tfParExp.setText(Double.toString(settings.getDecDistrP1()));
			tfParExp.setEditable(true);
		}
		FluxSimulatorGUI.repaintEverywhere(this, false);
	}
	
	private void uiUpdateStats() {
		
		if (settings== null|| settings.getProfiler()== null|| settings.getProfiler().getMolecules()== null) {
			tfTxTot.setText(FluxSimulatorGUI.emptyString);
			tfTxHi.setText(FluxSimulatorGUI.emptyString);
			tfTxMed.setText(FluxSimulatorGUI.emptyString);
			tfTxLo.setText(FluxSimulatorGUI.emptyString);
			tfTxNo.setText(FluxSimulatorGUI.emptyString);
			tfTxHiFrac.setText(FluxSimulatorGUI.emptyString);
			tfTxMedFrac.setText(FluxSimulatorGUI.emptyString);
			tfTxLoFrac.setText(FluxSimulatorGUI.emptyString);
			tfTxNoFrac.setText(FluxSimulatorGUI.emptyString);
			tfRNAtot.setText(FluxSimulatorGUI.emptyString);
			tfRNAhi.setText(FluxSimulatorGUI.emptyString);
			tfRNAmed.setText(FluxSimulatorGUI.emptyString);
			tfRNAlo.setText(FluxSimulatorGUI.emptyString);
			tfRNAhiFrac.setText(FluxSimulatorGUI.emptyString);
			tfRNAmedFrac.setText(FluxSimulatorGUI.emptyString);
			tfRNAloFrac.setText(FluxSimulatorGUI.emptyString);
			panGel.offScrImg= null;
			panRank1.offScrImg= null;
			panRank2.offScrImg= null;

		} else {
			int sfs= settings.getProfiler().getSfHi().size()+ settings.getProfiler().getSfMed().size()+ settings.getProfiler().getSfLo().size();
			tfTxTot.setText(Integer.toString(sfs));
			tfTxHi.setText(Integer.toString(settings.getProfiler().getSfHi().size()));
			tfTxMed.setText(Integer.toString(settings.getProfiler().getSfMed().size()));
			tfTxLo.setText(Integer.toString(settings.getProfiler().getSfLo().size()));
			tfTxNo.setText(Integer.toString(settings.getProfiler().getMolecules().length- sfs));
	
			tfTxHiFrac.setText(StringUtils.fprint(settings.getProfiler().getSfHi().size() * 100d / sfs, 2)+"%");
			tfTxMedFrac.setText(StringUtils.fprint(settings.getProfiler().getSfMed().size() * 100d / sfs, 2)+"%");
			tfTxLoFrac.setText(StringUtils.fprint(settings.getProfiler().getSfLo().size() * 100d / sfs, 2)+"%");
			tfTxNoFrac.setText(StringUtils.fprint((settings.getProfiler().getMolecules().length - sfs) * 100d / settings.getProfiler().getMolecules().length, 2)+"%");
	
			long molTot= 0, molHi= 0, molMed= 0, molLo= 0;
			for (int i = 0; i < settings.getProfiler().getMolecules().length; i++) {
				molTot+= settings.getProfiler().getMolecules()[i];
				if (settings.getProfiler().getSfHi().contains(settings.getProfiler().getIds()[i]))
					molHi+= settings.getProfiler().getMolecules()[i];
				else if (settings.getProfiler().getSfMed().contains(settings.getProfiler().getIds()[i]))
					molMed+= settings.getProfiler().getMolecules()[i];
				else if (settings.getProfiler().getSfLo().contains(settings.getProfiler().getIds()[i]))
					molLo+= settings.getProfiler().getMolecules()[i];
			}
			tfRNAtot.setText(Long.toString(molTot));
			tfRNAhi.setText(Long.toString(molHi));
			tfRNAmed.setText(Long.toString(molMed));
			tfRNAlo.setText(Long.toString(molLo));
			
			tfRNAhiFrac.setText(StringUtils.fprint(molHi * 100d / molTot, 2)+"%");
			tfRNAmedFrac.setText(StringUtils.fprint(molMed * 100d / molTot, 2)+"%");
			tfRNAloFrac.setText(StringUtils.fprint(molLo * 100d / molTot, 2)+"%");
			
			int[] allLen= new int[(int) molTot];	// TODO overflow int cast
			int cnt= 0;
			for (int i = 0; i < settings.getProfiler().getMolecules().length; i++) {
				for (int j = 0; j < settings.getProfiler().getMolecules()[i]; j++) {
					int len= settings.getProfiler().getLen()[i];
					// len= SimpleGelBinPlotterPanel.getEffectiveLength(len);	// f* slow
					allLen[cnt++]= len;
				}
			}
			panGel.paintOSI(allLen);
			allLen= null;
			System.gc();
			
			long[] mols= settings.getProfiler().getMolecules();
			Hashtable<CharSequence,Long> locMap= new Hashtable<CharSequence,Long>();
			long[] ranks= new long[settings.getProfiler().getSfHi().size()+ settings.getProfiler().getSfMed().size()+ settings.getProfiler().getSfLo().size()];
			cnt= 0;
			for (int i = 0; i < mols.length; i++) {
				if (mols[i]> 0) {
					if (cnt>= ranks.length)
						break;
					ranks[cnt++]= mols[i];
					if (locMap.containsKey(settings.getProfiler().getLocIDs()[i]))
						locMap.put(settings.getProfiler().getLocIDs()[i], locMap.get(settings.getProfiler().getLocIDs()[i])+ mols[i]);
					else
						locMap.put(settings.getProfiler().getLocIDs()[i], mols[i]);
				}
			}
			Arrays.sort(ranks);
			for (int i = 0; i < ranks.length/ 2; i++) {
				long h= ranks[i];
				ranks[i]= ranks[ranks.length- 1- i];
				ranks[ranks.length- 1- i]= h;
			}
			panRank1.paintOSI(ranks);

			HashMap<Long,Integer> mapRcountFreqLocus= new HashMap<Long,Integer>();
			Iterator<Long> iter= locMap.values().iterator();
			long max= 0;
			while (iter.hasNext()) {
				long nb= iter.next();
				if (mapRcountFreqLocus.containsKey(nb))
					mapRcountFreqLocus.put(nb, mapRcountFreqLocus.remove(nb)+1);
				else
					mapRcountFreqLocus.put(nb, 1);
				if (nb> max)
					max= nb;
			}			
			ranks= new long[(int) max+1];
			ranks[0]= settings.getProfiler().getIds().length- 
				(settings.getProfiler().getSfHi().size()+ 
						settings.getProfiler().getSfMed().size()+ 
						settings.getProfiler().getSfLo().size());
			iter= mapRcountFreqLocus.keySet().iterator();
			while(iter.hasNext()) {
				long nb= iter.next();
				ranks[(int) nb]= mapRcountFreqLocus.get(nb);
			}
			panRank2.paintOSI(ranks);
		}
		
		FluxSimulatorGUI.repaintEverywhere(this, false);
	}
	
	public boolean isFinished() {
		if (settings== null|| settings.getProfiler()== null)
			return false;
		return settings.getProfiler().isFinishedExpression();
	}

	boolean loadStats= false;
	public boolean getLoadStats() {
		return loadStats;
	}

	public void setLoadStats(boolean val) {
		loadStats= val;		
	}

	public void killResult() {
		if (settings!= null&& settings.getProfiler()!= null)
			settings.getProfiler().setMolecules(null);
		panRank1.offScrImg= null;
		panRank2.offScrImg= null;
		uiUpdateStats();
	}

	public boolean isStop() {		
		return settings== null|| settings.getProfiler()== null|| settings.getProfiler().isStop();
	}

	public boolean setStop(boolean stop) {
		if (stop) 
			return setStop();
		else {
			if (settings.getProfiler()== null)
				return false;
			return settings.getProfiler().setStop(stop);
		}
	}

}
