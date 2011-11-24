package fbi.genome.sequencing.rnaseq.simulation.gui;

import fbi.commons.ReadyOrNot;
import fbi.commons.StringUtils;
import fbi.commons.gui.SimpleBinPlotterPanel;
import fbi.commons.thread.StoppableRunnable;
import fbi.genome.model.commons.MyFile;
import fbi.genome.sequencing.rnaseq.simulation.FluxSimulatorSettings;
import fbi.genome.sequencing.rnaseq.simulation.Profiler;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.File;

public class AnnotationReaderGUI extends JPanel implements ReadyOrNot, StoppableRunnable {

	
	public static void main(String[] args) {
		AnnotationReaderGUI annotGUI= new AnnotationReaderGUI();
		annotGUI.set(FluxSimulatorSettings.createDefaults());
		
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
	FluxSimulatorSettings settings;
	JButton butLoadRef;
	public AnnotationReaderGUI() {
		setLayout(new BorderLayout());
		
		JPanel panInPan= new JPanel();
		panInPan.setBorder(BorderFactory.createTitledBorder("Input")); // or on lineborder black
		panInPan.setLayout(new BorderLayout(10,0));
		panInPan.add(new JLabel("Reference Annotation"), BorderLayout.WEST);
		refTextField= new JTextField();
		refTextField.setEditable(false);
		panInPan.add(refTextField, BorderLayout.CENTER);
		butLoadRef= new JButton("Change");
		butLoadRef.addActionListener(new ActionListener(){
			public void actionPerformed(ActionEvent e) {
			    JFileChooser chooser = new JFileChooser();
			    int returnVal = chooser.showOpenDialog(AnnotationReaderGUI.this);
			    if(returnVal == JFileChooser.APPROVE_OPTION) {
			    	String absPath= chooser.getSelectedFile().getAbsolutePath();
					File f= new File(absPath);
					AnnotationReaderGUI.this.settings.setRefFile(f);
					refTextField.setText(f.getName());
					refTextField.repaint();
			    }				
			}
		});
		butLoadRef.setEnabled(false);
		panInPan.add(butLoadRef, BorderLayout.EAST);
		add(panInPan, BorderLayout.NORTH);

		JPanel panLeft= new JPanel();
		panLeft.setLayout(new BoxLayout(panLeft,BoxLayout.Y_AXIS));
		panLeft.add(getParaPanel());
		panLeft.add(getStatsPanel());
		add(panLeft, BorderLayout.WEST);
		
		JPanel pan= new JPanel();
		pan.setLayout(new BoxLayout(pan, BoxLayout.Y_AXIS));
		pan.setBorder(BorderFactory.createTitledBorder("Histogram of Spliceform Lengths"));
		JPanel pan1= new JPanel();
		pan1.setLayout(new BorderLayout());
		lengthDistrPanel= new SimpleBinPlotterPanel();
		lengthDistrPanel.setTitleX("length [nt]");
		pan1.add(lengthDistrPanel, BorderLayout.CENTER);
		pan.add(pan1);
		JPanel panX= new JPanel();
		panX.setPreferredSize(new Dimension(100, 10));
		pan.add(panX);
		JPanel pan2= new JPanel();
		pan2.setLayout(new BorderLayout());
		lengthDistrPanel2= new SimpleBinPlotterPanel();
		lengthDistrPanel2.setTitleX("length [nt]");
		pan2.add(lengthDistrPanel2, BorderLayout.CENTER);
		pan.add(pan2);
		add(pan, BorderLayout.CENTER);		
	}

	JPanel paraPanel;
	private Component getParaPanel() {
		if (paraPanel == null) {
			paraPanel= new JPanel();
			paraPanel.setBorder(BorderFactory.createTitledBorder("Parameters"));
			paraPanel.setLayout(new GridBagLayout());
			GridBagConstraints c= new GridBagConstraints();
			c.anchor= GridBagConstraints.FIRST_LINE_START;
			c.fill = GridBagConstraints.NONE;
			c.gridwidth= 1; c.gridheight= 1;
			c.weightx= 0.1; c.weighty= 0.1;
			c.insets= new Insets(1,10,1,10);
			
			c.gridx= 0;
			c.gridy= 0;
			paraPanel.add(new JLabel("Load"), c);
			c.gridx= 1;
			c.gridy= 1;
			c.insets= new Insets(1,0,1,0);
			cbLoadCoding= new JCheckBox("Coding Transcripts");
			cbLoadCoding.setOpaque(false);
			cbLoadCoding.setEnabled(false);
			cbLoadCoding.addActionListener(new ActionListener(){
				public void actionPerformed(ActionEvent e) {
					settings.setLoadCoding(cbLoadCoding.isSelected());
				}
			});
			paraPanel.add(cbLoadCoding, c);
			c.gridy= 2;
			cbLoadNoncoding= new JCheckBox("Non-Coding Transcripts");
			cbLoadNoncoding.setOpaque(false);
			cbLoadNoncoding.setEnabled(false);
			cbLoadCoding.addActionListener(new ActionListener(){
				public void actionPerformed(ActionEvent e) {
					settings.setLoadNoncoding(cbLoadNoncoding.isSelected());
				}
			});
			paraPanel.add(cbLoadNoncoding, c);
		}
		return paraPanel;
	}
	
	JPanel statsPanel;
	JTextField tfFSize, tfTx, tfLoc, tfLenAvg, tfLenMed, tfASAvg, tfASMed, tfAS1Q, tfAS3Q, tfASSTD, tfLenMin, tfLenMax, tfLen1Q, tfLen3Q, tfLenSTD;
	JCheckBox cbLoadCoding, cbLoadNoncoding;
	private Component getStatsPanel() {
		if (statsPanel == null) {
			statsPanel = new JPanel();
			statsPanel.setBorder(BorderFactory.createTitledBorder("Reference Transcriptome"));
			statsPanel.setLayout(new GridBagLayout());
			GridBagConstraints c= new GridBagConstraints();
			c.anchor= GridBagConstraints.FIRST_LINE_START;
			c.fill = GridBagConstraints.NONE;
			c.insets= new Insets(1,10,1,10);
			c.gridwidth= 1; c.gridheight= 1;

			c.weighty= 0.01;
			c.gridx = 0; 
			c.gridy = 0;
			statsPanel.add(new JLabel("File Size"), c);
			c.gridy = 1;
			statsPanel.add(new JLabel("Spliceforms"), c);
			c.gridy = 2;
			statsPanel.add(new JLabel("Loci"), c);
			c.gridy = 3;
			statsPanel.add(new JLabel("Spliceforms per Locus"), c);
			c.gridy = 4;
			c.insets= new Insets(1,20,1,10);
			statsPanel.add(new JLabel("Average"), c);
			c.gridy = 5;
			c.insets= new Insets(1,20,1,10);
			statsPanel.add(new JLabel("Std Deviation"), c);
			c.gridy = 6;
			c.insets= new Insets(1,20,1,10);
			statsPanel.add(new JLabel("1st Quartile"), c);
			c.gridy = 7;
			c.insets= new Insets(1,20,1,10);
			statsPanel.add(new JLabel("Median"), c);
			c.gridy = 8;
			c.insets= new Insets(1,20,1,10);
			statsPanel.add(new JLabel("3rd Quartile"), c);
			c.gridy = 9;
			c.insets= new Insets(1,10,1,10);
			statsPanel.add(new JLabel("Spliceform Length"), c);
			c.gridy = 10;
			c.insets= new Insets(1,20,1,10);
			statsPanel.add(new JLabel("Minimum"), c);
			c.gridy = 11;
			c.insets= new Insets(1,20,1,10);
			statsPanel.add(new JLabel("Average"), c);
			c.gridy = 12;
			c.insets= new Insets(1,20,1,10);
			statsPanel.add(new JLabel("Std Deviation"), c);
			c.gridy = 13;
			c.insets= new Insets(1,20,1,10);
			statsPanel.add(new JLabel("Maximum"), c);
			c.gridy = 14;
			c.insets= new Insets(1,20,1,10);
			statsPanel.add(new JLabel("1st Quartile"), c);
			c.gridy = 15;
			c.insets= new Insets(1,20,1,10);
			statsPanel.add(new JLabel("Median"), c);
			c.gridy = 16;
			c.insets= new Insets(1,20,1,10);
			statsPanel.add(new JLabel("3rd Quartile"), c);
			
			int tfWidth= 8;
			c.gridx = 2; 
			c.insets= new Insets(1,10,1,10);
			tfFSize= new JTextField(tfWidth);
			tfFSize.setHorizontalAlignment(JTextField.RIGHT);
			tfFSize.setEditable(false);
			c.gridy = 0;
			statsPanel.add(tfFSize, c);
			tfTx= new JTextField(tfWidth);
			tfTx.setHorizontalAlignment(JTextField.RIGHT);
			tfTx.setEditable(false);
			c.gridy = 1;
			statsPanel.add(tfTx, c);
			tfLoc= new JTextField(tfWidth);
			tfLoc.setHorizontalAlignment(JTextField.RIGHT);
			tfLoc.setEditable(false);
			c.gridy = 2;
			statsPanel.add(tfLoc, c);
			
			tfASAvg= new JTextField(tfWidth);
			tfASAvg.setHorizontalAlignment(JTextField.RIGHT);
			tfASAvg.setEditable(false);
			c.gridy = 4;
			statsPanel.add(tfASAvg, c);
			tfASSTD= new JTextField(tfWidth);
			tfASSTD.setHorizontalAlignment(JTextField.RIGHT);
			tfASSTD.setEditable(false);
			c.gridy = 5;
			statsPanel.add(tfASSTD, c);
			tfAS1Q= new JTextField(tfWidth);
			tfAS1Q.setHorizontalAlignment(JTextField.RIGHT);
			tfAS1Q.setEditable(false);
			c.gridy = 6;
			statsPanel.add(tfAS1Q, c);
			tfASMed= new JTextField(tfWidth);
			tfASMed.setHorizontalAlignment(JTextField.RIGHT);
			tfASMed.setEditable(false);
			c.gridy = 7;
			statsPanel.add(tfASMed, c);
			tfAS3Q= new JTextField(tfWidth);
			tfAS3Q.setHorizontalAlignment(JTextField.RIGHT);
			tfAS3Q.setEditable(false);
			c.gridy = 8;
			statsPanel.add(tfAS3Q, c);
			
			tfLenMin= new JTextField(tfWidth);
			tfLenMin.setHorizontalAlignment(JTextField.RIGHT);
			tfLenMin.setEditable(false);
			c.gridy = 10;
			statsPanel.add(tfLenMin, c);
			c.gridy = 11;
			tfLenAvg= new JTextField(tfWidth);
			tfLenAvg.setHorizontalAlignment(JTextField.RIGHT);
			tfLenAvg.setEditable(false);
			statsPanel.add(tfLenAvg, c);
			tfLenSTD= new JTextField(tfWidth);
			tfLenSTD.setHorizontalAlignment(JTextField.RIGHT);
			tfLenSTD.setEditable(false);
			c.gridy = 12;
			statsPanel.add(tfLenSTD, c);
			tfLenMax= new JTextField(tfWidth);
			tfLenMax.setHorizontalAlignment(JTextField.RIGHT);
			tfLenMax.setEditable(false);
			c.gridy = 13;
			statsPanel.add(tfLenMax, c);
			tfLen1Q= new JTextField(tfWidth);
			tfLen1Q.setHorizontalAlignment(JTextField.RIGHT);
			tfLen1Q.setEditable(false);
			c.gridy = 14;
			statsPanel.add(tfLen1Q, c);
			tfLenMed= new JTextField(tfWidth);
			tfLenMed.setHorizontalAlignment(JTextField.RIGHT);
			tfLenMed.setEditable(false);
			c.gridy = 15;
			statsPanel.add(tfLenMed, c);
			tfLen3Q= new JTextField(tfWidth);
			tfLen3Q.setHorizontalAlignment(JTextField.RIGHT);
			tfLen3Q.setEditable(false);
			c.gridy = 16;
			statsPanel.add(tfLen3Q, c);
			
			c.insets= new Insets(1,10,1,10);
			c.fill = GridBagConstraints.VERTICAL;
			c.gridx = 0; c.gridy = 17;
			c.weighty= 1;
			c.gridwidth= 3;
			statsPanel.add(new JPanel(), c);

		}

		return statsPanel;
	}

	public boolean isReady() {
		//if (settings== null|| settings.getProfiler()== null) 
		if (settings== null|| settings.getParFile()== null|| (!settings.getParFile().exists())
				|| (!settings.getParFile().canRead()))
			return false;
		return true;	// settings.getProfiler().isReady()
	}
	
	public void set(Object o) {
		if (o== null|| !(o instanceof FluxSimulatorSettings))
			return;
		FluxSimulatorSettings settings= (FluxSimulatorSettings) o;
		this.settings= settings;
		if (settings== null) {
			uiUpdateStats();
			cbLoadCoding.setEnabled(false);
			cbLoadNoncoding.setEnabled(false);
		} else if (settings.getProfiler()== null) {
			settings.setProfiler(new Profiler(settings));
			cbLoadCoding.setEnabled(true);
			cbLoadNoncoding.setEnabled(true);
			cbLoadCoding.setSelected(settings.isLoadCoding());
			cbLoadNoncoding.setSelected(settings.isLoadNoncoding());
		} 
		uiUpdatePars();
		
	}

	public boolean setStop() {
		if (settings.getProfiler()== null)
			return false;
		return settings.getProfiler().setStop(true);
	}

	public void run() {		
				
		settings.getProfiler().setStop(false);
		if (loadStats) {
			loadStats();
			setLoadStats(false);
			return;
		}
		
		if (!settings.getProfiler().readAnnotation())
			return;
		
		if (!settings.getProfiler().isStop()) {
			settings.getProfiler().writeProfile();
			settings.save();
		}
		
		uiUpdateStats();
	}
	
	private void uiUpdatePars() {
		if (settings!= null&& settings.getRefFile()!= null) {
			refTextField.setText(settings.getRefFile().getAbsolutePath());
			//butLoadRef.setEnabled(true);
			
			cbLoadCoding.repaint();
//			cbLoadCoding.setEnabled(true);
			cbLoadCoding.revalidate();
			cbLoadNoncoding.revalidate();
		} else
			refTextField.setText(FluxSimulatorGUI.emptyString);
		
		FluxSimulatorGUI.repaintEverywhere(this, false);
	}
	
	private void uiUpdateStats() {
		
		if (settings== null|| settings.getProfiler()== null|| (!settings.getProfiler().isFinishedReadAnnotation())) {
			tfFSize.setText(FluxSimulatorGUI.emptyString);
			tfTx.setText(FluxSimulatorGUI.emptyString);
			tfLoc.setText(FluxSimulatorGUI.emptyString);
			tfASAvg.setText(FluxSimulatorGUI.emptyString);
			tfASMed.setText(FluxSimulatorGUI.emptyString);
			tfAS1Q.setText(FluxSimulatorGUI.emptyString);
			tfAS3Q.setText(FluxSimulatorGUI.emptyString);
			tfASSTD.setText(FluxSimulatorGUI.emptyString);
			tfLenMin.setText(FluxSimulatorGUI.emptyString);
			tfLenMax.setText(FluxSimulatorGUI.emptyString);
			tfLenMed.setText(FluxSimulatorGUI.emptyString);
			tfLenSTD.setText(FluxSimulatorGUI.emptyString);
			tfLen1Q.setText(FluxSimulatorGUI.emptyString);
			tfLen3Q.setText(FluxSimulatorGUI.emptyString);
			tfLenAvg.setText(FluxSimulatorGUI.emptyString);
			lengthDistrPanel.offScrImg= null;
			lengthDistrPanel2.offScrImg= null;

		} else {
			tfFSize.setText(MyFile.humanReadableSize(settings.getRefFile().length()));
			tfTx.setText(Integer.toString(settings.getProfiler().getLen().length));
			tfLoc.setText(Integer.toString(settings.getProfiler().getCntLoci()));	//MyFormatter.american(
			tfASAvg.setText(StringUtils.fprint(settings.getProfiler().getTxLocAvg(), 2));
	//		String s;
	//		int p;
	//		s= MyFormatter.fprint(profiler.getTxLocMed(), 1);
	//		p= s.indexOf(".");
	//		if (p>= 0)
	//			s= MyFormatter.american(Integer.parseInt(s.substring(0,p)))+ "."+ s.substring(p+1);
			tfASMed.setText(StringUtils.fprint(settings.getProfiler().getTxLocMed(), 2));
			tfAS1Q.setText(StringUtils.fprint(settings.getProfiler().getTxLoc1Q(), 2));
			tfAS3Q.setText(StringUtils.fprint(settings.getProfiler().getTxLoc3Q(), 2));
			tfASSTD.setText(StringUtils.fprint(settings.getProfiler().getTxLocSTD(), 2));
			tfLenMin.setText(StringUtils.fprint(settings.getProfiler().getLenMin(), 2));
			tfLenMax.setText(StringUtils.fprint(settings.getProfiler().getLenMax(), 2));
			tfLenMed.setText(StringUtils.fprint(settings.getProfiler().getLenMed(), 2));
			tfLenSTD.setText(StringUtils.fprint(settings.getProfiler().getLenSTD(), 2));
			tfLen1Q.setText(StringUtils.fprint(settings.getProfiler().getLen1Q(), 2));
			tfLen3Q.setText(StringUtils.fprint(settings.getProfiler().getLen3Q(), 2));
			tfLenAvg.setText(StringUtils.fprint(settings.getProfiler().getLenAvg(), 2));
			lengthDistrPanel.paintOSI(settings.getProfiler().getLen());
			lengthDistrPanel.repaint();
			lengthDistrPanel2.setThrUp(settings.getProfiler().getLen3Q());
			lengthDistrPanel2.paintOSI(settings.getProfiler().getLen());
			lengthDistrPanel2.repaint();
		}
		
		FluxSimulatorGUI.repaintEverywhere(this, false);
	}
	
	public boolean loadStats() {
		if (settings.getProfiler().getIds()== null|| settings.getProfiler().getLen()== null) {
			settings.getProfiler().setStop(false);
			if (!settings.getProfiler().loadStats())
				return false;
		}
		if (settings.getProfiler().isFinishedReadAnnotation()) {
			uiUpdateStats();
		}
		return true;
	}
	
	public boolean isFinished() {
		if (settings== null|| settings.getProfiler()== null)
			return false;
		return settings.getProfiler().isFinishedReadAnnotation();
	}

	public boolean getLoadStats() {
		return loadStats;
	}

	boolean loadStats= false;
	public void setLoadStats(boolean val) {
		loadStats= val;
	}

	public void killResult() {
		if (settings!= null&& settings.getProfiler()!= null) {
			settings.setProfiler(null);
		}
		if (settings!= null&& settings.getProFile()!= null&& settings.getProFile().exists()) {
			settings.getProFile().delete();
		}
		set(settings);
		lengthDistrPanel.offScrImg= null;
		lengthDistrPanel2.offScrImg= null;
		//uiUpdate();
	}

	public boolean isStop() {		
		return settings== null|| settings.getProfiler()== null|| settings.getProfiler().isStop();
	}

	public boolean setStop(boolean stop) {
		if (settings.getProfiler()== null)
			return false;
		return settings.getProfiler().setStop(stop);
	}

}
