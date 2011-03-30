package fbi.genome.sequencing.rnaseq.reconstruction.gui;

import java.awt.BorderLayout;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.io.File;
import java.util.Enumeration;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;

import fbi.genome.model.constants.Constants;
import fbi.genome.sequencing.rnaseq.reconstruction.FluxCapacitor;

import javax.swing.BorderFactory;
import javax.swing.DefaultListModel;
import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextField;
import javax.swing.UIManager;

import commons.ReadyOrNot;
import commons.thread.StoppableRunnable;

public class DecompositionGUI extends JPanel implements ReadyOrNot, StoppableRunnable {

	public static void main(String[] args) {
		
		FluxCapacitorGUI.setColors(UIManager.getLookAndFeel().getName());        
		
		JFrame frame= new JFrame("Test");
		DecompositionGUI gui= new DecompositionGUI();
		frame.getContentPane().add(gui);
		File f= new File("c:\\workspace\\Genome\\resources\\formats\\CME_W1_CI.8-PE_lp_mod.zip");
				//"c:\\Documents and Settings\\micha\\Local Settings\\Temp\\capacitor.1001101446450958.hg18_knownGenes_fromUCSC081019_sorted__seq-150_200_lp.null");
		gui.setFileLP(f);
		gui.updateFileList();
		
		frame.pack();
		frame.setVisible(true);
	}
	
	JList listLoci;
	LPVizualizer panLPviz;
	JTextField tfIn;
	JButton butChangeIn;
	String[] fileList= null;
	
	public DecompositionGUI() {
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
		butChangeIn.setEnabled(true);
		panInPan.add(butChangeIn, c);
		
		JScrollPane scroller= new JScrollPane(getListLoci());
		scroller.setBorder(BorderFactory.createTitledBorder("Loci"));
		add(scroller, BorderLayout.WEST);
		
		JPanel panHelp= new JPanel(new BorderLayout());
		add(panHelp, BorderLayout.CENTER);
		panHelp.add(panInPan, BorderLayout.NORTH);		
		panHelp.add(getLPvizPanel(), BorderLayout.CENTER);
	}
	
	private JList getListLoci() {
		if (listLoci == null) {
			DefaultListModel model= new DefaultListModel();
			for (int i = 0; fileList!= null&& i< fileList.length; i++) 
				model.addElement(fileList[i]);
			listLoci = new JList(model);
			
			listLoci.addMouseListener(new MouseAdapter(){
				@Override
				public void mouseClicked(MouseEvent e) {
					if(e.getClickCount() == 2){
					     int index = listLoci.locationToIndex(e.getPoint());
					     if (index< 0|| index>= listLoci.getModel().getSize())
					    	 return;
					     Object item = listLoci.getModel().getElementAt(index);;
					     //System.err.println(item);
					     
					     getLPvizPanel().parseEntry(item.toString());
					     getLPvizPanel().revalidate();
					}

					super.mouseClicked(e);
				}
			});
		}

		return listLoci;
	}
	
	private LPVizualizer getLPvizPanel() {
		if (panLPviz == null) {
			panLPviz = new LPVizualizer();
			panLPviz.setBorder(BorderFactory.createTitledBorder("Locus Geometry"));
		}

		return panLPviz;
	}
	
	public boolean getLoadStats() {
		// TODO Auto-generated method stub
		return false;
	}

	public boolean isFinished() {
		// TODO Auto-generated method stub
		return false;
	}

	public boolean isReady() {
		return FluxCapacitorGUI.singleton.capacitor.isInputReady();
				//&& FluxCapacitorGUI.singleton.capacitor.getFunc()!= null;
	}

	public void killResult() {
		// TODO Auto-generated method stub
		
	}

	public boolean loadStats() {
		// TODO Auto-generated method stub
		return false;
	}

	public void set(Object o) {
		uiUpdate();
	}

	public void setLoadStats(boolean val) {
		// TODO Auto-generated method stub
		
	}

	public boolean isStop() {
		// TODO Auto-generated method stub
		return false;
	}

	public boolean setStop() {
		// TODO Auto-generated method stub
		return false;
	}

	public boolean setStop(boolean stop) {
		// TODO Auto-generated method stub
		return false;
	}

	File lpFile;
	private void updateFileList() {
		
		if (lpFile== null)
			return;
		
		if (lpFile.isDirectory())
			fileList= lpFile.list();
		else {
			try {
				ZipFile zf= new ZipFile(lpFile);
				fileList= new String[zf.size()];
				Enumeration enu= zf.entries();
				for (int i = 0; i < fileList.length; i++) 
					fileList[i]= ((ZipEntry) enu.nextElement()).getName();
				zf.close();
			} catch (Exception e) { // IO.., Zip..
				;	// :)
			} 
			
		}
		
		if (fileList== null|| fileList.length== getListLoci().getModel().getSize())
			return;
		
		DefaultListModel model= ((DefaultListModel) getListLoci().getModel());
		for (int i = 0; i < fileList.length; i++) {
			if (!model.contains(fileList[i]))
				model.addElement(fileList[i]);
		}
		getListLoci().revalidate();
	}
	
	public void setFileLP(File filelp) {
		this.lpFile= filelp;
		getLPvizPanel().setBase(filelp);
	}
	
	private void uiUpdate() {
		File f= FluxCapacitorGUI.singleton.capacitor.getFileLP();
		if (f== null) {
			tfIn.setText("");
			tfIn.setToolTipText(null);
		} else {
			setFileLP(f);
			tfIn.setText(f.getName());
			tfIn.setToolTipText(f.getPath());
		}
		
		FluxCapacitorGUI.repaintEverywhere(this, false);	
	}
	public void run() {
		
		if (lpFile== null)
			return;
		
		final FluxCapacitor capacitor= FluxCapacitorGUI.singleton.capacitor;
		final Thread baseT= Thread.currentThread();
		
		// polling thread
		Thread poll= new Thread("polling"){
			@Override
			public void run() {
				while (baseT.isAlive()) {
					updateFileList();
					try {
						Thread.currentThread().sleep(5000);
					} catch (InterruptedException e) {
						; // :)
					}
				}
				updateFileList();
			}
		};
		poll.start();
		
		
		try {
			capacitor.explore(FluxCapacitor.MODE_RECONSTRUCT);
		} catch (Throwable e) {
			FluxCapacitorGUI.singleton.dialogError(e);
			if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
				System.err.println("[FATAL] Error occured during scanning\n\t"+ e.getMessage());
		}

	}

	
}
