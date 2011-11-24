package fbi.genome.sequencing.rnaseq.simulation.gui;

import fbi.commons.Dialogable;
import fbi.commons.Log;
import fbi.commons.ReadyOrNot;
import fbi.commons.file.FileHelper;
import fbi.commons.gui.MyProgressBar;
import fbi.commons.gui.SimpleBinPlotterPanel;
import fbi.commons.system.SystemInspector;
import fbi.commons.thread.StoppableRunnable;
import fbi.genome.model.constants.Constants;
import fbi.genome.sequencing.rnaseq.simulation.FluxSimulator;
import fbi.genome.sequencing.rnaseq.simulation.FluxSimulatorSettings;

import javax.swing.*;
import javax.swing.UIManager.LookAndFeelInfo;
import javax.swing.border.LineBorder;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.event.DocumentEvent;
import javax.swing.event.DocumentListener;
import javax.swing.filechooser.FileFilter;
import javax.swing.plaf.BorderUIResource;
import javax.swing.plaf.ColorUIResource;
import java.awt.*;
import java.awt.event.*;
import java.awt.geom.Rectangle2D;
import java.awt.image.FilteredImageSource;
import java.awt.image.ImageFilter;
import java.awt.image.ImageProducer;
import java.awt.image.RGBImageFilter;
import java.io.*;
import java.util.ArrayList;

public class FluxSimulatorGUI extends JFrame {

	public static FluxSimulatorGUI singleton; 
	
	public class GlassPane extends JPanel implements MouseListener {
	    //RootPaneContainer interface
	    private RootPaneContainer m_rootPane = null;
	    
	    //Stores the previous Glass Pane Component
	    private Component m_prevGlassPane = null;
	    private boolean m_handleMouseEvents = false;
	    private boolean m_drawing = false;    

	    /**
	    * Constructor
	    *
	    * @exception
	    */
	    public GlassPane()
	    {
	    }
	    
	    //Mouse Events
	    public void mousePressed(MouseEvent e) 
	    {
	        //sound beep
	        Toolkit.getDefaultToolkit().beep(); 
	    }
	    
	    public void mouseReleased(MouseEvent e) 
	    {
	    }
	    
	    public void mouseClicked(MouseEvent e) 
	    {
	    	setHandleMouseEvents(false);
	    	setDrawing(false);
	    }
	    
	    public void mouseEntered(MouseEvent e) 
	    {
	    }
	    
	    public void mouseExited(MouseEvent e) 
	    {
	    }
	    
	    /**
	    * Set the glassPane
	    */
	    public void setGlassPane(RootPaneContainer rootPane)
	    {
	        m_rootPane = rootPane;
	        
	        //store the current glass pane
	        m_prevGlassPane = m_rootPane.getGlassPane();
	        
	        //set this as new glass pane
	        m_rootPane.setGlassPane(this);
	        
	        //set opaque to false, i.e. make transparent
	        setOpaque(false);
	    }
	    
	    /**
	    * remove the glassPane
	    */
	    public void removeGlassPane()
	    {
	        //set the glass pane visible false
	        setVisible(false);
	        
	        //reset the previous glass pane
	        m_rootPane.setGlassPane(m_prevGlassPane);
	    }
	    
		@Override
		public Dimension getPreferredSize() {
			if (imgAbout!= null) 
				return new Dimension(imgAbout.getWidth(this), imgAbout.getHeight(this));
			return super.getPreferredSize();
		}
	    
	    /**
	    * Set the handleMoveEvents
	    * This will add the mouseListener and all the event will be
	    * trapped my the glass pane
	    */
	    public void setHandleMouseEvents(boolean handleMouseEvents)
	    {
	        if (m_handleMouseEvents == handleMouseEvents)
	        {
	            //ignore if the state is same
	            return;
	        }
	        
	        m_handleMouseEvents = handleMouseEvents;
	        
	        if (m_handleMouseEvents)
	        {
	            //add this as mouse event listener and trap all the
	            //mouse events
	            addMouseListener(this);
	        }
	        else
	        {
	            //remove mouse events listener
	            removeMouseListener(this);
	        }
	        
	        //This is important otherwise the glass pane will no catch events
	        setVisible(true);
	    } 
	    
	    /**
	     * setDrawing
	     *
	     */
	     public void setDrawing(boolean drawing)
	     {
	         m_drawing = drawing;
	         
	         //This is important otherwise the glass pane will not be visible
	         setVisible(true);
	         
	         //Call repaint
	         repaint();
	     }
	     
	     /**
	     * get the drawing state
	     */
	     public boolean getDrawing()
	     {
	         return m_drawing;
	     }
	     
	     
	     @Override
	    protected void paintComponent(Graphics g) {
	    	
	        if ((!m_drawing)|| imgAbout== null)
	        	return;
	        
	        int diaX= FluxSimulatorGUI.this.getWidth()/2- imgAbout.getWidth(this)/2,
	        	diaY= FluxSimulatorGUI.this.getHeight()/2- imgAbout.getHeight(this)/2;
	        g.drawImage(imgAbout, diaX, diaY, this);
	        
	        String msg= "build "+ FluxSimulator.FLUX_VERSION;
	        g.setFont(new Font("Arial", Font.BOLD, 18));
	        g.setColor(COL_YEL_DARK);
	        Rectangle2D rect= g.getFontMetrics().getStringBounds(msg, g);
	        int strX= (int) (diaX+ (imgAbout.getWidth(this)/ 2)- (rect.getWidth()/ 2)+ 18),
	        	strY= diaY+ imgAbout.getHeight(this)- 1; 
	        g.drawString(msg, strX, strY);
	        g.setColor(COL_BG_YEL);
	        g.drawString(msg, strX- 1, strY- 1);
	    }
	}
	
	public static void dialogError(Throwable e) {
		JOptionPane.showMessageDialog(singleton, e.getMessage(), "Error Occured", JOptionPane.ERROR_MESSAGE);
	}
	
	public static void repaintEverywhere(final Component c, boolean now) {
		Runnable runner= new Runnable() {
			public void run() {
				c.repaint();
			}
		};
		if (EventQueue.isDispatchThread())
			runner.run();
		else {
			if (now)
				try {
					SwingUtilities.invokeAndWait(runner);
				} catch (Exception e) {
					; // :)
				}
			else
				SwingUtilities.invokeLater(runner);
		}
	}
	
	public static final String emptyString= " ";
	
	
	class ParFilter extends FileFilter {
		
        public boolean accept(File f) {
            return f.isDirectory() || f.getName().toLowerCase().endsWith(".par");
        }
        
        public String getDescription() {
            return "Parameter Files (*.par)";
        }
    }

	class Dialoguer implements Dialogable {

		public boolean checkOverwrite(String s) {
			String[] options= new String[] {"Go ahead, kill it!", "Maybe I check, just in case.."};
			int n = JOptionPane.showOptionDialog(FluxSimulatorGUI.this,
					s,
	                "Confirm overwrite",
	                JOptionPane.YES_NO_OPTION,
	                JOptionPane.WARNING_MESSAGE,
	                null,
	                options,
	                options[1]);
			if (n== JOptionPane.YES_OPTION)
				return true;
			return false;
		}

		public void showError(String s) {
			int n = JOptionPane.showOptionDialog(FluxSimulatorGUI.this,
					s,
	                "Error",
	                JOptionPane.NO_OPTION,
	                JOptionPane.ERROR_MESSAGE,
	                null,
	                null,
	                null);
		}

		public void showInfo(String s) {
			int n = JOptionPane.showOptionDialog(FluxSimulatorGUI.this,
					s,
	                "Information",
	                JOptionPane.NO_OPTION,
	                JOptionPane.INFORMATION_MESSAGE,
	                null,
	                null,
	                null);
		}

		public void showWarning(String s) {
			int n = JOptionPane.showOptionDialog(FluxSimulatorGUI.this,
					s,
	                "Warning",
	                JOptionPane.NO_OPTION,
	                JOptionPane.WARNING_MESSAGE,
	                null,
	                null,
	                null);
		}
		
	}
	
	class BlockingThread extends Thread {
		StoppableRunnable target;
		BlockingThread q;
		int forceIdx= -1;
		
		public BlockingThread(StoppableRunnable target, BlockingThread queue) {
			super(target, "Blocking Thread");
			this.target= target;
			this.q= queue;
		}
		
		@Override
		public void run() {
			while (q!= null&& q.isAlive())
				try {
					q.join();
				} catch (InterruptedException e1) {
					if (target.isStop()) {
						q.cancel();
						q.interrupt();
					}
				}
			if (target.isStop()&& loadActive)
				return;
			if (q!= null)
				try {	// eye
					sleep(1000);
				} catch (Exception e) {
					; // :)
				}
				
			butNew.setEnabled(false);
			butLoad.setEnabled(false);
			butCopy.setEnabled(false);
			butRun.setEnabled(false);
			boolean tabWasActive= getTabPane().isEnabled();
			getTabPane().setEnabled(false);
			butCancel.setEnabled(true);
						
			if (forceIdx>= 0) {
//				for (int i = 0; i < tabPane.getTabCount(); i++) 
//					if (i!= forceIdx)
//						tabPane.setEnabledAt(i,false);
				getTabPane().setSelectedIndex(forceIdx);
			}
			// needed for resetOSI
			try {
				SwingUtilities.invokeLater(new Runnable(){
					public void run() {
						getTabPane().paintImmediately(tabPane.getBounds());
						labStatus.setString(emptyString);
						labStatus.paintImmediately(labStatus.getBounds());
					}
				});
			} catch (Exception e) {
				; // :)
			}
			
			super.run();
			
			if (tabWasActive) {
				getTabPane().setEnabled(true);
				butNew.setEnabled(true);
				butLoad.setEnabled(true);
				butCopy.setEnabled(true);
			}
			if (!loadActive) {
				butRun.setEnabled(true);
				butCancel.setEnabled(false);
			}
//			if (forceIdx>= 0) {
//				for (int i = 0; i < tabPane.getTabCount(); i++) 
//					tabPane.setEnabledAt(i,true);
//			}

			if (target.isStop())
				try {
					SwingUtilities.invokeLater(new Runnable(){
						public void run() {
							FluxSimulatorGUI.this.uiUpdate("aborted");
							getTabPane().paintImmediately(tabPane.getBounds());
						}
					});
				} catch (Exception e) {
					; // :)
				}
			else {
				for (int i = 0; i < getTabPane().getTabCount(); i++) {
					((ReadyOrNot) getTabPane().getComponentAt(i)).set(settings);	// update
				}
				try {
					SwingUtilities.invokeLater(new Runnable(){
						public void run() {
							FluxSimulatorGUI.this.uiUpdate("done");
							getTabPane().paintImmediately(tabPane.getBounds());
						}
					});
				} catch (Exception e) {
					; // :)
				}
			}
		}

		public final void cancel() {
			target.setStop();
		}

		public int getForceIdx() {
			return forceIdx;
		}

		public void setForceIdx(int forceIdx) {
			this.forceIdx = forceIdx;
		}
	}
	
	class GTFFilter extends FileFilter {
		
	    public boolean accept(File f) {
	        return f.isDirectory() || f.getName().toLowerCase().endsWith(".gtf") || f.getName().toLowerCase().endsWith(".gff");
	    }
	    
	    public String getDescription() {
	        return "Gene Transfer Format (*.gtf, *.gff)";
	    }
	}

	static class AsyncImageLoader implements Runnable {
	    String  imageFileName;
	    Thread  loaderThread;
	    Frame  parentFrame;
	
	    public AsyncImageLoader(Frame parent, String fileName){
	        parentFrame   = parent;
	        imageFileName = fileName;
	        loaderThread  = new Thread(this);
	        loaderThread.start();
	    }
	    
	    public void run() {
	        // start loading the image
	        splashImage = Toolkit.getDefaultToolkit().createImage(ClassLoader.getSystemResource(imageFileName));
	        
	        // wait for image to be loaded
	        MediaTracker tracker = new MediaTracker(parentFrame);
	        tracker.addImage(splashImage,0);
	        try {
	            tracker.waitForID(0);
	        }
	        catch(InterruptedException e){
	            e.printStackTrace();
	            System.exit(1);
	        }
	                
	        // check to ensure the image loaded okay. It would be nice to give
	        // a more specific error message here, but the Image load/observe
	        // API doesn't give us further details.
	        if(tracker.isErrorID(0)){
	            System.err.println("splashloader: error loading image \"" +
	                               imageFileName +
	                               "\"");
	
	            // this isn't a fatal error - the target class should be able
	            // to load.
	            return;
	        }
	
	        // resize frame to match size of image, and keep frame at centre of screen
	        Dimension screenSize = Toolkit.getDefaultToolkit().getScreenSize();
	        parentFrame.setBounds((screenSize.width-splashImage.getWidth(null))/2,
	                  (screenSize.height-splashImage.getHeight(null))/2,
	                  splashImage.getWidth(null),
	                  splashImage.getHeight(null));
	
	        // signal a redraw, so the image can be displayed
	        //imageLoaded = true;
	        parentFrame.setVisible(true);
	        parentFrame.repaint();
	    }
	} /* end of static inner class AsyncImageLoader */

	public final static byte MODE_NONE= -1, MODE_ANN= 0, MODE_EXP= 1, MODE_RTFRAG= 2, MODE_SEQ= 3;
	byte mode= MODE_NONE;
	
	public static void main(String[] args) {
		createGUI();
	}
	
	public static FluxSimulatorGUI createGUI() {
		Constants.verboseLevel= Constants.VERBOSE_SHUTUP;
		
		Frame f = new Frame() {
			Image image= null;
			@Override
			public void paint(Graphics g) {
				if (image== null) {
					Rectangle rect = new Rectangle(getBounds()); 
					Robot robot;
					Image image2= null;
					try {
						robot = new Robot();
					} catch (AWTException e) {
						return; // :)
					}  
					image2 = robot.createScreenCapture(rect);
					
					ImageFilter filter = new RGBImageFilter() {
					      // the color we are looking for... Alpha bits are set to opaque
					      public int markerRGB = Color.white.getRGB() | 0xFF000000;

					      public final int filterRGB(int x, int y, int rgb) {
					        if ( ( rgb | 0xFF000000 ) == markerRGB ) {
					          // Mark the alpha bits as zero - transparent
					          return 0x00FFFFFF & rgb;
					          }
					        else {
					          // nothing to do
					          return rgb;
					          }
					        }
					      }; 

					 ImageProducer ip = new FilteredImageSource(image2.getSource(), filter);
					 image= Toolkit.getDefaultToolkit().createImage(ip);

				}
				if (splashImage== null)
					return;
				if (image!= null)
					g.drawImage(image, 0, 0, null);
				g.drawImage(splashImage,0,0,null);
				
		       String msg= "build "+ FluxSimulator.FLUX_VERSION;
		        g.setFont(new Font("Arial", Font.BOLD, 18));
		        g.setColor(Color.black);
		        Rectangle2D rect= g.getFontMetrics().getStringBounds(msg, g);
		        int strX= (int) ((splashImage.getWidth(this)/ 2)- (rect.getWidth()/ 2)), 
		        	strY= splashImage.getHeight(this)- 3; 
		        g.drawString(msg, strX, strY);
			}
		};
        f.setUndecorated(true);
        if (SystemInspector.getOSgroup()== SystemInspector.OS_GROUP_MACOSX)
        	f.setBackground(new Color(1f,1f,1f,0.5f));	// transparent on mac
        AsyncImageLoader loader= new AsyncImageLoader(f, "pics/splash.png");

		LookAndFeelInfo[] lafs= UIManager.getInstalledLookAndFeels();			
		for (int i = 0; i < lafs.length; i++) {
			try {
				UIManager.setLookAndFeel(lafs[i].getClassName());
			} catch (Exception e) {
				; // :)
			}
		}
        if (!UIManager.getLookAndFeel().getName().contains("Metal")) {
			for (int i = 0; i < lafs.length; i++) {
				if (lafs[i].getName().contains("Metal")) {
					try {
						UIManager.setLookAndFeel(lafs[i].getClassName());
					} catch (Exception e) {
						; // :)
					}
					break;
				}
			}
        }
        setColors(UIManager.getLookAndFeel().getName());        
		UIManager.put("FileChooser.readOnly", Boolean.TRUE);
		
		
		final FluxSimulatorGUI gui= new FluxSimulatorGUI();
		try {
			Thread.sleep(2000);
		} catch (InterruptedException e) {
			; // :)
		}
		f.setVisible(false);
		gui.setVisible(true);

		// compromise, try here another instance for preload
		JFileChooser myChooser= new JFileChooser();	// not in the EDT !!!
		
		return gui;
	}
	
	public static final String TITLE= "Flux Simulator";
	public FluxSimulatorGUI() {
		super(TITLE+ " "+ FluxSimulator.FLUX_VERSION);

		singleton= this;
		setSize(new Dimension(1024,746));
		setLocation(getToolkit().getScreenSize().width/2- 512, getToolkit().getScreenSize().height/ 2- 373);
		setResizable(false);
		initUI();
		addWindowListener(new WindowAdapter() {
			@Override
			public void windowClosing(WindowEvent e) {
				if (settings!= null)
					settings.save();
				exit();
			}
//			@Override
//			public void windowStateChanged(WindowEvent e) {
//				super.windowStateChanged(e);
//				System.err.println("state changed");
//			}
		});
	}
	
	private void exit() {
		while (blockThread!= null&& blockThread.isAlive()) {
			killBlockThread();
		}
		System.exit(0);
	}
	
	public final static String NAME_ANN= "Annotation", NAME_EXP= "Expression", NAME_RTFRAG= "Library", NAME_SEQ= "Sequencing";
	JTabbedPane tabPane;
	MyProgressBar labStatus;
	public final static Color COL_BG_YEL= new Color(220, 199, 152), 
		COL_BG_GREEN= new Color(154, 158, 111), 
		COL_SEL_GREEN= new Color(110,121,87), 
		COL_YEL_SEMI_DARK= new Color(192,154,70),
		COL_YEL_DARK= new Color(152,122,54),
		COL_SHADOW_DARK= new Color(85,85,53),
		COL_SHADOW_SEMI_DARK= new Color(96,95,59);
	
	private static UIDefaults uiDefs;
	// http://home.tiscali.nl/~bmc88/java/sbook/061.html
	private static void setColors(String lafName) {
		if (uiDefs== null) {
			uiDefs= (UIDefaults) UIManager.getLookAndFeel().getDefaults().clone();
//			Object[] kk= uiDefs.keySet().toArray();
//			for (int i = 0; i < kk.length; i++) 
//				System.out.println(kk[i]+"\t->\t"+uiDefs.get(kk[i]));
		}
		
		if (lafName.contains("Metal")) { 
			SimpleBinPlotterPanel.defaultPanBg= COL_BG_YEL;
			
			UIManager.put("Panel.background", new ColorUIResource(COL_BG_GREEN));		
			UIManager.put("TabbedPane.selected", new ColorUIResource(COL_BG_GREEN));	// COL_SEL_GREEN, but cannot set text to white		
			UIManager.put("TabbedPane.contentAreaColor", new ColorUIResource(COL_BG_GREEN));	// tab-hover + border um components ++ non-inited bkg
			UIManager.put("TabbedPane.borderHightlightColor", new ColorUIResource(COL_SEL_GREEN));	// border um tab und component
			UIManager.put("TabbedPane.focus", new ColorUIResource(COL_SHADOW_SEMI_DARK));	// focus Kasten
			UIManager.put("TabbedPane.light", new ColorUIResource(COL_SHADOW_SEMI_DARK));	// li + oben light fx in tab			
			UIManager.put("TabbedPane.unselectedBackground", new ColorUIResource(COL_BG_YEL));
			
			//UIManager.put("TabbedPane.selectHighlight", COL_SHADOW_DARK);
			//UIManager.put("TabbedPane.darkShadow", COL_SHADOW_DARK);
			//UIManager.put("TabbedPane.shadow", COL_SHADOW_DARK);
			//UIManager.put("TabbedPane.borderHightlightColor", COL_SHADOW_DARK);
			//UIManager.put("TabbedPane.tabAreaBackground", COL_SHADOW_DARK);
			
			UIManager.put("ProgressBar.foreground", new ColorUIResource(COL_SEL_GREEN));
			UIManager.put("ProgressBar.background", new ColorUIResource(COL_BG_YEL));
			UIManager.put("ProgressBar.border", new ColorUIResource(COL_SHADOW_DARK));
			UIManager.put("ProgressBar.border", new BorderUIResource(new LineBorder(COL_SEL_GREEN)));
			//UIManager.put("ProgressBar.selectionBackground", new ColorUIResource(COL_SEL_GREEN));
			
			UIManager.put("TextField.background", new ColorUIResource(COL_YEL_SEMI_DARK));
			UIManager.put("TextField.inactiveBackground", new ColorUIResource(COL_YEL_DARK));
			UIManager.put("TextField.selectionBackground", new ColorUIResource(COL_SEL_GREEN));
			UIManager.put("TextField.selectionForeground", new ColorUIResource(Color.white));
			
//			UIManager.put("RadioButton.background", new ColorUIResource(COL_BG_YEL));	
			UIManager.put("RadioButton.disabledText", uiDefs.get("Label.foreground"));	
//			UIManager.put("CheckBox.background", new ColorUIResource(COL_BG_YEL));	
			UIManager.put("CheckBox.disabledText", uiDefs.get("Label.foreground"));	
			
			ArrayList a= new ArrayList();
			a.add(new Float(1.0));	// 1.0
			a.add(new Float(0.0));	// 0.0
			a.add(new ColorUIResource(COL_BG_YEL));	// [r=255,g=255,b=255]
			a.add(new ColorUIResource(COL_YEL_SEMI_DARK));	// [r=218,g=218,b=218]
			a.add(new ColorUIResource(COL_YEL_DARK));	// [r=218,g=218,b=218]
			UIManager.put("MenuBar.gradient", a);

			a= new ArrayList();
			a.add(new Float(0.3));	// 0.3
			a.add(new Float(0.0));	// 0.0
			a.add(new ColorUIResource(COL_BG_YEL));	// [r=221,g=232,b=243]
			a.add(new ColorUIResource(Color.white));	// [r=255,g=255,b=255]
			a.add(new ColorUIResource(COL_YEL_SEMI_DARK));	// [r=184,g=207,b=229]
			UIManager.put("Button.gradient", a);
			//UIManager.put("RadioButton.gradient", a);	// the same
			UIManager.put("Button.background", new ColorUIResource(COL_YEL_SEMI_DARK));
			//UIManager.put("Button.foreground", new ColorUIResource(COL_YEL_SEMI_DARK));
			UIManager.put("Button.disabledToolBarBorderBackground", new ColorUIResource(COL_YEL_DARK));			
			//UIManager.put("Button.border", new ColorUIResource(COL_YEL_DARK));
			UIManager.put("Button.disabledText", uiDefs.get("Button.foreground"));
			UIManager.put("Button.shadow", new ColorUIResource(Color.blue));			
			UIManager.put("Button.highlight", new ColorUIResource(Color.red));			
			UIManager.put("Button.darkShadow", new ColorUIResource(Color.green));			
			
			
			UIManager.put("ToolBar.dockingForeground", new ColorUIResource(COL_SEL_GREEN));	// higlight rect when hovering docking area
			UIManager.put("ToolBar.dockingBackground", new ColorUIResource(COL_BG_GREEN));
			UIManager.put("ToolBar.floatingForeground", new ColorUIResource(COL_YEL_DARK));
			UIManager.put("ToolBar.floatingBackground", new ColorUIResource(COL_BG_YEL));

			UIManager.put("MenuItem.background", new ColorUIResource(COL_BG_YEL));
			UIManager.put("MenuItem.selectionForeground", new ColorUIResource(Color.white));
			UIManager.put("MenuItem.selectionBackground", new ColorUIResource(COL_SEL_GREEN));
			
			//UIManager.put("RadioButton.background", new ColorUIResource(COL_YEL_DARK));

			//UIManager.put("ComboBox.foreground", new ColorUIResource(COL_BG_YEL));	// Schrift
			//UIManager.put("ComboBox.buttonHighlight", new ColorUIResource(COL_SEL_GREEN));
			//UIManager.put("ComboBox.buttonDarkShadow", new ColorUIResource(COL_SEL_GREEN));
			UIManager.put("ComboBox.buttonBackground", new ColorUIResource(COL_BG_YEL));
			UIManager.put("ComboBox.background", new ColorUIResource(COL_BG_YEL));
			UIManager.put("ComboBox.selectionBackground", new ColorUIResource(COL_SEL_GREEN));
			UIManager.put("ComboBox.selectionForeground", new ColorUIResource(Color.white));
			UIManager.put("ComboBox.disabledBackground", new ColorUIResource(COL_YEL_DARK));

			
			a= new ArrayList();
			a.add(new Float(0.3));	// 0.3
			a.add(new Float(0.0));	// 0.0
			a.add(new ColorUIResource(COL_BG_YEL));	// [r=221,g=232,b=243]
			a.add(new ColorUIResource(Color.white));	// [r=255,g=255,b=255]
			a.add(new ColorUIResource(COL_YEL_DARK));	// [r=184,g=207,b=229]
			UIManager.put("ScrollBar.gradient", a);
			//UIManager.put("ScrollBar.thumbDarkShadow", new ColorUIResource(COL_YEL_DARK));
			//UIManager.put("ScrollBar.darkShadow", new ColorUIResource(COL_YEL_DARK));
			UIManager.put("ScrollBar.background", new ColorUIResource(COL_BG_GREEN));
			UIManager.put("ScrollBar.thumb", new ColorUIResource(COL_BG_GREEN));
			UIManager.put("ScrollBar.shadow", new ColorUIResource(COL_SHADOW_SEMI_DARK));
			//UIManager.put("ScrollBar.highlight", new ColorUIResource(COL_YEL_DARK));
			//UIManager.put("ScrollBar.trackHighlight", new ColorUIResource(COL_YEL_DARK));
			UIManager.put("ScrollBar.track", new ColorUIResource(COL_SHADOW_SEMI_DARK));
			//UIManager.put("ScrollBar.foreground", new ColorUIResource(COL_YEL_DARK));
			UIManager.put("ScrollBar.thumbShadow", new ColorUIResource(COL_SEL_GREEN));
			//UIManager.put("ScrollBar.thumbHighlight", new ColorUIResource(COL_YEL_DARK));
			UIManager.put("scrollbar", new ColorUIResource(COL_BG_YEL));
			
			//UIManager.put("ScrollPane.foreground", new ColorUIResource(COL_BG_YEL));
			UIManager.put("ScrollPane.background", new ColorUIResource(COL_BG_GREEN));
			
			//UIManager.put("Table.foreground", new ColorUIResource(COL_BG_YEL));
			UIManager.put("Table.background", new ColorUIResource(COL_BG_YEL));
			UIManager.put("Table.selectionForeground", new ColorUIResource(Color.white));
			UIManager.put("Table.selectionBackground", new ColorUIResource(COL_SEL_GREEN));
			//UIManager.put("TableHeader.foreground", new ColorUIResource(COL_BG_YEL));
			UIManager.put("TableHeader.background", new ColorUIResource(COL_BG_GREEN));
			//UIManager.put("Table.gridColor", new ColorUIResource(COL_BG_YEL));
			//UIManager.put("Table.focusCellBackground", new ColorUIResource(COL_BG_YEL));
			//UIManager.put("Table.focusCellForeground", new ColorUIResource(COL_BG_YEL));
			UIManager.put("Table.scrollPaneBorder", new ColorUIResource(COL_YEL_DARK));
			
			UIManager.put("OptionPane.background", new ColorUIResource(COL_SEL_GREEN));			
			
			UIManager.put("List.background", new ColorUIResource(COL_BG_YEL));	// file chooser
			UIManager.put("List.selectionBackground", new ColorUIResource(COL_SEL_GREEN));	
			UIManager.put("List.selectionForeground", new ColorUIResource(Color.white));	
			
			UIManager.put("Menu.selectionForeground", new ColorUIResource(Color.white));	
			UIManager.put("Menu.selectionBackground", new ColorUIResource(COL_SEL_GREEN));	
			UIManager.put("PopupMenu.background", new ColorUIResource(COL_BG_YEL));	
			UIManager.put("PopupMenu.foreground", new ColorUIResource(COL_YEL_DARK));	
			
			UIManager.put("Menu.background", new ColorUIResource(COL_BG_YEL));	
			UIManager.put("RadioButtonMenuItem.selectionForeground", new ColorUIResource(Color.white));	
			UIManager.put("RadioButtonMenuItem.selectionBackground", new ColorUIResource(COL_SEL_GREEN));	
			UIManager.put("RadioButtonMenuItem.background", new ColorUIResource(COL_BG_YEL));	
						
			UIManager.put("ToolTip.background", new ColorUIResource(COL_BG_YEL));	
			
			// UIManager.put("Table.background", new ColorUIResource(COL_BG_GREEN)); // already defined
			
		} else {
			UIManager.put("Panel.background", uiDefs.get("Panel.background"));
			UIManager.put("TabbedPane.selected", uiDefs.get("TabbedPane.selected"));	// COL_SEL_GREEN, but cannot set text to white		
			UIManager.put("TabbedPane.contentAreaColor", uiDefs.get("TabbedPane.contentAreaColor"));
			UIManager.put("TabbedPane.borderHightlightColor", uiDefs.get("TabbedPane.borderHightlightColor"));
			UIManager.put("TabbedPane.focus", uiDefs.get("TabbedPane.focus"));
			UIManager.put("TabbedPane.light", uiDefs.get("TabbedPane.light"));
			UIManager.put("TabbedPane.unselectedBackground", uiDefs.get("TabbedPane.unselectedBackground"));
			SimpleBinPlotterPanel.defaultPanBg= (Color) uiDefs.get("TabbedPane.unselectedBackground");

			UIManager.put("ProgressBar.foreground", uiDefs.get("ProgressBar.foreground"));
			UIManager.put("ProgressBar.background", uiDefs.get("ProgressBar.background"));
			UIManager.put("ProgressBar.border", uiDefs.get("ProgressBar.border"));

			UIManager.put("TextField.background", uiDefs.get("TextField.background"));
			UIManager.put("TextField.inactiveBackground", uiDefs.get("TextField.inactiveBackground"));
			UIManager.put("TextField.selectionBackground", uiDefs.get("TextField.selectionBackground"));
			UIManager.put("TextField.selectionForeground", uiDefs.get("TextField.selectionForeground"));
			
			UIManager.put("RadioButton.background", uiDefs.get("RadioButton.background"));
			UIManager.put("RadioButton.disabledText", uiDefs.get("RadioButton.disabledText"));
			UIManager.put("CheckBox.background", uiDefs.get("CheckBox.background"));
			UIManager.put("CheckBox.disabledText", uiDefs.get("CheckBox.disabledText"));
			
			UIManager.put("Button.gradient", uiDefs.get("Button.gradient"));
			UIManager.put("Button.background", uiDefs.get("Button.background"));
			UIManager.put("Button.disabledToolBarBorderBackground", uiDefs.get("Button.disabledToolBarBorderBackground"));			
			UIManager.put("Button.disabledText", uiDefs.get("Button.disabledText"));
			//UIManager.put("Button.border", uiDefs.get("Button.border"));

			UIManager.put("RadioButton.background", uiDefs.get("RadioButton.background"));
			
			UIManager.put("ToolBar.dockingBackground", uiDefs.get("ToolBar.dockingBackground"));
			UIManager.put("ToolBar.dockingForeground", uiDefs.get("ToolBar.dockingForeground"));
			UIManager.put("ToolBar.floatingBackground", uiDefs.get("ToolBar.floatingBackground"));
			UIManager.put("ToolBar.floatingForeground", uiDefs.get("ToolBar.floatingForeground"));
			
			UIManager.put("MenuBar.gradient", uiDefs.get("MenuBar.gradient"));
			UIManager.put("MenuItem.background", uiDefs.get("MenuItem.background"));
			UIManager.put("MenuItem.selectionForeground", uiDefs.get("MenuItem.selectionForeground"));
			UIManager.put("MenuItem.selectionBackground", uiDefs.get("MenuItem.selectionBackground"));

			UIManager.put("ComboBox.buttonBackground", uiDefs.get("ComboBox.buttonBackground"));
			UIManager.put("ComboBox.selectionBackground", uiDefs.get("ComboBox.selectionBackground"));
			UIManager.put("ComboBox.disabledBackground", uiDefs.get("ComboBox.disabledBackground"));
			UIManager.put("ComboBox.background", uiDefs.get("ComboBox.background"));
			UIManager.put("ComboBox.selectionForeground", uiDefs.get("ComboBox.selectionForeground"));
			
			UIManager.put("ScrollBar.gradient", uiDefs.get("ScrollBar.gradient"));
			UIManager.put("ScrollBar.background", uiDefs.get("ScrollBar.background"));
			UIManager.put("ScrollBar.thumb", uiDefs.get("ScrollBar.thumb"));
			UIManager.put("ScrollBar.shadow", uiDefs.get("ScrollBar.shadow"));
			UIManager.put("ScrollBar.track", uiDefs.get("ScrollBar.track"));
			UIManager.put("ScrollBar.thumbShadow", uiDefs.get("ScrollBar.thumbShadow"));
			UIManager.put("scrollbar", uiDefs.get("scrollbar"));
			
			UIManager.put("ScrollPane.background", uiDefs.get("ScrollPane.background"));
			
			UIManager.put("Table.background", uiDefs.get("Table.background"));
			UIManager.put("Table.selectionForeground", uiDefs.get("Table.selectionForeground"));
			UIManager.put("Table.selectionBackground", uiDefs.get("Table.selectionBackground"));
			UIManager.put("TableHeader.background", uiDefs.get("TableHeader.background"));
			UIManager.put("Table.scrollPaneBorder", uiDefs.get("Table.scrollPaneBorder"));
			
			UIManager.put("OptionPane.background", uiDefs.get("OptionPane.background"));			
			
			UIManager.put("List.background", uiDefs.get("List.background"));	// file chooser
			UIManager.put("List.selectionBackground", uiDefs.get("List.selectionBackground"));	
			UIManager.put("List.selectionForeground", uiDefs.get("List.selectionForeground"));	
			
			UIManager.put("Menu.selectionForeground", uiDefs.get("Menu.selectionForeground"));	
			UIManager.put("Menu.selectionBackground", uiDefs.get("Menu.selectionBackground"));	
			UIManager.put("PopupMenu.background", uiDefs.get("PopupMenu.background"));	
			UIManager.put("PopupMenu.foreground", uiDefs.get("PopupMenu.foreground"));	

			UIManager.put("Menu.background", uiDefs.get("Menu.background"));	
			UIManager.put("RadioButtonMenuItem.selectionForeground", uiDefs.get("RadioButtonMenuItem.selectionForeground"));	
			UIManager.put("RadioButtonMenuItem.selectionBackground", uiDefs.get("RadioButtonMenuItem.selectionBackground"));	
			UIManager.put("RadioButtonMenuItem.background", uiDefs.get("RadioButtonMenuItem.background"));
			
			UIManager.put("ToolTip.background", uiDefs.get("ToolTip.background"));	

			// UIManager.put("Table.background", uiDefs.get("Table.background")); // already defined
		}
		try {
			UIManager.setLookAndFeel(lafName);
		} catch (Exception e) {
			; // :)
		}

	}

	Image favimage, imgAbout;
	public static ImageIcon favicon, imgGreen, imgRed, imgGrey;
	ImageIcon iconNew, iconNewDis, iconLoad, iconLoadDis, iconCopy, iconCopyDis, iconRun, iconRunDis, iconStop, iconStopDis;
	private void loadIcons() {

//	    try {
//    	MediaTracker m = new MediaTracker(this);
//     	m.addImage(faviImage, 0);
//     	m.waitForAll();
//    } catch (Exception e) {
//    	; // :)
//    }

		try {
			favimage = Toolkit.getDefaultToolkit().getImage(ClassLoader.getSystemResource("pics/favicon.png"));
		} catch (Exception e) {
			; // :)
		}
		if (favimage!= null) {
			setIconImage(favimage);
			favicon= new ImageIcon(favimage);
		}

	    // take classloader for jar file
	    // also possible: this.getClass().getClassLoader().getSystemResource
		try {
		    imgGreen= new ImageIcon(ClassLoader.getSystemResource("pics/icon_green.gif"), "finished");	
		} catch (Exception e) {
			; // :)
		}
		try {
		    imgRed= new ImageIcon(ClassLoader.getSystemResource("pics/icon_red.gif"), "ready");
		} catch (Exception e) {
			; // :)
		}
		try {
		    imgGrey= new ImageIcon(ClassLoader.getSystemResource("pics/icon_disabled.gif"), "not inited");
		} catch (Exception e) {
			; // :)
		}

		try {
		    iconNew= new ImageIcon(ClassLoader.getSystemResource("pics/new.gif"), "New");
		} catch (Exception e) {
			; // :)
		}
		try {
		    iconNewDis= new ImageIcon(ClassLoader.getSystemResource("pics/new_disabled.gif"), "New");
		} catch (Exception e) {
			; // :)
		}
		try {
		    iconLoad= new ImageIcon(ClassLoader.getSystemResource("pics/open.gif"), "Load");
		} catch (Exception e) {
			; // :)
		}
		try {
		    iconLoadDis= new ImageIcon(ClassLoader.getSystemResource("pics/open_disabled.gif"), "Load");
		} catch (Exception e) {
			; // :)
		}
		try {
		    iconCopy= new ImageIcon(ClassLoader.getSystemResource("pics/copy.gif"), "Copy");
		} catch (Exception e) {
			; // :)
		}
		try {
		    iconCopyDis= new ImageIcon(ClassLoader.getSystemResource("pics/copy_disabled.gif"), "Copy");
		} catch (Exception e) {
			; // :)
		}
		try {
		    iconRun= new ImageIcon(ClassLoader.getSystemResource("pics/run.gif"), "Run");
		} catch (Exception e) {
			; // :)
		}
		try {
		    iconRunDis= new ImageIcon(ClassLoader.getSystemResource("pics/run_disabled.gif"), "Run");
		} catch (Exception e) {
			; // :)
		}
		try {
		    iconStop= new ImageIcon(ClassLoader.getSystemResource("pics/stop.gif"), "Stop");
		} catch (Exception e) {
			; // :)
		}
		try {
		    iconStopDis= new ImageIcon(ClassLoader.getSystemResource("pics/stop_disabled.gif"), "Stop");
		} catch (Exception e) {
			; // :)
		}

		
		try {
		    imgAbout= Toolkit.getDefaultToolkit().getImage(ClassLoader.getSystemResource("pics/rnaseq_group_1.jpg"));
		} catch (Exception e) {
			; // :)
		}
	}
	
	private JTabbedPane getTabPane() {
		if (tabPane == null) {
			tabPane= new JTabbedPane();
			
			tabPane.addTab(NAME_ANN, getPanelAnnotation());	// getPanelAnnotation()
			
			tabPane.addTab(NAME_EXP, getPanExpression());	
			
			tabPane.addTab(NAME_RTFRAG, getRTFragPanel());	

			tabPane.addTab(NAME_SEQ, getSeqPanel());	
			
			readyOrNot= new ReadyOrNot[] {
					(ReadyOrNot) getPanelAnnotation(),
					(ReadyOrNot) getPanExpression(),
					(ReadyOrNot) getRTFragPanel(),
					(ReadyOrNot) getSeqPanel()
			};
			
			tabPane.addChangeListener(new ChangeListener() {
		        // This method is called whenever the selected tab changes
		        public void stateChanged(ChangeEvent evt) {
		        	if (blockThread!= null&& blockThread.isAlive())
		        		return;
		            JTabbedPane pane = (JTabbedPane) evt.getSource();
//		            if (UIManager.getLookAndFeel().getName().contains("Metal"))
//		            	pane.setForeground(Color.white);
		            int sel= pane.getSelectedIndex();	// removingAll throws (-1) then 0	            
		            if (pane.getComponentAt(sel) instanceof ReadyOrNot&& 
		            		((ReadyOrNot) pane.getComponentAt(sel)).isReady()) 
		            	getButRun().setEnabled(true);
		            else
		            	getButRun().setEnabled(false);
		            labStatus.setString(emptyString);
		            labStatus.repaint();
		        }
			});
			add(tabPane, BorderLayout.CENTER);

		}

		return tabPane;
	}
	
	
	ReadyOrNot[] readyOrNot;
	 
	private void initUI() {
		SwingUtilities.updateComponentTreeUI(FluxSimulatorGUI.this);		

		loadIcons();
		
		setLayout(new BorderLayout());

		add(createToolBar(), BorderLayout.PAGE_START);
		
		Dialoguer dia= new Dialoguer();
		Constants.dialog= dia;
		
		labStatus= new MyProgressBar();
		labStatus.setString(emptyString);
		labStatus.setMaximum(9);
		labStatus.setValue(-1);
		//labStatus.setBorder(BorderFactory.createLineBorder(Color.black));
		//Constants.progress= labStatus;
		add(labStatus, BorderLayout.SOUTH);
		
		setJMenuBar(createJMenuBar());
		//setMenuBar(createMenuBar());
		
		uiUpdate(null);
	}
	
	FragmenterGUI rtFragPanel;
	private Component getRTFragPanel() {
		if (rtFragPanel == null) {
			rtFragPanel = new FragmenterGUI();
			rtFragPanel.setName(NAME_RTFRAG);
		}

		return rtFragPanel;
	}

	SequencerGUI seqPanel;
	private Component getSeqPanel() {
		if (seqPanel == null) {
			seqPanel = new SequencerGUI();
			seqPanel.setName(NAME_RTFRAG);
		}

		return seqPanel;	
	}
	private void swap(JComponent comp, String name1, String name2) {
		Component[] c= comp.getComponents();
		int posRT= -1, posFrag= -1;
		for (int i = 0; i < c.length; i++) {
			if (c[i].getName()!= null&& c[i].getName().equals(name1))
				posRT= i;
			else if (c[i].getName()!= null&& c[i].getName().equals(name2))
				posFrag= i;
		}
		assert(posRT>= 0&& posFrag>= 0);
		
		int hi= Math.max(posRT,posFrag), lo= Math.min(posRT, posFrag);
		Component hiComp= comp.getComponent(hi), loComp= comp.getComponent(lo);
		comp.remove(loComp);
		comp.remove(hiComp);
		comp.add(hiComp, lo);
		comp.add(loComp, hi);
		comp.revalidate();				
	}
	
	/**
	 * @param name1
	 * @param name2
	 */
	private void swapTabs(String name1, String name2) {
		//swap(tabPane, NAME_RT, NAME_FRAG);
		Component[] c= tabPane.getComponents();
		int posRT= -1, posFrag= -1;
		for (int i = 0; i < c.length; i++) {
			if (c[i].getName()!= null&& c[i].getName().equals(name1))
				posRT= i;
			else if (c[i].getName()!= null&& c[i].getName().equals(name2))
				posFrag= i;
		}
		assert(posRT>= 0&& posFrag>= 0);
		
		int hi= Math.max(posRT,posFrag), lo= Math.min(posRT, posFrag);
		Component hiComp= tabPane.getComponent(hi), loComp= tabPane.getComponent(lo);
		
		tabPane.removeAll();
		for (int i = 0; i < c.length; i++) {
			Component comp= null;
			if (i== lo)
				comp= hiComp;
			else if (i== hi)
				comp= loComp;
			else
				comp= c[i];
			tabPane.addTab(comp.getName(), comp);
			if (i>= lo)
				tabPane.setEnabledAt(i, false);
		}
		tabPane.revalidate();				
	}
	
	JButton butRun, butCancel, butNew, butLoad, butCopy;
	JToolBar toolBar;
	BlockingThread blockThread;
	private JToolBar createToolBar() {
		if (toolBar == null) {
			toolBar = new JToolBar("Control bar", JToolBar.HORIZONTAL);
			
			butNew= new JButton();
			butNew.setAction(getActionNew());
			if (iconNew== null) 
				butNew.setText("New");
			else {
				butNew.setIcon(iconNew);
				butNew.setDisabledIcon(iconNewDis);
				butNew.setOpaque(false);
			}
			butNew.setToolTipText("Create New Run");
			toolBar.add(butNew);
			
			butLoad= new JButton();
			butLoad.setAction(getActionLoad());
			if (iconLoad== null)
				butLoad.setText("Open");
			else {
				butLoad.setIcon(iconLoad);
				butLoad.setDisabledIcon(iconLoadDis);
				butLoad.setOpaque(false);
			}
			butLoad.setToolTipText("Load Run");
			toolBar.add(butLoad);

			butCopy= new JButton();
			butCopy.setAction(getActionCopy());
			if (iconCopy== null)
				butCopy.setText("Copy");
			else {
				butCopy.setIcon(iconCopy);
				butCopy.setDisabledIcon(iconCopyDis);
				butCopy.setOpaque(false);
			}
			butCopy.setToolTipText("Copy Run");
			toolBar.add(butCopy);
			
			toolBar.add(new JToolBar.Separator());
			
			toolBar.add(getButRun());
			toolBar.add(getButCancel());
			
			toolBar.setFloatable(false);
		}

		return toolBar;
	}
	
	
	private JButton getButRun() {
		if (butRun == null) {
			butRun = new JButton();
			if (iconRun== null)
				butRun.setText("Do");
			else {
				butRun.setIcon(iconRun);
				butRun.setDisabledIcon(iconRunDis);
				butRun.setOpaque(false);
			}
			butRun.setToolTipText("(Re-)Run Current Step");
			butRun.addActionListener(new ActionListener(){
				public void actionPerformed(ActionEvent e) {
					if (settings!= null)
						settings.save();
						
					Component c= tabPane.getSelectedComponent();
					if (!(c instanceof StoppableRunnable))
						return;
					if (((ReadyOrNot) c).isFinished()) {
						if (Constants.dialog.checkOverwrite("You are about to overwrite data from the project, proceed?")) {
							int x= tabPane.getSelectedIndex();
							for (int i = x; i < tabPane.getTabCount(); i++) {
								if (tabPane.getComponentAt(i) instanceof ReadyOrNot)
									((ReadyOrNot) tabPane.getComponentAt(i)).killResult();
							}
							uiUpdate(null);
						} else 
							return;
					}
					blockThread= new BlockingThread((StoppableRunnable) c, blockThread);
					blockThread.start();
				}

			});

			butRun.setEnabled(false);
		}

		return butRun;
	}

	private void killBlockThread() {
		if (blockThread== null|| (!blockThread.isAlive()))
			return;
		blockThread.cancel();
		blockThread.interrupt();
		// deadlock when updating progress bar, in sequencing -> gtfreader
//		while (blockThread.isAlive())
//			try {
//				blockThread.interrupt();
//				Thread.sleep(100);
//			} catch (InterruptedException e1) {
//				; // :)
//			}
	}
	
	private JButton getButCancel() {
		if (butCancel == null) {
			butCancel = new JButton();
			if (iconStop== null)
				butCancel.setText("Stop");
			else {
				butCancel.setIcon(iconStop);
				butCancel.setDisabledIcon(iconStopDis);
				butCancel.setOpaque(false);
			}
			butCancel.setToolTipText("Cancel Run");
			butCancel.addActionListener(new ActionListener(){
				public void actionPerformed(ActionEvent e) {
					killBlockThread();
				}
			});
			butCancel.setEnabled(false);
		}

		return butCancel;
	}
	
	MenuBar menuBar;
	private MenuBar createMenuBar() {
		if (menuBar == null) {
			menuBar = new MenuBar();
			
			Menu menu;
			MenuItem item;
			menu= new Menu("Project");
			
			item= new MenuItem("New");
			item.addActionListener(getActionNew());
			menu.add(item);
			
			item= new MenuItem("Open");
			item.addActionListener(getActionLoad());
			menu.add(item);

			MenuItem saveItem;
			saveItem= new MenuItem("Save");
			saveItem.addActionListener(getActionCopy());
			//saveItem.setEnabled(getActionSave().isEnabled());
			menu.add(saveItem);

			item= new MenuItem("-");
			menu.add(item);
			
			item= new MenuItem("Exit");
			item.addActionListener(new AbstractAction(){
				public void actionPerformed(ActionEvent e) {
					System.exit(0);
				}
			});
			menu.add(item);
			menuBar.add(menu);
			
			menu= new Menu("View");
			item= new MenuItem("Parameters");
			menu.add(item);
			menuBar.add(menu);
		}

		return menuBar;
	}

	Action actionNew;
	private Action getActionNew() {
		if (actionNew == null) {
			actionNew = new AbstractAction() {
				public void actionPerformed(ActionEvent e) {

					if (checkProjectContainsData()&& (!askUserOverwriteProject())) 
						return;					

					FluxSimulatorGUI.this.settings= settings.createDefaults();
					Dialog dia= getProjectDialog();
					tfAnn.setText("");
					tfName.setText("");
					tfPDir.setText("");
					tfPar.setText("");
					tfPro.setText("");
					tfFrg.setText("");
					tfSeq.setText("");
					tfTmp.setText(System.getProperty("java.io.tmpdir"));
					dia.setTitle(TITLE_PROJECT_NEW);
					dia.setVisible(true);	// waits for the dia to close (win)
				}
			};
		}
		return actionNew;
	}

	private void updateTitle(File parFile) {
		String title= parFile.getName();
		title= title.substring(0,title.lastIndexOf("."));
		setTitle(title+" - "+ FluxSimulator.FLUX_VERSION);
	}
	
	
	JDialog newDialog;
	File parFile;
	public static final String TITLE_PROJECT_NEW= "Create New Run", TITLE_PROJECT_COPY= "Copy Current Run";
	JTextField tfAnn, tfName, tfPDir, tfPar, tfPro, tfFrg, tfSeq, tfTmp;
	private JTextField getTFtmp() {
		if (tfTmp == null) {
			tfTmp = new JTextField(30);
		}

		return tfTmp;
	}
	
	private JDialog getProjectDialog() {
		if (newDialog == null) {
			newDialog = new JDialog(FluxSimulatorGUI.this, true);
			newDialog.setLayout(new BorderLayout());
			
			JPanel grid= new JPanel();
			grid.setLayout(new GridBagLayout());
			JButton but;
			int tfWidth= 30;
			final JButton butOKnew= new JButton("OK");
			butOKnew.setEnabled(false);
			KeyAdapter keyCR= new KeyAdapter(){
				@Override
				public void keyPressed(KeyEvent e) {
					int key= e.getKeyCode();
				    if (key!= KeyEvent.VK_ENTER) 
				    	return;
				    if (butOKnew.isEnabled())
				    	butOKnew.doClick();
				}
			};

			tfAnn= new JTextField(tfWidth);
			tfName= new JTextField(tfWidth);
			tfPDir= new JTextField(tfWidth);
			tfPar= new JTextField(tfWidth);
			tfPro= new JTextField(tfWidth);
			tfFrg= new JTextField(tfWidth);
			tfSeq= new JTextField(tfWidth);
			getTFtmp();
			tfGenome= new JTextField(tfWidth);

			DocumentListener docChecker= new DocumentListener() {
				public void changedUpdate(DocumentEvent e) {
				}
				public void insertUpdate(DocumentEvent e) {
					doIt();
				}
				public void removeUpdate(DocumentEvent e) {
					doIt();
				}
				private void doIt() {
					if (tfAnn.getText()!= null&& tfAnn.getText().length()> 0
							&& tfPar.getText()!= null&& tfPar.getText().length()> 0
							&& tfPro.getText()!= null&& tfPro.getText().length()> 0
							&& tfFrg.getText()!= null&& tfFrg.getText().length()> 0
							&& tfSeq.getText()!= null&& tfSeq.getText().length()> 0)
						butOKnew.setEnabled(true);
					else
						butOKnew.setEnabled(false);
				}
			};
			
			tfAnn.addKeyListener(keyCR);
			tfName.addKeyListener(keyCR);
			tfPDir.addKeyListener(keyCR);
			tfPro.addKeyListener(keyCR);
			tfFrg.addKeyListener(keyCR);
			tfSeq.addKeyListener(keyCR);
			getTFtmp().addKeyListener(keyCR);
			tfGenome.addKeyListener(keyCR);
			
			tfAnn.getDocument().addDocumentListener(docChecker); 
			tfName.getDocument().addDocumentListener(docChecker);
			tfPDir.getDocument().addDocumentListener(docChecker);
			tfPro.getDocument().addDocumentListener(docChecker);
			tfFrg.getDocument().addDocumentListener(docChecker);
			tfSeq.getDocument().addDocumentListener(docChecker);
			getTFtmp().getDocument().addDocumentListener(docChecker);

			/*tfAnn.getDocument().addDocumentListener(new DocumentListener(){
				public void changedUpdate(DocumentEvent e) {
					downcast(e);
				}
				public void insertUpdate(DocumentEvent e) {
					downcast(e);
				}
				public void removeUpdate(DocumentEvent e) {
					downcast(e);
				}
				void downcast(DocumentEvent e) {
					int p= tfAnn.getText().lastIndexOf(File.separator);
					if (p== -1)
						return;
					String spath= tfAnn.getText().substring(0, p+1);
					if (tfPDir.getText().length()== 0) {
						tfPDir.setText(spath);
						if (tfPar.getText().length()== 0)
							tfPar.setText(spath);
						if (tfPro.getText().length()== 0)
							tfPro.setText(spath);
						if (tfFrg.getText().length()== 0)
							tfFrg.setText(spath);
						if (tfSeq.getText().length()== 0)
							tfSeq.setText(spath);
					}
				}
			}); */
			
			tfPDir.getDocument().addDocumentListener(new DocumentListener(){
				public void changedUpdate(DocumentEvent e) {
					downcast(e);
				}

				public void insertUpdate(DocumentEvent e) {
					downcast(e);
				}

				public void removeUpdate(DocumentEvent e) {
					downcast(e);
				}				
				void downcast(DocumentEvent e) {
					String spath= tfPDir.getText();
					if (spath.length()== 0)
						return;
					if (spath.charAt(spath.length()-1)!= File.separatorChar)
						spath+= File.separator;
					if (tfPar.getText().length()== 0
							|| tfPar.getText().charAt(tfPar.getText().length()-1)== File.separatorChar)
						tfPar.setText(spath);
					else {
						int p= tfPar.getText().lastIndexOf(File.separator)+ 1;
						tfPar.setText(spath+tfPar.getText().substring(p));
					}
					if (tfPro.getText().length()== 0
							|| tfPro.getText().charAt(tfPro.getText().length()-1)== File.separatorChar)
						tfPro.setText(spath);
					else {
						int p= tfPro.getText().lastIndexOf(File.separator)+ 1;
						tfPro.setText(spath+tfPro.getText().substring(p));
					}
					if (tfFrg.getText().length()== 0
							|| tfFrg.getText().charAt(tfFrg.getText().length()-1)== File.separatorChar)
						tfFrg.setText(spath);
					else {
						int p= tfFrg.getText().lastIndexOf(File.separator)+ 1;
						tfFrg.setText(spath+tfFrg.getText().substring(p));
					}
					if (tfSeq.getText().length()== 0
							|| tfSeq.getText().charAt(tfSeq.getText().length()-1)== File.separatorChar)
						tfSeq.setText(spath);
					else {
						int p= tfSeq.getText().lastIndexOf(File.separator)+ 1;
						tfSeq.setText(spath+tfSeq.getText().substring(p));
					}
				}
			});
			
			tfName.getDocument().addDocumentListener(new DocumentListener(){

				public void changedUpdate(DocumentEvent e) {
					downcast(e);
				}

				public void insertUpdate(DocumentEvent e) {
					downcast(e);
				}

				public void removeUpdate(DocumentEvent e) {
					downcast(e);
				}
				void downcast(DocumentEvent e) {
					int p= tfPar.getText().lastIndexOf(File.separator);					
					String s= tfPar.getText().substring(0,Math.max(0,p));
					tfPar.setText(s+(s.length()==0?"":File.separator)+tfName.getText()+FluxSimulatorSettings.DEF_SFX_PAR);
					p= tfPro.getText().lastIndexOf(File.separator);
					s= tfPro.getText().substring(0,Math.max(0,p));
					tfPro.setText(s+(s.length()==0?"":File.separator)+tfName.getText()+FluxSimulatorSettings.DEF_SFX_PRO);
					p= tfFrg.getText().lastIndexOf(File.separator);
					s= tfFrg.getText().substring(0,Math.max(0,p));
					tfFrg.setText(s+(s.length()==0?"":File.separator)+tfName.getText()+FluxSimulatorSettings.DEF_SFX_LIB);
					p= tfSeq.getText().lastIndexOf(File.separator);
					s= tfSeq.getText().substring(0,Math.max(0,p));
					tfSeq.setText(s+(s.length()==0?"":File.separator)+tfName.getText()+FluxSimulatorSettings.DEF_SFX_SEQ);
				}
			});
			
			GridBagConstraints c= new GridBagConstraints();
			c.fill = GridBagConstraints.HORIZONTAL;
			c.insets= new Insets(1,10,1,10);
			c.gridheight= 1;
			
			c.gridwidth= 2; c.gridx = 0; c.gridy = 0;
			grid.add(new JLabel("Reference Annotation"), c);
			c.gridwidth= 1; c.gridx = 2; c.gridy = 0;
			grid.add(tfAnn, c);
			but= new JButton("Browse..");
			but.addActionListener(new ActionListener(){
				public void actionPerformed(ActionEvent e) {
					JFileChooser dia= getChooser();	// FileDialog dia= getLoadDialog();
					if (settings.getRefFile()!= null)
						try {
							dia.setCurrentDirectory(settings.getRefFile().getParentFile());
						} catch (Exception e2) {
							; // :)
						}
					else if (parFile!= null) 
						try {
							dia.setCurrentDirectory(parFile.getParentFile());
						} catch (Exception e2) {
							; // :)
						}
					dia.setDialogTitle("Choose Reference Annotation");	// dia.setTitle("Choose Reference Annotation");
					dia.setFileSelectionMode(JFileChooser.FILES_ONLY);
					dia.removeChoosableFileFilter(getFilterPar());
					dia.addChoosableFileFilter(getFilterGTF());
					int rval= dia.showOpenDialog(FluxSimulatorGUI.this);					
					//dia.setVisible(true);
					
//					String fname= dia.getFile();					
//					if (fname== null)
//						return;
//					File f= new File(dia.getDirectory()+ File.separator+ fname);
					File f= dia.getSelectedFile();
					if (f== null|| f.isDirectory()|| (!f.canRead())|| rval!= JFileChooser.APPROVE_OPTION)
						return;
					tfAnn.setText(f.getAbsolutePath());
				}
			});
			c.gridwidth= 1; c.gridx = 3; c.gridy = 0;
			grid.add(but, c);
			c.gridwidth= 2; c.gridx = 0; c.gridy = 1;
			grid.add(new JLabel("Run Name"), c);
			c.gridwidth= 1; c.gridx = 2; c.gridy = 1;
			grid.add(tfName, c);
			c.gridwidth= 1; c.gridx = 0; c.gridy = 2;
			grid.add(new JLabel("Run Folder"), c);
			c.gridwidth= 1; c.gridx = 2; c.gridy = 2;
			grid.add(tfPDir, c);
			but= new JButton("Browse..");
			but.addActionListener(new ActionListener(){
				public void actionPerformed(ActionEvent e) {
					JFileChooser dia= getChooser();
					dia.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
					dia.setDialogTitle("Choose Run Folder");
					dia.removeChoosableFileFilter(getFilterGTF());
					dia.removeChoosableFileFilter(getFilterPar());
					
					//dia.setVisible(true);
					int rval= dia.showOpenDialog(FluxSimulatorGUI.this);
					File f= dia.getSelectedFile(); // should be a dir
					if (f== null)
						f= dia.getCurrentDirectory();
					if (f== null|| rval!= JFileChooser.APPROVE_OPTION)
						return;
					
					tfPDir.setText(f.getAbsolutePath());
				}
			});
			c.gridwidth= 1; c.gridx = 3; c.gridy = 2;
			grid.add(but, c);
			c.gridwidth= 1; c.gridx = 1; c.gridy = 3;
			grid.add(new JLabel("Parameter File"), c);
			c.gridwidth= 1; c.gridx = 2; c.gridy = 3;
			grid.add(tfPar, c);
			c.gridwidth= 1; c.gridx = 1; c.gridy = 4;
			grid.add(new JLabel("Profile File"), c);
			c.gridwidth= 1; c.gridx = 2; c.gridy = 4;
			grid.add(tfPro, c);
			c.gridwidth= 1; c.gridx = 1; c.gridy = 5;
			grid.add(new JLabel("Fragment File"), c);
			c.gridwidth= 1; c.gridx = 2; c.gridy = 5;
			grid.add(tfFrg, c);
			c.gridwidth= 1; c.gridx = 1; c.gridy = 6;
			grid.add(new JLabel("Sequencing File"), c);
			c.gridwidth= 1; c.gridx = 2; c.gridy = 6;
			grid.add(tfSeq, c);
			c.gridwidth= 2; c.gridx = 0; c.gridy = 7;
			grid.add(new JLabel("Temporary Folder"), c);
			c.gridwidth= 1; c.gridx = 2; c.gridy = 7;
			grid.add(tfTmp, c);
			but= new JButton("Browse..");
			but.addActionListener(getTmpAction());
			c.gridwidth= 1; c.gridx = 3; c.gridy = 7;
			grid.add(but, c);
			if (settings!= null&& settings.getGenDir()!= null)
				tfGenome.setText(settings.getGenDir().getAbsolutePath());
			tfGenome.getDocument().addDocumentListener(new DocumentListener(){
				public void changedUpdate(DocumentEvent e) {
					File f= new File(tfGenome.getText());
					if ((!f.exists())|| (!f.isDirectory()))
						return;
					settings.setGenDir(f);
				}

				public void insertUpdate(DocumentEvent e) {
				}

				public void removeUpdate(DocumentEvent e) {
				}
			});
			c.gridwidth= 2; c.gridx = 0; c.gridy = 8;
			grid.add(new JLabel("Genome Folder"), c);
			c.gridwidth= 1; c.gridx = 2; c.gridy = 8;
			grid.add(tfGenome, c);
			but= new JButton("Browse..");
			but.addActionListener(getActionGen());
			c.gridwidth= 1; c.gridx = 3; c.gridy = 8;
			grid.add(but, c);

			newDialog.add(grid, BorderLayout.CENTER);
			
			JPanel buts= new JPanel();
			but= new JButton("Clear");
			but.addActionListener(new ActionListener(){
				public void actionPerformed(ActionEvent e) {
					tfAnn.setText("");
					tfName.setText("");
					tfPDir.setText("");
					tfPar.setText("");
					tfPro.setText("");
					tfFrg.setText("");
					tfSeq.setText("");
					tfTmp.setText(System.getProperty("java.io.tmpdir"));
				}
			});
			buts.add(but);
			buts.add(new JSeparator(JSeparator.HORIZONTAL));
			but= new JButton("Cancel");
			but.addActionListener(new ActionListener(){
				public void actionPerformed(ActionEvent e) {
					newDialog.setVisible(false);
				}
			});
			buts.add(but);
			buts.add(new JSeparator(JSeparator.HORIZONTAL));
			butOKnew.addActionListener(new ActionListener(){
				public void actionPerformed(ActionEvent e) {

					File f= new File(tfPDir.getText());
					File f1= new File(f.getAbsolutePath());
					f= f1;
					if (!f.exists()) {
						int ret= JOptionPane.showOptionDialog(FluxSimulatorGUI.this, 
								"Hey, the specified directory does not exist!\n" +
								"Do you want to create the folder\n"+f.getAbsolutePath(), 
								"Create New Folder", 
								JOptionPane.YES_NO_OPTION, 
								JOptionPane.QUESTION_MESSAGE, 
								null, 
								null, 
								JOptionPane.YES_OPTION);
						if (ret== JOptionPane.NO_OPTION)
							return;
						if (!f.mkdir()) {
							JOptionPane.showMessageDialog(FluxSimulatorGUI.this, "Oh no, could not create folder\n"+ f.getAbsolutePath());
							return;
						}
						if ((!f.canRead())|| (!f.canWrite())) {
							JOptionPane.showMessageDialog(FluxSimulatorGUI.this, "Oh no, cannot read/write in folder\n"+ f.getAbsolutePath());
							return;
						}
					}
					
					if (!FluxSimulatorSettings.checkNewDialogOK(true, tfAnn.getText(), tfTmp.getText(), 
							new String[] {tfPar.getText(), tfPro.getText(), tfFrg.getText(), tfSeq.getText()}, true)) {
						return;
					}
					
					FluxSimulatorSettings set= settings;
					if (set== null) 
						set= FluxSimulatorSettings.createDefaults();
					
					set.setRefFile(new File(tfAnn.getText()));
					set.setName(tfName.getText());
					set.setParFile(new File(tfPar.getText()));
					set.setProFile(new File(tfPro.getText()));
					set.setFrgFile(new File(tfFrg.getText()));
					set.setSeqFile(new File(tfSeq.getText()));
					set.setTmpDir(new File(tfTmp.getText()));
//					if (!set.checkComplete()) {
//						showParfileIncomplete(set.getParFile());
//						return;
//					}					
					updateTitle(set.getParFile());
					newDialog.setVisible(false);
					if (newDialog.getTitle().equals(TITLE_PROJECT_NEW)) {
						panAnn= null;
						panXpr= null;
						rtFragPanel= null; 
						seqPanel= null;
						FluxSimulatorGUI.this.remove(tabPane);
						tabPane= null;
						getTabPane();
						//init(set); 
						uiUpdate("Initialized with default values.");
						settings.save();
						
					} else if (newDialog.getTitle().equals(TITLE_PROJECT_COPY)) {
						final File origProFile= (settings== null)?null:settings.getProFile(),
								origFrgFile= (settings== null)?null:settings.getFrgFile(),
								origSeqFile= (settings== null)?null:settings.getSeqFile(),
								origErrFile= (settings== null)?null:settings.getErrFile();
							
						blockThread= new BlockingThread(new StoppableRunnable(){
						boolean stop= false;
						public boolean isStop() {
							return stop;
						}
						public boolean setStop() {
							stop= true;
							return true;
						}
						
						public boolean setStop(boolean stop) {
							if (stop)
								return setStop();
							this.stop= stop;
							return true;
						}
						
						public void run() {
							boolean copied= false;
							if (origProFile.exists()&& !settings.getProFile().getAbsolutePath().equals(origProFile.getAbsolutePath())) { 
								try {
									labStatus.setString("Copying .pro file");
									FileHelper.fastChannelCopy(origProFile, settings.getProFile(), false);
								} catch (IOException e) {
									; // :)
								}					
								copied= true;
							}
							if (origFrgFile.exists()&& !settings.getFrgFile().getAbsolutePath().equals(origFrgFile.getAbsolutePath())) { 
								try {
									labStatus.setString("Copying .frg file");
									FileHelper.fastChannelCopy(origFrgFile, settings.getFrgFile(), false);
								} catch (IOException e) {
									; // :)
								}									
								copied= true;
							}
							if (origSeqFile.exists()&& !settings.getSeqFile().getAbsolutePath().equals(origSeqFile.getAbsolutePath())) { 
								try {
									labStatus.setString("Copying .seq file");
									FileHelper.fastChannelCopy(origSeqFile, settings.getSeqFile(), false);
								} catch (IOException e) {
									; // :)
								}									
								copied= true;
							}
							if (origErrFile!= null&& origErrFile.exists()&& !settings.getErrFile().getAbsolutePath().equals(origErrFile.getAbsolutePath())) { 
								try {
									labStatus.setString("Copying .err file");
									FileHelper.fastChannelCopy(origErrFile, settings.getErrFile(), false);
								} catch (IOException e) {
									; // :)
								}									
								copied= true;
							}
							if (copied)
								settings.save();
						}
					}, blockThread);
					blockThread.start();
					}
					init(set);
					updateTitle(settings.getParFile());
				}
			});
			buts.add(butOKnew);
			newDialog.add(buts, BorderLayout.SOUTH);
			
			newDialog.pack(); // Bug ID: 4221414
			newDialog.setSize(newDialog.getPreferredSize());
			newDialog.setResizable(false);
		}
		if (settings!= null) {
			if (settings.getParFile()!= null) {
				String path= settings.getParFile().getAbsolutePath();
				tfPar.setText(path);
				String name= settings.getParFile().getName();
				name= name.substring(0, name.lastIndexOf("."));
				tfName.setText(name);
				path= path.substring(0,path.lastIndexOf(File.separator));
				tfPDir.setText(path);
			}
			if (settings.getRefFile()!= null)
				tfAnn.setText(settings.getRefFile().getAbsolutePath());
			if (settings.getProFile()!= null)
				tfPro.setText(settings.getProFile().getAbsolutePath());
			if (settings.getFrgFile()!= null)
				tfFrg.setText(settings.getFrgFile().getAbsolutePath());
			if (settings.getSeqFile()!= null)
				tfSeq.setText(settings.getSeqFile().getAbsolutePath());
			if (settings.getTmpDir()!= null)
				tfTmp.setText(settings.getTmpDir().getAbsolutePath());
			else
				tfTmp.setText(System.getProperty("java.io.tmpdir"));
		}

		Dimension d= newDialog.getPreferredSize();
		newDialog.setLocation(FluxSimulatorGUI.this.getLocation().x+ FluxSimulatorGUI.this.getSize().width/ 2- d.width/ 2,
				FluxSimulatorGUI.this.getLocation().y+ FluxSimulatorGUI.this.getSize().height/ 2- d.height/ 2);

		return newDialog;
	}

	private Action actionTmp;	
	private Action getTmpAction() {		
		if (actionTmp == null) {
			actionTmp = new AbstractAction() {
				public void actionPerformed(ActionEvent e) {
					JFileChooser dia= getChooser();
					if (settings.getTmpDir()!= null)
						try {
							dia.setCurrentDirectory(settings.getTmpDir());
							dia.setSelectedFile(settings.getTmpDir());
						} catch (Exception e2) {
							; // :)
						}
					else 
						try {
							File f= new File(System.getProperty("java.io.tmpdir"));
							dia.setCurrentDirectory(f);
							dia.setSelectedFile(f);
						} catch (Exception e2) {
							; // :)
						}

					dia.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
					dia.setDialogTitle("Choose Temporary Folder");
					dia.removeChoosableFileFilter(getFilterPar());
					dia.removeChoosableFileFilter(getFilterGTF());
					//dia.setVisible(true);
					int rval= dia.showOpenDialog(FluxSimulatorGUI.this);
					File f= dia.getSelectedFile();	// getCurrentDirectory()
					if (f== null|| (!f.isDirectory())|| (!f.canRead())|| (!f.canWrite())|| rval!= JFileChooser.APPROVE_OPTION)
						return;
					getTFtmp().setText(f.getAbsolutePath());
				}
			};
		}
		return actionTmp;
	}
	
	private Action actionGen;
	JTextField tfGenome;
	private Action getActionGen() {
		if (actionGen == null) {
			actionGen = new AbstractAction(){
				public void actionPerformed(ActionEvent e) {
					JFileChooser dia= getChooser();
					dia.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
					dia.setDialogTitle("Choose Genome Folder");
					if (settings.getGenDir()!= null) {
						dia.setCurrentDirectory(settings.getGenDir());
						dia.setSelectedFile(settings.getGenDir());
					}
					dia.removeChoosableFileFilter(getFilterPar());
					dia.removeChoosableFileFilter(getFilterGTF());
					//dia.setVisible(true);
					int rval= dia.showOpenDialog(FluxSimulatorGUI.this);
					File f= dia.getSelectedFile();
					if (f== null|| (!f.isDirectory())|| (!f.canRead())|| rval!= JFileChooser.APPROVE_OPTION)
						return;
					if (tfGenome!= null)
						tfGenome.setText(f.getAbsolutePath());
					settings.setGenDir(f);
					if (parFile!= null){
						settings.save();
                    }
					init(settings);	// dangerous, can have null values in new dialog, Issue 25
				}
			};
			
		}

		return actionGen;
	}

	protected void uiUpdate(String msg) {
		
		Component c= getTabPane().getSelectedComponent();
		if (c instanceof ReadyOrNot&& ((ReadyOrNot) c).isReady()) {
			butRun.setEnabled(true);
		} else {
			butRun.setEnabled(false);
		}
		if (msg!= null) {
			labStatus.setString(msg);
			labStatus.setValue(-1);
			labStatus.repaint();
		}
		for (int i = 0; i < getTabPane().getComponentCount(); i++) {
			if (getTabPane().getComponentAt(i) instanceof ReadyOrNot) {
				ReadyOrNot ready= (ReadyOrNot) getTabPane().getComponentAt(i);
				if (ready.isFinished()) 
					getTabPane().setIconAt(i, imgGreen);
				else if (ready.isReady())
					getTabPane().setIconAt(i, imgRed);
				else
					getTabPane().setIconAt(i, imgGrey);
				getTabPane().getComponentAt(i).repaint();
			}
		}

		butCopy.setEnabled(settings!= null&& settings.checkComplete().length()== 0);
			
	}

	private boolean askUserOverwriteProject() {
		JOptionPane dia= new JOptionPane("You are about to replace the existing run, make sure you have saved its data or it will be lost.",
				JOptionPane.WARNING_MESSAGE);
		String[] options= new String[] {"Go ahead, kill it!", "Maybe I check, just in case.."};
		int n = JOptionPane.showOptionDialog(this,
				"You are about to replace the existing run, make sure you have saved its data or it will be lost.",
                "Warning: overwriting existing project",
                JOptionPane.YES_NO_OPTION,
                JOptionPane.WARNING_MESSAGE,
                null,
                options,
                options[1]);
		if (n== JOptionPane.YES_OPTION) {
			if (settings!= null)
				settings.save();
			return true;
		}
		return false;
	}
	
	private boolean checkProjectContainsData() {
		boolean hihaData= false;
		for (int i = 0; i < readyOrNot.length; i++) 
			hihaData|= FluxSimulatorGUI.this.readyOrNot[i].isFinished();
		return hihaData;
	}
	
	private ParFilter filterPar;
	private ParFilter getFilterPar() {
		if (filterPar == null) {
			filterPar = new ParFilter();
		}

		return filterPar;
	}
	
	Action actionLoad;
	boolean loadActive= false;
	private Action getActionLoad() {
		if (actionLoad == null) {
			actionLoad = new AbstractAction() {
				public void actionPerformed(ActionEvent e) {
					
					if (checkProjectContainsData()&& (!askUserOverwriteProject())) 
						return;					
					
					JFileChooser dia= getChooser();	// getLoadDialog()
					if (parFile!= null) 
						try {
							dia.setCurrentDirectory(parFile.getParentFile());
							dia.setSelectedFile(parFile);
						} catch (Exception e2) {
							; // :)
						}

					dia.setDialogTitle("Open Run");	//dia.setTitle()
					dia.setFileSelectionMode(JFileChooser.FILES_ONLY);
					dia.removeChoosableFileFilter(getFilterGTF());
					dia.addChoosableFileFilter(getFilterPar());
					int rval= dia.showOpenDialog(FluxSimulatorGUI.this);
					File f= dia.getSelectedFile(); // dia.getFile()
					//dia.setVisible(false);	// auf gar keinen Fall
					if (f== null|| rval!= JFileChooser.APPROVE_OPTION)
						return;
					
					settings= load(f);
					if (settings!= null)
						loadInit();
				}
			};
		}
		return actionLoad;
	}

	public void loadInit() {
		if (settings!= null)
			init(settings);
		
		StoppableRunnable[] comp= new StoppableRunnable[] {
			(StoppableRunnable) getPanelAnnotation(),
			(StoppableRunnable) getPanExpression(),
			(StoppableRunnable) getRTFragPanel(),
			(StoppableRunnable) getSeqPanel()
		};
		getTabPane().setEnabled(false);
		butNew.setEnabled(false);
		butLoad.setEnabled(false);
		butCopy.setEnabled(false);
		butCancel.setEnabled(true);
		loadActive= true;
		for (int i = 0; i < comp.length; i++) {						
			((ReadyOrNot) comp[i]).setLoadStats(true);
            Log.progressStart("Loading " + ((Component) comp[i]).getName());
			blockThread= new BlockingThread(comp[i], blockThread);
			blockThread.setForceIdx(i);
			blockThread.start();
		}
		Thread tEnable= new Thread() {
			@Override
			public void run() {
				while (blockThread.isAlive())
					try {
						blockThread.join();
					} catch (InterruptedException e) {
						; // :)
					}
				getTabPane().setEnabled(true);
				butNew.setEnabled(true);
				butLoad.setEnabled(true);
				butCopy.setEnabled(true);
				butCancel.setEnabled(false);
				loadActive= false;
				if (blockThread.target.isStop())
					init(null);
			}
		};
		tEnable.start();
		
		uiUpdate(null);
	}
	
	void init(FluxSimulatorSettings settings) {
		this.settings= settings;
		
		Component[] c= getTabPane().getComponents();
		for (int i = 0; i < c.length; i++) {
			if (c[i] instanceof ReadyOrNot)
				((ReadyOrNot) c[i]).set(settings);
		}
		getActionCopy().setEnabled(settings!= null);	// the toolbar but
		copyItem.setEnabled(settings!= null);
		
		uiUpdate(null);
	}
	
	FileDialog saveDia, loadDia;
	private FileDialog getSaveDialog() {
		if (saveDia == null) {
			JFrame invis= new JFrame();
			invis.setVisible(false);
			invis.setLocation(FluxSimulatorGUI.this.getLocation().x+ 100, FluxSimulatorGUI.this.getLocation().y+ 100);
			saveDia= new FileDialog(invis, "Save", FileDialog.SAVE);
			
		}

		return saveDia;
	}
	
	Frame invisLoad;
	private FileDialog getLoadDialog() {
		if (invisLoad == null) {
			invisLoad= new Frame();
		}
		invisLoad.setLocation(FluxSimulatorGUI.this.getLocation().x+ 100, FluxSimulatorGUI.this.getLocation().y+ 100);
		
		if (loadDia != null) { 
			loadDia.dispose();
		}
		loadDia= new FileDialog(invisLoad, "Load", FileDialog.LOAD);

		return loadDia;
	}
	
	Action actionSave;
	private Action getActionCopy() {
		if (actionSave == null) {
			actionSave = new AbstractAction() {
				public void actionPerformed(ActionEvent e) {
					JDialog dia= getProjectDialog();	// FileDialog dia= getSaveDialog();
					dia.setTitle(TITLE_PROJECT_COPY);
					dia.setVisible(true);
				}
			};
			actionSave.setEnabled(false);
		}
		return actionSave;
	}
	
	FluxSimulatorSettings settings= null;
	
	public FluxSimulatorSettings load(File f) {
		if (f== null|| f.isDirectory()|| (!f.exists())|| (!f.canRead()))
			return settings;
		settings= FluxSimulatorSettings.createSettings(f);
		String ss= "";
		if (settings!= null&& (ss= settings.checkComplete()).length()== 0) {
			FluxSimulatorSettings settingsNew= null;
			try {
				settingsNew= (FluxSimulatorSettings) settings.clone();
			} catch (CloneNotSupportedException e) {
				; //:)
			}
			String s= FluxSimulatorSettings.fillDefaults(settingsNew);
			if (s.length()> 0) {
				int res= JOptionPane.showOptionDialog(
						FluxSimulatorGUI.this,
					    "The parameter file "+ f.getName()
					    + "\nseems to be from an older version of the Flux Simulator."
					    + "\nDo you want me to add the following parameters for compatibility"+ s,
					    "Incomplete Parameter File",
					    JOptionPane.YES_NO_OPTION,
					    JOptionPane.QUESTION_MESSAGE,
					    null,
					    null,
					    null
				);
				if (res== JOptionPane.YES_OPTION) {
					settings= settingsNew;
					settings.save();
				} else
					return settings;
			}
			updateTitle(settings.getParFile());
			return settings;
		}		
		//else 
		showParfileIncomplete(f, ss);
		return null;
	}
	
	private void showParfileIncomplete(File f, String ss) {
		JOptionPane.showMessageDialog(
				this,
				"Incomplete set of parameters in file\n"+ f.getName()+ "\n"+ ss,
				"Incomplete Parameter File",
				JOptionPane.ERROR_MESSAGE,
				null
		);
	}
	
	
	JPanel panAnn;
	JTextField annTFRef;
	private JPanel getPanelAnnotation() {
		if (panAnn == null) {
			panAnn = new AnnotationReaderGUI();
		}

		return panAnn;
	}
	
	JPanel panXpr;
	private JPanel getPanExpression() {
		if (panXpr == null) {
			panXpr = new ExpressionSimulatorGUI();
			panXpr.setName(NAME_EXP);
		}

		return panXpr;
	}
	
	boolean save(File f) {		
		if ((!f.getParentFile().canWrite()))
			return false;

        OutputStream out = null;
        try {
            out = new FileOutputStream(f);
            settings.write(out);
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }finally {
            try {out.close();} catch (IOException e) {}
        }

		return  true;
	}

	private JFileChooser chooser;
	private JFileChooser getCommonChooser() {
		if (chooser == null) {
			chooser= new JFileChooser();	// not in the EDT !!!
			chooser.setMultiSelectionEnabled(false);
			chooser.setPreferredSize(new Dimension(chooser.getPreferredSize().width, 600));

			//chooser.putClientProperty("FileChooser.useShellFolder", Boolean.FALSE);
		
//			while (chooser== null) // main still busy
//				try {
//					Thread.sleep(500);
//				} catch (InterruptedException e1) {
//					; // :)
//				} 

			
//			AbstractButton button = SwingUtils.getDescendantOfType(AbstractButton.class,
//				      chooser, "Icon", UIManager.getIcon("FileChooser.detailsViewIcon"));
//			button.doClick();
			
//			File[] cbFolders = (File[]) ShellFolder.get("fileChooserComboBoxFolders");
//			for (int i = 0; i < cbFolders.length; i++) {
//				cbFolders[i].list();
//			}
//			chooser.doLayout();
		}
		return chooser;
	}
		
	private JFileChooser getChooser() {
		JFileChooser chooser= getCommonChooser();
		chooser.setLocation(FluxSimulatorGUI.this.getLocation().x+ 100, FluxSimulatorGUI.this.getLocation().y+ 100);
		chooser.setVisible(true);
		return chooser;
	}

	JMenuBar jmenuBar;
	JMenuItem copyItem, tmpItem;
	private JMenuBar createJMenuBar() {
		if (jmenuBar == null) {
			jmenuBar = new JMenuBar();
			
			JMenu menu;
			JMenuItem item;
			menu= new JMenu("Project");
			
			item= new JMenuItem("New");
			item.addActionListener(getActionNew());
			menu.add(item);
			
			item= new JMenuItem("Open");
			item.addActionListener(getActionLoad());
			menu.add(item);
			
			copyItem= new JMenuItem("Copy");
			copyItem.addActionListener(getActionCopy());
			copyItem.setEnabled(settings!= null);
			menu.add(copyItem);

			menu.add(new JSeparator());

			tmpItem= new JMenuItem("Set Temporary Folder");
			//tmpItem.setEnabled(settings!= null);
			tmpItem.addActionListener(getTmpAction());
			menu.add(tmpItem);
			
			item= new JMenuItem("Set Genome folder");
			//tmpItem.setEnabled(settings!= null);
			item.addActionListener(getActionGen());
			menu.add(item);
			
			menu.add(new JSeparator());
			
			item= new JMenuItem("Exit");
			item.addActionListener(new AbstractAction(){
				public void actionPerformed(ActionEvent e) {
					exit();
				}
			});
			menu.add(item);
			jmenuBar.add(menu);
			
			menu= new JMenu("View");
			LookAndFeelInfo[] lafs= UIManager.getInstalledLookAndFeels();			
			for (int i = 0; i < lafs.length; i++) {
				final String lafname= lafs[i].getClassName();
				item= new JMenuItem(lafs[i].getName());
				item.addActionListener(new ActionListener(){
					public void actionPerformed(ActionEvent e) {
						setColors(lafname);
						SwingUtilities.updateComponentTreeUI(FluxSimulatorGUI.this);
						SwingUtilities.updateComponentTreeUI(chooser);
						SwingUtilities.invokeLater(new Runnable() {
							public void run() {
								newDialog= null;
								tfTmp= null;
								getProjectDialog();
							}
						});						
					}
				});
				menu.add(item);
			}

//			item= new JMenuItem("Parameters");
//			menu.add(item);
			jmenuBar.add(menu);
			
			menu= new JMenu("Help");
			item= new JMenuItem("About");
			item.addActionListener(new AbstractAction(){
				public void actionPerformed(ActionEvent e) {
					getAboutScreen().setDrawing(true);
					getAboutScreen().setHandleMouseEvents(true);
				}
			});
			menu.add(item);
			jmenuBar.add(menu);

		}
	
		return jmenuBar;
	}
	
	GlassPane aboutScreen;
	private GlassPane getAboutScreen() {
		if (aboutScreen == null) {
			aboutScreen= new GlassPane();
		    aboutScreen.setGlassPane(this); 			
		}

		return aboutScreen;
	}

	GTFFilter filterGTF;
	static Image   splashImage;
	private GTFFilter getFilterGTF() {
		if (filterGTF == null) {
			filterGTF = new GTFFilter();
		}
	
		return filterGTF;
	}

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

	static long parseLong(JTextField tf, long defVal) {
		try {
			return Long.parseLong(tf.getText());
		} catch (NumberFormatException xxx) {
			tf.setText(Long.toString(defVal));
		}
		return defVal;
	}

}
