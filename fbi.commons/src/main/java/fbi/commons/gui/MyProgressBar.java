package fbi.commons.gui;

import fbi.commons.Progressable;
import fbi.commons.StringConstants;

import javax.swing.*;
import java.awt.*;

public class MyProgressBar extends JPanel implements Progressable {

	JProgressBar bar;
	JLabel lab;
	
	public MyProgressBar() {
		bar= new JProgressBar();
		lab= new JLabel(" "); // size
		setLayout(new FlowLayout(FlowLayout.LEFT));
		add(lab);
		add(bar);
		bar.setVisible(false);
	}
	
	public void progress() {
		Runnable runner= new Runnable() {
			public void run() {
				bar.setValue(bar.getValue()+ 1);
				bar.setVisible(true);
				bar.paintImmediately(bar.getVisibleRect());	// getBounds()
			}
		};
		runFromEverywhere(runner);
	}

	public void clear() {
		setString(""); // was constant but hardcoded should have the same effect
		setValue(-1);
	}
	
	@Override
	public Dimension getPreferredSize() {
		int w= Math.max(lab.getPreferredSize().width, 250);
		int h= Math.max(80, bar.getPreferredSize().height);
		return super.getPreferredSize();	//new Dimension(w,h);
	}
	
	public void setMaximum(int newValue) {
		bar.setMaximum(newValue);
	}

	private void runFromEverywhere(Runnable runner) {
		if (EventQueue.isDispatchThread())
			runner.run();
		else
			try {
				SwingUtilities.invokeAndWait(runner);
			} catch (Exception e) {
				//System.err.println("everywhere disrupt");
				System.currentTimeMillis(); // :)
			}
	}
	
	public void setString(String value) {
		//System.out.println("want set text "+value);
		final String val= value;
		Runnable runner= new Runnable() {
			public void run() {
				//System.out.println("set text "+val);
				lab.setText(val);
				lab.revalidate();
				lab.paintImmediately(lab.getVisibleRect());
				bar.setVisible(false);
				paintImmediately(getBounds());
			}
		};
		runFromEverywhere(runner);
	}

	public void setValue(final int newValue) {
		//System.out.println("want set val "+newValue);
		Runnable runner= new Runnable() {
			public void run() {
				//System.out.println("set val "+newValue);
				if (newValue< 0)
					bar.setVisible(false);
				else {
					bar.setValue(newValue);
					bar.setVisible(true);
					bar.paintImmediately(bar.getVisibleRect());
				}
			}
		};
		
		runFromEverywhere(runner);
	}

	public void finish() {
		bar.setVisible(false);
		bar.setValue(bar.getMinimum());
		lab.setText(StringConstants.SPACE);
	}

    public void finish(String msg, boolean time) {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    public void failed(String msg) {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    public int steps() {
        return 0;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public int currentStep() {
        return 0;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public void message(String value) {
		final String val= value;
		Runnable runner= new Runnable() {
			public void run() {
				lab.setText(val);
				lab.paintImmediately(lab.getBounds());
			}
		};
		runFromEverywhere(runner);
	}

	public void finish(String msg, long time) {
		// TODO Auto-generated method stub
	}

    public void start(String message) {
        throw new UnsupportedOperationException("Not implemented yet");
    }
}
