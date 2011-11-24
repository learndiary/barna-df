package fbi.genome.sequencing.rnaseq.simulation.gui;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Container;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;

import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;
import javax.swing.event.ChangeListener;

// TODO improve preferred size
public class IncDecField extends JPanel {

	public static void main(String[] args) {
		IncDecField plusMinus= new IncDecField();
		JFrame aFrame= new JFrame();
		aFrame.getContentPane().setLayout(new BorderLayout());
		aFrame.getContentPane().add(plusMinus, BorderLayout.CENTER);
		aFrame.getContentPane().add(new JPanel(), BorderLayout.NORTH);
		aFrame.getContentPane().add(new JPanel(), BorderLayout.SOUTH);
		aFrame.getContentPane().add(new JPanel(), BorderLayout.WEST);
		aFrame.getContentPane().add(new JPanel(), BorderLayout.EAST);
		aFrame.addWindowListener(new WindowAdapter(){
			@Override
			public void windowClosing(WindowEvent e) {
				System.exit(0);
			}
		});
		aFrame.setSize(new Dimension(200,100));
		aFrame.setVisible(true);
	}
	
	JLabel label;
	JTextField valField;
	JButton decBut, incBut; 
	int minVal= 0, maxVal= 100, val= 5, delta= 1;

	UpdateRunnable runner;
	Object target= null;
	Method updateMethod= null;
	
	class UpdateRunnable implements Runnable {
		public boolean active= true;
		public void run() {
			//System.out.println("run "+active);
			int nr= Integer.parseInt(valField.getText());
			while (active) {
				nr+= delta;
				//System.out.println(minVal+"<"+nr+"<"+maxVal);
				if (nr< minVal|| nr> maxVal) {
					nr-= delta;
					continue;
				}
				
				valField.setText(Integer.toString(nr));
				if (!valField.getSize().equals(valField.getPreferredSize())) {
					valField.setSize(valField.getPreferredSize());
					IncDecField.this.doLayout();
				}
				try {
					Thread.sleep(10);
				} catch (InterruptedException e) {
					e.printStackTrace();
				}
			}
			//System.out.println("run "+active);
			IncDecField.this.setVal(nr);
			Container parent= IncDecField.this;
			while(parent!= null) {
				parent.doLayout();
				parent= parent.getParent();
			}
		};
	}

	public IncDecField() {
		
		super();
		setBorder(BorderFactory.createLineBorder(Color.black));
		setLayout(new FlowLayout(FlowLayout.CENTER));
		add(getLabel());
		add(getDecBut());
		add(getValField());
		add(getIncBut());
		
		setEnabled(false);
	}

//	@Override
//	public Dimension getPreferredSize() {
//		int x= getDecBut().getPreferredSize().width
//			+getValField().getPreferredSize().width
//			+getIncBut().getPreferredSize().width;
//		int y= Math.max(Math.max(getDecBut().getPreferredSize().height, 
//				getIncBut().getPreferredSize().height), 
//				getValField().getPreferredSize().height);
//		return new Dimension(x,y);
//	}
	
	public void setEnabled(boolean enabled) {
		getDecBut().setEnabled(enabled);
		getIncBut().setEnabled(enabled);
		getValField().setEnabled(enabled);
	}
	
	public synchronized int getMinVal() {
		return minVal;
	}

	public synchronized void setMinVal(int minVal) {
		this.minVal = minVal;
	}

	public synchronized int getMaxVal() {
		return maxVal;
	}

	public synchronized void setMaxVal(int maxVal) {
		this.maxVal = maxVal;
	}

	public synchronized Object getTarget() {
		return target;
	}

	public synchronized void setTarget(Object target) {
		this.target = target;
	}

	public synchronized Method getUpdateMethod() {
		return updateMethod;
	}

	public synchronized void setUpdateMethod(Method updateMethod) {
		this.updateMethod = updateMethod;
	}

	public synchronized int getNr() {
		return val;
	}

	public synchronized void setNr(int nr) {
		this.val = nr;
	}

	public synchronized int getVal() {
		return val;
	}

	public synchronized void setVal(int val) {
		this.val = val;
		getValField().setText(Integer.toString(val));
		getValField().setSize(getValField().getPreferredSize());
		doLayout();
		if (getParent()!= null)
			getParent().doLayout(); 
	}

	private synchronized JTextField getValField() {
		if (valField == null) {
			valField= new JTextField(Integer.toString(val));
			valField.addActionListener(new ActionListener(){
				//@Override
				public void actionPerformed(ActionEvent e) {
					//System.out.println("action");
					val= maxVal;
					try {
						val= Integer.parseInt(valField.getText());
					} catch (Exception ex) {
						; // :)
					}
					if (val== maxVal)
						valField.setText(Integer.toString(val));
					if (updateMethod!= null&& target!= null)
						try {
							updateMethod.invoke(target, new Object[]{new Integer(val)});
						} catch (IllegalArgumentException e1) {
							e1.printStackTrace();
						} catch (IllegalAccessException e1) {
							e1.printStackTrace();
						} catch (InvocationTargetException e1) {
							e1.printStackTrace();
						}
					if (!valField.getSize().equals(valField.getPreferredSize()))
						doLayout();
				}
			});
			if (!valField.getSize().equals(valField.getPreferredSize()))
				doLayout();
		}

		return valField;
	}

	private JButton getIncBut() {
		if (incBut == null) {
			incBut= new JButton(">>");
			incBut.setFont(new Font(incBut.getFont().getName(), incBut.getFont().getStyle(), incBut.getFont().getSize()/ 2));
			incBut.addChangeListener(new ChangeListener() {
				//@Override
				public void stateChanged(javax.swing.event.ChangeEvent e) {
					if (incBut.getModel().isPressed()) {
						if (runner== null) {
							runner= new UpdateRunnable();
							delta= Math.abs(1);
							Thread t= new Thread(runner);
							t.start();
						}
					} else {
						if (runner!= null) 
							runner.active= false;					
						runner= null;
						if (updateMethod!= null&& target!= null)
							try {
								updateMethod.invoke(target, new Object[]{new Integer(val)});
							} catch (IllegalArgumentException e1) {
								e1.printStackTrace();
							} catch (IllegalAccessException e1) {
								e1.printStackTrace();
							} catch (InvocationTargetException e1) {
								e1.printStackTrace();
							}
					}
				}
			});		
		}

		return incBut;
	}
	
	private JButton getDecBut() {
		
		if (decBut == null) {
			decBut= new JButton("<<");
			//downNr.setSize(new Dimension(downNr.getPreferredSize().width/ 2, downNr.getPreferredSize().height/2));
			decBut.setFont(new Font(decBut.getFont().getName(), decBut.getFont().getStyle(), decBut.getFont().getSize()/ 2));
			decBut.addChangeListener(new ChangeListener() {
				//@Override
				public void stateChanged(javax.swing.event.ChangeEvent e) {
					if (decBut.getModel().isPressed()) {
						//System.out.println("on");
						if (runner== null) {
							runner= new UpdateRunnable();
							delta= -Math.abs(delta);
							Thread t= new Thread(runner);
							t.start();
						}
					} else {
						//System.out.println("off");
						if (runner!= null) 
							runner.active= false;					
						runner= null;
						if (updateMethod!= null&& target!= null)
							try {
								updateMethod.invoke(target, new Object[]{new Integer(val)});
							} catch (IllegalArgumentException e1) {
								e1.printStackTrace();
							} catch (IllegalAccessException e1) {
								e1.printStackTrace();
							} catch (InvocationTargetException e1) {
								e1.printStackTrace();
							}
					}
				}
			});			
		}

		return decBut;
	}

	private JLabel getLabel() {
		if (label == null) {
			label = new JLabel();
		}

		return label;
	}
	
	public void setText(String labelText) {
		getLabel().setText(labelText);
	}
}
