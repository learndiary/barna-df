package fbi.genome.sequencing.rnaseq.reconstruction.gui;

import fbi.genome.sequencing.rnaseq.reconstruction.LParser;
import fbi.genome.sequencing.rnaseq.reconstruction.LParser.LP;

import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Rectangle;
import java.io.File;
import java.util.Arrays;

import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.SwingConstants;
import javax.swing.UIManager;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import javax.swing.plaf.ColorUIResource;
import javax.swing.table.DefaultTableCellRenderer;
import javax.swing.table.DefaultTableModel;
import javax.swing.table.JTableHeader;
import javax.swing.table.TableCellRenderer;
import javax.swing.table.TableColumn;
import javax.swing.table.TableModel;

public class LPVizualizer extends JPanel {

	public static final Dimension LTABLE_CELL_DIM= new Dimension(20,20);
	
	public static void main(String[] args) {
		
		File myZip= null;		
		try {
			myZip= new File("c:\\workspace\\Genome\\resources\\formats\\CME_W1_CI.8-PE_lp_mod.zip");
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		LParser myParser= new LParser(myZip);
		String name= myParser.getNames()[0];
		System.err.println("Reading "+ name);
		LP myLP= myParser.parseEntry(name);
		
		System.err.println(myLP);
		LPVizualizer tableLP= new LPVizualizer();
		JFrame myFrame= new JFrame();
		myFrame.getContentPane().add(tableLP);
		myFrame.pack();
		myFrame.setVisible(true);
	}

	public static class LineLabel extends JLabel {
		public LineLabel() {
			super();
			setOpaque(true);
		}
		
		@Override
		protected void paintComponent(Graphics g) {
			super.paintComponent(g);
			
			g.setColor(getBackground());
			g.fillRect(0,0,getWidth(),getHeight());
			int y= getHeight()/ 2;
			g.setColor(Color.black);
			g.drawRect(0, y, getWidth(), 0);
		}
	}
	
	
		public class SelectionListener implements ListSelectionListener {
	        JTable table;
	    
	        // It is necessary to keep the table since it is not possible
	        // to determine the table from the event's source
	        SelectionListener(JTable table) {
	            this.table = table;
	        }
	        public void valueChanged(ListSelectionEvent e) {
	            // If cell selection is enabled, both row and column change events are fired
	            if (e.getSource() == table.getSelectionModel()
	                  && table.getRowSelectionAllowed()) {
	                // Column selection changed
	                int first = e.getFirstIndex();
	                int last = e.getLastIndex();
	            } else if (e.getSource() == table.getColumnModel().getSelectionModel()
	                   && table.getColumnSelectionAllowed() ){
	                // Row selection changed
	                int first = e.getFirstIndex();
	                int last = e.getLastIndex();
	            }
	    
	            if (e.getValueIsAdjusting()) {
	                // The mouse button has not yet been released
	            }
	        }
	    }

	class RTableRenderer extends DefaultTableCellRenderer implements TableCellRenderer {
		
		JLabel labHeader;
		
		RTableRenderer(JTable table) {
			labHeader= new JLabel();
		    JTableHeader header = table.getTableHeader();
		    labHeader.setOpaque(true);
		    labHeader.setBorder(UIManager.getBorder("TableHeader.cellBorder"));
		    labHeader.setHorizontalAlignment(CENTER);
		    labHeader.setForeground(header.getForeground());
		    labHeader.setBackground(header.getBackground());
		    labHeader.setFont(header.getFont());
		}
		
		@Override
		public Component getTableCellRendererComponent(JTable table,
				Object value, boolean isSelected, boolean hasFocus, int row,
				int column) {
			
			if (column== 0|| row== -1) {
				labHeader.setToolTipText(null);
				if (row== -1) {
					if (column== 0)
						labHeader.setText(value.toString());
					else
						labHeader.setText(Integer.toString(column));
					labHeader.setToolTipText((String) value);
				} else {
					labHeader.setText((value == null)?"" :value.toString());
					labHeader.setToolTipText(Double.toString(sums[row]));
				}
				return labHeader;
			}
			
			Component supa= super.getTableCellRendererComponent(table, value, isSelected, hasFocus, row, column);
			if (value instanceof Double) {
				double val= (Double) value;
				double ratio= Math.abs(val)/ sums[row];
				Color c1= Color.white, c2= Color.red;
				if (val< 0)
					c2= Color.blue;
				
				int red = (int)(c2.getRed() * ratio + c1.getRed() * (1 - ratio));
				int green = (int)(c2.getGreen() * ratio +
				                    c1.getGreen() * (1 - ratio));
				int blue = (int)(c2.getBlue() * ratio +
				                   c1.getBlue() * (1 - ratio));
				Color c = new Color(red, green, blue);
				supa.setBackground(c);
				JLabel label= (JLabel) supa;
				label.setText(" ");
				label.setToolTipText(value.toString());
				if (column== table.getSelectedColumn()) {
					if (row== 0)
						label.setBorder(BorderFactory.createMatteBorder(1, 1, 0, 1, Color.black));
					else if (row== table.getRowCount()- 1)
						label.setBorder(BorderFactory.createMatteBorder(0, 1, 1, 1, Color.black));
					else
						label.setBorder(BorderFactory.createMatteBorder(0, 1, 0, 1, Color.black));
				}
				supa.setPreferredSize(new Dimension(20, 20));
			} else {
				supa.setBackground(null);
				supa.setPreferredSize(new Dimension(50, 20));
			}
			
			//supa.setText((value == null)?"" :value.toString());
			return supa;
		}
	}
	class LTableRenderer extends DefaultTableCellRenderer implements TableCellRenderer {
		
		JLabel labHeader= null;
		
		LTableRenderer(JTable table) {
			labHeader= new JLabel();
		    JTableHeader header = table.getTableHeader();
		    labHeader.setOpaque(true);
		    labHeader.setBorder(UIManager.getBorder("TableHeader.cellBorder"));
		    labHeader.setHorizontalAlignment(CENTER);
		    labHeader.setForeground(header.getForeground());
		    labHeader.setBackground(header.getBackground());
		    labHeader.setFont(header.getFont());
		  }
		
		
		public Component getTableCellRendererComponent(JTable table,
				Object value, boolean isSelected, boolean hasFocus, int row,
				int column) {
			
			//Component supa= super.getTableCellRendererComponent(table, value, isSelected, hasFocus, row, column);
			JLabel supaLabel= new JLabel();
			supaLabel.setOpaque(true);
//			if (supa instanceof JLabel)
//				supaLabel= (JLabel) supa;
			if (value instanceof String)
				supaLabel.setText(value.toString());
			
			supaLabel.setBackground((ColorUIResource) UIManager.get("Table.background"));
			if (row== table.getSelectedRow())
				supaLabel.setBackground((ColorUIResource) UIManager.get("Table.selectionBackground"));
			//supa.setSize(LTABLE_CELL_DIM);
			//((JLabel) supa).setBorder(BorderFactory.createEmptyBorder());
			
			if (column== 0|| row== -1) {
				labHeader.setToolTipText(null);
				if (row== -1) {
					if (column== 0)
						labHeader.setText(value.toString());
					else
						labHeader.setText(Integer.toString(column));
					labHeader.setToolTipText((String) value);
					return labHeader;
				} else {
					//labHeader.setText((value == null)?"" :value.toString());
					supaLabel.setHorizontalAlignment(SwingConstants.RIGHT);
					supaLabel.setToolTipText(Double.toString(lp.getTranscriptExp(value.toString())));
					return supaLabel;
				}
				
			}
			
			
			TableModel model= table.getModel();
			
			if (value instanceof Integer) {
				supaLabel.setText("");
				Integer val= (Integer) value;
				if (val== 1) {				
					int left= 0, right= 0;
					if (column== 1) 
						left= 1;
					else {
						Integer leftVal= (Integer) model.getValueAt(row, column- 1);
						if (leftVal!= 1)
							left= 1;
					}
					if (column== model.getColumnCount()- 1) 
						right= 1;
					else {
						Integer rightVal= (Integer) model.getValueAt(row, column+ 1);
						if (rightVal!= 1)
							right= 1;
					}
					supaLabel.setBorder(BorderFactory.createMatteBorder(1, left, 1, right, Color.black));
					
					if (ltabSelCol!= null) {
						Color c= supaLabel.getBackground();
						if (ltabSelCol[column]> 0) 
							c= (ColorUIResource) UIManager.get("Table.selectionForeground");
						if (ltabSelCol[column]> 1)
							c= (ColorUIResource) UIManager.get("Table.focusCellBackground");
						supaLabel.setBackground(c);
					}

				} else if (val== 2) {
					LineLabel lab= new LineLabel();
					lab.setBackground(supaLabel.getBackground());
					return lab;
				}
			}
			
			
			//supa.setText((value == null)?"" :value.toString());
			return supaLabel;
		}
	}

	private JTable tabL= null;
	double[] sums;
	private DefaultTableModel getLTableModel(LP lp) {
		
		Object[] head= new Object[lp.coords.length];
		head[0]= "Segment";
		for (int i = 1; i < head.length; i++) 
			head[i]= lp.coords[i- 1]+ ".."+ lp.coords[i];
		
		Object[][] oo= new Object[lp.getTranscriptIDs().length][lp.coords.length];
		for (int i = 0; i < oo.length; i++) {
			oo[i][0]= lp.getTranscriptIDs()[i];
			int cnt= 0, val= -1;
			for (int j = 1; j < oo[i].length; j++) {
				char type= lp.types[j-1];
				int lo= lp.coords[j- 1];
				int idx= Arrays.binarySearch(lp.tstruct[i], lo);
				if (idx< 0) {
					oo[i][j]= val;
					continue;
				}
				val= -1;
				if (type== '['|| type== '-')
					val= 1;
				else if (type== '^')
					val= 2;
				oo[i][j]= val;
			}
		}
		
/*			for (int i = 0; i < lp.getTranscriptIDs().length; i++) { 
			
			oo[i][0]= lp.getTranscriptIDs()[i];
			for (int j = 1; j < oo[i].length; j++) {
				boolean b= lp.transcriptPresent(lp.getTranscriptIDs()[i], new int[] {lp.coords[j- 1], lp.coords[j]});
				oo[i][j]=(double) (b?1:2);	// 0..1 exon, 2 out/intron
			}
			for (int j = 1; j < oo[i].length; j++) {
				if (((Double) oo[i][j])>=0&& ((Double) oo[i][j])<= 1)
					break;
				oo[i][j]= -1; // out
			}
			for (int j = oo[i].length- 1; j > 0; --j) {
				if (((Double) oo[i][j])>=0&& ((Double) oo[i][j])<= 1)
					break;
				oo[i][j]= -1; // out
			}
		}
*/			
		DefaultTableModel model= new DefaultTableModel(oo, head){
			@Override
			public boolean isCellEditable(int row, int column) {
				return false;
			}
		};
		
		return model;
	}
	
	
	LTableRenderer defaultLTableCellRenderer;
	RTableRenderer defaultRTableCellRenderer;
	JTable getLTable() {
	
		if (tabL == null) {
			
			DefaultTableModel model= null;
			if (lp!= null)
				model= getLTableModel(lp);
			tabL= new JTable(model);
			tabL.setOpaque(true);
			tabL.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
			tabL.setRowSelectionAllowed(false);
			tabL.setColumnSelectionAllowed(false);
			tabL.setShowGrid(false);
			tabL.setIntercellSpacing(new Dimension(0, 2)); 
			defaultLTableCellRenderer= new LTableRenderer(tabL);
			tabL.setDefaultRenderer(Object.class, defaultLTableCellRenderer);
			
		}

		return tabL;
		
	}

	
	private JTable tabR= null;
	int[][] chains;
	private DefaultTableModel getRTableModel(LP lp) {
		
		Object[] head= new Object[lp.m.length+ 1];
		head[0]= "Restriction";
		for (int i = 1; i < head.length; i++) { 
			StringBuffer sb= new StringBuffer();
			for (int j = 0; j < lp.getChains()[i-1].length; j++) 
				sb.append(Integer.toString(lp.getChains()[i-1][j])
						+ lp.getSiteSymbol(lp.getChains()[i-1][j]));
			head[i]= sb.toString();
		}
		
		Object[][] oo= new Object[lp.getTranscriptIDs().length+ 2][lp.m.length+ 1];
		chains= lp.getChains();
		if (chains.length!= lp.m.length)
			System.currentTimeMillis();

		oo[0][0]= "observation";
		oo[1][0]= "correction";
		for (int i = 2; i < oo.length; i++) 
			oo[i][0]= lp.getTranscriptIDs()[i- 2];
		
		for (int i = 0; i < chains.length; i++) {
			double[] restr= lp.getRestriction(chains[i]);
			if (restr== null) {
				for (int j = 1; j < oo.length; j++) 
					oo[j][i+1]= 0;
			} else {
				oo[0][i+1]= restr[0];	// obs
				double corr= lp.getCorrection(chains[i]);
				oo[1][i+1]= corr;
				for (int j = 0; j < lp.getTranscriptIDs().length; ++j) 
					oo[oo.length- 1- j][i+1]= restr[restr.length- 1- j];
			}
		}
		
		sums= new double[oo.length];
		for (int i = 0; i < oo.length; i++) 
			for (int j = 1; j < oo[i].length; j++) 
				sums[i]+= Math.abs((Double) oo[i][j]);

		DefaultTableModel model= new DefaultTableModel(oo, head){
			@Override
			public boolean isCellEditable(int row, int column) {
				return false;
			}
		};
		
		return model;
	}
	
	JTable getRTable() {
		if (tabR == null) {
			
			DefaultTableModel model= null;
			if (this.lp!= null)
				model= getRTableModel(this.lp);
			tabR= new JTable(model){
				@Override
				public Dimension getPreferredSize() {
					// TODO Auto-generated method stub
					return super.getPreferredSize();
				}
				
				
			};
			tabR.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
			tabR.setRowSelectionAllowed(false);
			tabR.setColumnSelectionAllowed(true);
//			JTextArea area1 = (JTextArea) tabR.getEditorComponent();
//			area1.setEditable(false);
			defaultRTableCellRenderer= new RTableRenderer(tabR);
			tabR.setDefaultRenderer(Object.class, defaultRTableCellRenderer);

			ListSelectionListener listener= new ListSelectionListener() {
		        public void valueChanged(ListSelectionEvent e) {
		            if (e.getValueIsAdjusting()) 
		            	return;
		            if (e.getSource() == tabR.getColumnModel().getSelectionModel() 
		                   && tabR.getColumnSelectionAllowed() ){
		            	int idx= tabR.getSelectedColumn();
						if (idx> 0&& idx< tabR.getColumnCount())
							setLocusRegions(lp.getChains()[idx-1]);
						else if (idx== 0)
							getLTable().scrollRectToVisible(new Rectangle(0,0,1,1));
		            }
		    
		        }
		    };
			tabR.getSelectionModel().addListSelectionListener(listener);
		    tabR.getColumnModel().getSelectionModel().addListSelectionListener(listener);
		}

		return tabR;
	}
	
	int[] ltabSelCol;
	void setLocusRegions(int[] is) {
//		getLTable().removeColumnSelectionInterval(0, getLTable().getColumnCount()- 1);
//		getLTable().repaint();
			
		for (int i = 0; i < ltabSelCol.length; i++) 
			ltabSelCol[i]= 0;
		for (int i = 0; i < is.length; i+= 2) {
			int lo= Arrays.binarySearch(lp.coords, is[i]);
			int hi= Arrays.binarySearch(lp.coords, is[i+1]);
			//getLTable().addColumnSelectionInterval(lo+ 1, hi);
			for (int j = lo+ 1; j <= hi; j++) 
				++ltabSelCol[j];
			int row= getLTable().getSelectedRow();
			Rectangle rectLo= getLTable().getCellRect(row, lo+ 1, true),
					rectHi= getLTable().getCellRect(row, hi+ 1, true);
			Rectangle rectot= new Rectangle(rectLo.x, rectLo.y,
					rectHi.x+ rectHi.width- rectLo.x,
					rectLo.height);
			//getLTable().scrollRectToVisible(rectHi);	// funny, works incredibly fine
			getLTable().scrollRectToVisible(rectot);

		}
		
		getLTable().repaint();
		
	}
	
	
	LP lp;
	LParser parser= null;
	JScrollPane scrollerLocus, scrollerRestrictions;
	public LPVizualizer() {
		
		BoxLayout boxLayout= new BoxLayout(this, BoxLayout.Y_AXIS);
		setLayout(boxLayout);

		scrollerLocus= new JScrollPane(getLTable());
		add(scrollerLocus);
		
		scrollerRestrictions= new JScrollPane(getRTable());
		add(scrollerRestrictions);
		
		add(Box.createVerticalGlue());
	}
	
	public void setBase(File base) {
		this.parser= new LParser(base);
	}
	
	public LP parseEntry(String name) {
		
		this.lp= parser.parseEntry(name);
		
		getLTable().setModel(getLTableModel(this.lp));
		ltabSelCol= new int[getLTable().getColumnCount()];
		TableColumn column = null;
		for (int i = 0; i < getLTable().getColumnCount(); i++) {
		    column = getLTable().getColumnModel().getColumn(i);
		    column.setResizable(false);
		    column.setHeaderRenderer(defaultLTableCellRenderer);
		    if (i == 0) {
		        column.setPreferredWidth(100);
		    } else {
		        column.setPreferredWidth(20);
		    }
		}
		int cols= getLTable().getColumnModel().getTotalColumnWidth();
		int rows = getLTable().getRowHeight()* (getLTable().getModel().getRowCount()+ 1)+ 5;	// 1 4hdr, 5 extra
		Dimension d = new Dimension(cols, rows);
		scrollerLocus.setMaximumSize(d);
		getLTable().setPreferredScrollableViewportSize(d);		
		//getLTable().doLayout();
		
		
		getRTable().setModel(getRTableModel(this.lp));		
		column = null;
		for (int i = 0; i < getRTable().getColumnCount(); i++) {
		    column = getRTable().getColumnModel().getColumn(i);
		    column.setResizable(false);
		    column.setHeaderRenderer(defaultRTableCellRenderer);
		    if (i == 0) {
		        column.setPreferredWidth(100);
		    } else {
		        column.setPreferredWidth(20);
		    }
		}
		cols= getRTable().getColumnModel().getTotalColumnWidth();
		rows = getRTable().getRowHeight()* (getRTable().getModel().getRowCount()+ 1)+ 5;	// 1 4hdr, 5 extra
		d = new Dimension(cols, rows);
		scrollerRestrictions.setMaximumSize(d);
		getRTable().setPreferredScrollableViewportSize(d);		
		getRTable().doLayout();
		
		doLayout();
		getLayout().layoutContainer(this);
		
		return this.lp;
	}
	
	
}
