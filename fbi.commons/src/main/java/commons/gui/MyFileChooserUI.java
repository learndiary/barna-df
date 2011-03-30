package commons.gui;

import java.awt.BorderLayout;

import javax.swing.JFileChooser;
import javax.swing.JList;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.ListSelectionModel;
import javax.swing.plaf.metal.MetalFileChooserUI;

/**
 * @deprecated never used
 * @author micha
 *
 */
public class MyFileChooserUI extends MetalFileChooserUI
{
    protected class MyFileRenderer extends FileRenderer {
       ;
    }
	MyFileChooserUI(JFileChooser fileChooser)
	  {
	    super(fileChooser);
	  }
	
  protected JPanel createList(JFileChooser fc) {
    JPanel p = new JPanel(new BorderLayout());
    JList list = new JList();
    list.setCellRenderer(new MyFileRenderer());
 
    if(fc.isMultiSelectionEnabled()) 
    {
      list.setSelectionMode(ListSelectionModel.MULTIPLE_INTERVAL_SELECTION);
    }
    else
    {
      list.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
    }
 
    list.setModel(getModel());
    list.addListSelectionListener(createListSelectionListener(fc));
    list.addMouseListener(createDoubleClickListener(fc, list));
    // Just not add the single click listener to the list
    //list.addMouseListener(createSingleClickListener(fc, list));
 
    JScrollPane scrollpane = new JScrollPane(list);
    p.add(scrollpane, BorderLayout.CENTER);
    return p;
  }
}

