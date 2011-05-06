package fbi.commons.gui;

import javax.swing.*;
import javax.swing.plaf.metal.MetalFileChooserUI;
import java.awt.*;

/**
 * @author micha
 * @deprecated never used
 */
public class MyFileChooserUI extends MetalFileChooserUI {
    protected class MyFileRenderer extends FileRenderer {
        ;
    }

    MyFileChooserUI(JFileChooser fileChooser) {
        super(fileChooser);
    }

    protected JPanel createList(JFileChooser fc) {
        JPanel p = new JPanel(new BorderLayout());
        JList list = new JList();
        list.setCellRenderer(new MyFileRenderer());

        if (fc.isMultiSelectionEnabled()) {
            list.setSelectionMode(ListSelectionModel.MULTIPLE_INTERVAL_SELECTION);
        } else {
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

