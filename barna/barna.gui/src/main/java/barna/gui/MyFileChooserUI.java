/*
 * This file is part of the Flux Library.
 *
 * The code of the Flux Library may be freely distributed and modified under the terms of the
 * European Union Public Licence (EUPL) published on the web site <http://www.osor.eu/eupl/european-union-public-licence-eupl-v.1.1>.
 * Copyright for the code is held jointly by the individual authors, who should be listed
 * in @author doc comments. According to Article 5 and Article 11 of the EUPL, publications that
 * include results produced by the Flux Library are liable to reference the Work,
 * see the Flux Library homepage <http://flux.sammeth.net> for more information.
 */

package barna.gui;

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

