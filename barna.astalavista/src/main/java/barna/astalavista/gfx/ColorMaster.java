package barna.astalavista.gfx;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.File;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Random;

import javax.swing.JColorChooser;
import javax.swing.JFrame;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.ListSelectionModel;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import javax.swing.table.DefaultTableCellRenderer;
import javax.swing.table.TableCellRenderer;

/**
 * Created with IntelliJ IDEA.
 * User: micha
 * Date: 10/3/12
 * Time: 10:19 AM
 * To change this template use File | Settings | File Templates.
 */
public class ColorMaster {

        public static final String DEFAULT_COLOR_TABLE_FNAME= "colTable.txt";	// colorTable.txt

        public static class ColorComparator implements Comparator {
            public int compare(Object o1, Object o2) {
                int rgb1= ((Circle) ((Object[]) o1)[0]).getColor().getRGB();
                int rgb2= ((Circle) ((Object[]) o2)[0]).getColor().getRGB();

                if (rgb1< rgb2)
                    return -1;
                if (rgb1> rgb2)
                    return 1;
                return 0;
            }
        }

        static class CellComponentRenderer implements TableCellRenderer {

            public static final DefaultTableCellRenderer DEFAULT_RENDERER = new DefaultTableCellRenderer();

            public Component getTableCellRendererComponent(JTable table, Object value,
                                                           boolean isSelected, boolean hasFocus, int row, int column) {
                Component renderer = DEFAULT_RENDERER.getTableCellRendererComponent(
                        table, value, isSelected, hasFocus, row, column);
                if (value instanceof Component)
                    return ((Component) value);
                return renderer;
            }

        }

        public static Color getRandomColor() {
            Random random= new Random();
            return new Color(random.nextInt(255),random.nextInt(255),random.nextInt(255));
        }

        class SelectionListener implements ListSelectionListener {

            // It is necessary to keep the table since it is not possible
            // to determine the table from the event's source
            SelectionListener() {
            }
            public void valueChanged(ListSelectionEvent e) {
                // If cell selection is enabled, both row and column change events are fired
                if (e.getSource() == ColorMaster.this.getTable().getSelectionModel()
                        && ColorMaster.this.getTable().getRowSelectionAllowed()) {
                    // Column selection changed
                    int first = e.getFirstIndex();
                    int last = e.getLastIndex();	// works better :)
                    ColorMaster.this.getColChooser().setColor(((Circle) ColorMaster.this.getTable().getValueAt(getTable().getSelectedRow(), 0)).getColor());
                } else if (e.getSource() == ColorMaster.this.getTable().getColumnModel().getSelectionModel()
                        && ColorMaster.this.getTable().getColumnSelectionAllowed() ){
                    // Row selection changed
                    int first = e.getFirstIndex();
                    int last = e.getLastIndex();
                }

                if (e.getValueIsAdjusting()) {
                }
            }


        }

        public static void main(String[] args) {
            ColorMaster myCM= new ColorMaster(SpliceOSigner.DEFAULT_DOMAIN_COLOR_MAP);
        }

        JFrame window= null;
        JColorChooser colChooser= null;
        String fName= null;
        JTable table= null;

        public ColorMaster(String newFName) {
            this.fName= newFName;
            window= new JFrame();
            window.getContentPane().setLayout(new BorderLayout());
            window.getContentPane().add(getColChooser(), BorderLayout.NORTH);
            JScrollPane scroller= new JScrollPane(getTable());
            window.getContentPane().add(scroller, BorderLayout.SOUTH);

            window.pack();
            window.setVisible(true);
            window.addWindowListener(new WindowAdapter() {

                public void windowClosing(WindowEvent e) {
                    super.windowClosing(e);

                    HashMap<String, Color> map= new HashMap<String,Color>();
                    for (int i = 0; i < ColorMaster.this.getTable().getRowCount(); i++)
                        map.put((String) ColorMaster.this.getTable().getValueAt(i, 1),
                                ((Circle) ColorMaster.this.getTable().getValueAt(i, 0)).getColor());
                    SpliceOSigner.setColorMap(map);
                    SpliceOSigner.writeOutDomainColorMap();
                    System.out.println("wrote color map.");
                    System.exit(0);
                }
            });
        }

        JTable getTable() {
            if (table== null) {
                HashMap<String, Color> colMap=
                        SpliceOSigner.readInDomainColorMap(new File(fName));
                Object[] keys= colMap.keySet().toArray();
                Object[][] rowData= new Object[keys.length][];
                for (int i = 0; i < keys.length; i++) {
                    rowData[i]= new Object[2];
                    rowData[i][0]= new Circle(colMap.get(keys[i]), 10);
                    rowData[i][1]= keys[i];
                }
                Arrays.sort((Object[]) rowData, new ColorComparator());
                Object[] colNames= new String[] {"color", "ID"};
                table= new JTable(rowData, colNames) {
                    public boolean isCellEditable(int row, int column) {
                        return false;
                    }
                };
                table.setDefaultRenderer(Object.class, new CellComponentRenderer());

                // listener
                table.getSelectionModel().addListSelectionListener(new SelectionListener());
                table.setSelectionMode(ListSelectionModel.SINGLE_INTERVAL_SELECTION);
            }
            return table;
        }

        public JColorChooser getColChooser() {
            if (colChooser == null) {
                colChooser= new JColorChooser();
                colChooser.getSelectionModel().addChangeListener(new ChangeListener() {
                    public void stateChanged(ChangeEvent e) {
                        Color newColor = ColorMaster.this.getColChooser().getColor();
                        for (int i = 0; i < getTable().getSelectedRows().length; i++)
                            getTable().setValueAt(new Circle(newColor, 10),
                                    getTable().getSelectedRows()[i], 0);
                        getTable().repaint();
                    }
                });

            }

            return colChooser;
        }

        static HashMap<String,Color> defaultColorMap= null;

        public static HashMap<String,Color> getDefaultColorMap() {
            if (defaultColorMap == null) {
                defaultColorMap = readInColorMap();
            }

            return defaultColorMap;
        }

        public static HashMap<String,Color> readInColorMap() {
            return readInColorMap(DEFAULT_COLOR_TABLE_FNAME);
        }

        public static HashMap<String,Color> readInColorMap(String fName) {
            File inFile= new File(fName);
            HashMap<String, Color> colorMap= new HashMap<String, Color>();
            if (!inFile.exists()) {
                System.out.println("No color configuration found for domains.");
                return colorMap;
            }
            try {
                TabDelimitedFormatWrapper reader= new TabDelimitedFormatWrapper(inFile.getAbsolutePath());
                reader.read();
                String[][] tab= reader.getTable();
                for (int i = 0; tab!= null&& i < tab.length; i++) {
                    if (tab[i].length!= 4)
                        continue;
                    Color c= new Color(Integer.parseInt(tab[i][1]),
                            Integer.parseInt(tab[i][2]), Integer.parseInt(tab[i][3]));
                    colorMap.put(tab[i][0], c);
                }
            } catch (Exception e) {
                e.printStackTrace();
            }
            return colorMap;
        }


    }