/*
 * Copyright (c) 2010, Micha Sammeth
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * The names of its contributors may be not used to endorse or promote
 *       products derived from this software without specific prior written
 *       permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL MICHA SAMMETH BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

/**
 * 
 */

import javax.swing.*;
import javax.swing.text.JTextComponent;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.*;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.ResultSet;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Vector;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * @author Colin Kingswood
 *
 */
public class DdfGui implements ActionListener {

	
	String technology; // the technology used by this ddf (only gonna be one) 
	String outputDir = "/users/rg/colin/Encode/workdir" ; // the string relative to the tar path
	String fileStartDirectory  = "/home/colin/disk2" ; 
	
	String dafFn ;  
	String ddfFn ; 
	String daf_dir = "/users/rg/colin/Encode/DDF_DAFs" ; 
	
	List views ; // a list of the views contained in the daf
	HashMap<String, String[]> hm ; // used to store database tablenames, to a list of contents 
	HashMap<File , JComponent> buttonRows ; 
	
	GridLayout gridLayout ; 
	// some components that will not change
	JFrame frame ; 
	JButton startDirButton ; 
	JButton chooseFilesButton ; 
	JButton getDafButton ;
	JButton createDdfButton ; 
	JCheckBox updateAllBox ; 
	JLabel statusLabel; 
	// a file chooser
	final JFileChooser fc = new JFileChooser();
	
	Vector componentList = new Vector() ; 
	int rowWidgetsStart ; 
	
	String ddfLocation ; 
	
	private String[] ddfColumns ; // this will store the columns used in the ddf file
	// two dimensional_array, containing the field name at the top of the daf, and the corresponding database table if available
	public static String[][] ddf_fields = {	
									{"files","FILE"} , 
									{"view", "DAF"} ,  
									{"cell", "cell_type"}, 
									{"localization","localization"},
									{"rnaExtract", "rna_type"},
									{"accession" ,  "TEXTFIELD"}, 
									{"labVersion", "TEXTFIELD"} , 
									{"softwareVersion" , "TEXTFIELD"},
									{"fragSize" , "TEXTFIELD"} , 
									{"replicate" , "TEXTFIELD"} , 
									{"DELETE" , "DELETE"}
								}; 
	
	
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		 //Schedule a job for the event-dispatching thread:
        //creating and showing this application's GUI.
        javax.swing.SwingUtilities.invokeLater(new Runnable() {
            public void run() {
                DdfGui ddf = new DdfGui() ; 
            	ddf.createAndShowGUI();
            }
        });		
	}

	
	
	
	/**
	 *  Start the application off
	 */
	private void createAndShowGUI() {
		
		connectToDatabase(); 
		gridLayout = new GridLayout(0, 8 , 12, 12);
		
		//Create and set up the window.
        frame = new JFrame("Colly's dff creator");
        frame.setMaximizedBounds(new Rectangle(1500,2500)); 
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.getContentPane().setLayout(gridLayout) ;
        
        JLabel label = new JLabel("ddf generator");
        addToGui(label);

        startDirButton = new JButton("output directory") ; 
        startDirButton.addActionListener(this); 
        addToGui(startDirButton) ; 
        
        getDafButton = new JButton("get daf / ddf") ;
        getDafButton.addActionListener(this) ; 
        addToGui(getDafButton) ; 
        
        chooseFilesButton = new JButton("Choose files") ;
        chooseFilesButton.addActionListener(this) ; 
        chooseFilesButton.setEnabled(false) ; 
        addToGui(chooseFilesButton) ; 
        
        createDdfButton  =new JButton("Generate ddf") ; 
        createDdfButton.addActionListener(this) ; 
        createDdfButton.setEnabled(false) ; 
        addToGui(createDdfButton) ;
        
        updateAllBox = new JCheckBox("Update all rows"); 
        updateAllBox.setSelected(true);
        updateAllBox.addActionListener(this);
        addToGui(updateAllBox) ; 
        
        statusLabel = new JLabel("Choose a daf"); 
        statusLabel.setOpaque(true);
        statusLabel.setHorizontalAlignment(JLabel.CENTER) ; 
        statusLabel.setBackground(Color.WHITE) ;
        statusLabel.setBorder(BorderFactory.createLineBorder(Color.BLACK)) ; 
        addToGui(statusLabel) ; 
        
        //Display the window.
        frame.pack();
        frame.setVisible(true);
        System.out.println(gridLayout.getColumns() + " x " + gridLayout.getRows() ) ;  
    }
	
	// adds a row to the gui for each file added
	private void addRows(File[] files){

		for (File file : files){
           	System.out.println(file.getPath());
            addWidget("DELETE" , file) ; 
           	for(String column : ddfColumns){           		
           		addWidget(column, file) ; 
           	}
         }
		System.out.println(gridLayout.getColumns() + " x " + gridLayout.getRows() ) ;  
		createDdfButton.setEnabled(true) ; 
		frame.pack(); 
		
		Component[] components = frame.getContentPane().getComponents(); 
		for(Component comp : components){
			System.out.println() ; 
		}
	}
	
	Object currentSource= null;
	private void addWidget(String column, File file){
		String type = "" ;
		for(String[] pair : ddf_fields){
			if( pair[0].equals(column) ){
				type = pair[1] ; 
			}
		}
//		if (column.equals("DELETE")) type = "DELETE" ; 
		if (type.length() == 0 )	System.exit(1) ; // shouldn't exit, but incase we do then it shows I have messed up this bit
		
		if (type.equals("DELETE")){
			JButton deleteButton = new JButton("remove") ; 
			deleteButton.addActionListener(this); 
			addToGui(deleteButton) ; 
		}
		else if (type.equals("TEXTFIELD")){
			JTextField tf = new JTextField(""); 
   			tf.addActionListener(this) ; 
			addToGui(tf) ; 
   		}else if(type.equals("DAF")){
   			// poulate a combo box with contents of views from daf
   			String[] contents = (String[])views.toArray(new String[0]); 
   			JComboBox menu = new JComboBox(contents){
   				@Override
   				protected void fireActionEvent() {
   					if (currentSource!= null)
   						super.fireActionEvent();
   				}
   			};
   			addToGui(menu) ; 
   		}else if(type.equals("FILE")){
   			JTextField fileField = new JTextField(file.getPath()); 
   			fileField.setHorizontalAlignment(JTextField.RIGHT) ;
   			fileField.setScrollOffset(50) ; 
   			addToGui(fileField) ; 
   		}else{
   			// this is where we will populate the dropdown menus 
   			String[] contents = (String[])hm.get(type); 
   			JComboBox menu = new JComboBox(contents) ;
   			menu.addActionListener(this) ; 
   			addToGui(menu) ; 
   		}
	}
  
	// the user should select a ddf or daf file. This will work out what one we have, and grab the other
	private void processDafDdf(File file){
		String filename = file.getAbsolutePath() ; 
		
		Pattern daf = Pattern.compile("daf$");
		Pattern ddf = Pattern.compile("ddf$") ; 
		
		ddfFn = ""; 
		dafFn = ""; 
		
		Matcher daf_m = daf.matcher(filename) ; 
		Matcher ddf_m = ddf.matcher(filename) ; 
		if (daf_m.find()){
			System.out.println("We have the daf file") ;
			dafFn = filename ; 
			ddfFn = daf_m.replaceFirst("ddf");
		}
		else if (ddf_m.find()){
			System.out.println("We have the ddf file") ;
			ddfFn = filename ; 
			dafFn = ddf_m.replaceFirst("daf");	
		}
		
		System.out.println("DAf: " + dafFn); 
		System.out.println("DDF: " + ddfFn); 
		
		FileInputStream fin;		
		views = new ArrayList()  ;

		try
		{
		    // Open and read the daf
		    fin = new FileInputStream (dafFn); 
		    BufferedReader br = new BufferedReader(new InputStreamReader(fin));
		    String line ; 
			line =  br.readLine() ;  

			while (line != null)  {
		    	System.out.println(line);	
		    	if (line.startsWith("view") ){
		    		String[] parts = line.split("\\s+"); 
		    		views.add(parts[1]);
		    	} 
		    	line =  br.readLine() ;
		    } 
		    fin.close();	

		    // get the daf fields
		    fin = new FileInputStream (ddfFn); 
		    br = new BufferedReader(new InputStreamReader(fin));
		    line =  br.readLine() ;
    		ddfColumns = line.split("\\s+"); 
    		
    		
    		gridLayout.setColumns(ddfColumns.length + 1); // add 1 for the delete button 
    		System.out.println("We have " + ddfColumns.length + " columns"); 
		    
    		
            // pad out the layout
            for (int i = 0 ; i < (ddfColumns.length + 1 - frame.getContentPane().getComponentCount() ) ; i ++  ){
            	addToGui(new JLabel("")) ; 
            }
            // do a row of column headings
            addToGui(new JLabel("remove")) ; 
            for (String col  : ddfColumns ){
            	addToGui(new JLabel(col));
            }

            //rowWidgetsStart = frame.getContentPane().getComponentCount() ;
            rowWidgetsStart = 2 ; 
            statusLabel.setText("add files") ; 
    		frame.pack();
		}
		// Catches any error conditions
		catch (IOException e)
		
		{
			System.err.println ("Unable to read from file");
			e.printStackTrace();
			System.exit(-1);
		}
		
		chooseFilesButton.setEnabled(true) ; 
	}
	
	// so I can add an index to each component, and do stuff with that
	private void addToGui(Component component){
		frame.getContentPane().add(component) ;
		componentList.add(component) ; 
	}
	
	// this will actually generate the dff file
	private void createDdf(){
		
		// Stream to write file
		FileOutputStream fout;		
		
		int columns = gridLayout.getColumns(); 
		JTextField  label = (JTextField)componentList.get((columns * 3) + 1); // this should take us to the first filename (3rd row, 1st column, vector starts at index 0)
		System.out.println( "label : " + label.getText()) ; 
		
		String  exptName= label.getText().replaceAll( "^.*/", ""); //should remove all leading directories 
		exptName = exptName.replaceAll("\\..*$", "") ; 				// should remove the file extension
		System.out.println( "Expt_name : " + exptName) ; 
		if (exptName.length() == 0) exptName = "ucsc_submission_" + System.nanoTime() ;  
		
		System.out.println( "Expt_name : " + exptName) ; 
		
		try
		{
		    // Open an output stream
		    fout = new FileOutputStream (outputDir + "/" + exptName + ".ddf");

			for(Object obj :componentList){
				int prevRow = -1 ; 
				int row = getRow(obj) ; 
				int col = getColumn(obj) ; 
				// first row is only buttons
				PrintStream ps = new PrintStream(fout);  
				if (row > 1 && col > 1){
				    // Print a line of text
					String value = "" ; 
					if (obj.getClass().equals( new JTextField().getClass()) ){
						value = ((JTextComponent)obj).getText() ;
					}else if(obj.getClass().equals(new JComboBox().getClass())){
						value = (String)((JComboBox)obj).getSelectedItem() ; 
					}else if(obj.getClass().equals(new JLabel().getClass())){
						value = (String)((JLabel)obj).getText() ; 
					}
					String seperator = "\t"; 
					if ( col %  gridLayout.getColumns() == 0) seperator = "\n" ; 
					ps.print(value + seperator) ; 
				}
			}
		    // Close our output stream
		    fout.close();	
		    System.out.println("ddf created in " + outputDir + "/" + exptName + ".ddf") ; 
		    ddfLocation = outputDir + "/" + exptName + ".ddf" ; 
		
		    createTar(outputDir + "/" + exptName) ; 
		}
		// Catches any error conditions
		catch (IOException e)
		{
			System.err.println ("Unable to write to file\n" + e.getMessage());
			System.exit(-1);
		}		
	}
	
	private void createTar(String tarName){
		
		statusLabel.setText("Creating tar") ; 
		
		try {
			// copy the daf file over
			String dafOnlyFile = dafFn.replaceAll("^.*/" , "" ) ; 
			fileCopy(dafFn , outputDir + dafOnlyFile ) ; 
			
			String tarCmd = "tar -czvf " + tarName + ".tar.gz " + dafFn + " " + ddfFn + "  "  ;
			// add each of the files 
			for(Object o : componentList){
				if (getColumn(o) ==  2  && getRow(o) > 2 ){
					tarCmd = tarCmd   + ((JTextField)o).getText() + " " ;
				}
			}
			
			System.out.println("TAR:" + tarCmd) ;
			startDirButton.setEnabled(false) ; 
			chooseFilesButton.setEnabled(false) ; 
			getDafButton.setEnabled(false) ; 
			createDdfButton.setEnabled(false) ; 
			Process createTar = Runtime.getRuntime().exec(tarCmd) ; 
			JOptionPane.showMessageDialog(frame , "Doing daf, please wait") ; 
			createTar.waitFor() ; 
			JOptionPane.showMessageDialog(frame , "Done daf:\n" + tarName) ; 
			startDirButton.setEnabled(true) ; 
			chooseFilesButton.setEnabled(true) ; 
			getDafButton.setEnabled(true) ; 
			createDdfButton.setEnabled(true) ; 


			statusLabel.setText("done") ;
		}catch(Exception e){
			statusLabel.setText("done") ; 
			e.printStackTrace() ; 
			System.exit(1) ; 
		}
	} 
	
	/**
	 * @param callingComponent
	 * This will take a component from the GUI (textbox / dropdown) and update all relevant components (ones in the same row) with the same value
	 */
	private void updateOtherComponents(Object callingComponent){
		// convert to original object
		int col = getColumn(callingComponent) ; // getColumnRow(callingComponent)[0] ;
		int row = getRow(callingComponent) ;//getColumnRow(callingComponent)[1] ; 
		
		System.out.println("Column" + col + "  row " + row ) ; 
		if (callingComponent.getClass().equals(new JComboBox().getClass())){
			JComboBox currentMenu = (JComboBox)callingComponent ;  
			for (Component c : frame.getContentPane().getComponents()){
				if (getColumn(c) == col && getRow(c) > 2  ){
					if (c!= callingComponent)
						try{
							//c.setEnabled(false);
							JComboBox jc= (JComboBox) c;
							((ComboBoxModel) jc.getModel()).setSelectedItem(currentMenu.getSelectedItem());
//							if (jc.getSelectedItem()!= currentMenu.getSelectedItem())
//								((JComboBox)c).setSelectedItem(currentMenu.getSelectedItem()) ; // update any component in the same row with the callers selection 
							//c.setEnabled(true);
						}
						catch (ClassCastException e){
							System.out.println("Problem with " + componentList.indexOf(c)  + "\n" + e.getMessage()); 
						}
				}
			}
		}
		else if (callingComponent.getClass().equals(new JTextField().getClass())){
			JTextField currentText = (JTextField)callingComponent ;  
			for (Component c : frame.getContentPane().getComponents()){
				if (getColumn(c) == col && getRow(c) > 2  ){
					try{
						((JTextField)c).setText(  currentText.getText()) ; // update any component in the same row with the callers selection 				
					}
					catch(ClassCastException e){
						System.out.println(c.getClass()); 
						//System.out.println("error on: col" + getColumn(c) + " row " + getRow(c) + " index:" + componentList.indexOf(c) ) ; 
						System.out.println(e.getMessage()); 
					}
				}
			}
		}

		System.out.println("----") ; 
	}
	
	
	
	// takes a component from the gui, and returns the column and row it is from (as an array of 2 integers)
	private int[] getColumnRow(Object component){
		int index = componentList.indexOf(component) ; 
		int guiColumns = gridLayout.getColumns() ;  
	
//		System.out.print("index: " + index + " cols: " + guiColumns + "\t\t") ; 
		
		int[] returnInts = {(index % guiColumns) + 1  ,  ((int)(index / guiColumns)+ 1) } ; 
		System.out.println(returnInts[0] + ":" + returnInts[1]); 
		return returnInts ; 
	}
	
	private int getColumn(Object component){
		return getColumnRow(component)[0] ; 
	}
	private int getRow(Object component){
			return getColumnRow(component)[1]  ; 
	}
	
	
	
	/**
	 *  connects to the database and populates the hm hashmap - of tablenames to arrays of contents
	 */
	private void connectToDatabase(){
		hm = new HashMap<String, String[]>() ; 

		// first connect to the database, and load the hashmap with the relevant tables
		try{
			Class.forName("com.mysql.jdbc.Driver");
		
			String url = "jdbc:mysql://monstre1.crg.es:3306/encode_transcriptome_v1";
			Connection con = DriverManager.getConnection( url, "encode", "3ncod3"); //"colin", "encode");

			System.out.println("URL: " + url);
			System.out.println("Connection: " + con);

		     //Get a Statement object
			Statement stmt = con.createStatement();
			 
			 
			for (String[] field_n_type: ddf_fields){
				// get the ucsc names from database (unless the field is the wrong type)
				if ( ! (field_n_type[1].equals("DAF")  || field_n_type[1].equals("TEXTFIELD") || field_n_type[1].equals("FILE")) ){
					String tablename = field_n_type[1] ;
					ResultSet rs = stmt.executeQuery("SELECT ucsc_name FROM " + tablename) ; 
			 		Vector ucscnames = new Vector() ; 

					while (rs.next()){
			 			String str = rs.getString("ucsc_name");
			 			System.out.println(str);
			 			
			 			ucscnames.add(rs.getString("ucsc_name")) ; 
				 	} 
					String[] strings = (String[]) ucscnames.toArray(new String[0]);  // 'cos java is crap	
					hm.put( tablename ,  strings  ); 
				}
			}
		
		}catch(Exception e){
			System.out.println(e.toString()); 
		}
		
	}
	
	
	
	public void fileCopy(String inputPath , String outputPath) throws IOException {
	    File inputFile = new File(inputPath);
	    File outputFile = new File(outputPath);

	    FileReader in = new FileReader(inputFile);
	    FileWriter out = new FileWriter(outputFile);
	    int c;

	    while ((c = in.read()) != -1)
	      out.write(c);

	    in.close();
	    out.close();
	}

///-------------------------------------- event handlers ----------------------------------------------------------------------------
	
	  public void actionPerformed(ActionEvent e) {
	        System.out.println("ACTION: "+ e.getActionCommand() + " "+ componentList.indexOf(e.getSource())) ; 
	        currentSource= e.getSource();
	        
	        if (e.getSource()  == chooseFilesButton) {
	        	//Create a file chooser
		        final JFileChooser fc = new JFileChooser(new File(fileStartDirectory));
		        fc.setMultiSelectionEnabled(true) ; 
	        	int returnVal = fc.showOpenDialog(frame);
	            if (returnVal == JFileChooser.APPROVE_OPTION) {
	                File[] files = fc.getSelectedFiles();
	                this.addRows(files) ; 
	            }
	        }
	        else if(e.getSource() == startDirButton){
	        	//Create a file chooser
		        final JFileChooser fc = new JFileChooser();
	        	fc.setDialogTitle("Start Directory");
	            fc.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);

		        int returnVal = fc.showOpenDialog(frame); 

	        	if (returnVal == JFileChooser.APPROVE_OPTION) {
	                File file = fc.getSelectedFile();
	                //This is where a real application would open the file
	                System.out.println(file.getPath()); 
	                outputDir = file.getPath() ; 
	        	}	   
	       }
	       else if(e.getSource() == getDafButton){
	        	//Create a file chooser
		        final JFileChooser fc = new JFileChooser(daf_dir);
	        	fc.setDialogTitle("Select a daf / ddf");

		        int returnVal = fc.showOpenDialog(frame); 

	        	if (returnVal == JFileChooser.APPROVE_OPTION) {
	                File file = fc.getSelectedFile();
	                processDafDdf(file) ;  
	            }
	       }
	       else if (e.getSource() == createDdfButton){
	    	    createDdf() ; 
	       }
	       else if(componentList.indexOf(e.getSource()) > rowWidgetsStart && getColumn(e.getSource()) > 2 ){
	    	   // an event from one of the data rows. 
	    	   if (updateAllBox.isSelected()){
	    		   updateOtherComponents(e.getSource()) ; 
	    	   }
	       }else if (getColumn(e.getSource()) == 1 && getRow(e.getSource()) > 1){
	    	   System.out.println("Removing line: " + getRow(e.getSource()) ) ; 
	    	   Vector indexs = new Vector() ;  
	    	   for(Object o: componentList){   
	    		   if (getRow(o) == getRow(e.getSource()) && o != e.getSource()){  
	    			   frame.getContentPane().remove((JComponent)o) ;
	    			   indexs.add(o); 
	    		   }
	    	   }
	    	   frame.remove((JComponent)e.getSource()) ; 
	    	   for(Object o: indexs){
	    		   componentList.remove(o) ; 
 	    	   }
	    	   componentList.remove(e.getSource()) ; 
	    	   frame.pack() ;  
	       }
	       currentSource= null;
	    }

}
