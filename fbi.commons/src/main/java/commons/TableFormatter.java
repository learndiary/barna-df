package commons;

import java.util.Vector;

public class TableFormatter {
	Vector<String[]> v;
	int col= 0;
	boolean header= false;
	int[] max= null;
	char charFill= ' ', charSpace= ' ', charHbar= '=';
	boolean tabRow= false;
	
	public TableFormatter(int col) {
		this.col= col;
		v= new Vector<String[]>();
		max= new int[col];
		for (int i = 0; i < max.length; i++) 
			max[i]= 0;
	}
	public void add(String[] ss) {
		if (ss.length!= col)
			return;
		for (int i = 0; i < ss.length; i++) 
			max[i]= ss[i].length()> max[i]? ss[i].length(): max[i];
		v.add(ss);
	}
	@Override
	public String toString() {
		if (v.size()== 0)
			return "";
		
		StringBuilder s= new StringBuilder();
		for (int i = 0; i < v.size(); i++) {
			s.append("\t");
			for (int j = 0; j < v.elementAt(i).length; j++) {
				s.append(v.elementAt(i)[j]);
				for (int m = v.elementAt(i)[j].length(); m < max[j]; m++) 
					s.append(charFill);
				s.append(charSpace);
			}
			s.append("\n");
			if (i== 0&& header) {
				for (int j = 0; j < max.length; j++) {
					for (int m = 0; m < max[j]+ 1; m++)	// +1 for spacer 
						s.append(charHbar);
				}
				s.append("\n");
			}
		}
		return s.toString();
	}
	public boolean isTabRow() {
		return tabRow;
	}
	public void setTabRow(boolean tabRow) {
		this.tabRow = tabRow;
	}
	
}
