package fbi.genome.model.constants;

import java.io.PrintStream;

import commons.Progressable;



public class PrintstreamProgressable implements Progressable {

	public static final char DEFAULT_PROG_CHAR= '*', DEFAULT_PROG_MILESTONE= '+';
	public static char progressChar= DEFAULT_PROG_CHAR;
	
	PrintStream p;
	int min= 0, max= 9, val= 0;
	public PrintstreamProgressable(PrintStream stream) {
		p= stream;
	}
	
	public void progress() {
		p.print(progressChar);
		p.flush();
		++val;
	}
	
	public void finish() {
		for (int i = val; i < max; i++) 
			p.print("*");
		p.println();
	}

	public void setMaximum(int newValue) {
		max= newValue;
	}

	public void setMinimum(int newValue) {
		min= newValue;
	}

	public void setString(String value) {
		p.print("\t"+value+" ");
		p.flush();
	}

	public void setValue(int newValue) {
		val= newValue;
	}

	public void message(String value) {
		p.println(value);
	}

	public void finish(String msg, long time) {
		
		for (int i = val; i < max; i++) 
			p.print("*");
		
		if (msg!= null) 
			p.print(Constants.SPACE+ msg);
		
		if (time>= 0) {
			int hh= (int) (time/ 3600000);
			if (hh> 0)
				time%= (hh* 3600000);
			int mm= (int) (time/ 60000);
			if (mm> 0)
				time%= (mm* 60000);
			int ss= (int) (time/ 1000);
			if (ss> 0)
				time%= (ss* 1000);
			String s= hh< 10? Constants.NULL+ Integer.toString(hh): Integer.toString(hh);
			s+= Constants.COLON;
			s+= mm< 10? Constants.NULL+ Integer.toString(mm): Integer.toString(mm);
			s+= Constants.COLON;
			s+= ss< 10? Constants.NULL+ Long.toString(ss): Long.toString(ss);
			
			p.print(Constants.SPACE+ Constants.PAROPEN+ s+ Constants.PARCLOSE);
		}
		
		
		p.println();
	}
}
