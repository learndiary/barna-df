package fbi.genome.model.constants;

import fbi.commons.Dialogable;

import java.io.*;
import java.util.HashSet;

public class PrintstreamDialogable implements Dialogable {
	PrintStream p;
	InputStream in;
	public static HashSet<String> setPos= new HashSet<String>(), setNeg= new HashSet<String>(); 
	static {
		setPos.add("YES");
		setPos.add("Y");
		setPos.add("YO");
		setPos.add("YEP");
		setNeg.add("NO");
		setNeg.add("N");
		setNeg.add("NOPE");
		setNeg.add("NONE");
	}
	public PrintstreamDialogable(PrintStream p, InputStream in) {
		this.p= p;
	}
	public boolean checkOverwrite(String str) {
		if (str!= null)
			p.println(str);
		BufferedReader r= new BufferedReader(new InputStreamReader(System.in));
		String s= "";
		while ((!setPos.contains(s))&& (!setNeg.contains(s))) 
			try {
				s= r.readLine().trim().toUpperCase();
			} catch (IOException e) {
				;	// :)
			}
		try {
			r.close();
		} catch (IOException e) {
			;	// :)
		}
		if (setPos.contains(s))
			return true;
		return false;
	}
	public void showError(String s) {
		if (Constants.verboseLevel>= Constants.VERBOSE_ERRORS)
			p.println(s);
	}
	public void showInfo(String s) {
		if (Constants.verboseLevel> Constants.VERBOSE_SHUTUP)
			p.println(s);
	}
	public void showWarning(String s) {
		if (Constants.verboseLevel>= Constants.VERBOSE_ERRORS)
			p.println(s);
	}
}
