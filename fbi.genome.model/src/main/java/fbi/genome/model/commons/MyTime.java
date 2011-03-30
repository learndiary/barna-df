package fbi.genome.model.commons;

import java.util.Calendar;

public class MyTime {

	long time= -1l;

	public MyTime(long newTime) {
		this.time= newTime;
	}
	
	public int getHours() {
		return (int) (time/ 3600000);
	}
	
	private long getRestHours() {
		return time% 3600000;
	}
	
	public int getMinutes() {
		return (int) (getRestHours()/ 60000);
	}

	private long getRestMinutes() {
		return getRestHours()% 60000;
	}
	
	public int getSeconds() {
		return (int) (getRestMinutes()/ 1000);
	}

	public long getMillis() {
		return getRestMinutes()% 1000;
	}
	
	public String toString() {
		return getHours()+ ":"+ getMinutes()+":"+ getSeconds()+"."+getMillis();
	}
	
	public static String toTime(long temps) {
		MyTime t= new MyTime(temps);
		return t.toString();
	}
	
	public static String getHexDate() {
		Calendar rightNow= Calendar.getInstance();
		StringBuffer sb= new StringBuffer(6);
		int x= rightNow.get(Calendar.YEAR);
		String s= java.lang.Integer.toString(x);
		sb.append(s.substring(2));		
		x= rightNow.get(Calendar.MONTH);
		++x;
		s= java.lang.Integer.toString(x);
		if (s.length()== 1)
			sb.append("0");
		sb.append(s);
		x= rightNow.get(Calendar.DAY_OF_MONTH);
		s= java.lang.Integer.toString(x);
		if (s.length()== 1)
			sb.append("0");
		sb.append(s);
		return sb.toString();
	}
}
