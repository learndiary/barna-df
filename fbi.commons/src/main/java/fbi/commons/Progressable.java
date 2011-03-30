package fbi.commons;

public interface Progressable {

	public void setValue(int newValue);
	public void setMaximum(int newValue);
	public void setMinimum(int newValue);
	public void setString(String value);
	public void progress();
	public void finish();
	public void finish(String msg, long time);
	public void message(String value);
}
