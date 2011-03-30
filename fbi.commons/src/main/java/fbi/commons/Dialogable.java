package fbi.commons;

public interface Dialogable {
	public boolean checkOverwrite(String s);
	public void showWarning(String s);
	public void showError(String s);
	public void showInfo(String s);
}
