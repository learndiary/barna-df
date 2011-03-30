package commons;


public interface ReadyOrNot {

	public boolean isReady();

	public void set(Object settings);
	
	public boolean loadStats();
	
	public boolean isFinished();
	
	public void setLoadStats(boolean val);
	
	public boolean getLoadStats();
	
	public void killResult();
}
