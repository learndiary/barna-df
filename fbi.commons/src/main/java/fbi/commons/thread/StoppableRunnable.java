package fbi.commons.thread;

public interface StoppableRunnable extends Runnable {
	public boolean setStop();
	public boolean setStop(boolean stop);
	public boolean isStop();
}
