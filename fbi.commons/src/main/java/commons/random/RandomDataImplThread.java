package commons.random;

import org.apache.commons.math.random.RandomDataImpl;

public class RandomDataImplThread extends Thread {

	public static final byte MODE_UNIF= 0, MODE_GAUSS= 1, MODE_POISS= 2;
	double[] data;
	int pos;
	byte mode;
	double p1, p2; 	// parameters
	RandomDataImpl rnd;
	boolean stop= false;
	
	public RandomDataImplThread(int size, byte mode, double p1, double p2) {
		rnd= new RandomDataImpl();
		data= new double[size];
		this.mode= mode;
		this.p1= p1;
		this.p2= p2;
	}
	
	void fill() {
		synchronized (data) {
			while (pos< data.length)
				data[pos++]= nextVal();
		}
	}
	
	double nextVal() {
		if (mode== MODE_UNIF)
			return rnd.nextUniform(p1, p2);
		if (mode== MODE_GAUSS)
			return rnd.nextGaussian(p1, p2);
		if (mode== MODE_POISS)
			return rnd.nextPoisson(p1);
		return 0d;
	}
	
	@Override
	public void run() {
		while (!stop) {
			if (pos< data.length)
				fill();
			try {
				sleep(10);
			} catch (InterruptedException e) {
				; //:)
			}
		}
	}
	
	public void setStop() {
		this.stop= true;
	}
}
