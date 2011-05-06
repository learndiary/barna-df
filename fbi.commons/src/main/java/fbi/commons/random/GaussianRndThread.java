package fbi.commons.random;

import java.util.Random;

public class GaussianRndThread extends Thread {
    double[] rnds;
    int p1 = -1, p2 = -1, cap = -1;
    ;
    Random rnd;
    boolean stop = false;
    Thread thisThread;

    public GaussianRndThread(int capacity) {
        rnd = new Random();
        this.cap = capacity;
        rnds = new double[capacity * 10];
        p1 = 0;
        p2 = 0;
    }

    @Override
    public void run() {
        thisThread = Thread.currentThread();
        while (!stop) {
            if (p1 >= rnds.length - cap) {
                synchronized (this) {
                    System.arraycopy(rnds, p1, rnds, 0, p2 - p1);
                    p1 = 0;
                    p2 -= p1;
                }
            }
            while (p2 - p1 < cap && p2 < rnds.length)
                rnds[p2++] = rnd.nextGaussian();
            if (p1 == 0)
                try {
                    sleep(10);
                } catch (InterruptedException e) {
                    ; //;
                }
        }
    }

    public void setStop() {
        stop = true;
    }

    public double nextGaussian() {
        while (p1 >= p2) {
            thisThread.interrupt();
        }
        return rnds[p1++];
    }

}
