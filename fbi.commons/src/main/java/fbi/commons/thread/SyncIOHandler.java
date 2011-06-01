package fbi.commons.thread;

import fbi.commons.ByteArrayCharSequence;

import java.io.*;

public class SyncIOHandler extends Thread {

    public static void main(String[] args) {
        try {
            // H_sapiens-HepG2_10879_42KJPAAXX_75_sorted.bed
            // H_sapiens-HepG2-r2.nextgeneid.AS_sorted.gtf
            // big.gtf
            InputStream fis = new FileInputStream("N:\\tmp\\H_sapiens-HepG2-r2.nextgeneid.AS_sorted.gtf");
            OutputStream fos = new FileOutputStream("N:\\tmp\\H_sapiens-HepG2-r2.nextgeneid.AS_sorted.gtf_3");
            SyncIOHandler tst = new SyncIOHandler(new InputStream[]{fis}, new OutputStream[]{fos}, 10 * 1024);
            tst.start();

            ByteArrayCharSequence cs = new ByteArrayCharSequence(128);
            while ((cs = tst.readLine(0)) != null) {
                tst.writeLine(cs, 1);
            }
            tst.close();

        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public static byte BYTE_NL = '\n', BYTE_CR = '\r';
    InputStream[] readers;
    OutputStream[] writers;
    byte[][] buf;
    int[] pos;
    long t0, totIn, totOut;
    boolean stop = false;
    int minVol = 1024;
    Monitor mon;
    boolean monitor = false;

    class Monitor extends Thread {

        long delta = 5000;

        public Monitor() {
            super("RequestMonitor");
        }

        @Override
        public void run() {
            while (!stop) {
                try {
                    sleep(delta);
                } catch (InterruptedException e) {
                    ; // :)
                }
                float dt = ((System.currentTimeMillis() - t0) / 1000f);
                System.err.println(
                        "sec: " + dt + "\ttot: " + ((totIn + totOut) / dt) + "\t in:" + (totIn / dt) + "\tout: " + (totOut / dt));
            }
        }

    }


    public SyncIOHandler(InputStream[] readers, OutputStream[] writers, int bufSize) {
        super("SyncDiskThread");
        this.readers = readers;
        this.writers = writers;
        buf = new byte[readers.length + writers.length][];
        for (int i = 0; i < buf.length; i++) {
            buf[i] = new byte[bufSize];
        }
        pos = new int[buf.length];
        for (int i = 0; i < pos.length; i++) {
            pos[i] = 0;
        }
    }

    @Override
    public void run() {

        t0 = System.currentTimeMillis();
        if (monitor) {
            getMonitor().start();
        }
        int x;
        byte[] b;
        float thr;
        while (!stop) {
            x = -1;
/*			while (x< 0)
				synchronized (prio) {
					thr= 0.0f;
					for (int i = 0; i < prio.length; i++) {
						if (prio[i]> thr) {
							thr= prio[i];
							x= i;
						}
					}
					if (x< 0)
						try {
							prio.wait();
						} catch (InterruptedException e) {							
							;	// :)
						}
				}
			System.err.println("prio "+x+" "+prio[x]);
*/
/*			for (int i = 0; i < readers.length; i++) {
                pos[i]= fill(buf[i], pos[i], readers[i]);
            }
            for (int i = 0; i < writers.length; i++) {
                x= i+ readers.length;
                pos[x]= flush(buf[x], pos[x], writers[i]);
            }
*/
            //x= mIdx< readers.length? mIdx: mIdx- readers.length;

            int max = -1;
            for (int i = 0; i < pos.length; i++) {
                int curr;
                if (i < readers.length) {
                    if (readers[i] == null) {
                        continue;
                    }
                    curr = buf[i].length - pos[i];
                } else {
                    curr = pos[i];
                }

                if (curr > minVol && curr > max) {
                    max = curr;
                    x = i;
                }
            }
            if (x < 0) {
                continue;
            }
            b = buf[x];
            synchronized (b) {
                if (x < readers.length) {
                    fill(x);
                } else {
                    flush(x);
                }

            }

        }
        if (monitor) {
            getMonitor().interrupt();
            try {
                getMonitor().join();
            } catch (InterruptedException e) {
                ; // :)
            }
        }
    }

    public void close() {
        stop = true;
        while (isAlive()) {
            try {
                this.join();
            } catch (InterruptedException e) {
                ; // :)
            }
        }
        for (int i = 0; i < readers.length; i++) {
            if (readers[i] == null) {
                continue;
            }
            try {
                readers[i].close();
            } catch (IOException e) {
                e.printStackTrace();
            }
        }

        for (int i = 0; i < writers.length; i++) {
            flush(readers.length + i);
            try {
                writers[i].close();
            } catch (IOException e) {
                e.printStackTrace();
            }
        }


    }

    public ByteArrayCharSequence readLine(int idx, ByteArrayCharSequence cs) {
        byte n = BYTE_NL, r = BYTE_CR;
        byte[] b = buf[idx]; // , bb= null;
        synchronized (b) {
//			while (pos[idx]== 0) 
//				fill(idx);
            while (pos[idx] < minVol && readers[idx] != null) {
                try {
                    b.wait();
                } catch (InterruptedException e) {
                    ; // :)
                }
            }
            int p = pos[idx];
            if (p == 0) {
                assert (readers[idx] == null);
                return null;
            }
            int i = 0;
            for (; i < p && b[i] != n; ++i) {
                ;    // find lsep
            }

            assert (i != p || readers[idx] == null);
            if (i > 0 && b[i - 1] == r) {
                --i;
            }
            //bb= new byte[i];
            System.arraycopy(b, 0, cs.chars, 0, i);
            cs.start = 0;
            cs.end = i;

            while (++i < p && b[i] == r || b[i] == n) {
                ;
            }
            System.arraycopy(b, i, b, 0, p - i);
            pos[idx] -= i;
        }
//		if (bb== null) 
//			return null;
        //System.err.println("got "+bb.length);
        //ByteArrayCharSequence cs= new ByteArrayCharSequence(bb);
        return cs;
    }


    private int fill(int x) {
        byte[] b = buf[x];
        int len = -1;
        synchronized (b) {
            int p = pos[x];
            len = b.length - p;
            //if (c< x)
            //	x= c;
            try {
                len = readers[x].read(b, p, len);
                if (len < 0) {
                    readers[x].close();
                    readers[x] = null;
                    b.notifyAll();
                } else {
                    pos[x] += len;
                    totIn += len;
                    b.notifyAll();
                }
            } catch (IOException e) {
                e.printStackTrace();
                return -1;
            }
        }
        //System.err.println("read "+ len);
        return len;
    }

    public boolean writeLine(ByteArrayCharSequence cs, int idx) {
        int len = cs.length();
        byte[] b = buf[idx];
        synchronized (b) {
            while (pos[idx] + len + 1 > b.length)
            //flush(idx);
            {
                try {
                    b.wait();
                } catch (InterruptedException e) {
                    ; // :)
                }
            }
            System.arraycopy(cs.chars, cs.start, b, pos[idx], len);
            pos[idx] += len;
            b[pos[idx]++] = BYTE_NL;
        }
        //System.out.println("put "+len);
        return true;
    }

    private int flush(int idx) {
        //int y, z;
        OutputStream out = writers[idx - readers.length];
        byte[] b = buf[idx];
        int p = -1;
        synchronized (b) {
            p = pos[idx];
            //y= c> p? p: c;
            try {
                out.write(b, 0, p);
                out.flush();
                pos[idx] = 0;
                totOut += p;
                b.notifyAll();
            } catch (IOException e) {
                e.printStackTrace();
                return -1;
            }
            //z= p- y;
            //System.arraycopy(b, p- 1, b, 0, p);

        }
        //System.out.println("wrote "+ p);
        return 0;
    }

    public Monitor getMonitor() {
        if (mon == null) {
            mon = new Monitor();
        }

        return mon;
    }

    public ByteArrayCharSequence readLine(int idx) {
        byte n = BYTE_NL, r = BYTE_CR;
        byte[] b = buf[idx], bb = null;
        synchronized (b) {
            //			while (pos[idx]== 0)
            //				fill(idx);
            while (pos[idx] < minVol && readers[idx] != null) {
                try {
                    b.wait();
                } catch (InterruptedException e) {
                    ; // :)
                }
            }
            int p = pos[idx];
            if (p == 0) {
                assert (readers[idx] == null);
                return null;
            }
            int i = 0;
            for (; i < p && b[i] != n; ++i) {
                ;    // find lsep
            }

            assert (i != p || readers[idx] == null);
            if (i > 0 && b[i - 1] == r) {
                --i;
            }
            if (i == 0) {
                System.currentTimeMillis();
            }
            bb = new byte[i];
            System.arraycopy(b, 0, bb, 0, i);

            while (++i < p && (b[i] == r || b[i] == n)) {
                ;
            }
            System.arraycopy(b, i, b, 0, p - i);
            pos[idx] -= i;
        }
        if (bb == null) {
            return null;
        }
        //System.err.println("got "+bb.length);
        ByteArrayCharSequence cs = new ByteArrayCharSequence(bb);
        if (cs.length() == 0) {
            System.currentTimeMillis();
        }
        return cs;
    }

    public void setMonitor(boolean monitor) {
        this.monitor = monitor;
    }
}
