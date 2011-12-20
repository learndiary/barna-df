/*
 * This file is part of the Flux Library.
 *
 * The code of the Flux Library may be freely distributed and modified under the terms of the
 * European Union Public Licence (EUPL) published on the web site <http://www.osor.eu/eupl/european-union-public-licence-eupl-v.1.1>.
 * Copyright for the code is held jointly by the individual authors, who should be listed
 * in @author doc comments. According to Article 5 and Article 11 of the EUPL, publications that
 * include results produced by the Flux Library are liable to reference the Work,
 * see the Flux Library homepage <http://flux.sammeth.net> for more information.
 */

package barna.commons.io;

import barna.commons.ByteArrayCharSequence;

import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.HashSet;
import java.util.Hashtable;

class ThreadedIOHandler extends Thread implements IOHandler {
    Hashtable<Object, byte[]> bufHash;
    Hashtable<Object, Integer> posHash;
    HashSet<Object> closedHash;
    byte[][] reuseVec;
    int reusePos = 0;
    Object[] a;

    int maxStreamCap = -1;
    long t0, totIn, totOut;
    boolean stop = false;
    int minVol = 1024;
    Monitor mon;
    boolean monitor = false;
    private ByteArrayCharSequence bufferSequence;

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


    public ThreadedIOHandler(int maxStreamCapacity) {
        super("SyncDiskThread");
        this.maxStreamCap = maxStreamCapacity;
        bufHash = new Hashtable<Object, byte[]>(maxStreamCapacity * 2);
        posHash = new Hashtable<Object, Integer>(bufHash.size());
        closedHash = new HashSet<Object>(bufHash.size());
        reuseVec = new byte[maxStreamCapacity][];
        a = new Object[maxStreamCapacity];
    }

    public void addStream(Object stream) {
        addStream(stream, IOHandler.DEFAULT_BUFFER_SIZE);
    }

    public void addStream(Object stream, int bufSize) {
        byte[] buf = null;
//		System.err.print("add() locking rv..");
//		System.err.flush();
        synchronized (reuseVec) {
//			System.err.print("locked..");
//			System.err.flush();
            if (reusePos > 0 && reuseVec[reusePos - 1].length >= bufSize) {
                buf = reuseVec[--reusePos];
            }
        }
//		System.err.println("free.");

        if (buf == null) {
            buf = new byte[bufSize];    // no init
        }
//		System.err.print("add() locking buf#..");
//		System.err.flush();
        synchronized (bufHash) {
//			System.err.print("locked..");
//			System.err.flush();
            while (bufHash.size() >= maxStreamCap) {
                try {
                    bufHash.wait();
                } catch (InterruptedException e) {
                    ; // :)
                }
            }
            bufHash.put(stream, buf);
            posHash.put(stream, 0);
        }
//		System.err.println("free");

    }

    /**
     * for InputStreams, blocks until finished reading
     * for OutputStreams, blocks untif finished writing
     *
     * @param stream
     * @return
     */
    public void removeStream(Object stream) {
        byte[] b = null;
        if (stream instanceof InputStream) {
            synchronized (closedHash) {
                while (!closedHash.contains(stream)) {
                    try {
                        closedHash.wait();
                    } catch (InterruptedException e) {
                        ; // :)
                    }
                }
            }
        }

        synchronized (bufHash) {
            b = bufHash.get(stream);
            int p = posHash.get(stream);
            while (p > 0) {
                synchronized (b) {
                    try {
                        b.wait();
                    } catch (InterruptedException e) {
                        ; // :)
                    }
                }
            }
            bufHash.remove(stream);
            posHash.remove(stream);
            bufHash.notifyAll();
        }
        synchronized (reuseVec) {
            if (reusePos < reuseVec.length) {
                reuseVec[reusePos++] = b;
            }
            reuseVec.notify();
        }
        synchronized (closedHash) {
            closedHash.remove(stream);
        }

        //return posHash.get(stream);
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
            int max = -1;
            int n = -1;
//			System.err.print("run() locking bHash..");
//			System.err.flush();
            synchronized (bufHash) {
//				System.err.print("locked..");
//				System.err.flush();
                a = bufHash.keySet().toArray(a);
                n = bufHash.size();
                for (int i = 0; i < n; i++) {
                    int curr = -1;
                    if (a[i] instanceof InputStream) {
                        if (!closedHash.contains(a[i]))    // not closed
                        {
                            curr = bufHash.get(a[i]).length
                                    - posHash.get(a[i]);
                        }
                    } else //if (a[i] instanceof OutputStream)
                    {
                        curr = posHash.get(a[i]);
                    }

                    if (curr > minVol && curr > max) {
                        max = curr;
                        x = i;
                    }
                }
            }
//			System.err.println("free");

            if (x < 0) {
                try {
                    sleep(10);
                } catch (InterruptedException e) {
                    ; // :)
                }
                continue;
            }
            b = bufHash.get(a[x]);
//			System.err.print("run() locking b..");
//			System.err.flush();
            synchronized (b) {
//				System.err.print("locked..");
//				System.err.flush();
                if (a[x] instanceof InputStream) {
                    //System.err.println("filling "+ bufHash.get(a[x]).length+ " - "+ posHash.get(a[x]));
                    fill((InputStream) a[x]);
                } else if (a[x] instanceof OutputStream) {
                    try {
                        flush((OutputStream) a[x]);
                    } catch (IOException e) {
                        throw new RuntimeException(e);
                    }
                }

            }
//			System.err.println("free.");

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

    public boolean close() {
        stop = true;
        while (isAlive()) {
            try {
                this.join();
            } catch (InterruptedException e) {
                ; // :)
            }
        }

        int n = -1;
        synchronized (posHash) {
            a = posHash.keySet().toArray(a);
            n = posHash.size();
        }
        for (int i = 0; i < n; ++i) {
            if (a[i] instanceof OutputStream) {
                OutputStream out = (OutputStream) a[i];
                try {
                    flush(out);
                    out.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            } else {
                InputStream in = (InputStream) a[i];
                try {
                    in.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        }
        return true;
    }

    public int readLine(InputStream in, ByteArrayCharSequence cs) {
        byte n = IOHandler.BYTE_NL, r = IOHandler.BYTE_CR;
        byte[] b = bufHash.get(in); // , bb= null;
//		System.err.print("rLine() locking b..");
//		System.err.flush();
        int i;
        synchronized (b) {
//			System.err.print("locked..");
//			System.err.flush();
            int pos = posHash.get(in);
            while (pos < minVol && !closedHash.contains(in)) {
                if (isAlive()) {
                    try {
                        interrupt();
                        b.wait(10);
                    } catch (InterruptedException e) {
                        ; // :)
                    }
                } else {
                    fill(in);
                }
                pos = posHash.get(in);
            }
            int p = pos;
            if (p == 0) {
                assert (closedHash.contains(in));
                return -1;
            }
            i = 0;
            for (; i < p && b[i] != n; ++i) {
                ;    // find lsep
            }

            if (i > 0 && b[i - 1] == r) {
                --i;
            }
            //bb= new byte[i];
            System.arraycopy(b, 0, cs.chars, 0, i);
            cs.start = 0;
            cs.end = i;

            while (++i < p && (b[i] == r || b[i] == n)) {
                ;
            }
//			if (i< 0|| p> b.length)
//				System.currentTimeMillis();
            if (p - i > 0) {
                System.arraycopy(b, i, b, 0, p - i);
                posHash.put(in, pos - i);
            } else {
                posHash.put(in, 0); // last line
            }
            //System.err.println("read "+ i+ " from "+ idx);
        }
//		System.err.println("free.");
        return i;
    }


    private int fill(InputStream x) {

        byte[] b = bufHash.get(x);
        int len = -1;
//		System.err.println("fill() locking b..");
//		System.err.flush();
        synchronized (b) {
//			System.err.print("locked..");
//			System.err.flush();
            int p = posHash.get(x);
            len = b.length - p;
            //if (c< x)
            //	x= c;
            try {
                len = x.read(b, p, len);
                if (len < 0) {
                    x.close();
                    closedHash.add(x);
                    b.notifyAll();
                    //System.err.println("closed "+ x);
                } else {
                    posHash.put(x, posHash.get(x) + len);
                    totIn += len;
                    b.notifyAll();
                    //System.err.println("filled "+ len+ " from "+ x);
                }
                //System.err.println("->notify");
            } catch (IOException e) {
                e.printStackTrace();
                return -1;
            }
        }
        //System.err.println("read "+ len);
//		System.err.println("free.");
        return len;
    }

    public void writeLine(Object object, OutputStream out) throws IOException {
        if (object == null) {
            throw new NullPointerException();
        }
        if (bufferSequence == null) {
            bufferSequence = new ByteArrayCharSequence(object.toString());
        }
        bufferSequence.clear();
        bufferSequence.append(object.toString());
        writeLine(bufferSequence, out);
    }


    public void writeLine(ByteArrayCharSequence cs, OutputStream out) throws IOException {
        int len = cs.length();
        byte[] b = bufHash.get(out);
//		System.err.print("wLine() locking b..");
//		System.err.flush();
        synchronized (b) {
//			System.err.print("locked..");
//			System.err.flush();
            int pos = posHash.get(out);
            while (pos + len + 1 > b.length) {
                //flush(idx);
                if (isAlive()) {
                    try {
                        interrupt();
                        b.wait(10);
                    } catch (InterruptedException e) {
                        ; // :)
                    }
                } else {
                    flush(out);
                }
                pos = posHash.get(out);
            }
            System.arraycopy(cs.chars, cs.start, b, pos, len);
            pos += len;
            b[pos++] = IOHandler.BYTE_NL;
            posHash.put(out, pos);
        }
        //System.out.println("put "+len);
//		System.err.println("free.");


    }

    public void write(byte[] bb, int x, int len, OutputStream idx) throws IOException {
        byte[] b = bufHash.get(idx);
//		System.err.print("write() locking b..");
//		System.err.flush();
        synchronized (b) {
//			System.err.print("locked..");
            System.err.flush();
            int pos = posHash.get(idx);
            while (pos + len + 1 > b.length) {
                //flush(idx);
                try {
                    interrupt();
                    b.wait();
                } catch (InterruptedException e) {
                    ; // :)
                }
                pos = posHash.get(idx);
            }
            System.arraycopy(bb, x, b, pos, len);
            pos += len;
            posHash.put(idx, pos);
        }
        //System.out.println("put "+len);
//		System.err.println("free.");

    }

    public static final Integer INTEGER_0 = new Integer(0);

    private void flush(OutputStream out) throws IOException {
        byte[] b = bufHash.get(out);
        int p = -1;
        synchronized (b) {
            p = posHash.get(out);
            out.write(b, 0, p);
            out.flush();
            posHash.put(out, INTEGER_0);
            totOut += p;
            b.notify();


        }
    }


    public Monitor getMonitor() {
        if (mon == null) {
            mon = new Monitor();
        }

        return mon;
    }

    /**
     * @param idx
     * @return
     * @deprecated
     */
    public ByteArrayCharSequence readLine(InputStream idx) throws IOException {
        if (1 == 1) {
            ByteArrayCharSequence cs = new ByteArrayCharSequence(IOHandler.DEFAULT_BUFFER_SIZE);
            int i = readLine(idx, cs);
            if (i < 0) {
                return null;
            }
            return cs;
        }


        byte n = IOHandler.BYTE_NL, r = IOHandler.BYTE_CR;
        byte[] b = bufHash.get(idx), bb = null;
        synchronized (b) {
            //			while (pos[idx]== 0)
            //				fill(idx);
            int pos = posHash.get(idx);
            while (pos < minVol && idx != null) {
                if (isAlive()) {
                    try {
                        interrupt();
                        b.wait();
                    } catch (InterruptedException e) {
                        ; // :)
                    }
                } else {
                    fill(idx);
                }
                pos = posHash.get(idx);
            }
            int p = posHash.get(idx);
            if (p == 0) {
                assert (idx == null);
                return null;
            }
            int i = 0;
            for (; i < p && b[i] != n; ++i) {
                ;    // find lsep
            }

            assert (i != p || idx == null);
            if (i > 0 && b[i - 1] == r) {
                --i;
            }
            bb = new byte[i];
            System.arraycopy(b, 0, bb, 0, i);

            while (++i < p && (b[i] == r || b[i] == n)) {
                ;
            }
            System.arraycopy(b, i, b, 0, p - i);
            posHash.put(idx, pos - i);
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
