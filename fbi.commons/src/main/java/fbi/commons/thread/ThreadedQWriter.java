package fbi.commons.thread;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.concurrent.ConcurrentLinkedQueue;

public class ThreadedQWriter extends Thread {

    private static int id = 0;

    private boolean closed = false;
    protected File file = null;
    protected BufferedWriter writer = null;
    protected String lineSep = "\n";
    protected ConcurrentLinkedQueue q;
    protected int maxSize = 100;    // 1024, 8192
    private boolean silent = true, append = false, stop = false;

    public ThreadedQWriter(File file) {
        super("QWriter-" + (id++));
        this.q = new ConcurrentLinkedQueue();
        this.file = file;
    }

    public ThreadedQWriter(BufferedWriter writer) {
        super("QWriter-" + (id++));
        this.q = new ConcurrentLinkedQueue();
        this.writer = writer;
    }

    public boolean init() {
        try {
            this.writer = new BufferedWriter(new FileWriter(file, append));
        } catch (IOException e) {
            if (!silent) {
                e.printStackTrace();
            }
            return false;
        }

        return true;
    }

    @Override
    public void run() {

        if (writer == null) {
            init();
        }
        while (!closed) {
            if (stop) {
                break;
            } else {
                writeAll();
                try {
                    sleep(100);
                } catch (InterruptedException e) {
                    ; // :)
                }
            }
        }
        if (stop || closed) {
            while (!q.isEmpty()) {
                writeAll();     // q.poll();
            }
        }
        try {
            writer.flush();
            if (file != null) {
                writer.close();
            }
        } catch (IOException e) {
            ; // :)
        }
    }

    public void flush() {
        interrupt();    // let only the writer thread modify the queue !!
        while (!q.isEmpty()) {
            try {
                Thread.sleep(100);
            } catch (Exception e) {
                ; // :)
            }
        }
    }

    public void writeAll() {
        try {
            while (!q.isEmpty()) {
                String s = q.poll().toString();
                if (bytes > 0) {
                    bytes -= s.length();
                    writer.write(s);
                } else {
                    String t = s + lineSep;
                    writer.write(t);
                }
            }
            writer.flush();
        } catch (Exception e) {
            ; // :)
        }
    }

    int added = 0;

    public void add(Object o) {
        if (o == null) {
            return;
        }
        q.add(o);
        ++added;
        if (added > maxSize) {    // q.size()> ... f* slow
            interrupt();
            added = 0;
        }
    }

    protected long bytes = 0, limitBytes = -1;

    public void addString(Object o) {
        if (o == null) {
            return;
        }
        String s = o.toString();
        q.add(s);
        bytes += s.length();
        if (limitBytes > 0) {
            while (bytes >= limitBytes) {
                interrupt();
                try {
                    Thread.currentThread().sleep(100);
                } catch (InterruptedException e) {
                    ; // :)
                }
            }
        } else if (q.size() > maxSize) {
            interrupt();
        }
    }

    public void close() {
        closed = true;
        interrupt();
    }

    public boolean isSilent() {
        return silent;
    }

    public void setSilent(boolean silent) {
        this.silent = silent;
    }

    public boolean isAppend() {
        return append;
    }

    public void setAppend(boolean append) {
        this.append = append;
    }

    public boolean isStop() {
        return stop;
    }

    public boolean setStop() {
        if (isAlive()) {
            stop = true;
            interrupt();
            return true;
        }
        return false;
    }

    public long getLimitBytes() {
        return limitBytes;
    }

    public void setLimitBytes(long limitBytes) {
        this.limitBytes = limitBytes;
    }

    public boolean setStop(boolean stop) {
        if (stop) {
            return setStop();
        } else {
            this.stop = stop;
            return true;
        }
    }
}
