package barna.astalavista;

import barna.model.ASEvent;
import barna.model.commons.MyFile;

import java.io.BufferedOutputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.Queue;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.zip.GZIPOutputStream;

/**
 * Created with IntelliJ IDEA.
 * User: micha
 * Date: 9/16/12
 * Time: 5:14 PM
 * To change this template use File | Settings | File Templates.
 */
public class WriterThread extends Thread {
    Queue<ASEvent> queue = new ConcurrentLinkedQueue<ASEvent>();
    Thread writingThread;
    boolean kill = false, outputASTA = false, outputGTF = true;
    int maxQevents = 10000, minQevents = 1000;
    static boolean writeHeader = true;
    public String outputFname;

    public static boolean outputSeq = false;



    public WriterThread() {
        super();
        setName("event_writing_thread");
    }

    public void addEvent(ASEvent ev) {
        queue.add(ev);
        writingThread.interrupt();
    }

    public void setKill(boolean newKill) {
        kill = newKill;
    }


    public void run() {
        writingThread = Thread.currentThread();

        //BufferedWriter buffy= new BufferedWriter(new OutputStreamWriter(
        OutputStreamWriter globalWriter = null;

        try {
            if (outputFname != null) {
                MyFile f = new MyFile(outputFname);
                if (!f.exists())
                    writeHeader = true;
                GZIPOutputStream zipperStream = new GZIPOutputStream(new BufferedOutputStream(new FileOutputStream(f, false)));    // WAS true, for appending
                globalWriter = new OutputStreamWriter(zipperStream);
            } else
                globalWriter = new OutputStreamWriter(System.out);

            //				if (writeHeader) {
            //					//Writer writer= new OutputStreamWriter(zipperStream);
            //					//printSettings(writer);
            //					printSettings(globalWriter);
            //					writeHeader= false;
            //					//writer.close();
            //				}
        } catch (IOException e) {
            e.printStackTrace();
        }

        while (true) {
            try {
                while (!queue.isEmpty()) {
                    ASEvent ev = queue.poll();
                    ev.init();
                    String s = null;
                    if (outputGTF)
                        s = ev.toStringGTF(outputSeq);    // true
                    if (outputASTA)
                        s = ev.toStringASTA();
                    //for (int i = 0; i < s.length(); i++)
                    //zipperStream.write((byte) s.charAt(i));
                    globalWriter.write(s.toCharArray(), 0, s.length());
                    //zipperStream.write('\n');
                    globalWriter.write('\n');
                }
                //zipperStream.flush(); zipperStream.close();

                // better compression if not flushed?
//					globalWriter.flush();

                // close only at the end
//					if (outputFname!= null)
//						globalWriter.close();

                //buffy.flush(); buffy.close();

            } catch (IOException e) {
                e.printStackTrace();
            }

            if (kill) {
                try {
                    globalWriter.flush();
                    if (outputFname != null)
                        globalWriter.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
                return;
            }

            try {
                writingThread.sleep(0);
            } catch (InterruptedException e) {
                ; // :)
            }
        }
    }

    public int getMaxQevents() {
        return maxQevents;
    }

    public void setMaxQevents(int maxQevents) {
        this.maxQevents = maxQevents;
    }

    public int getMinQevents() {
        return minQevents;
    }

    public void setMinQevents(int minQevents) {
        this.minQevents = minQevents;
    }

}
