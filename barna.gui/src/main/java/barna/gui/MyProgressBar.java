/*
 * Copyright (c) 2010, Micha Sammeth
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * The names of its contributors may be not used to endorse or promote
 *       products derived from this software without specific prior written
 *       permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL MICHA SAMMETH BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

package barna.gui;

import barna.commons.Progressable;

import javax.swing.*;
import java.awt.*;

public class MyProgressBar extends JPanel implements Progressable {

    JProgressBar bar;
    JLabel lab;

    public MyProgressBar() {
        bar = new JProgressBar();
        lab = new JLabel(" "); // size
        setLayout(new FlowLayout(FlowLayout.LEFT));
        add(lab);
        add(bar);
        bar.setVisible(false);
    }

    public void progress() {
        Runnable runner = new Runnable() {
            public void run() {
                bar.setValue(bar.getValue() + 1);
                bar.setVisible(true);
                bar.paintImmediately(bar.getVisibleRect());    // getBounds()
            }
        };
        runFromEverywhere(runner);
    }

    public void clear() {
        setString(""); // was constant but hardcoded should have the same effect
        setValue(-1);
    }

    @Override
    public Dimension getPreferredSize() {
        int w = Math.max(lab.getPreferredSize().width, 250);
        int h = Math.max(80, bar.getPreferredSize().height);
        return super.getPreferredSize();    //new Dimension(w,h);
    }

    public void setMaximum(int newValue) {
        bar.setMaximum(newValue);
    }

    private void runFromEverywhere(Runnable runner) {
        if (EventQueue.isDispatchThread()) {
            runner.run();
        } else {
            try {
                SwingUtilities.invokeAndWait(runner);
            } catch (Exception e) {
                //System.err.println("everywhere disrupt");
                System.currentTimeMillis(); // :)
            }
        }
    }

    public void setString(String value) {
        //System.out.println("want set text "+value);
        final String val = value;
        Runnable runner = new Runnable() {
            public void run() {
                //System.out.println("set text "+val);
                lab.setText(val);
                lab.revalidate();
                lab.paintImmediately(lab.getVisibleRect());
                bar.setVisible(false);
                paintImmediately(getBounds());
            }
        };
        runFromEverywhere(runner);
    }

    public void setValue(final int newValue) {
        //System.out.println("want set val "+newValue);
        Runnable runner = new Runnable() {
            public void run() {
                //System.out.println("set val "+newValue);
                if (newValue < 0) {
                    bar.setVisible(false);
                } else {
                    bar.setValue(newValue);
                    bar.setVisible(true);
                    bar.paintImmediately(bar.getVisibleRect());
                }
            }
        };

        runFromEverywhere(runner);
    }

    public void finish() {
        bar.setVisible(false);
        bar.setValue(bar.getMinimum());
        lab.setText(" ");
    }

    public void finish(String msg, boolean time) {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    public void failed(String msg) {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    public int steps() {
        return 0;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public int currentStep() {
        return 0;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public void message(String value) {
        final String val = value;
        Runnable runner = new Runnable() {
            public void run() {
                lab.setText(val);
                lab.paintImmediately(lab.getBounds());
            }
        };
        runFromEverywhere(runner);
    }

    public void finish(String msg, long time) {
        // TODO Auto-generated method stub
    }

    public void start(String message) {
        throw new UnsupportedOperationException("Not implemented yet");
    }
}
