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

package barna.model.commons;

import java.io.IOException;
import java.io.InputStream;
/**  
 * Ein Thread, der eine InputStream ausliest und die erhaltenen Daten ins Nirvana schickt   
 * (bzw. nach /dev/null f??r den Linux-user ;-)). Eigentlich werden sie nur dem GarbageCollector  
 * ?berlassen.  
 * <br><br>  
 * Wird ein Name f??r den Stream ??bergeben (z.B. "standard_out"), so wird der Inhalt des Streams  
 * unter Angabe des Namens nach System.out geschrieben.<br>  
 * Der Thread beendet sich selbst, sobald der stream geschlossen wird oder ein .interrupt()   
 * aufgerufen wird.<br>  
 * Erstellungsdatum: (03.04.01 21:55:05)  
 * @author Tobias Vogele   
 */  
public class DevNullReaderThread extends Thread {  

// needed for run()
// inherited from interface Constants
protected boolean DEBUG= false;

/**  
 * Der Inputstream, der ausgelesen wird.  
 */  
private java.io.InputStream in;  
/**  
 * Falls ein Name f?r den Stream angegeben wird, wird er hier gespeichert und mit dem Inhalt  
 * des Streams ausgegeben.  
 */  
private java.lang.String streamName = null;  
/**  
 * DevNullThread - leerer Konstruktor.  
 */  
public DevNullReaderThread() {  
        super("DevNullReader");  
}  
/**  
 * DevNullThread - erstellt einen neunen DevNullReader mit dem ?bergebenen Stream ohne Name.  
 * @param in input stream
 */
public DevNullReaderThread(InputStream in) {  
        this(in, null);  
}  
/**  
 * DevNullThread - erstellt einen neuen DevNullReader mit den ?bergebenen Stream und Name.
 * @param in input stream
 * @param name name
 */  
public DevNullReaderThread(InputStream in, String name) {  
        super("DevNullReader");  
        this.in = in;  
        streamName = name;  
}  
/**  
 * Gibt den Inputstream zur??ck, von dem gelesen werden soll.<br>  
 * Erstellungsdatum: (03.04.01 21:55:55)  
 * @return java.io.InputStream  
 */  
public java.io.InputStream getIn() {  
        return in;  
}  
/**  
 * Gibt den Name des Inputstreams zur??ck.<br>  
 * Erstellungsdatum: (03.04.01 21:57:24)  
 * @return java.lang.String  
 */  
public java.lang.String getStreamName() {  
        return streamName;  
}  
/**  
 * Die Hauptmethode.  
 * Hier wird immerfort der Stream ausgelesen, solange was da ist zum lesen und solange  
 * keine interrupt() aufgerufen wurde.<br>  
 * Erstellungsdatum: (03.04.01 21:56:13)  
 */  
public void run() {  
                if (DEBUG && streamName != null)
					System.out.println(streamName+"+ >");   

// evtl. std-out und std-err auslesen.  
                int w;  
                do {  
                        try {  
                                w = in.read();  
                        } catch (IOException e) {  
                                w = -1;  
                        }  
                        if (w != -1 && DEBUG && streamName != null)  
System.out.print((char)   
w);  
                } 
                while (w != -1 && ! isInterrupted());  
                if (DEBUG && streamName != null) System.out.println("\n<"); 

  
}  
/**  
 * Setzt den Inputstream.<br>  
 * Erstellungsdatum: (03.04.01 21:55:55)  
 * @param newIn java.io.InputStream  
 */  
public void setIn(java.io.InputStream newIn) {  
        in = newIn;  
}  
/**  
 * Setzt den Name des Streams.<br>  
 * Erstellungsdatum: (03.04.01 21:57:24)  
 * @param newStreamName java.lang.String  
 */  
public void setStreamName(java.lang.String newStreamName) {  
        streamName = newStreamName;  
}  
}  

